##############################################################################
# Soybean seed longevity
# Stage 4c: sample connectivity QC (Z.K) and body part x BioProject confounding
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: (1) flag curated samples with abnormally low global connectivity using
# the logic of BioNERO::ZKfiltering() (Oldham et al., 2012):
#   A = WGCNA adjacency of type "distance" on Spearman correlation,
#   k = colSums(A) - 1, Z.K = scale(k), outlier if Z.K < -2;
# (2) check whether a global outlier is also an outlier within its own body
# part, so that a rare tissue is not removed just for being different;
# (3) quantify how much of the PCA structure goes with body part, with
# BioProject, or with both (they are partly confounded). Nothing is removed
# here; the sample sets are defined in script 04.
#
# Run from the repository root, after script 02 of this stage:
#   source("scripts/04_coexpression_network/03_connectivity_qc_confounding.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2, matrixStats, WGCNA, BioNERO
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2", "matrixStats"))
for (p in c("WGCNA", "BioNERO")) if (!requireNamespace(p, quietly = TRUE)) stop("Install first: ", p)

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_expr <- file.path("data", "processed", "SEA_seedrelated_gene_TPM.tsv.gz")
f_keep <- file.path(dirs$res, "05_curated_samples.txt")
f_dec <- file.path(dirs$res, "03_sample_curation_decisions.csv")
check_inputs(c(f_expr, f_keep, f_dec))

keep <- readr::read_lines(f_keep)
meta <- readr::read_csv(f_dec, show_col_types = FALSE) %>%
  dplyr::filter(sample_id %in% keep) %>%
  dplyr::mutate(BioProject = as.character(BioProject))
mat <- as_expr_matrix(read_table_df(f_expr))[, keep, drop = FALSE]

# Prevalence filter used for the QC: TPM >= 1 in >= 25% of the curated samples
mat_qc <- mat[rowSums(mat >= NET$min_tpm) >= NET$min_prop * ncol(mat), , drop = FALSE]
nrow(mat_qc)
log_qc <- log2(mat_qc + 1)

# ---------------------------------------------------------------------------
# 1) Z.K, global and within body part
# ---------------------------------------------------------------------------

# Standardized sample connectivity, as in BioNERO::ZKfiltering()
zk_scores <- function(m) {
  A <- WGCNA::adjacency(as.data.frame(m), type = "distance", corOptions = list(use = "p", method = "spearman"))
  k <- colSums(A) - 1
  setNames(as.numeric(scale(k)), colnames(m))
}
zk_global <- zk_scores(mat_qc)

# Within body part only when the group has >= 15 samples (endosperm and
# suspensor are too small to be judged this way)
zk_within <- unlist(lapply(split(meta$sample_id, meta$body_part), function(ids) {
  if (length(ids) < 15) setNames(rep(NA_real_, length(ids)), ids) else zk_scores(mat_qc[, ids, drop = FALSE])
}), use.names = FALSE)
names(zk_within) <- unlist(split(meta$sample_id, meta$body_part), use.names = FALSE)

zk <- meta %>%
  dplyr::select(sample_id, body_part, BioProject, cultivar_raw, Treatment, Study_title) %>%
  dplyr::mutate(ZK_global = zk_global[sample_id], ZK_within_bodypart = zk_within[sample_id],
                outlier_global = ZK_global < NET$zk_cut,
                ZK_interpretation = dplyr::case_when(
                  outlier_global & !is.na(ZK_within_bodypart) & ZK_within_bodypart < NET$zk_cut ~
                    "global_and_within_bodypart_outlier",
                  outlier_global & !is.na(ZK_within_bodypart) ~ "global_only_possible_biological_structure",
                  outlier_global ~ "global_outlier_small_bodypart_review_manually",
                  TRUE ~ "not_global_outlier")) %>%
  dplyr::arrange(ZK_global)
readr::write_csv(zk, file.path(dirs$res, "09_sample_connectivity_ZK.csv"))

# Validation against BioNERO itself and sensitivity to the cutoff
removed_bionero <- setdiff(colnames(mat_qc), colnames(BioNERO::ZKfiltering(as.data.frame(mat_qc), zk = NET$zk_cut,
                                                                             cor_method = "spearman")))
(zk_check <- tibble::tibble(
  metric = c("samples", "outliers_ZK_lt_-2", "outliers_ZK_lt_-2.5", "outliers_ZK_lt_-3", "identical_to_BioNERO"),
  value = c(ncol(mat_qc), sum(zk$ZK_global < -2), sum(zk$ZK_global < -2.5), sum(zk$ZK_global < -3),
            setequal(removed_bionero, zk$sample_id[zk$outlier_global]))))
readr::write_csv(zk_check, file.path(dirs$res, "10_ZK_threshold_check.csv"))
# expected: 235 samples, 4 outliers at -2, 1 at -3, identical = TRUE (1)

p_zk <- ggplot2::ggplot(zk, ggplot2::aes(body_part, ZK_global)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.75, size = 1.8) +
  ggplot2::geom_hline(yintercept = NET$zk_cut, linetype = "dashed") +
  ggplot2::labs(x = "Body part", y = "Global Z.K") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(dirs$fig, "04_ZK_by_bodypart.png"), p_zk, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------------------------
# 2) BioProject x body part
# ---------------------------------------------------------------------------
cross <- meta %>%
  dplyr::count(BioProject, body_part, name = "n_samples") %>%
  tidyr::complete(BioProject, body_part, fill = list(n_samples = 0))
purity <- cross %>%
  dplyr::group_by(BioProject) %>%
  # counts per body part are summarised BEFORE the total, so that the new
  # column does not overwrite n_samples inside summarise
  dplyr::summarise(n_bodyparts = sum(n_samples > 0), dominant_bodypart = body_part[which.max(n_samples)],
                   dominant_fraction = max(n_samples) / sum(n_samples), n_total = sum(n_samples), .groups = "drop") %>%
  dplyr::rename(n_samples = n_total) %>%
  dplyr::arrange(dplyr::desc(n_samples))
readr::write_csv(purity, file.path(dirs$res, "11_BioProject_bodypart_purity.csv"))
(purity_summary <- tibble::tibble(
  metric = c("n_BioProjects", "n_single_bodypart", "prop_samples_in_single_bodypart_projects"),
  value = c(nrow(purity), sum(purity$n_bodyparts == 1), sum(purity$n_samples[purity$n_bodyparts == 1]) / nrow(meta))))
readr::write_csv(purity_summary, file.path(dirs$res, "12_BioProject_bodypart_confounding_summary.csv"))

# ---------------------------------------------------------------------------
# 3) PC1-PC10 against body part and BioProject (commonality decomposition)
#   unique_body    = R2(full) - R2(project only)
#   unique_project = R2(full) - R2(body only)
#   shared         = R2(body only) + R2(project only) - R2(full)
# A large shared part is direct evidence of confounding; the two R2 are not
# independent causal effects.
# ---------------------------------------------------------------------------
top <- head(names(sort(matrixStats::rowVars(log_qc, useNames = TRUE), decreasing = TRUE)), 2000)
pca <- prcomp(t(log_qc[top, meta$sample_id]), center = TRUE, scale. = FALSE)
var_frac <- (pca$sdev^2 / sum(pca$sdev^2))[1:10]

# R2 of body part, BioProject and both for one PC
pc_r2 <- function(i) {
  d <- data.frame(y = pca$x[, i], body = factor(meta$body_part), proj = factor(meta$BioProject))
  r2 <- c(body = summary(lm(y ~ body, d))$r.squared, proj = summary(lm(y ~ proj, d))$r.squared,
          full = summary(lm(y ~ body + proj, d))$r.squared)
  tibble::tibble(PC = paste0("PC", i), variance_fraction = var_frac[i], R2_bodypart = r2[["body"]],
                 R2_BioProject = r2[["proj"]], R2_full = r2[["full"]], unique_bodypart = r2[["full"]] - r2[["proj"]],
                 unique_BioProject = r2[["full"]] - r2[["body"]], shared = r2[["body"]] + r2[["proj"]] - r2[["full"]])
}
pc_assoc <- dplyr::bind_rows(lapply(1:10, pc_r2))
readr::write_csv(pc_assoc, file.path(dirs$res, "13_PC_bodypart_BioProject_R2.csv"))
(pc_weighted <- pc_assoc %>%
    dplyr::summarise(dplyr::across(R2_bodypart:shared, ~ weighted.mean(.x, variance_fraction))) %>%
    dplyr::mutate(variance_PC1_PC10 = sum(pc_assoc$variance_fraction)))
readr::write_csv(pc_weighted, file.path(dirs$res, "14_PC_confounding_weighted_summary.csv"))

# Reading: body part is the biological signal we want the modules to carry.
# When BioProject explains more than body part AND the shared part is large,
# most projects sample a single tissue: regressing BioProject out would also
# remove anatomy. In that case the network is built without batch correction
# and every module-tissue association is re-checked within BioProject
# (stage 6, script 05). Expected here: R2 BioProject ~0.68, body part ~0.40,
# shared ~0.27.

close_stage("04_coexpression_network_zk", dirs)
