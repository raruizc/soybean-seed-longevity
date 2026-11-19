##############################################################################
# Soybean seed longevity
# Stage 4d: network samples, gene universe and soft-thresholding power
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: define the samples (curated, Z.K < -2 removed) and the genes of the
# network, and choose the soft-thresholding power beta. Genes must pass the
# prevalence filter (TPM >= 1 in >= 25% of samples); then a grid keeps the
# 100, 90, 80, 75, 70, 60 or 50% most variable genes (variance of log2(TPM+1))
# and BioNERO::SFT_fit() finds, for each set, the first beta with scale-free
# fit R2 >= 0.80 in a signed hybrid network (Pearson).
# The decision weighs three things: keeping genes and candidates, reaching
# R2 >= 0.80, and avoiding a very high beta with very low mean connectivity.
# The chosen set and beta are frozen in R/config.R (NET$var_keep, NET$power).
# A second check repeats the SFT keeping the samples with -3 < Z.K < -2, to
# see whether beta depends on how outliers were treated.
#
# Each SFT run over ~20-30 thousand genes takes several minutes; diagnostics
# already in results/ are reused instead of being recomputed.
#
# Run from the repository root, after script 03 of this stage:
#   source("scripts/04_coexpression_network/04_gene_filter_soft_threshold.R")
#
# Packages: readr, dplyr, tibble, ggplot2, matrixStats, WGCNA, BioNERO
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tibble", "ggplot2", "matrixStats"))
for (p in c("WGCNA", "BioNERO")) if (!requireNamespace(p, quietly = TRUE)) stop("Install first: ", p)

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_expr <- file.path("data", "processed", "SEA_seedrelated_gene_TPM.tsv.gz")
f_zk <- file.path(dirs$res, "09_sample_connectivity_ZK.csv")
f_cand <- from_stage("01_candidates", "01_candidates_39.csv")
check_inputs(c(f_expr, f_zk, f_cand))

cand_ids <- readr::read_csv(f_cand, show_col_types = FALSE)$gene_id
zk <- readr::read_csv(f_zk, show_col_types = FALSE)
mat <- as_expr_matrix(read_table_df(f_expr))

# ---------------------------------------------------------------------------
# 1) Sample sets: primary (Z.K < -2 removed) and sensitivity (only Z.K < -3)
# ---------------------------------------------------------------------------
samples_primary <- zk$sample_id[zk$ZK_global >= NET$zk_cut]
samples_sens <- zk$sample_id[zk$ZK_global >= -3]
c(curated = nrow(zk), primary = length(samples_primary), sensitivity = length(samples_sens))   # expected: 235, 231, 234
readr::write_lines(samples_primary, file.path(dirs$res, "15_network_samples_primary.txt"))

# ---------------------------------------------------------------------------
# 2) Prevalence filter and variance of each gene (primary samples)
# ---------------------------------------------------------------------------
m_prim <- mat[, samples_primary, drop = FALSE]
prev <- rowMeans(m_prim >= NET$min_tpm)
log_prev <- log2(m_prim[prev >= NET$min_prop, , drop = FALSE] + 1)
nrow(log_prev)                                 # expected: 31,531 genes
gene_var <- matrixStats::rowVars(log_prev, useNames = TRUE)
gene_diag <- tibble::tibble(gene_id = rownames(log_prev), prop_TPM_ge_1 = prev[rownames(log_prev)],
                            variance_log2TPM = gene_var, MAD_log2TPM = matrixStats::rowMads(log_prev),
                            variance_percentile = rank(gene_var) / length(gene_var),
                            is_candidate = rownames(log_prev) %in% cand_ids) %>%
  dplyr::arrange(dplyr::desc(variance_log2TPM))
readr::write_csv(gene_diag, file.path(dirs$res, "16_gene_variance_diagnostics.csv"))
saveRDS(log_prev, file.path("data", "processed", "network_prevalence_log2TPM.rds"))

# ---------------------------------------------------------------------------
# 3) Variance grid and scale-free fit
# ---------------------------------------------------------------------------
FRACTIONS <- c(all_prevalence = 1, top90 = 0.90, top80 = 0.80, top75 = 0.75, top70 = 0.70, top60 = 0.60, top50 = 0.50)

# Genes kept when retaining the given fraction of most variable genes
variance_set <- function(frac) {
  if (frac >= 1) return(names(gene_var))
  names(gene_var)[gene_var >= quantile(gene_var, probs = 1 - frac, names = FALSE, type = 7)]
}

# SFT diagnostics (power, fit, meank) for one gene set; reuses a saved file
sft_diagnostics <- function(m, label) {
  f <- file.path(dirs$res, paste0("SFT_", label, ".csv"))
  if (file.exists(f)) return(readr::read_csv(f, show_col_types = FALSE))
  stopifnot(WGCNA::goodSamplesGenes(t(m), verbose = 0)$allOK)
  sft <- BioNERO::SFT_fit(as.data.frame(m), net_type = NET$type, rsquared = 0.8, cor_method = NET$cor)
  d <- dplyr::select(as.data.frame(sft$plot$data), power, fit, meank)
  readr::write_csv(d, f)
  d
}

# First power with R2 >= 0.80; if never reached, the largest power tested
pick_power <- function(d) {
  d <- dplyr::arrange(d, power)
  hit <- d[d$fit >= 0.80, ]
  row <- if (nrow(hit)) hit[1, ] else d[nrow(d), ]
  tibble::tibble(target_reached = nrow(hit) > 0, power = row$power, R2 = row$fit, mean_connectivity = row$meank,
                 max_R2 = max(d$fit))
}

grid <- dplyr::bind_rows(lapply(names(FRACTIONS), function(nm) {
  ids <- variance_set(FRACTIONS[[nm]])
  d <- sft_diagnostics(log_prev[ids, , drop = FALSE], nm)
  tibble::tibble(gene_set = nm, fraction = FRACTIONS[[nm]], n_genes = length(ids),
                 n_candidates = sum(cand_ids %in% ids)) %>%
    dplyr::bind_cols(pick_power(d))
}))
readr::write_csv(grid, file.path(dirs$res, "17_soft_threshold_grid.csv"))
print(grid)

# Reading: expected top75 = 23,648 genes, 18 candidates, beta 16 (R2 0.80, mean
# k ~ 45). Keeping all genes needs beta 20 with mean k ~ 22; top50 reaches the
# fit at beta 8 but keeps only 10 of the 32 prevalence-eligible candidates.

# Outlier sensitivity: same genes, samples with Z.K >= -3
sens <- log2(mat[rownames(log_prev), samples_sens, drop = FALSE] + 1)
(sft_sens <- pick_power(sft_diagnostics(sens, "sensitivity_all_prevalence")) %>%
    dplyr::mutate(primary_power = grid$power[grid$gene_set == "all_prevalence"]))
readr::write_csv(sft_sens, file.path(dirs$res, "18_soft_threshold_outlier_sensitivity.csv"))
# Same beta in both = the choice of beta does not depend on the Z.K cutoff.

# ---------------------------------------------------------------------------
# 4) Candidates along the grid and the frozen network input
# ---------------------------------------------------------------------------
cand_grid <- tibble::tibble(gene_id = cand_ids) %>%
  dplyr::left_join(dplyr::select(gene_diag, gene_id, prop_TPM_ge_1, variance_log2TPM, variance_percentile), by = "gene_id") %>%
  dplyr::mutate(pass_prevalence = gene_id %in% rownames(log_prev))
for (nm in names(FRACTIONS)) cand_grid[[paste0("retained_", nm)]] <- cand_grid$gene_id %in% variance_set(FRACTIONS[[nm]])
readr::write_csv(cand_grid, file.path(dirs$res, "19_candidate_retention_grid.csv"))

net_genes <- variance_set(NET$var_keep)
readr::write_lines(net_genes, file.path(dirs$res, "20_network_genes.txt"))
p_var <- ggplot2::ggplot(gene_diag, ggplot2::aes(variance_log2TPM)) +
  ggplot2::geom_histogram(bins = 60) +
  ggplot2::geom_vline(xintercept = quantile(gene_var, 1 - NET$var_keep), linetype = "dashed") +
  ggplot2::geom_rug(data = dplyr::filter(gene_diag, is_candidate)) +
  ggplot2::labs(x = "Variance of log2(TPM + 1)", y = "Genes", subtitle = "Rug: candidates; dashed: network cutoff") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "05_gene_variance_candidates.png"), p_var, width = 10, height = 7, dpi = 300)

close_stage("04_coexpression_network_sft", dirs)
