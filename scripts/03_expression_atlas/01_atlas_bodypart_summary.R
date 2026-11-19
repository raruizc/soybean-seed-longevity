##############################################################################
# Soybean seed longevity
# Stage 3a: Soybean Expression Atlas v2, expression of the 39 candidates by body part
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: add an expression-context layer that is independent of GSE98199. The
# Soybean Expression Atlas v2 (Almeida-Silva et al., 2023) gathers public
# RNA-seq samples quantified on Wm82.a4.v1 IDs. For each candidate and body
# part we summarize median, mean, quartiles and prevalence of TPM, compute Tau
# on the medians and reproduce the atlas classification (null, weak, broad,
# part-specific).
# The median per body part is a robust summary but the atlas is heterogeneous;
# sample- and study-level spread is examined in script 02.
#
# Manual download (once): run section 0, then in the atlas "Search gene list"
# paste the IDs of results/03_expression_atlas/00_candidate_ids_Wm82a4.txt,
# select ALL body parts and TPM, and save the expression matrix and the sample
# metadata as data/raw/atlas/SEA_39genes_TPM.tsv and SEA_39genes_metadata.tsv.
# All body parts are kept so that "seed expression" has a comparison group.
#
# Run from the repository root:
#   source("scripts/03_expression_atlas/01_atlas_bodypart_summary.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))

stage <- "03_expression_atlas"
dirs <- setup_stage(stage)
f_cand <- from_stage("01_candidates", "01_candidates_39.csv")
check_inputs(f_cand)

# ---------------------------------------------------------------------------
# 0) Candidate IDs for the manual atlas download
# ---------------------------------------------------------------------------
cand_ids <- sort(readr::read_csv(f_cand, show_col_types = FALSE)$gene_id)
stopifnot(length(cand_ids) == 39)
readr::write_lines(cand_ids, file.path(dirs$res, "00_candidate_ids_Wm82a4.txt"))
check_inputs(c(RAW$atlas_tpm, RAW$atlas_meta))

# ---------------------------------------------------------------------------
# 1) Expression matrix, metadata and audit
# Sample ID in the metadata = BioSample; body part = Part.
# ---------------------------------------------------------------------------
expr <- read.delim(RAW$atlas_tpm, check.names = FALSE, stringsAsFactors = FALSE)
names(expr)[1] <- "gene_id"
sample_cols <- setdiff(names(expr), "gene_id")
expr[sample_cols] <- lapply(expr[sample_cols], as.numeric)
stopifnot(!anyNA(expr[sample_cols]))

meta <- read.delim(RAW$atlas_meta, check.names = FALSE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(sample_id = as.character(BioSample), body_part = as.character(Part)) %>%
  dplyr::distinct()

(audit <- tibble::tibble(
  metric = c("candidate_ids_expected", "candidate_ids_found", "expression_samples",
             "samples_with_metadata", "samples_without_metadata", "body_parts"),
  n = c(length(cand_ids), sum(cand_ids %in% expr$gene_id), length(sample_cols),
        sum(sample_cols %in% meta$sample_id), sum(!sample_cols %in% meta$sample_id),
        dplyr::n_distinct(meta$body_part[meta$sample_id %in% sample_cols]))))
readr::write_csv(audit, file.path(dirs$res, "01_atlas_audit.csv"))
if (any(!sample_cols %in% meta$sample_id)) stop("Matrix samples without metadata")

# Long format, kept in data/processed for script 02
expr_long <- expr %>%
  dplyr::filter(gene_id %in% cand_ids) %>%
  tidyr::pivot_longer(dplyr::all_of(sample_cols), names_to = "sample_id", values_to = "TPM") %>%
  dplyr::inner_join(meta, by = "sample_id")
readr::write_csv(expr_long, file.path("data", "processed", "SEA_39genes_TPM_long.csv"))

(part_counts <- expr_long %>%
    dplyr::distinct(sample_id, body_part) %>%
    dplyr::count(body_part, name = "n_samples") %>%
    dplyr::arrange(dplyr::desc(n_samples)))
readr::write_csv(part_counts, file.path(dirs$res, "02_bodypart_sample_counts.csv"))

# ---------------------------------------------------------------------------
# 2) Summary per gene x body part
# ---------------------------------------------------------------------------
part_summary <- expr_long %>%
  dplyr::group_by(gene_id, body_part) %>%
  dplyr::summarise(n_samples = dplyr::n(), median_TPM = median(TPM), mean_TPM = mean(TPM),
                   q25_TPM = quantile(TPM, 0.25), q75_TPM = quantile(TPM, 0.75),
                   prop_TPM_ge_1 = mean(TPM >= 1), prop_TPM_ge_5 = mean(TPM >= 5), .groups = "drop")
readr::write_csv(part_summary, file.path(dirs$res, "03_candidate_bodypart_summary.csv"))

# ---------------------------------------------------------------------------
# 3) Tau and atlas classification
# Tau (Yanai et al., 2005) on log2(median TPM + 1) across body parts:
#   tau = sum(1 - x_i / max(x)) / (n - 1), 0 = uniform, 1 = one part only.
# Classes as in the atlas: max median < 1 null; < 5 weak; otherwise broad
# (tau < 0.85) or part-specific (tau >= 0.85).
# ---------------------------------------------------------------------------
tau_index <- function(x) {
  x <- log2(as.numeric(x) + 1)
  xmax <- max(x)
  if (!is.finite(xmax) || xmax <= 0) return(NA_real_)
  sum(1 - x / xmax) / (length(x) - 1)
}
tau_index(c(10, 0, 0, 0))                      # check: 1

med_wide <- part_summary %>%
  dplyr::select(gene_id, body_part, median_TPM) %>%
  tidyr::pivot_wider(names_from = body_part, values_from = median_TPM, values_fill = 0)
med_mat <- as.matrix(med_wide[, -1])

# Body parts with median > 5, used to name where part-specific genes are expressed
specific_parts <- part_summary %>%
  dplyr::filter(median_TPM > 5) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(parts_median_gt5 = paste(body_part[order(median_TPM, decreasing = TRUE)], collapse = "; "))

classification <- tibble::tibble(gene_id = med_wide$gene_id, max_median_TPM = apply(med_mat, 1, max),
                                 tau = apply(med_mat, 1, tau_index)) %>%
  dplyr::mutate(expression_group = dplyr::case_when(
    max_median_TPM < 1 ~ "null",
    max_median_TPM < 5 ~ "weak",
    !is.na(tau) & tau < 0.85 ~ "broad",
    !is.na(tau) & tau >= 0.85 ~ "part-specific",
    TRUE ~ "unclassified")) %>%
  dplyr::left_join(specific_parts, by = "gene_id") %>%
  dplyr::mutate(specific_parts = dplyr::if_else(expression_group == "part-specific", parts_median_gt5, NA_character_)) %>%
  dplyr::select(-parts_median_gt5) %>%
  dplyr::arrange(dplyr::desc(tau))
readr::write_csv(classification, file.path(dirs$res, "04_gene_classification_tau.csv"))
print(dplyr::count(classification, expression_group))   # expected: 29 broad, 4 weak, 3 null, 3 part-specific

p_heat <- ggplot2::ggplot(part_summary, ggplot2::aes(body_part, gene_id, fill = log2(median_TPM + 1))) +
  ggplot2::geom_tile() +
  ggplot2::labs(x = "Body part", y = "Candidate gene (Wm82.a4.v1)", fill = "log2(median TPM + 1)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 6))
ggplot2::ggsave(file.path(dirs$fig, "01_heatmap_median_TPM_bodypart.png"), p_heat, width = 12, height = 11, dpi = 300)

close_stage("03_expression_atlas_summary", dirs)
