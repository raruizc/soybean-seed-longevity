##############################################################################
# Soybean seed longevity
# Stage 4b: audit and curation of the seed-related atlas samples
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: avoid a network dominated by cultivar, treatment or study design. The
# primary network keeps only Williams 82 samples from studies of normal seed
# development, and excludes transgenic/mutant, stress/pathogen, germination
# and somatic embryogenesis studies, and explicit non-neutral treatments.
# Nothing is deleted from the raw files: every decision and its reason is
# written to 03_sample_curation_decisions.csv. Expression filters and PCA are
# recomputed after curation.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/04_coexpression_network/02_sample_audit_curation.R")
#
# Packages: readr, dplyr, tibble, stringr, ggplot2, matrixStats
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tibble", "stringr", "ggplot2", "matrixStats"))

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_expr <- file.path("data", "processed", "SEA_seedrelated_gene_TPM.tsv.gz")
f_meta <- file.path("data", "processed", "SEA_seedrelated_metadata.tsv.gz")
f_cand <- from_stage("01_candidates", "01_candidates_39.csv")
check_inputs(c(f_expr, f_meta, f_cand))

cand_ids <- readr::read_csv(f_cand, show_col_types = FALSE)$gene_id
mat <- as_expr_matrix(read_table_df(f_expr))
dim(mat)                                       # expected: 52,837 genes x 933 samples
meta <- read_table_df(f_meta) %>%
  dplyr::mutate(sample_id = as.character(BioSample), body_part = as.character(download_body_part),
                cultivar_raw = as.character(Cultivar),
                study_title_lc = stringr::str_to_lower(dplyr::coalesce(as.character(Study_title), "")),
                treatment_lc = stringr::str_to_lower(dplyr::coalesce(as.character(Treatment), ""))) %>%
  dplyr::filter(sample_id %in% colnames(mat)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)
stopifnot(setequal(meta$sample_id, colnames(mat)))

# ---------------------------------------------------------------------------
# 1) Composition before curation
# ---------------------------------------------------------------------------
composition <- dplyr::bind_rows(lapply(c("body_part", "Cultivar", "BioProject", "Treatment"), function(v) {
  meta %>%
    dplyr::count(level = as.character(.data[[v]]), name = "n_samples") %>%
    dplyr::mutate(variable = v) %>%
    dplyr::arrange(dplyr::desc(n_samples))
})) %>%
  dplyr::select(variable, level, n_samples)
readr::write_csv(composition, file.path(dirs$res, "01_sample_composition_before_curation.csv"))

# Keyword scan of titles, abstracts and treatments: for MANUAL REVIEW only
keywords <- paste(c("stress", "drought", "salt", "salinity", "heat", "cold", "pathogen", "infect", "fung",
                    "bacter", "virus", "mutant", "transgenic", "overexpress", "rna[- ]?i", "crispr",
                    "knockout", "knock[- ]?down", "treatment", "treated"), collapse = "|")
scan_text <- tolower(paste(meta$Treatment, meta$Cultivar, meta$Study_title, meta$Study_abstract, sep = " | "))
readr::write_csv(meta[stringr::str_detect(scan_text, keywords), c("sample_id", "BioProject", "body_part", "Cultivar",
                                                                  "Treatment", "Study_title")],
                 file.path(dirs$res, "02_keyword_flagged_samples_review_only.csv"))

# ---------------------------------------------------------------------------
# 2) Curation rules for the PRIMARY network
# "Williams" and "Williams 43" are deliberately NOT treated as Williams 82.
# ---------------------------------------------------------------------------
W82_ALIASES <- c("Williams 82", "Williams82", "William 82", "strain Williams 82",
                 "/genotypeWilliams (PI548631)", "Wm 82")
NEUTRAL_TREATMENTS <- c("", "fresh", "dry", "dried", "control", "control1", "control3", "control5")
title_has <- function(terms) stringr::str_detect(meta$study_title_lc, paste(terms, collapse = "|"))

meta <- meta %>%
  dplyr::mutate(
    is_Williams82 = cultivar_raw %in% W82_ALIASES,
    flag_transgenic_mutant = title_has(c("transgenic", "genetically modified", "gm and non-gm", "mutant",
                                         "overexpression", "over-expression", "low phytic acid",
                                         "defective seed coat mutation")),
    flag_stress_pathogen = title_has(c("drought", "cold stress", "heat stress", "field mold stress",
                                       "phytophthora", "biostimulant")),
    flag_germination = title_has(c("germinat", "imbib", "paclobutrazol", "seedling")),
    flag_somatic_embryogenesis = title_has("somatic embryogenesis"),
    flag_nonstandard_treatment = !treatment_lc %in% NEUTRAL_TREATMENTS,
    exclusion_reason = dplyr::case_when(
      !is_Williams82 ~ "non_Williams82_or_unknown_cultivar",
      flag_transgenic_mutant ~ "transgenic_mutant_or_low_phytate_study",
      flag_stress_pathogen ~ "stress_pathogen_or_biostimulant_study",
      flag_germination ~ "germination_imbibition_or_seedling_study",
      flag_somatic_embryogenesis ~ "somatic_embryogenesis_study",
      flag_nonstandard_treatment ~ "explicit_nonstandard_treatment",
      TRUE ~ "KEEP_PRIMARY_NETWORK"),
    keep_primary_network = exclusion_reason == "KEEP_PRIMARY_NETWORK")

readr::write_csv(dplyr::select(meta, sample_id, BioProject, body_part, cultivar_raw, Treatment, Study_title,
                               keep_primary_network, exclusion_reason, dplyr::starts_with("flag_")),
                 file.path(dirs$res, "03_sample_curation_decisions.csv"))
(curation <- dplyr::count(meta, keep_primary_network, exclusion_reason, name = "n_samples") %>%
    dplyr::arrange(dplyr::desc(keep_primary_network), dplyr::desc(n_samples)))
readr::write_csv(curation, file.path(dirs$res, "04_sample_curation_summary.csv"))
# expected: 235 kept; 636 non-Williams 82, 36 transgenic/mutant, 26 germination

keep <- meta$sample_id[meta$keep_primary_network]
readr::write_lines(keep, file.path(dirs$res, "05_curated_samples.txt"))
curated_counts <- dplyr::bind_rows(
  dplyr::count(dplyr::filter(meta, keep_primary_network), level = body_part, name = "n_samples") %>%
    dplyr::mutate(variable = "body_part"),
  dplyr::count(dplyr::filter(meta, keep_primary_network), level = as.character(BioProject), name = "n_samples") %>%
    dplyr::mutate(variable = "BioProject")) %>%
  dplyr::select(variable, level, n_samples) %>%
  dplyr::arrange(variable, dplyr::desc(n_samples))
readr::write_csv(curated_counts, file.path(dirs$res, "06_curated_sample_counts.csv"))

# ---------------------------------------------------------------------------
# 3) Gene expression filters after curation
# ---------------------------------------------------------------------------
mat_cur <- mat[, keep, drop = FALSE]
log_cur <- log2(mat_cur + 1)
gene_diag <- tibble::tibble(gene_id = rownames(mat_cur), median_TPM = matrixStats::rowMedians(mat_cur),
                            mean_TPM = rowMeans(mat_cur), prop_TPM_ge_1 = rowMeans(mat_cur >= 1),
                            prop_TPM_ge_5 = rowMeans(mat_cur >= 5), variance_log2TPM = matrixStats::rowVars(log_cur),
                            is_candidate = rownames(mat_cur) %in% cand_ids) %>%
  dplyr::mutate(pass_median_TPM_ge_1 = median_TPM >= 1, pass_prevalence = prop_TPM_ge_1 >= NET$min_prop)
readr::write_csv(dplyr::filter(gene_diag, is_candidate) %>% dplyr::arrange(dplyr::desc(variance_log2TPM)),
                 file.path(dirs$res, "07_candidate_eligibility_curated.csv"))
(filter_diag <- tibble::tibble(
  criterion = c("median_TPM_ge_1", "TPM_ge_1_in_25pct_samples", "candidates_median_TPM_ge_1",
                "candidates_TPM_ge_1_in_25pct_samples"),
  n_genes = c(sum(gene_diag$pass_median_TPM_ge_1), sum(gene_diag$pass_prevalence),
              sum(gene_diag$is_candidate & gene_diag$pass_median_TPM_ge_1),
              sum(gene_diag$is_candidate & gene_diag$pass_prevalence))))
readr::write_csv(filter_diag, file.path(dirs$res, "08_expression_filter_diagnostics.csv"))

# ---------------------------------------------------------------------------
# 4) PCA of the 2,000 most variable genes (unscaled; scaled as sensitivity)
# Reading: body parts should separate. If samples group by BioProject instead,
# study design competes with anatomy (quantified in script 03).
# ---------------------------------------------------------------------------
top <- head(gene_diag$gene_id[order(-gene_diag$variance_log2TPM)], 2000)
pca_df <- function(scale) {
  x <- log_cur[top, , drop = FALSE]
  if (scale) x <- x[matrixStats::rowVars(x) > 0, , drop = FALSE]
  pca <- prcomp(t(x), center = TRUE, scale. = scale)
  pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
  tibble::as_tibble(pca$x[, 1:2], rownames = "sample_id") %>%
    dplyr::left_join(dplyr::select(meta, sample_id, body_part, BioProject), by = "sample_id") %>%
    dplyr::mutate(xlab = paste0("PC1: ", pct[1], "%"), ylab = paste0("PC2: ", pct[2], "%"))
}

# Scatter of PC1 x PC2 colored by one variable
plot_pca <- function(d, color, file) {
  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = .data[[color]])) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::labs(x = d$xlab[1], y = d$ylab[1], color = color) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(dirs$fig, file), p, width = 9, height = 7, dpi = 300)
}
d_pca <- pca_df(scale = FALSE)
plot_pca(d_pca, "body_part", "01_PCA_curated_by_bodypart.png")
plot_pca(pca_df(scale = TRUE), "body_part", "02_PCA_curated_scaled_by_bodypart.png")
top_proj <- head(dplyr::filter(curated_counts, variable == "BioProject")$level, 12)
plot_pca(dplyr::mutate(d_pca, BioProject = ifelse(BioProject %in% top_proj, BioProject, "Other")),
         "BioProject", "03_PCA_curated_by_BioProject.png")

close_stage("04_coexpression_network_curation", dirs)
