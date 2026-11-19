##############################################################################
# Soybean seed longevity
# Stage 3b: seed-related context, support across studies and final atlas table
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: describe, for each candidate, how it is expressed in seed-related body
# parts (seed, embryo, cotyledon, endosperm, seed coat, suspensor) compared with
# the other parts, and whether that signal is repeated across independent
# studies of the atlas. Body parts enter through their MEDIANS, so that leaf
# or root do not dominate just because they have thousands of samples.
# The temporal evidence of stage 2 is attached as a descriptive cross-layer
# pattern, without any score.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/03_expression_atlas/02_seed_context_study_support.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2, stringr
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2", "stringr"))

stage <- "03_expression_atlas"
dirs <- setup_stage(stage)
f_long <- file.path("data", "processed", "SEA_39genes_TPM_long.csv")
f_parts <- file.path(dirs$res, "03_candidate_bodypart_summary.csv")
f_tau <- file.path(dirs$res, "04_gene_classification_tau.csv")
f_temporal <- from_stage("02_maturation_timecourse", "09_candidate_temporal_evidence.csv")
check_inputs(c(f_long, f_parts, f_tau, f_temporal, RAW$atlas_meta))

parts <- readr::read_csv(f_parts, show_col_types = FALSE)
tau <- readr::read_csv(f_tau, show_col_types = FALSE)
temporal <- readr::read_csv(f_temporal, show_col_types = FALSE)
stopifnot(dplyr::n_distinct(parts$gene_id) == 39, all(SEED_PARTS %in% parts$body_part))

# ---------------------------------------------------------------------------
# 1) Seed-related context of each gene
# ---------------------------------------------------------------------------

# One row per gene: best seed-related part, best non-seed part, seed rank etc.
seed_context_gene <- function(x) {
  seed_x <- dplyr::filter(x, body_part %in% SEED_PARTS)
  non_x <- dplyr::filter(x, !body_part %in% SEED_PARTS)
  ranked <- dplyr::arrange(x, dplyr::desc(median_TPM), body_part)
  top_all <- ranked[1, ]
  top_seed <- dplyr::slice_max(seed_x, median_TPM, n = 1, with_ties = FALSE)
  top_non <- dplyr::slice_max(non_x, median_TPM, n = 1, with_ties = FALSE)
  seed_row <- dplyr::filter(x, body_part == "seed")
  tibble::tibble(
    top_body_part = top_all$body_part, top_bodypart_median_TPM = top_all$median_TPM,
    seed_median_TPM = seed_row$median_TPM, seed_mean_TPM = seed_row$mean_TPM,
    seed_q25_TPM = seed_row$q25_TPM, seed_q75_TPM = seed_row$q75_TPM,
    seed_prop_TPM_ge_1 = seed_row$prop_TPM_ge_1, seed_prop_TPM_ge_5 = seed_row$prop_TPM_ge_5,
    seed_bodypart_rank = which(ranked$body_part == "seed"),
    max_seed_related_part = top_seed$body_part, max_seed_related_median_TPM = top_seed$median_TPM,
    max_nonseed_part = top_non$body_part, max_nonseed_median_TPM = top_non$median_TPM,
    log2_max_seed_vs_nonseed = log2((top_seed$median_TPM + 1) / (top_non$median_TPM + 1)),
    n_seed_parts_median_ge_1 = sum(seed_x$median_TPM >= 1),
    n_seed_parts_median_ge_5 = sum(seed_x$median_TPM >= 5))
}
seed_context_gene(dplyr::filter(parts, gene_id == "Glyma.08G228400"))   # check: endosperm on top

seed_ctx <- parts %>%
  dplyr::group_by(gene_id) %>%
  dplyr::group_modify(~ seed_context_gene(.x)) %>%
  dplyr::ungroup()

# Project class (descriptive, not an official atlas class)
atlas_ctx <- seed_ctx %>%
  dplyr::left_join(dplyr::select(tau, gene_id, max_median_TPM, tau, expression_group, specific_parts), by = "gene_id") %>%
  dplyr::mutate(project_seed_context = dplyr::case_when(
    expression_group == "part-specific" &
      stringr::str_detect(dplyr::coalesce(specific_parts, ""), "seed|endosperm|embryo|cotyledon|suspensor") ~
      "part-specific_seed_related",
    max_seed_related_median_TPM >= 5 & top_body_part %in% SEED_PARTS ~ "seed_related_is_top_context",
    max_seed_related_median_TPM >= 5 ~ "expressed_in_seed_related_but_broad",
    max_seed_related_median_TPM >= 1 ~ "weak_seed_related_context",
    TRUE ~ "little_or_no_seed_related_expression")) %>%
  dplyr::arrange(dplyr::desc(max_seed_related_median_TPM))
readr::write_csv(atlas_ctx, file.path(dirs$res, "05_seed_context_39.csv"))

# ---------------------------------------------------------------------------
# 2) Descriptive cross-layer pattern with the temporal evidence (stage 2)
# Not a score or ranking: it only names the combination of evidence.
# ---------------------------------------------------------------------------
crosslayer <- temporal %>%
  dplyr::select(Wm82_a4_id, tested_DESeq2, filtered_low_count, shrunk_log2FoldChange_dry, padj_dry,
                significant_dry, shrunk_log2FoldChange_maturation, padj_maturation, significant_maturation,
                evidence_pattern) %>%
  dplyr::left_join(atlas_ctx, by = c("Wm82_a4_id" = "gene_id")) %>%
  dplyr::mutate(cross_layer_pattern = dplyr::case_when(
    evidence_pattern == "maturation_only" & project_seed_context == "part-specific_seed_related" ~
      "maturation_dominant + seed-part-specific",
    evidence_pattern == "maturation_only" & max_seed_related_median_TPM >= 5 ~
      "maturation_dominant + seed-related-expression",
    evidence_pattern == "both" & max_seed_related_median_TPM >= 5 ~
      "maturation+drying + seed-related-expression",
    evidence_pattern == "low_count_not_tested" & project_seed_context == "part-specific_seed_related" ~
      "low_in_GSE98199 + seed-part-specific_in_Atlas",
    evidence_pattern == "low_count_not_tested" ~ "low_in_GSE98199 + Atlas_context_available",
    TRUE ~ paste0(evidence_pattern, " + ", project_seed_context)))

# ---------------------------------------------------------------------------
# 3) Sample heterogeneity and support across studies
# Studies are identified by the Study_abstract field of the atlas metadata
# (one abstract per submitted study). A seed-related signal counts as
# replicated when at least 2 studies have median TPM >= 5 in that body part.
# ---------------------------------------------------------------------------
meta <- read.delim(RAW$atlas_meta, check.names = FALSE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(sample_id = as.character(BioSample), study_id = as.character(Study_abstract),
                   cultivar = as.character(Cultivar)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)
seed_expr <- readr::read_csv(f_long, show_col_types = FALSE) %>%
  dplyr::filter(body_part %in% SEED_PARTS) %>%
  dplyr::left_join(meta, by = "sample_id") %>%
  dplyr::mutate(study_id = dplyr::coalesce(study_id, "study_metadata_missing"))

heterogeneity <- seed_expr %>%
  dplyr::group_by(gene_id, body_part) %>%
  dplyr::summarise(n_samples = dplyr::n(), median_TPM = median(TPM), mean_TPM = mean(TPM),
                   q25_TPM = quantile(TPM, 0.25), q75_TPM = quantile(TPM, 0.75), IQR_TPM = IQR(TPM),
                   prop_TPM_ge_1 = mean(TPM >= 1), prop_TPM_ge_5 = mean(TPM >= 5), prop_zero = mean(TPM == 0),
                   .groups = "drop")
readr::write_csv(heterogeneity, file.path(dirs$res, "06_seed_sample_heterogeneity.csv"))

by_study <- seed_expr %>%
  dplyr::filter(!is.na(study_id), study_id != "") %>%
  dplyr::group_by(gene_id, body_part, study_id) %>%
  dplyr::summarise(n_samples = dplyr::n(), median_TPM = median(TPM), mean_TPM = mean(TPM),
                   prop_TPM_ge_5 = mean(TPM >= 5), .groups = "drop")
readr::write_csv(by_study, file.path(dirs$res, "07_gene_bodypart_by_study.csv"))

support <- by_study %>%
  dplyr::group_by(gene_id, body_part) %>%
  dplyr::summarise(n_studies = dplyr::n_distinct(study_id), n_studies_median_TPM_ge_1 = sum(median_TPM >= 1),
                   n_studies_median_TPM_ge_5 = sum(median_TPM >= 5),
                   prop_studies_median_TPM_ge_5 = n_studies_median_TPM_ge_5 / n_studies,
                   median_of_study_medians = median(median_TPM), max_study_median_TPM = max(median_TPM),
                   .groups = "drop")
readr::write_csv(support, file.path(dirs$res, "08_study_support_by_bodypart.csv"))

# Support in the best seed-related part of each gene
best_support <- atlas_ctx %>%
  dplyr::select(gene_id, max_seed_related_part) %>%
  dplyr::left_join(support, by = c("gene_id", "max_seed_related_part" = "body_part")) %>%
  dplyr::mutate(study_support_class = dplyr::case_when(
    is.na(n_studies) ~ "no_study_information",
    n_studies >= 2 & n_studies_median_TPM_ge_5 >= 2 ~ "replicated_high_expression_across_studies",
    n_studies >= 2 & n_studies_median_TPM_ge_5 == 1 ~ "high_expression_driven_by_one_study",
    n_studies == 1 & n_studies_median_TPM_ge_5 == 1 ~ "high_expression_single_study_only",
    n_studies_median_TPM_ge_1 >= 2 ~ "replicated_low_to_moderate_expression",
    TRUE ~ "little_cross_study_support"))
readr::write_csv(best_support, file.path(dirs$res, "09_best_seedpart_study_support.csv"))

# ---------------------------------------------------------------------------
# 4) Final atlas table (input of stages 4 to 6)
# ---------------------------------------------------------------------------
atlas_final <- crosslayer %>%
  dplyr::left_join(best_support, by = c("Wm82_a4_id" = "gene_id", "max_seed_related_part")) %>%
  dplyr::mutate(atlas_final_pattern = dplyr::case_when(
    project_seed_context == "part-specific_seed_related" &
      study_support_class == "replicated_high_expression_across_studies" ~ "seed-part-specific_with_cross-study_support",
    project_seed_context == "part-specific_seed_related" ~ paste0("seed-part-specific; ", study_support_class),
    project_seed_context == "seed_related_is_top_context" &
      study_support_class == "replicated_high_expression_across_studies" ~ "seed-related-top-context_with_cross-study-support",
    project_seed_context == "expressed_in_seed_related_but_broad" ~ "broad_but_consistently_expressed_in_seed-related-parts",
    project_seed_context == "weak_seed_related_context" ~ "weak_seed-related-context",
    project_seed_context == "little_or_no_seed_related_expression" ~ "little_or_no_seed-related-expression",
    TRUE ~ project_seed_context))
readr::write_csv(atlas_final, file.path(dirs$res, "10_atlas_evidence_39.csv"))

(atlas_summary <- dplyr::bind_rows(
  dplyr::count(atlas_final, level = expression_group) %>% dplyr::mutate(category = "Atlas expression group"),
  dplyr::count(atlas_final, level = project_seed_context) %>% dplyr::mutate(category = "Seed context"),
  dplyr::count(atlas_final, level = study_support_class) %>% dplyr::mutate(category = "Study support"),
  dplyr::count(atlas_final, level = atlas_final_pattern) %>% dplyr::mutate(category = "Final atlas pattern")) %>%
  dplyr::select(category, level, n))
readr::write_csv(atlas_summary, file.path(dirs$res, "11_atlas_summary.csv"))

# Reading: 30 replicated high + 4 replicated low/moderate + 5 little support is
# the expected split. A gene that is "part-specific" but supported by one study
# only would be a warning sign (protocol effect), not seed biology.

# ---------------------------------------------------------------------------
# 5) Figures
# ---------------------------------------------------------------------------
seed_plot <- parts %>%
  dplyr::filter(body_part %in% SEED_PARTS) %>%
  dplyr::mutate(body_part = factor(body_part, SEED_PARTS))
gene_order <- seed_plot %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(m = max(median_TPM)) %>%
  dplyr::arrange(m) %>%
  dplyr::pull(gene_id)
p_bubble <- ggplot2::ggplot(dplyr::mutate(seed_plot, gene_id = factor(gene_id, gene_order)),
                            ggplot2::aes(body_part, gene_id, size = prop_TPM_ge_5, fill = log2(median_TPM + 1))) +
  ggplot2::geom_point(shape = 21) +
  ggplot2::scale_size_continuous(range = c(0.5, 6), limits = c(0, 1), name = "Fraction of samples\nwith TPM >= 5") +
  ggplot2::labs(x = "Seed-related body part", y = "Candidate gene (Wm82.a4.v1)", fill = "log2(median TPM + 1)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.text.y = ggplot2::element_text(size = 6))
ggplot2::ggsave(file.path(dirs$fig, "02_seed_related_bubbleplot.png"), p_bubble, width = 10, height = 12, dpi = 300)

focal <- crosslayer %>%
  dplyr::filter(evidence_pattern %in% c("maturation_only", "both") | project_seed_context == "part-specific_seed_related") %>%
  dplyr::pull(Wm82_a4_id)
p_cmp <- atlas_ctx %>%
  dplyr::mutate(label = ifelse(gene_id %in% focal, gene_id, NA_character_)) %>%
  ggplot2::ggplot(ggplot2::aes(log2(max_nonseed_median_TPM + 1), log2(max_seed_related_median_TPM + 1))) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(ggplot2::aes(label = label), na.rm = TRUE, nudge_y = 0.12, size = 2.7, check_overlap = TRUE) +
  ggplot2::labs(x = "log2(max median TPM in non-seed parts + 1)", y = "log2(max median TPM in seed-related parts + 1)") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "03_seed_vs_nonseed_context.png"), p_cmp, width = 8, height = 7, dpi = 300)

p_dist <- seed_expr %>%
  dplyr::filter(gene_id %in% focal) %>%
  dplyr::mutate(body_part = factor(body_part, SEED_PARTS)) %>%
  ggplot2::ggplot(ggplot2::aes(body_part, log2(TPM + 1))) +
  ggplot2::geom_boxplot(outlier.size = 0.45) +
  ggplot2::facet_wrap(~ gene_id, scales = "free_y", ncol = 3) +
  ggplot2::labs(x = "Seed-related body part", y = "log2(TPM + 1)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), strip.text = ggplot2::element_text(size = 8))
ggplot2::ggsave(file.path(dirs$fig, "04_focal_candidates_sample_distributions.png"), p_dist, width = 13, height = 12, dpi = 300)

close_stage("03_expression_atlas_context", dirs)
