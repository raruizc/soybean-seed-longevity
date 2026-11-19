##############################################################################
# Soybean seed longevity
# Stage 6a: integrated evidence matrix and primary ranking
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: join, for the 39 Wm82.a4 loci, the evidence of stages 1 to 5, give each
# layer an internal rank and combine the layers into a primary ranking
# (temporal domain + SEA domain).
# Principle: genetics is the entry criterion (no weight); GSE98199 is the
# temporal domain; atlas and network come from the same Soybean Expression
# Atlas and are collapsed into one SEA domain; GO derives from the network and
# is kept only as interpretive support.
#
# Run from the repository root, after stages 1 to 5:
#   source("scripts/06_evidence_integration/01_integrated_matrix_primary_ranking.R")
#
# Packages: readr, readxl, dplyr, tidyr, tibble, stringr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "readxl", "dplyr", "tidyr", "tibble", "stringr", "ggplot2"))

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "01_integrated_matrix"))
f_temporal <- from_stage("02_maturation_timecourse", "09_candidate_temporal_evidence.csv")
f_atlas <- from_stage("03_expression_atlas", "10_atlas_evidence_39.csv")
f_network <- from_stage("04_coexpression_network", "33_candidate_network_evidence_39.csv")
f_function <- from_stage("05_functional_enrichment", "12_candidate_functional_profile_39.csv")
check_inputs(c(RAW$gwas_xlsx, f_temporal, f_atlas, f_network, f_function))

# Rank fraction in (0, 1], lower = better; NA stays NA and is not in the denominator
rank_fraction <- function(x, higher_is_better = TRUE) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  if (sum(ok) == 0) return(out)
  out[ok] <- rank(if (higher_is_better) -x[ok] else x[ok], ties.method = "average") / sum(ok)
  out
}
rank_pos <- function(x) rank(x, na.last = "keep", ties.method = "average")   # 1 = best, NA kept

# ---------------------------------------------------------------------------
# 1) Genetics: provenance only
# Wu et al. (2024) and Thenveettil et al. (2025) have different designs, so
# genetics does not become a weight.
# ---------------------------------------------------------------------------
genetics <- readxl::read_excel(RAW$gwas_xlsx, sheet = "Final_Loci") %>%
  dplyr::filter(!is.na(Wm82_a4_id)) %>%
  dplyr::distinct(Wm82_a4_id, .keep_all = TRUE) %>%
  dplyr::transmute(gene_id = Wm82_a4_id, genetic_study = studies, genetic_traits = traits,
                   genetic_evidence_class = evidence_class, genetic_annotation_stage1 = functional_annotation,
                   harmonization_confidence = confidence,
                   n_genetic_traits = dplyr::if_else(is.na(traits), NA_integer_, stringr::str_count(traits, ";") + 1L))
stopifnot(nrow(genetics) == 39)

# ---------------------------------------------------------------------------
# 2) Temporal domain (GSE98199)
# R_T = sqrt(rank padj_maturation x rank |shrunken LFC maturation|). Drying is
# kept as annotation only, and the sign of the LFC does not enter.
# Note: padj and |LFC| are correlated (partial double counting), and ordering
# padj among non-significant genes carries little information.
# ---------------------------------------------------------------------------
temporal <- readr::read_csv(f_temporal, show_col_types = FALSE) %>%
  dplyr::transmute(gene_id = Wm82_a4_id, tested_DESeq2, filtered_low_count, evidence_pattern,
                   significant_maturation, significant_dry, padj_maturation,
                   shrunk_lfc_maturation = shrunk_log2FoldChange_maturation,
                   abs_shrunk_lfc_maturation = abs(shrunk_log2FoldChange_maturation),
                   padj_dry, shrunk_lfc_dry = shrunk_log2FoldChange_dry) %>%
  dplyr::mutate(temporal_rank_padj = rank_fraction(padj_maturation, higher_is_better = FALSE),
                temporal_rank_effect = rank_fraction(abs_shrunk_lfc_maturation),
                temporal_domain_score = sqrt(temporal_rank_padj * temporal_rank_effect))

# ---------------------------------------------------------------------------
# 3) Atlas: anatomical context
# Geometric mean of three ranks: seed vs non-seed (context), support across
# studies, and maximum level in a seed-related part. A barely expressed gene
# then does not look "specific" from the ratio alone, and a very abundant
# ubiquitous gene does not dominate on its own.
# ---------------------------------------------------------------------------
atlas <- readr::read_csv(f_atlas, show_col_types = FALSE) %>%
  dplyr::transmute(gene_id = Wm82_a4_id, project_seed_context, atlas_final_pattern, expression_group, tau,
                   specific_parts, seed_median_TPM, max_seed_related_part, max_seed_related_median_TPM,
                   max_nonseed_part, max_nonseed_median_TPM, log2_max_seed_vs_nonseed,
                   prop_studies_median_TPM_ge_5, study_support_class, n_studies) %>%
  dplyr::mutate(atlas_rank_seed_vs_nonseed = rank_fraction(log2_max_seed_vs_nonseed),
                atlas_rank_cross_study = rank_fraction(prop_studies_median_TPM_ge_5),
                atlas_rank_seed_expression = rank_fraction(max_seed_related_median_TPM),
                atlas_context_score = (atlas_rank_seed_vs_nonseed * atlas_rank_cross_study *
                                         atlas_rank_seed_expression)^(1/3))

# ---------------------------------------------------------------------------
# 4) Network: gene-eigengene coupling
# Direct members use kME in their own module; genes outside the topology use
# the strongest positive r of the post hoc projection. Hub and kWithin are
# secondary support. Note: kME has no anchor in longevity, and the maximum r
# over 21 modules tends to exceed an own-module kME (selection of the maximum).
# ---------------------------------------------------------------------------
network <- readr::read_csv(f_network, show_col_types = FALSE) %>%
  dplyr::mutate(module_coupling = dplyr::case_when(
                  network_status == "direct_member_of_network" & !is.na(direct_kME) ~ direct_kME,
                  !is.na(strongest_positive_r) ~ strongest_positive_r),
                module_for_integration = dplyr::case_when(
                  network_status == "direct_member_of_network" & !is.na(direct_module) &
                    direct_module != "grey" ~ direct_module,
                  !is.na(strongest_positive_module) ~ strongest_positive_module,
                  TRUE ~ direct_module)) %>%
  dplyr::transmute(gene_id, network_status, network_evidence_type, module_for_integration, direct_module,
                   direct_kME, direct_kWithin_percentile, direct_hub, strongest_positive_module,
                   strongest_positive_r, strongest_positive_padj, projected_module_bodypart,
                   projected_module_bodypart_cor, projected_module_bodypart_padj, module_coupling,
                   network_rank_coupling = rank_fraction(module_coupling))

# ---------------------------------------------------------------------------
# 5) Function (GO): interpretive support, outside the ranking
# The number of GO terms depends on annotation density and ontology
# redundancy, so it only becomes a descriptive class. The class depends on the
# BH family used for the neighbourhoods (stage 5, script 02).
# ---------------------------------------------------------------------------
function_profile <- readr::read_csv(f_function, show_col_types = FALSE) %>%
  dplyr::transmute(gene_id, functional_module, module_BP_themes, module_MF_themes, module_CC_themes,
                   neighborhood_BP_themes, neighborhood_MF_themes, neighborhood_CC_themes, direct_GO_n,
                   direct_GO_terms, n_GO_robust_both, n_GO_within_module_only, n_GO_network_only, robust_GO_terms,
                   functional_support_class = dplyr::case_when(
                     dplyr::coalesce(n_GO_robust_both, 0) > 0 ~ "robust_local_GO_both_backgrounds",
                     dplyr::coalesce(n_GO_within_module_only, 0) > 0 ~ "local_GO_within_module_only",
                     dplyr::coalesce(n_GO_network_only, 0) > 0 ~ "module_level_GO_only",
                     !is.na(module_BP_themes) | !is.na(module_MF_themes) | !is.na(module_CC_themes) ~
                       "module_function_available_no_local_GO",
                     dplyr::coalesce(direct_GO_n, 0) > 0 ~ "direct_GO_only",
                     TRUE ~ "limited_or_no_GO_support"))

# ---------------------------------------------------------------------------
# 6) Integrated matrix and combined scores
# Strict SEA requires atlas + network. The expanded score (3 layers) exists only
# for complete cases, since means of 2 and 3 layers are not comparable.
# Primary ranking: R_primary = sqrt(R_T x R_SEA) and Borda = (R_T + R_SEA) / 2.
# As R_SEA = sqrt(R_A x R_N), the implicit weights are 50% temporal, 25% atlas,
# 25% network; this is stated in the methods.
# ---------------------------------------------------------------------------
integrated <- genetics %>%
  dplyr::left_join(temporal, by = "gene_id") %>%
  dplyr::left_join(atlas, by = "gene_id") %>%
  dplyr::left_join(network, by = "gene_id") %>%
  dplyr::left_join(function_profile, by = "gene_id") %>%
  dplyr::mutate(
    sea_domain_score_strict = sqrt(atlas_context_score * network_rank_coupling),
    n_quant_layers_available = rowSums(!is.na(cbind(temporal_domain_score, atlas_context_score, network_rank_coupling))),
    expanded_score_complete = (temporal_domain_score * atlas_context_score * network_rank_coupling)^(1/3),
    primary_rank_product_score = sqrt(temporal_domain_score * sea_domain_score_strict),
    primary_borda_score = (temporal_domain_score + sea_domain_score_strict) / 2,
    evidence_coverage_class = dplyr::case_when(
      !is.na(temporal_domain_score) & !is.na(sea_domain_score_strict) ~ "two_independent_domains_complete",
      is.na(temporal_domain_score) & !is.na(sea_domain_score_strict) ~ "SEA_complete_temporal_missing",
      !is.na(temporal_domain_score) & is.na(sea_domain_score_strict) ~ "temporal_available_SEA_incomplete",
      TRUE ~ "limited_quantitative_coverage"),
    # positions (1 = best); the "_available_all" ranks are within-layer references only
    primary_rank_RP = rank_pos(primary_rank_product_score),
    primary_rank_Borda = rank_pos(primary_borda_score),
    primary_consensus_rank = (primary_rank_RP + primary_rank_Borda) / 2,
    expanded_rank_complete = rank_pos(expanded_score_complete),
    temporal_rank_available_all = rank_pos(temporal_domain_score),
    atlas_rank_all_39 = rank_pos(atlas_context_score),
    network_rank_available_all = rank_pos(network_rank_coupling),
    sea_rank_available_all = rank_pos(sea_domain_score_strict),
    # descriptive labels (not a new score); low count in GSE98199 is information
    temporal_support_label = dplyr::case_when(
      evidence_pattern == "maturation_only" ~ "maturation_dominant",
      evidence_pattern == "both" ~ "maturation_and_drying",
      evidence_pattern == "drying_only" ~ "drying_dominant",
      evidence_pattern == "low_count_not_tested" ~ "not_tested_low_count",
      TRUE ~ "no_threshold_level_DE"),
    network_support_label = dplyr::case_when(
      direct_hub %in% TRUE ~ "direct_hub",
      network_status == "direct_member_of_network" & !is.na(direct_kME) ~ "direct_module_member",
      dplyr::coalesce(strongest_positive_r, -1) >= 0.80 ~ "strong_posthoc_projection",
      dplyr::coalesce(strongest_positive_r, -1) >= 0.60 ~ "moderate_posthoc_projection",
      !is.na(strongest_positive_r) ~ "weak_posthoc_projection",
      TRUE ~ "not_evaluated_or_no_projection"))
stopifnot(!anyDuplicated(integrated$gene_id))

(audit <- tibble::tibble(
  metric = c("candidate_genes", "temporal_scores", "atlas_scores", "network_scores", "strict_SEA_scores",
             "primary_two_domain_scores", "expanded_complete_case_scores", "robust_local_GO_support"),
  value = c(nrow(integrated), sum(!is.na(integrated$temporal_domain_score)), sum(!is.na(integrated$atlas_context_score)),
            sum(!is.na(integrated$network_rank_coupling)), sum(!is.na(integrated$sea_domain_score_strict)),
            sum(!is.na(integrated$primary_rank_product_score)), sum(!is.na(integrated$expanded_score_complete)),
            sum(integrated$functional_support_class == "robust_local_GO_both_backgrounds"))))
readr::write_csv(audit, file.path(dirs$res, "01_integration_audit.csv"))
# expected: 39, 31, 39, 32, 32, 30, 30, 12

# ---------------------------------------------------------------------------
# 7) Output tables
# ---------------------------------------------------------------------------
integrated_out <- dplyr::arrange(integrated, primary_consensus_rank, expanded_rank_complete)
readr::write_csv(integrated_out, file.path(dirs$res, "02_integrated_evidence_matrix_39.csv"))
primary_table <- integrated_out %>%
  dplyr::select(gene_id, genetic_study, genetic_traits, temporal_support_label, evidence_pattern,
                temporal_domain_score, project_seed_context, atlas_final_pattern, atlas_context_score,
                network_support_label, module_for_integration, module_coupling, network_rank_coupling,
                sea_domain_score_strict, functional_support_class, neighborhood_BP_themes, neighborhood_MF_themes,
                neighborhood_CC_themes, primary_rank_product_score, primary_borda_score, primary_rank_RP,
                primary_rank_Borda, primary_consensus_rank, n_quant_layers_available, evidence_coverage_class,
                expanded_score_complete, expanded_rank_complete, temporal_rank_available_all, atlas_rank_all_39,
                network_rank_available_all, sea_rank_available_all)
readr::write_csv(primary_table, file.path(dirs$res, "03_primary_consensus_ranking.csv"))
incomplete <- dplyr::filter(integrated_out, is.na(primary_consensus_rank))
readr::write_csv(incomplete, file.path(dirs$res, "04_incomplete_coverage_candidates.csv"))
print(head(dplyr::select(dplyr::filter(primary_table, !is.na(primary_consensus_rank)), gene_id, primary_consensus_rank,
                         temporal_support_label, project_seed_context, network_support_label), 10))
# expected top: 04G192900 and 02G149600 (tied at 2), then 04G192100

# ---------------------------------------------------------------------------
# 8) Figure: three layers with NA marked (missing is not negative evidence)
# Complete cases ordered by consensus; incomplete ones after them, by the rank
# of the layer they have (+1000 only pushes them to the end of the axis).
# ---------------------------------------------------------------------------
plot_order <- integrated_out %>%
  dplyr::mutate(ref = dplyr::case_when(
                  evidence_coverage_class == "SEA_complete_temporal_missing" ~ sea_rank_available_all,
                  evidence_coverage_class == "temporal_available_SEA_incomplete" ~ temporal_rank_available_all,
                  TRUE ~ atlas_rank_all_39),
                ord = dplyr::coalesce(primary_consensus_rank, 1000 + ref)) %>%
  dplyr::arrange(ord) %>%
  dplyr::pull(gene_id)
heat <- integrated_out %>%
  dplyr::select(gene_id, temporal_domain_score, atlas_context_score, network_rank_coupling) %>%
  tidyr::pivot_longer(-gene_id, names_to = "layer", values_to = "rank_fraction") %>%
  dplyr::mutate(gene_id = factor(gene_id, rev(plot_order)),
                layer = factor(layer, c("temporal_domain_score", "atlas_context_score", "network_rank_coupling"),
                               c("Temporal maturation", "Atlas seed context", "Network module coupling")))
p_heat <- ggplot2::ggplot(heat, ggplot2::aes(layer, gene_id, fill = rank_fraction)) +
  ggplot2::geom_tile(color = "grey85", linewidth = 0.25) +
  ggplot2::geom_text(data = dplyr::filter(heat, is.na(rank_fraction)), label = "NA", size = 2.5) +
  ggplot2::scale_fill_gradient(low = "grey95", high = "grey20", na.value = "white", limits = c(0, 1)) +
  ggplot2::labs(x = NULL, y = "Candidate gene", fill = "Rank fraction\n(lower = stronger)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(dirs$fig, "01_layer_rank_heatmap.png"), p_heat, width = 10, height = 12, dpi = 300)

close_stage("06_evidence_integration_matrix", dirs)
