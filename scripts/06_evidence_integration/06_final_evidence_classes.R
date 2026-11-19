##############################################################################
# Soybean seed longevity
# Stage 6f: final evidence table per gene and per locus
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: gather everything produced in stage 6 and give each gene an evidence
# CLASS by rules fixed before looking at the result (section 1). It is not a
# new score: the class says in which domain the gene has support. Maturation
# program (script 04) and tissue (script 05) enter as biological annotation and
# do not change the class.
#
# Run from the repository root, after scripts 01 to 05 of this stage:
#   source("scripts/06_evidence_integration/06_final_evidence_classes.R")
#
# Packages: readr, dplyr, tidyr, tibble, stringr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "stringr", "ggplot2"))

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "06_final_classes"))
f_int <- from_stage(stage, "01_integrated_matrix/02_integrated_evidence_matrix_39.csv")
f_weak <- from_stage(stage, "03_convergence_loci/04_weakest_link_2domains.csv")
f_w3 <- from_stage(stage, "03_convergence_loci/05_weight_sensitivity_3layers.csv")
f_w2 <- from_stage(stage, "03_convergence_loci/06_weight_sensitivity_2domains.csv")
f_anchor <- from_stage(stage, "04_maturation_anchor/05_candidate_module_anchor_39.csv")
f_tissue <- from_stage(stage, "05_tissue_within_project/02_module_tissue_within_vs_original.csv")
f_locus <- from_stage("01_candidates", "02_locus_map.csv")
check_inputs(c(f_int, f_weak, f_w3, f_w2, f_anchor, f_tissue, f_locus))

# ---------------------------------------------------------------------------
# 1) Class rules (in order; the first that applies wins)
#   Complete cases (temporal + SEA), K_TOP = 10 (1/3 of 30), K_MOD = 15:
#     A_convergent        rank <= 10 in BOTH domains
#     B_balanced_moderate both ranks <= 15 (not A)
#     C_temporal_dominant temporal <= 10 and SEA > 15
#     D_SEA_dominant      SEA <= 10 and temporal > 15
#     E_low_support       none of the above
#   Incomplete coverage (no global position, never sent to the bottom):
#     F_SEA_only_low_count  SEA complete, low count in GSE98199 (= little
#                           expression at R7-R9, which is information)
#     G_temporal_only       temporal available, SEA incomplete
#     H_insufficient_coverage
# ---------------------------------------------------------------------------
K_TOP <- 10
K_MOD <- 15
STABLE_FREQ <- 0.80                            # min fraction in top 10 under random weights

integrated <- readr::read_csv(f_int, show_col_types = FALSE) %>%
  dplyr::select(gene_id, genetic_study, genetic_traits, evidence_coverage_class, temporal_support_label,
                significant_maturation, shrunk_lfc_maturation, padj_maturation, significant_dry, shrunk_lfc_dry,
                project_seed_context, network_support_label, functional_support_class, neighborhood_BP_themes,
                temporal_rank_available_all, atlas_rank_all_39, network_rank_available_all, sea_rank_available_all)
weak <- readr::read_csv(f_weak, show_col_types = FALSE) %>%
  dplyr::select(gene_id, rank_T, rank_SEA, rank_A, rank_N, worst_rank_2dom, rank_weakest_link)
w3 <- readr::read_csv(f_w3, show_col_types = FALSE) %>%
  dplyr::select(gene_id, freq_top10_3layer = freq_top10, median_rank_3layer = median_rank)
w2 <- dplyr::select(readr::read_csv(f_w2, show_col_types = FALSE), gene_id, freq_top10_2dom, median_rank_2dom)
locus <- dplyr::select(readr::read_csv(f_locus, show_col_types = FALSE), gene_id, locus_id, chromosome, n_genes_locus)
anchor <- readr::read_csv(f_anchor, show_col_types = FALSE) %>%
  dplyr::select(gene_id, module = module_for_integration, module_trend = trend_class, module_mean_rho = mean_rho,
                module_robust_no71 = robust_without_7.1, maturation_program = anchor_class)
# Project-robust tissues per module: only "survives" associations, by |t| within project
tissue_mod <- readr::read_csv(f_tissue, show_col_types = FALSE) %>%
  dplyr::filter(class == "survives_within_project") %>%
  dplyr::arrange(ME, dplyr::desc(abs(t_within))) %>%
  dplyr::group_by(module = sub("^ME", "", ME)) %>%
  dplyr::summarise(robust_tissue_top = tissue[1],
                   robust_tissues_all = paste(sprintf("%s (%s)", tissue, ifelse(t_within > 0, "+", "-")), collapse = "; "))

# ---------------------------------------------------------------------------
# 2) Table per gene
# ---------------------------------------------------------------------------
genes <- integrated %>%
  dplyr::left_join(locus, by = "gene_id") %>%
  dplyr::left_join(weak, by = "gene_id") %>%
  dplyr::left_join(w3, by = "gene_id") %>%
  dplyr::left_join(w2, by = "gene_id") %>%
  dplyr::left_join(anchor, by = "gene_id") %>%
  dplyr::left_join(tissue_mod, by = "module") %>%
  dplyr::mutate(
    complete_case = !is.na(worst_rank_2dom),
    maturation_direction = dplyr::case_when(
      temporal_support_label == "not_tested_low_count" ~ "low_count",
      significant_maturation %in% TRUE & shrunk_lfc_maturation > 0 ~ "up",
      significant_maturation %in% TRUE & shrunk_lfc_maturation < 0 ~ "down",
      TRUE ~ "no_DE"),
    gene_module_concordance = dplyr::case_when(           # same direction as its module? (DE genes only)
      !maturation_direction %in% c("up", "down") | is.na(module_trend) ~ NA_character_,
      maturation_direction == "up" & module_trend %in% c("rises_in_maturation", "weak_rise") ~ "concordant",
      maturation_direction == "down" & module_trend %in% c("falls_in_maturation", "weak_fall") ~ "concordant",
      TRUE ~ "discordant"),
    weight_stable_top10 = freq_top10_3layer >= STABLE_FREQ & freq_top10_2dom >= STABLE_FREQ,
    evidence_class = dplyr::case_when(
      complete_case & rank_T <= K_TOP & rank_SEA <= K_TOP ~ "A_convergent",
      complete_case & rank_T <= K_MOD & rank_SEA <= K_MOD ~ "B_balanced_moderate",
      complete_case & rank_T <= K_TOP & rank_SEA > K_MOD ~ "C_temporal_dominant",
      complete_case & rank_SEA <= K_TOP & rank_T > K_MOD ~ "D_SEA_dominant",
      complete_case ~ "E_low_support",
      evidence_coverage_class == "SEA_complete_temporal_missing" ~ "F_SEA_only_low_count",
      evidence_coverage_class == "temporal_available_SEA_incomplete" ~ "G_temporal_only",
      TRUE ~ "H_insufficient_coverage")) %>%
  dplyr::group_by(locus_id) %>%
  dplyr::mutate(rank_in_locus = rank(rank_weakest_link, na.last = "keep", ties.method = "min"),
                best_in_locus = rank_in_locus %in% 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(evidence_class, rank_weakest_link)
stopifnot(nrow(genes) == 39, !anyDuplicated(genes$gene_id))
readr::write_csv(genes, file.path(dirs$res, "01_final_evidence_table_39_genes.csv"))
(class_counts <- dplyr::count(genes, evidence_class))
readr::write_csv(class_counts, file.path(dirs$res, "02_evidence_class_counts.csv"))
# expected: A 2 (04G192100, 04G192900); B 4; C 7; D 7; E 10; F 2; G 1; H 6
print(dplyr::count(genes, maturation_direction, gene_module_concordance))   # expected: all DE genes concordant

# ---------------------------------------------------------------------------
# 3) Table per locus
#   convergent_locus = some gene A; moderate_locus = best gene B;
#   single_domain_locus = best gene C or D; weak_locus = best gene E;
#   no_comparable_gene = no complete case
# ---------------------------------------------------------------------------
loci <- genes %>%
  dplyr::group_by(locus_id, genetic_study, chromosome) %>%
  dplyr::summarise(n_genes = dplyr::n(), n_complete = sum(complete_case), genes = paste(gene_id, collapse = "; "),
                   classes = paste(sort(unique(evidence_class)), collapse = "; "),
                   best_gene = if (any(best_in_locus)) gene_id[best_in_locus][1] else NA_character_,
                   best_class = if (any(best_in_locus)) evidence_class[best_in_locus][1] else NA_character_,
                   best_worst_rank = if (any(complete_case)) min(worst_rank_2dom, na.rm = TRUE) else NA_real_,
                   best_program = if (any(best_in_locus)) maturation_program[best_in_locus][1] else NA_character_,
                   best_robust_tissue = if (any(best_in_locus)) robust_tissue_top[best_in_locus][1] else NA_character_,
                   .groups = "drop") %>%
  dplyr::mutate(locus_status = dplyr::case_when(
    stringr::str_detect(classes, "A_convergent") ~ "convergent_locus",
    best_class == "B_balanced_moderate" ~ "moderate_locus",
    best_class %in% c("C_temporal_dominant", "D_SEA_dominant") ~ "single_domain_locus",
    best_class == "E_low_support" ~ "weak_locus",
    TRUE ~ "no_comparable_gene")) %>%
  dplyr::arrange(factor(locus_status, c("convergent_locus", "moderate_locus", "single_domain_locus", "weak_locus",
                                        "no_comparable_gene")), best_worst_rank)
readr::write_csv(loci, file.path(dirs$res, "03_final_locus_table.csv"))
print(dplyr::count(loci, locus_status))        # expected: 1 convergent, 4 moderate, 10 single, 6 weak/none

# ---------------------------------------------------------------------------
# 4) Compact table for the manuscript (classes A to D and F)
# ---------------------------------------------------------------------------
pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.0f%%", 100 * x))
manuscript <- genes %>%
  dplyr::filter(evidence_class %in% c("A_convergent", "B_balanced_moderate", "C_temporal_dominant",
                                      "D_SEA_dominant", "F_SEA_only_low_count")) %>%
  dplyr::transmute(Class = evidence_class, Locus = locus_id, Gene = gene_id, Study = genetic_study,
                   `Temporal rank` = rank_T, `SEA rank` = rank_SEA,
                   `Top-10 under weights (3 layers / 2 domains)` = paste(pct(freq_top10_3layer), "/", pct(freq_top10_2dom)),
                   `Maturation DE` = maturation_direction, Module = module, `Module program` = maturation_program,
                   `Robust tissue (within BioProject)` = robust_tissues_all) %>%
  dplyr::arrange(Class, `Temporal rank`)
readr::write_csv(manuscript, file.path(dirs$res, "04_manuscript_table_priority_genes.csv"))

# ---------------------------------------------------------------------------
# 5) Figure: temporal vs SEA rank with the final class
# ---------------------------------------------------------------------------
PAL <- c(A_convergent = "#1b7837", B_balanced_moderate = "#7fbf7b", C_temporal_dominant = "#2166ac",
         D_SEA_dominant = "#b2182b", E_low_support = "grey70")
cc <- dplyr::filter(genes, complete_case)
p_cls <- ggplot2::ggplot(cc, ggplot2::aes(rank_T, rank_SEA, colour = evidence_class)) +
  ggplot2::geom_hline(yintercept = c(K_TOP, K_MOD) + 0.5, linetype = c(2, 3), colour = "grey60") +
  ggplot2::geom_vline(xintercept = c(K_TOP, K_MOD) + 0.5, linetype = c(2, 3), colour = "grey60") +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::geom_text(data = dplyr::filter(cc, evidence_class != "E_low_support"), ggplot2::aes(label = gene_id),
                     size = 2.6, nudge_y = 0.8, check_overlap = TRUE, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = PAL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.05))) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = "Temporal rank (1 = strongest)", y = "SEA rank (atlas + network)", colour = "Evidence class",
                subtitle = paste0("Dashed: top ", K_TOP, "; dotted: top ", K_MOD, " (n = ", nrow(cc), ")")) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "01_final_classes_temporal_vs_SEA.png"), p_cls, width = 9, height = 8, dpi = 300)

close_stage("06_evidence_integration_classes", dirs)
