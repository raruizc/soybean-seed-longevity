##############################################################################
# Soybean seed longevity
# Stage 6e: module x tissue association controlling for BioProject
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: the module-tissue associations of stage 4 used all 231 samples without
# controlling for study, and BioProject is partly confounded with anatomy
# (stage 4, script 03). About 2/3 of the samples are in BioProjects with more
# than one tissue, so tissues can be compared WITHIN a project. In the LCM
# series (Goldberg/Harada) each project is one stage, so the within-project
# comparison also compares compartments at the same stage.
# Model, for each module and tissue (one vs the others), only in multi-tissue
# projects containing that tissue (>= 2 samples of it and >= 2 of others):
#   ME ~ tissue + BioProject     (BioProject as a fixed block)
# The same model without BioProject, on the same subset, shows how much the
# control changes the result.
#
# Run from the repository root, after stage 4 and script 04 of this stage:
#   source("scripts/06_evidence_integration/05_module_tissue_within_project.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "05_tissue_within_project"))
s4 <- "04_coexpression_network"
f_me <- from_stage(s4, "23_module_eigengenes.csv")
f_dec <- from_stage(s4, "03_sample_curation_decisions.csv")
f_assoc <- from_stage(s4, "25_module_bodypart_correlations.csv")
f_anchor <- from_stage(stage, "04_maturation_anchor/05_candidate_module_anchor_39.csv")
check_inputs(c(f_me, f_dec, f_assoc, f_anchor))
MIN_PER_TISSUE <- 2

me_long <- readr::read_csv(f_me, show_col_types = FALSE) %>%
  tidyr::pivot_longer(-sample_id, names_to = "ME", values_to = "eigengene") %>%
  dplyr::filter(ME != "MEgrey")
meta <- readr::read_csv(f_dec, show_col_types = FALSE) %>%
  dplyr::filter(sample_id %in% me_long$sample_id) %>%
  dplyr::distinct(sample_id, BioProject, body_part)
me_long <- dplyr::left_join(me_long, meta, by = "sample_id")
(n_samples <- nrow(meta))                      # expected: 231

(projects <- meta %>%
    dplyr::group_by(BioProject) %>%
    dplyr::summarise(n_samples = dplyr::n(), n_parts = dplyr::n_distinct(body_part),
                     parts = paste(sort(unique(body_part)), collapse = ", ")) %>%
    dplyr::arrange(dplyr::desc(n_samples)))
readr::write_csv(projects, file.path(dirs$res, "01_bioproject_tissue_composition.csv"))

# ---------------------------------------------------------------------------
# 1) Within-project model for each module and tissue
# estimate = mean eigengene difference (tissue - others) adjusted for project
# ---------------------------------------------------------------------------

# Projects that carry an internal "tissue t vs others" contrast
projects_for <- function(t) {
  meta %>%
    dplyr::group_by(BioProject) %>%
    dplyr::summarise(n_t = sum(body_part == t), n_o = sum(body_part != t)) %>%
    dplyr::filter(n_t >= MIN_PER_TISSUE, n_o >= MIN_PER_TISSUE) %>%
    dplyr::pull(BioProject)
}

# Fit ME ~ tissue (+ BioProject) for one module and one tissue
fit_within <- function(me, t) {
  pj <- projects_for(t)
  d <- dplyr::filter(me_long, ME == me, BioProject %in% pj) %>%
    dplyr::mutate(is_t = as.numeric(body_part == t))
  out <- tibble::tibble(ME = me, tissue = t, n_projects = length(pj), n_samples = nrow(d), n_tissue = sum(d$is_t))
  if (length(pj) == 0 || length(unique(d$is_t)) < 2) return(out)
  cw <- summary(lm(if (length(pj) > 1) eigengene ~ is_t + BioProject else eigengene ~ is_t, d))$coefficients["is_t", ]
  cn <- summary(lm(eigengene ~ is_t, d))$coefficients["is_t", ]
  diff_proj <- d %>%
    dplyr::group_by(BioProject) %>%
    dplyr::summarise(dd = mean(eigengene[is_t == 1]) - mean(eigengene[is_t == 0]))
  dplyr::mutate(out, estimate_within = cw[["Estimate"]], t_within = cw[["t value"]], p_within = cw[["Pr(>|t|)"]],
                estimate_naive = cn[["Estimate"]], p_naive = cn[["Pr(>|t|)"]],
                n_proj_same_sign = sum(sign(diff_proj$dd) == sign(cw[["Estimate"]])))
}
fit_within("MEpurple", "seed coat")            # check with one pair (stage 4: cor = -0.71)

grid <- tidyr::expand_grid(ME = sort(unique(me_long$ME)), tissue = sort(unique(meta$body_part)))
within <- dplyr::bind_rows(Map(fit_within, grid$ME, grid$tissue)) %>%
  dplyr::mutate(padj_within = p.adjust(p_within, "BH"), padj_naive = p.adjust(p_naive, "BH"))

# ---------------------------------------------------------------------------
# 2) Comparison with the original association (all 231 samples)
#   survives = significant in both, same sign; reversed = significant in both,
#   opposite sign; lost = only original; new = only within project.
# "survives" supports calling the module tissue-associated; "lost" means the
# original association came from differences between studies (stage, protocol,
# laboratory) and is described as exploratory.
# ---------------------------------------------------------------------------
original <- readr::read_csv(f_assoc, show_col_types = FALSE) %>%
  dplyr::transmute(ME, tissue = trait, cor_original = cor, padj_original = padj_BH_within_trait)
comparison <- within %>%
  dplyr::left_join(original, by = c("ME", "tissue")) %>%
  dplyr::mutate(sig_orig = padj_original < ALPHA, sig_within = padj_within < ALPHA,
                class = dplyr::case_when(
                  is.na(p_within) ~ "not_testable",
                  sig_orig & sig_within & sign(cor_original) == sign(estimate_within) ~ "survives_within_project",
                  sig_orig & sig_within ~ "reversed_within_project",
                  sig_orig & !sig_within ~ "lost_within_project",
                  !sig_orig & sig_within ~ "new_within_project",
                  TRUE ~ "ns_both")) %>%
  dplyr::arrange(ME, tissue)
readr::write_csv(comparison, file.path(dirs$res, "02_module_tissue_within_vs_original.csv"))
(class_summary <- dplyr::count(comparison, class))
readr::write_csv(class_summary, file.path(dirs$res, "03_class_summary.csv"))
# expected: of 72 original associations, 37 survive, 30 lost, 5 reversed; 24 new.
# Seed coat and embryo survive; whole "seed" is lost (study/protocol, not tissue).

# ---------------------------------------------------------------------------
# 3) Tissue context of the candidate modules
# ---------------------------------------------------------------------------
cand_context <- readr::read_csv(f_anchor, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(module_for_integration), module_for_integration != "grey") %>%
  dplyr::transmute(gene_id, module_for_integration, ME = paste0("ME", module_for_integration), anchor_class) %>%
  dplyr::left_join(comparison %>%
                     dplyr::filter(class %in% c("survives_within_project", "new_within_project")) %>%
                     dplyr::group_by(ME) %>%
                     dplyr::summarise(tissues_within_project = paste(sprintf("%s (%+.3f)", tissue, estimate_within),
                                                                     collapse = "; ")),
                   by = "ME")
readr::write_csv(cand_context, file.path(dirs$res, "04_candidate_module_tissue_context.csv"))

lab <- c(survives_within_project = "S", lost_within_project = "L", reversed_within_project = "R",
         new_within_project = "N", ns_both = "", not_testable = "")
p_hm <- comparison %>%
  dplyr::mutate(mark = unname(lab[class]), ME = factor(ME, rev(sort(unique(ME))))) %>%
  ggplot2::ggplot(ggplot2::aes(tissue, ME, fill = t_within)) +
  ggplot2::geom_tile(colour = "grey85") +
  ggplot2::geom_text(ggplot2::aes(label = mark), size = 3) +
  ggplot2::scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue", midpoint = 0, na.value = "grey95") +
  ggplot2::labs(x = "Body part", y = "Module eigengene", fill = "t within\nBioProject",
                subtitle = "S = survives, L = lost, R = reversed, N = new (vs original association)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(dirs$fig, "01_within_project_heatmap.png"), p_hm, width = 9, height = 9, dpi = 300)

close_stage("06_evidence_integration_tissue", dirs)
