##############################################################################
# Soybean seed longevity
# Stage 5c: GO redundancy reduction, background concordance and functional
# profile of the 39 candidates
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: summarize the enrichment results without re-running them.
# Redundancy: many significant GO terms are parent/child terms supported by
# almost the same genes. Within each gene set x ontology, terms are sorted by
# padj; the first unassigned term becomes the representative, and every term
# whose supporting genes have Jaccard >= 0.50 with it is collapsed into it:
#   Jaccard(A, B) = |A n B| / |A u B|
# The representative is still a real GO term (no new categories).
# Concordance of the neighbourhood backgrounds: robust = significant against
# the module AND the network; within_module_only = local specialization;
# network_only = unusual in the network but not special inside the module.
#
# Run from the repository root, after scripts 01 and 02 of this stage:
#   source("scripts/05_functional_enrichment/03_redundancy_functional_profile.R")
#
# Packages: readr, dplyr, tidyr, tibble, purrr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "purrr", "ggplot2"))

stage <- "05_functional_enrichment"
dirs <- setup_stage(stage)
f_modgo <- file.path(dirs$res, "03_GO_candidate_relevant_modules_significant.csv")
f_within <- file.path(dirs$res, "05_GO_neighborhood_within_module.csv")
f_network <- file.path(dirs$res, "06_GO_neighborhood_network_background.csv")
f_cgo <- file.path(dirs$res, "04_candidate_direct_GO_39.csv")
f_net <- from_stage("04_coexpression_network", "33_candidate_network_evidence_39.csv")
check_inputs(c(f_modgo, f_within, f_network, f_cgo, f_net))

JACCARD_MIN <- 0.50
TOP_THEMES <- 5
module_go <- readr::read_csv(f_modgo, show_col_types = FALSE)
nb_within <- readr::read_csv(f_within, show_col_types = FALSE)
nb_network <- readr::read_csv(f_network, show_col_types = FALSE)

# ---------------------------------------------------------------------------
# 1) Greedy reduction by supporting-gene Jaccard
# ---------------------------------------------------------------------------
jaccard <- function(a, b) if (!length(a) && !length(b)) 0 else length(intersect(a, b)) / length(union(a, b))
jaccard(c("g1", "g2"), c("g2", "g3"))          # check: 1/3

# Cluster the significant terms of one gene set x ontology
reduce_terms <- function(df) {
  df <- dplyr::arrange(df, padj_BH, dplyr::desc(fold_enrichment), dplyr::desc(overlap_n))
  genes <- strsplit(dplyr::coalesce(df$overlap_genes, ""), ";", fixed = TRUE)
  cluster <- integer(nrow(df))
  for (i in seq_len(nrow(df))) {
    if (cluster[i] > 0) next
    cluster[i] <- max(cluster) + 1L
    free <- which(cluster == 0)
    hit <- free[vapply(free, function(j) jaccard(genes[[i]], genes[[j]]) >= JACCARD_MIN, logical(1))]
    cluster[hit] <- cluster[i]
  }
  df$redundancy_cluster <- cluster
  df$is_representative <- !duplicated(cluster)
  reps <- df[df$is_representative, ]
  df %>%
    dplyr::left_join(dplyr::transmute(reps, redundancy_cluster, representative_GO_ID = GO_ID,
                                      representative_GO_term = GO_term, representative_padj = padj_BH,
                                      representative_overlap_n = overlap_n), by = "redundancy_cluster") %>%
    dplyr::group_by(redundancy_cluster) %>%
    dplyr::mutate(n_collapsed_terms = dplyr::n(), collapsed_GO_terms = paste(GO_term, collapse = " | ")) %>%
    dplyr::ungroup()
}

# Apply reduce_terms within each group of the given columns
reduce_by <- function(df, cols) {
  df %>%
    dplyr::filter(padj_BH < ALPHA) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(cols))) %>%
    dplyr::group_modify(~ reduce_terms(.x)) %>%
    dplyr::ungroup()
}
module_red <- reduce_by(module_go, c("gene_set", "ontology"))
nb_red <- reduce_by(nb_within, c("candidate_gene", "ontology"))
readr::write_csv(module_red, file.path(dirs$res, "09_GO_module_terms_reduced.csv"))
readr::write_csv(nb_red, file.path(dirs$res, "10_GO_neighborhood_terms_reduced.csv"))
module_rep <- dplyr::filter(module_red, is_representative)
nb_rep <- dplyr::filter(nb_red, is_representative)

# ---------------------------------------------------------------------------
# 2) Concordance of the two neighbourhood backgrounds (exact GO_ID)
# ---------------------------------------------------------------------------
pick <- function(df, prefix) {
  x <- dplyr::select(df, candidate_gene, module, ontology, GO_ID, GO_term, padj_BH, fold_enrichment, overlap_n)
  names(x)[6:8] <- paste0(prefix, c("_padj", "_fold_enrichment", "_overlap_n"))
  x
}
concordance <- dplyr::full_join(pick(nb_within, "within"), pick(nb_network, "network"),
                                by = c("candidate_gene", "module", "ontology", "GO_ID", "GO_term")) %>%
  dplyr::mutate(sig_within = dplyr::coalesce(within_padj < ALPHA, FALSE),
                sig_network = dplyr::coalesce(network_padj < ALPHA, FALSE),
                support_class = dplyr::case_when(sig_within & sig_network ~ "robust_both_backgrounds",
                                                 sig_within ~ "within_module_only",
                                                 sig_network ~ "network_only",
                                                 TRUE ~ "not_significant")) %>%
  dplyr::arrange(candidate_gene, ontology, support_class, GO_ID)
readr::write_csv(concordance, file.path(dirs$res, "11_GO_neighborhood_background_concordance.csv"))

# ---------------------------------------------------------------------------
# 3) Functional profile of the 39 candidates (input of stage 6)
# ---------------------------------------------------------------------------

# Top representative themes per group and ontology as one text column each
themes_wide <- function(rep, group_col, prefix) {
  rep %>%
    dplyr::arrange(representative_padj) %>%
    dplyr::group_by(key = .data[[group_col]], ontology) %>%
    dplyr::summarise(themes = paste(head(representative_GO_term, TOP_THEMES), collapse = " | "),
                     n_themes = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = ontology, values_from = c(themes, n_themes),
                       names_glue = paste0(prefix, "_{ontology}_{.value}"))
}
module_themes <- themes_wide(module_rep, "gene_set", "module")
nb_themes <- themes_wide(nb_rep, "candidate_gene", "neighborhood")

direct_go <- readr::read_csv(f_cgo, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(GO_ID)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(direct_GO_n = dplyr::n_distinct(GO_ID), direct_GO_terms = paste(unique(GO_term), collapse = " | "))
conc_by_cand <- concordance %>%
  dplyr::filter(support_class != "not_significant") %>%
  dplyr::group_by(gene_id = candidate_gene) %>%
  dplyr::summarise(n_GO_robust_both = sum(support_class == "robust_both_backgrounds"),
                   n_GO_within_module_only = sum(support_class == "within_module_only"),
                   n_GO_network_only = sum(support_class == "network_only"),
                   robust_GO_terms = paste(unique(GO_term[support_class == "robust_both_backgrounds"]), collapse = " | "))

profile <- readr::read_csv(f_net, show_col_types = FALSE) %>%
  dplyr::mutate(functional_module = dplyr::case_when(!is.na(direct_module) & direct_module != "grey" ~ direct_module,
                                                     !is.na(strongest_positive_module) ~ strongest_positive_module)) %>%
  dplyr::left_join(module_themes, by = c("functional_module" = "key")) %>%
  dplyr::left_join(nb_themes, by = c("gene_id" = "key")) %>%
  dplyr::left_join(direct_go, by = "gene_id") %>%
  dplyr::left_join(conc_by_cand, by = "gene_id")
readr::write_csv(profile, file.path(dirs$res, "12_candidate_functional_profile_39.csv"))

(summary_5 <- tibble::tibble(
  metric = c("significant_module_terms", "module_themes_after_reduction", "significant_neighborhood_terms",
             "neighborhood_themes_after_reduction", "candidates_robust_both_backgrounds",
             "candidates_within_module_only", "candidates_network_only"),
  value = c(nrow(module_red), nrow(module_rep), nrow(nb_red), nrow(nb_rep),
            dplyr::n_distinct(concordance$candidate_gene[concordance$support_class == "robust_both_backgrounds"]),
            dplyr::n_distinct(concordance$candidate_gene[concordance$support_class == "within_module_only"]),
            dplyr::n_distinct(concordance$candidate_gene[concordance$support_class == "network_only"]))))
readr::write_csv(summary_5, file.path(dirs$res, "13_functional_summary.csv"))
# expected: 438 -> 162 module themes; 145 -> 40 neighbourhood themes; 12 / 7 / 8 candidates

p_mod <- module_rep %>%
  dplyr::filter(ontology == "BP") %>%
  dplyr::group_by(gene_set) %>%
  dplyr::slice_min(representative_padj, n = TOP_THEMES, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(score = -log10(pmax(representative_padj, .Machine$double.xmin)),
                term = paste0(representative_GO_term, "___", gene_set)) %>%
  ggplot2::ggplot(ggplot2::aes(score, reorder(term, score), size = representative_overlap_n)) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ gene_set, scales = "free_y") +
  ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  ggplot2::labs(x = expression(-log[10]("BH-adjusted P")), y = "Representative GO biological process", size = "Genes") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "01_module_BP_themes.png"), p_mod, width = 16, height = 12, dpi = 300)

close_stage("05_functional_enrichment_profile", dirs)
