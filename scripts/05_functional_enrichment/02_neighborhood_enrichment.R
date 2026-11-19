##############################################################################
# Soybean seed longevity
# Stage 5b: GO enrichment of candidate neighbourhoods and BH-family sensitivity
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: for each candidate inside an informative module, test whether its local
# neighbourhood (candidate + 30 genes of highest adjacency, stage 4) is
# functionally specialized.
#   A) primary: background = the candidate's own module (is the neighbourhood
#      special WITHIN the module?);
#   B) complementary: background = whole network (is the function unusual in
#      the network?).
# Terms with 5 to 500 genes in the background and overlap >= 2; BH per
# ontology within each neighbourhood.
# Sensitivity: filtering terms by overlap >= 2 BEFORE BH makes the BH family
# depend on the outcome. The within-module test is repeated with three
# families: x2 (overlap >= 2, as above), x1 (overlap >= 1, clusterProfiler
# convention) and all (every testable term). A term counts as significant only
# with padj < 0.05 AND overlap >= 2 in all three.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/05_functional_enrichment/02_neighborhood_enrichment.R")
#
# Packages: readr, dplyr, tidyr, tibble
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble"))

stage <- "05_functional_enrichment"
dirs <- setup_stage(stage)
f_go <- file.path("data", "processed", "GO_gene_map_expanded.rds")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_nb <- from_stage("04_coexpression_network", "31_candidate_top30_neighbors.csv")
check_inputs(c(f_go, f_mod, f_nb))

go_map <- readRDS(f_go)
modules <- readr::read_csv(f_mod, show_col_types = FALSE)
network_genes <- unique(modules$gene_id)
sets <- readr::read_csv(f_nb, show_col_types = FALSE) %>%
  dplyr::group_by(candidate_gene, module) %>%
  dplyr::summarise(genes = list(unique(c(candidate_gene[1], neighbor_gene))), .groups = "drop")
nrow(sets)                                     # expected: 16 neighbourhoods

# ---------------------------------------------------------------------------
# 1) Neighbourhood ORA against the module and against the network
# ---------------------------------------------------------------------------

# Both ORAs for neighbourhood i of `sets`
neighborhood_ora <- function(i) {
  bg_mod <- modules$gene_id[modules$module == sets$module[i]]
  tag <- function(x, bg) dplyr::mutate(x, candidate_gene = sets$candidate_gene[i], module = sets$module[i],
                                        background_type = bg)
  list(within = tag(ora_go(sets$genes[[i]], bg_mod, sets$candidate_gene[i], go_map, 5, min(500, length(bg_mod))),
                    "within_module"),
       network = tag(ora_go(sets$genes[[i]], network_genes, sets$candidate_gene[i], go_map, 5, 500),
                     "whole_network"))
}
res <- lapply(seq_len(nrow(sets)), neighborhood_ora)
nb_within <- dplyr::bind_rows(lapply(res, `[[`, "within"))
nb_network <- dplyr::bind_rows(lapply(res, `[[`, "network"))
readr::write_csv(nb_within, file.path(dirs$res, "05_GO_neighborhood_within_module.csv"))
readr::write_csv(nb_network, file.path(dirs$res, "06_GO_neighborhood_network_background.csv"))

# ---------------------------------------------------------------------------
# 2) BH-family sensitivity (within-module background)
# ---------------------------------------------------------------------------

# All testable terms of one neighbourhood, with BH in three families
ora_three_families <- function(genes, bg) {
  bg <- intersect(unique(bg), go_map$gene_id)
  genes <- intersect(unique(genes), bg)
  go_map %>%
    dplyr::filter(gene_id %in% bg) %>%
    dplyr::group_by(ontology, GO_ID, GO_term) %>%
    dplyr::summarise(K = dplyr::n_distinct(gene_id), x = dplyr::n_distinct(gene_id[gene_id %in% genes]),
                     .groups = "drop") %>%
    dplyr::filter(K >= 5, K <= min(500, length(bg))) %>%
    dplyr::mutate(p = phyper(x - 1, K, length(bg) - K, length(genes), lower.tail = FALSE)) %>%
    dplyr::group_by(ontology) %>%
    dplyr::mutate(padj_x2 = ifelse(x >= 2, p.adjust(ifelse(x >= 2, p, NA), "BH"), NA),
                  padj_x1 = ifelse(x >= 1, p.adjust(ifelse(x >= 1, p, NA), "BH"), NA),
                  padj_all = p.adjust(p, "BH")) %>%
    dplyr::ungroup()
}
families <- dplyr::bind_rows(lapply(seq_len(nrow(sets)), function(i) {
  dplyr::mutate(ora_three_families(sets$genes[[i]], modules$gene_id[modules$module == sets$module[i]]),
                candidate_gene = sets$candidate_gene[i], module = sets$module[i], .before = 1)
})) %>%
  dplyr::mutate(sig_x2 = dplyr::coalesce(padj_x2 < ALPHA & x >= 2, FALSE),
                sig_x1 = dplyr::coalesce(padj_x1 < ALPHA & x >= 2, FALSE),
                sig_all = dplyr::coalesce(padj_all < ALPHA & x >= 2, FALSE))
readr::write_csv(dplyr::filter(families, sig_x2 | sig_x1 | sig_all),
                 file.path(dirs$res, "07_neighborhood_terms_three_BH_families.csv"))

# The top term is taken BEFORE the counts: inside summarise, a recreated
# column (n_sig_all = sum(sig_all)) would otherwise replace the original
(family_summary <- families %>%
    dplyr::group_by(candidate_gene, module) %>%
    dplyr::summarise(n_terms_tested = dplyr::n(),
                     top_term_all = if (any(sig_all)) GO_term[sig_all][which.min(padj_all[sig_all])] else NA_character_,
                     n_sig_x2 = sum(sig_x2), n_sig_x1 = sum(sig_x1), n_sig_all = sum(sig_all), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig_all), dplyr::desc(n_sig_x1)))
readr::write_csv(family_summary, file.path(dirs$res, "08_BH_family_summary_by_candidate.csv"))
cat("Neighbourhoods with >= 1 term: x2 =", sum(family_summary$n_sig_x2 > 0), "| x1 =",
    sum(family_summary$n_sig_x1 > 0), "| all =", sum(family_summary$n_sig_all > 0), "of", nrow(family_summary), "\n")
# expected: 12 / 10 / 6 of 16. Terms that survive "all" are the most robust;
# a neighbourhood with terms only in x2 depends on the correction convention
# and is reported as exploratory.

close_stage("05_functional_enrichment_neighborhoods", dirs)
