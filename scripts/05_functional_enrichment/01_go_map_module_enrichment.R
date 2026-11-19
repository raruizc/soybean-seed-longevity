##############################################################################
# Soybean seed longevity
# Stage 5a: GO annotation map and GO enrichment of the network modules
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: describe which biological functions each module carries, using the GO
# terms of the Wm82.a4 annotation (SoyBase info_annot, column GO). Each gene
# annotated to a term also counts for all ancestors of that term (GO.db), so
# that general terms are not missed.
# Over-representation analysis (ORA, one-sided hypergeometric test):
#   P = P(X >= x), X ~ Hypergeometric(K term genes, N - K, n set genes)
#   fold enrichment = (x / n) / (K / N)
# Background for modules = all genes that entered the network (23,648), not
# the genome. BP, MF and CC are tested separately, with BH within each gene
# set x ontology. Terms with 10 to 500 genes in the background and overlap >= 2.
# The direct GO of each candidate (no ancestors) is listed as annotation, not
# as enrichment. Enrichment describes; it does not show causality.
#
# Run from the repository root, after stage 4:
#   source("scripts/05_functional_enrichment/01_go_map_module_enrichment.R")
#
# Packages: readr, dplyr, tidyr, tibble, stringr, GO.db, AnnotationDbi
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "stringr"))
for (p in c("GO.db", "AnnotationDbi")) if (!requireNamespace(p, quietly = TRUE)) stop("Install first: ", p)

stage <- "05_functional_enrichment"
dirs <- setup_stage(stage)
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_cand <- from_stage("04_coexpression_network", "33_candidate_network_evidence_39.csv")
check_inputs(c(RAW$info_annot, f_mod, f_cand))

modules <- readr::read_csv(f_mod, show_col_types = FALSE)
cand <- readr::read_csv(f_cand, show_col_types = FALSE)
network_genes <- unique(modules$gene_id)

# ---------------------------------------------------------------------------
# 1) Gene-GO map with ancestors
# ---------------------------------------------------------------------------
annot <- readr::read_tsv(RAW$info_annot, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
  dplyr::rename(locusName = dplyr::any_of(c("locusName", "#locusName")))
direct_go <- annot %>%
  dplyr::transmute(gene_id = stringr::str_extract(locusName, "Glyma\\.\\d{2}G\\d+"), GO_ID = GO) %>%
  dplyr::filter(!is.na(gene_id), !is.na(GO_ID), GO_ID != "") %>%
  tidyr::separate_rows(GO_ID, sep = "[,;\\s]+") %>%
  dplyr::filter(stringr::str_detect(GO_ID, "^GO:[0-9]{7}$")) %>%
  dplyr::distinct()

# Term name and ontology of a set of GO IDs (obsolete IDs dropped)
go_info <- function(ids) {
  AnnotationDbi::select(GO.db::GO.db, keys = unique(ids), keytype = "GOID", columns = c("TERM", "ONTOLOGY")) %>%
    dplyr::transmute(GO_ID = GOID, GO_term = TERM, ontology = ONTOLOGY) %>%
    dplyr::filter(!is.na(GO_term), !is.na(ontology)) %>%
    dplyr::distinct()
}
direct_go <- dplyr::inner_join(direct_go, go_info(direct_go$GO_ID), by = "GO_ID")

# Ancestors of each term, looked up once per unique term and ontology
ANCESTOR_DB <- list(BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
terms <- dplyr::distinct(direct_go, GO_ID, ontology)
anc <-dplyr::bind_rows(lapply(names(ANCESTOR_DB), function(o) {
  ids <- terms$GO_ID[terms$ontology == o]
  a <- AnnotationDbi::mget(ids, ANCESTOR_DB[[o]], ifnotfound = NA)
  tibble::tibble(GO_ID = rep(ids, lengths(a)), ancestor = unlist(a, use.names = FALSE))
})) %>%
  dplyr::filter(!is.na(ancestor), ancestor != "all")
expanded <- dplyr::bind_rows(
  dplyr::select(direct_go, gene_id, GO_ID),
  direct_go %>% dplyr::inner_join(anc, by = "GO_ID", relationship = "many-to-many") %>%
    dplyr::transmute(gene_id, GO_ID = ancestor)) %>%
  dplyr::distinct()
go_map <- dplyr::inner_join(expanded, go_info(expanded$GO_ID), by = "GO_ID")
saveRDS(go_map, file.path("data", "processed", "GO_gene_map_expanded.rds"))

(coverage <- tibble::tibble(
  gene_set = c("annotation_genes", "network_genes", "candidates_39", "candidates_in_network"),
  n_total = c(dplyr::n_distinct(annot$locusName), length(network_genes), nrow(cand),
              sum(cand$network_status == "direct_member_of_network")),
  n_with_GO = c(dplyr::n_distinct(direct_go$gene_id), sum(network_genes %in% go_map$gene_id),
                sum(cand$gene_id %in% go_map$gene_id),
                sum(cand$gene_id[cand$network_status == "direct_member_of_network"] %in% go_map$gene_id))) %>%
    dplyr::mutate(fraction_with_GO = n_with_GO / n_total))
readr::write_csv(coverage, file.path(dirs$res, "01_GO_annotation_coverage.csv"))

# ---------------------------------------------------------------------------
# 2) Enrichment of every informative module (background = network)
# ---------------------------------------------------------------------------
bio_modules <- setdiff(unique(modules$module), "grey")
module_go <- dplyr::bind_rows(lapply(bio_modules, function(m) {
  ora_go(modules$gene_id[modules$module == m], network_genes, m, go_map, min_size = 10, max_size = 500)
}))
readr::write_csv(module_go, file.path(dirs$res, "02_GO_module_enrichment_all.csv"))

# Modules relevant to candidates: holding a candidate, or receiving a post hoc
# projection with r >= 0.60 and padj < 0.05
relevant <- sort(unique(c(
  cand$direct_module[cand$network_status == "direct_member_of_network" & cand$direct_module != "grey"],
  cand$strongest_positive_module[!is.na(cand$strongest_positive_r) & cand$strongest_positive_r >= 0.60 &
                                   cand$strongest_positive_padj < ALPHA])))
relevant
module_go_relevant <- dplyr::filter(module_go, gene_set %in% relevant, padj_BH < ALPHA) %>%
  dplyr::arrange(gene_set, ontology, padj_BH)
readr::write_csv(module_go_relevant, file.path(dirs$res, "03_GO_candidate_relevant_modules_significant.csv"))

# ---------------------------------------------------------------------------
# 3) Direct GO of the 39 candidates (annotation, no ancestors)
# ---------------------------------------------------------------------------
cand_go <- tibble::tibble(gene_id = cand$gene_id) %>%
  dplyr::left_join(dplyr::select(direct_go, gene_id, GO_ID, GO_term, ontology), by = "gene_id") %>%
  dplyr::arrange(gene_id, ontology, GO_ID)
readr::write_csv(cand_go, file.path(dirs$res, "04_candidate_direct_GO_39.csv"))

cat("Modules tested:", length(bio_modules), "| with significant GO:",
    dplyr::n_distinct(module_go$gene_set[module_go$padj_BH < ALPHA]),
    "| candidate-relevant:", length(relevant), "\n")

close_stage("05_functional_enrichment_modules", dirs)
