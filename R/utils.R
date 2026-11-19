##############################################################################
# Soybean seed longevity
# Helper functions used by more than one stage
# Author: Ricardo Antonio Ruiz Cardozo
##############################################################################

# Stop with a clear message if a package is missing, then attach all of them
load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install first: ", paste(missing, collapse = ", "))
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# Create results/<stage>, figures/<stage>, logs/ and data/processed/ and return the paths
setup_stage <- function(stage) {
  dirs <- list(res = file.path("results", stage), fig = file.path("figures", stage), log = "logs")
  for (d in c(dirs, file.path("data", "processed"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# Stop if any input file is missing, listing which ones
check_inputs <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))
  invisible(TRUE)
}

# Path of an output written by an earlier stage: results/<stage>/<file>
from_stage <- function(stage, file) file.path("results", stage, file)

# Glyma.01G000100.1 -> Glyma.01G000100
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# Rank fraction with lower = better: position / n available, ties averaged
rank_frac <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  out[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
  out
}

# Wm82.a4 annotation: best Arabidopsis hit and GO terms, one row per gene
read_info_annot <- function(f = RAW$info_annot) {
  readr::read_tsv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
    dplyr::rename(locusName = dplyr::any_of(c("locusName", "#locusName"))) %>%
    dplyr::distinct(locusName, .keep_all = TRUE) %>%
    dplyr::rename(gene_id = locusName)
}

# Gene coordinates from the Wm82.a4 GFF (chromosome, start, end)
read_gff_genes <- function(f = RAW$gff) {
  readr::read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE,
                  col_types = readr::cols(.default = "c")) %>%
    dplyr::filter(X3 == "gene") %>%
    dplyr::transmute(gene_id = stringr::str_extract(X9, "Glyma\\.(\\d{2}G\\d+|U\\d+)"),
                     chromosome = as.integer(stringr::str_match(X1, "Gm(\\d+)")[, 2]),
                     start = as.numeric(X4), end = as.numeric(X5)) %>%
    dplyr::filter(!is.na(gene_id), !is.na(chromosome)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

# Distance (kb) from a gene [start, end] to a position; 0 if the position is inside the gene
dist_kb <- function(pos, start, end) pmax(start - pos, pos - end, 0) / 1000

# Read a TSV/CSV (plain or .gz) as a plain data.frame, keeping column names as they are
read_table_df <- function(f) {
  readr::read_delim(f, delim = if (grepl("\\.csv", f)) "," else "\t", show_col_types = FALSE,
                    progress = FALSE, name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE)
}

# Expression table (gene column + samples) -> numeric matrix genes x samples
as_expr_matrix <- function(df, gene_col = 1) {
  m <- as.matrix(df[, -gene_col, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- df[[gene_col]]
  if (anyNA(m)) stop("NA values after numeric conversion")
  m
}

# Save sessionInfo for the stage and list what was written
close_stage <- function(stage, dirs) {
  writeLines(capture.output(sessionInfo()), file.path(dirs$log, paste0("sessionInfo_", stage, ".txt")))
  written <- c(list.files(dirs$res, full.names = TRUE), list.files(dirs$fig, full.names = TRUE))
  cat("\n== ", stage, " done ==\n", sep = "")
  cat(paste0("- ", written, "\n"), sep = "")
}

# GO over-representation of one gene set against a background (stage 5).
# One-sided hypergeometric P = P(X >= overlap); BH within each ontology.
# go_map: one row per gene x GO term (with ancestors), columns gene_id,
# GO_ID, GO_term, ontology. Terms outside [min_size, max_size] genes in the
# background or with overlap < min_overlap are not tested.
ora_go <- function(genes, background, set_name, go_map, min_size, max_size, min_overlap = 2) {
  background <- intersect(unique(background), go_map$gene_id)
  genes <- intersect(unique(genes), background)
  if (length(genes) < 2) return(tibble::tibble())
  go_map %>%
    dplyr::filter(gene_id %in% background) %>%
    dplyr::group_by(ontology, GO_ID, GO_term) %>%
    dplyr::summarise(term_background_n = dplyr::n_distinct(gene_id),
                     overlap_n = dplyr::n_distinct(gene_id[gene_id %in% genes]),
                     overlap_genes = paste(sort(unique(gene_id[gene_id %in% genes])), collapse = ";"),
                     .groups = "drop") %>%
    dplyr::filter(term_background_n >= min_size, term_background_n <= max_size, overlap_n >= min_overlap) %>%
    dplyr::mutate(background_annotated_n = length(background), geneset_annotated_n = length(genes),
                  pvalue = phyper(overlap_n - 1, term_background_n, background_annotated_n - term_background_n,
                                  geneset_annotated_n, lower.tail = FALSE),
                  fold_enrichment = (overlap_n / geneset_annotated_n) / (term_background_n / background_annotated_n),
                  gene_set = set_name) %>%
    dplyr::group_by(ontology) %>%
    dplyr::mutate(padj_BH = p.adjust(pvalue, "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(ontology, padj_BH, dplyr::desc(fold_enrichment))
}
