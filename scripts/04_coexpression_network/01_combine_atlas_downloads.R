##############################################################################
# Soybean seed longevity
# Stage 4a: combine the six seed-related downloads of the expression atlas
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: the network is built on genome-wide expression of seed-related samples,
# not on the 39 candidates alone (a 39-gene network would only say how the
# candidates correlate with each other). In the Soybean Expression Atlas v2,
# "Download by body part" with TPM at gene level gives one matrix and one
# metadata file per body part. Save them in data/raw/network/ as
#   SEA_<part>_TPM.tsv.gz and SEA_<part>_metadata.tsv(.gz)
# with <part> in seed, embryo, cotyledon, endosperm, seed_coat, suspensor.
# This script joins them into one matrix and one metadata table, keeping the
# body part of each download. TPM is used because this is a coexpression
# analysis, not a differential test; BioNERO accepts normalized expression.
#
# Run from the repository root:
#   source("scripts/04_coexpression_network/01_combine_atlas_downloads.R")
#
# Packages: readr, dplyr, purrr
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "purrr"))

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
PARTS <- c("seed", "embryo", "cotyledon", "endosperm", "seed_coat", "suspensor")

# First existing file among <stem>.tsv.gz, <stem>.gz and <stem>.tsv
find_input <- function(stem) {
  f <- file.path(RAW$network_dir, paste0(stem, c(".tsv.gz", ".gz", ".tsv")))
  if (!any(file.exists(f))) stop("Not found: ", stem, " (.tsv.gz, .gz or .tsv) in ", RAW$network_dir)
  f[file.exists(f)][1]
}
expr_files <- vapply(PARTS, function(p) find_input(paste0("SEA_", p, "_TPM")), character(1))
meta_files <- vapply(PARTS, function(p) find_input(paste0("SEA_", p, "_metadata")), character(1))

# Read one TPM download: first column = gene, empty columns dropped, numeric check
read_part <- function(f) {
  x <- read_table_df(f)
  x <- x[, !vapply(x, function(z) all(is.na(z) | trimws(as.character(z)) == ""), logical(1)), drop = FALSE]
  names(x)[1] <- "gene_id"
  x[-1] <- lapply(x[-1], as.numeric)
  stopifnot(!anyNA(x[-1]), !anyDuplicated(x$gene_id))
  message(basename(f), ": ", nrow(x), " genes x ", ncol(x) - 1, " samples")
  x
}

# ---------------------------------------------------------------------------
# 1) TPM matrices: full join by gene; a gene absent from one download = 0
# ---------------------------------------------------------------------------
expr <- purrr::reduce(lapply(expr_files, read_part), dplyr::full_join, by = "gene_id")
expr[is.na(expr)] <- 0
stopifnot(!anyDuplicated(names(expr)))         # a sample cannot appear in two body parts
readr::write_tsv(expr, file.path("data", "processed", "SEA_seedrelated_gene_TPM.tsv.gz"))

# ---------------------------------------------------------------------------
# 2) Metadata, with the body part of the download recorded explicitly
# ---------------------------------------------------------------------------
meta <- purrr::map2_dfr(meta_files, PARTS, function(f, p) {
  dplyr::mutate(read_table_df(f), dplyr::across(dplyr::everything(), as.character),
                download_body_part = gsub("_", " ", p))
}) %>%
  dplyr::distinct()
readr::write_tsv(meta, file.path("data", "processed", "SEA_seedrelated_metadata.tsv.gz"))

cat("Combined matrix:", nrow(expr), "genes x", ncol(expr) - 1, "samples\n")   # expected: 52,837 x 933
cat("Metadata:", nrow(meta), "rows\n")
