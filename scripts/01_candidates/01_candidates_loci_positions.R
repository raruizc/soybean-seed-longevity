##############################################################################
# Soybean seed longevity
# Stage 1: candidate genes, GWAS loci and physical position in Wm82.a4
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: start from the curated candidate list (gwas_candidates_curated.xlsx,
# sheet Master_Candidates, IDs already converted to Wm82.a4) and build:
#   1) the table of the 39 candidates with study, source IDs and ID system;
#   2) the GWAS locus map (Thenveettil: SNPs < 250 kb apart on the same
#      chromosome; Wu: one locus per QTN; gene without position = own locus);
#   3) the position of each gene in Wm82.a4 (GFF) and the span of each locus;
#   4) the Wm82.a4 genes actually within +-10 kb of each Thenveettil SNP and
#      the distance from each listed gene to its SNP.
# Item 4 exists because Thenveettil et al. (2025) called the SNPs on Wm82.a4
# (Glycine_max_v4) but looked up the genes in the window on Wm82.a2, a1 and
# ZH13 using the same position number. Where the assemblies differ, the listed
# gene ends up far from the SNP.
#
# Run from the repository root:
#   source("scripts/01_candidates/01_candidates_loci_positions.R")
#
# Packages: readr, readxl, dplyr, tidyr, stringr
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "readxl", "dplyr", "tidyr", "stringr"))

stage <- "01_candidates"
dirs <- setup_stage(stage)
check_inputs(c(RAW$gwas_xlsx, RAW$gff, RAW$info_annot))

# ---------------------------------------------------------------------------
# 1) Harmonized candidates
# ---------------------------------------------------------------------------

# Naming system inferred from the format of the source ID
id_system <- function(id) dplyr::case_when(
  stringr::str_detect(id, "^Glyma\\d{2}g") ~ "Wm82.a1",
  stringr::str_detect(id, "^SoyZH13") ~ "ZH13",
  stringr::str_detect(id, "^Glyma\\.\\d{2}G") ~ "Wm82.a2/a4",
  TRUE ~ "other")

# Smallest numeric position in strings such as "18738006; 18738006, 18737804"
min_pos <- function(x) {
  v <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x), "[;, ]+"))))
  if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
}
min_pos("18738006; 18738006, 18737804")   # check: 18737804

master <- readxl::read_excel(RAW$gwas_xlsx, sheet = "Master_Candidates") %>%
  dplyr::filter(!is.na(Wm82_a4_id)) %>%
  dplyr::mutate(chromosome = as.integer(readr::parse_number(as.character(chromosome))),
                study = dplyr::if_else(stringr::str_detect(study, "Wu"), "Wu 2024", "Thenveettil 2025"),
                snp_txt = as.character(snp_or_positions),
                pos = dplyr::if_else(study == "Wu 2024", NA_real_, vapply(snp_txt, min_pos, numeric(1))))

candidates <- master %>%
  dplyr::group_by(gene_id = Wm82_a4_id) %>%
  dplyr::summarise(study = paste(unique(study), collapse = "; "), chromosome = dplyr::first(chromosome),
                   source_ids = paste(unique(source_gene_id), collapse = "; "),
                   id_systems = paste(unique(id_system(source_gene_id)), collapse = "; "),
                   snp_pos = if (all(is.na(pos))) NA_real_ else min(pos, na.rm = TRUE),
                   .groups = "drop")
(n_genes <- nrow(candidates))                  # expected: 39
readr::write_csv(candidates, file.path(dirs$res, "01_candidates_39.csv"))

# ---------------------------------------------------------------------------
# 2) GWAS locus map (coordinates of the source study, used only for grouping)
# ---------------------------------------------------------------------------

# Cluster label: a new locus starts when the chromosome changes or the gap exceeds the window
cluster_pos <- function(chr, pos, window) {
  o <- order(chr, pos)
  new <- c(TRUE, diff(pos[o]) > window | diff(chr[o]) != 0)
  out <- integer(length(chr))
  out[o] <- cumsum(new)
  out
}

th <- dplyr::filter(candidates, study == "Thenveettil 2025")
th_pos <- th %>%
  dplyr::filter(!is.na(snp_pos)) %>%
  dplyr::mutate(cl = cluster_pos(chromosome, snp_pos, LOCUS_WIN_BP)) %>%
  dplyr::group_by(cl) %>%
  dplyr::mutate(locus_pos_bp = min(snp_pos),
                locus_id = sprintf("Th_chr%02d_%.2fMb", chromosome[1], locus_pos_bp / 1e6)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-cl)
th_no_pos <- th %>%
  dplyr::filter(is.na(snp_pos)) %>%
  dplyr::mutate(locus_pos_bp = NA_real_, locus_id = sprintf("Th_chr%02d_noPosition_%s", chromosome, gene_id))
wu <- master %>%
  dplyr::filter(study == "Wu 2024") %>%
  dplyr::distinct(gene_id = Wm82_a4_id, snp_txt) %>%
  dplyr::inner_join(dplyr::filter(candidates, study == "Wu 2024"), by = "gene_id") %>%
  dplyr::mutate(locus_pos_bp = NA_real_, locus_id = paste0("Wu_chr", sprintf("%02d", chromosome), "_", snp_txt)) %>%
  dplyr::select(-snp_txt)

locus_map <- dplyr::bind_rows(th_pos, th_no_pos, wu) %>%
  dplyr::group_by(locus_id) %>%
  dplyr::mutate(n_genes_locus = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_id, study, chromosome, snp_pos, locus_pos_bp, locus_id, n_genes_locus) %>%
  dplyr::arrange(chromosome, locus_id, gene_id)
stopifnot(!anyDuplicated(locus_map$gene_id), setequal(locus_map$gene_id, candidates$gene_id))
(n_loci <- dplyr::n_distinct(locus_map$locus_id))   # expected: 21
readr::write_csv(locus_map, file.path(dirs$res, "02_locus_map.csv"))

# ---------------------------------------------------------------------------
# 3) Candidate positions in Wm82.a4 and locus span
# ---------------------------------------------------------------------------
gff <- read_gff_genes()
annot <- read_info_annot() %>%
  dplyr::select(gene_id, arabi_hit = `Best-hit-arabi-name`, arabi_defline = `arabi-defline`)

positions <- locus_map %>%
  dplyr::left_join(dplyr::select(candidates, gene_id, source_ids, id_systems), by = "gene_id") %>%
  dplyr::left_join(dplyr::rename(gff, gff_chr = chromosome), by = "gene_id") %>%
  dplyr::mutate(mid_mb = (start + end) / 2 / 1e6, dist_snp_kb = dist_kb(snp_pos, start, end))
cat("Genes with GFF coordinates:", sum(!is.na(positions$start)), "of", nrow(positions), "\n")
readr::write_csv(positions, file.path(dirs$res, "03_candidate_positions_a4.csv"))

# Reading: a small span means the genes of a locus are also neighbours in
# Wm82.a4. A large span with IDs from different systems indicates that the
# source list mixed assemblies (see section 4).
(span <- positions %>%
   dplyr::filter(n_genes_locus >= 2) %>%
   dplyr::group_by(locus_id, study) %>%
   dplyr::summarise(n_genes = dplyr::n(), genes = paste(gene_id, collapse = "; "),
                    id_systems = paste(unique(unlist(strsplit(id_systems, "; "))), collapse = "; "),
                    span_mb = max(mid_mb, na.rm = TRUE) - min(mid_mb, na.rm = TRUE), .groups = "drop") %>%
   dplyr::arrange(dplyr::desc(span_mb)))
readr::write_csv(span, file.path(dirs$res, "04_locus_span.csv"))

# ---------------------------------------------------------------------------
# 4) Wm82.a4 genes in the window of each Thenveettil SNP and verdict per SNP
# ---------------------------------------------------------------------------

# a4 genes overlapping the +-SNP_WIN_BP window around a position
genes_in_window <- function(chr, pos) {
  gff %>%
    dplyr::filter(chromosome == chr, end >= pos - SNP_WIN_BP, start <= pos + SNP_WIN_BP) %>%
    dplyr::mutate(dist_kb = dist_kb(pos, start, end)) %>%
    dplyr::select(-chromosome) %>%
    dplyr::arrange(dist_kb)
}
genes_in_window(13, 43132863)                 # check: chr13 SNP

snps <- positions %>%
  dplyr::filter(study == "Thenveettil 2025", !is.na(snp_pos)) %>%
  dplyr::group_by(locus_id, chromosome, snp_pos) %>%
  dplyr::summarise(listed_genes = paste(gene_id, collapse = "; "),
                   min_dist_listed_kb = min(dist_snp_kb, na.rm = TRUE), .groups = "drop")

window <- snps %>%
  dplyr::mutate(a4 = mapply(genes_in_window, chromosome, snp_pos, SIMPLIFY = FALSE)) %>%
  dplyr::select(locus_id, chromosome, snp_pos, a4) %>%
  tidyr::unnest(a4) %>%
  dplyr::left_join(annot, by = "gene_id")
readr::write_csv(window, file.path(dirs$res, "05_a4_genes_in_snp_window.csv"))

# Reading: "consistent" = some listed gene within 50 kb of the SNP in a4;
# "borderline" = 50 to 200 kb; "displaced" = > 200 kb. For displaced SNPs, the
# genes in 05_ are the ones the authors' window would return on Wm82.a4.
(verdict <- snps %>%
   dplyr::left_join(window %>% dplyr::group_by(locus_id, snp_pos) %>%
                      dplyr::summarise(a4_genes_in_window = paste(gene_id, collapse = "; "), .groups = "drop"),
                    by = c("locus_id", "snp_pos")) %>%
   dplyr::mutate(verdict = dplyr::case_when(min_dist_listed_kb <= NEAR_SNP_KB ~ "consistent",
                                            min_dist_listed_kb <= FAR_SNP_KB ~ "borderline",
                                            TRUE ~ "displaced")) %>%
   dplyr::arrange(dplyr::desc(min_dist_listed_kb)))
readr::write_csv(verdict, file.path(dirs$res, "06_snp_verdict.csv"))
print(dplyr::count(verdict, verdict))         # expected: 8 consistent, 2 borderline, 8 displaced

close_stage(stage, dirs)
