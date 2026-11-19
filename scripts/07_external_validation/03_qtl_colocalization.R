##############################################################################
# Soybean seed longevity
# Stage 7c: colocalization of the GWAS loci with published QTL for seed
# storability, longevity, viability and vigor
# Author: Ricardo Antonio Ruiz Cardozo
#
# Question: do the 21 loci (and, above all, the convergent chromosome 4
# locus) fall inside QTL mapped by other studies, more often than expected
# for intervals of the same size placed at random?
# Input curated by hand: data/raw/external/qtl/storability_qtl_curated.csv,
# one row per QTL, coordinates already in Wm82.a4 (template and column rules in
# docs/qtl_table_template.csv). Rows on another assembly are reported and left
# out; they have to be converted first (e.g. flanking marker sequences or the
# genes at the interval ends searched in Wm82.a4).
# Locus interval: from the first to the last base of its candidate genes in
# Wm82.a4, extended to the GWAS SNP when the SNP has an a4 position; a gene
# displaced more than FAR_SNP_KB from its SNP is replaced by the SNP window
# (SNP_WIN_BP), so a misplaced gene does not stretch the locus. A locus
# overlaps a QTL when the two intervals share at least 1 bp after the QTL is
# extended by QTL_FLANK_BP on each side (0 = interval as published).
# Null model: each locus keeps its span and chromosome and is moved to a
# uniform random start (N_PERM_QTL draws). P = (1 + draws with as many or more
# overlapping loci) / (1 + N_PERM_QTL). A genome-wide version, chromosome drawn
# in proportion to its length, is reported as well.
#
# Run from the repository root, after stages 1 and 6:
#   source("scripts/07_external_validation/03_qtl_colocalization.R")
#
# Packages: readr, dplyr, tidyr, stringr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "stringr", "ggplot2"))

stage <- "07_external_validation"
dirs <- setup_stage(file.path(stage, "03_qtl"))
f_pos <- from_stage("01_candidates", "03_candidate_positions_a4.csv")
f_cls <- from_stage("06_evidence_integration", "06_final_classes/03_final_locus_table.csv")
check_inputs(c(RAW$qtl_table, f_pos, f_cls, RAW$gff))

# ---------------------------------------------------------------------------
# 1) QTL table, loci intervals and chromosome lengths
# ---------------------------------------------------------------------------
qtl_raw <- readr::read_csv(RAW$qtl_table, show_col_types = FALSE, comment = "#")
need <- c("qtl_id", "source", "trait_group", "chromosome", "start_bp", "end_bp", "assembly")
if (length(setdiff(need, names(qtl_raw)))) stop("Missing columns in the QTL table: ", paste(setdiff(need, names(qtl_raw)), collapse = ", "))
(qtl_off <- dplyr::filter(qtl_raw, assembly != "Wm82.a4" | is.na(start_bp) | is.na(end_bp)))
qtl <- qtl_raw %>%
  dplyr::filter(assembly == "Wm82.a4", !is.na(start_bp), !is.na(end_bp)) %>%
  dplyr::mutate(chromosome = as.integer(chromosome), start = pmin(start_bp, end_bp), end = pmax(start_bp, end_bp),
                span_mb = (end - start) / 1e6)
readr::write_csv(qtl_off, file.path(dirs$res, "01_qtl_left_out.csv"))
qtl %>% dplyr::count(trait_group, source)
# If many rows land in 01_qtl_left_out.csv, the test below runs on a partial
# map and a "no colocalization" reading is weaker than it looks.

# Wm82.a4 chromosome length approximated by the end of its last gene
chr_len <- read_gff_genes() %>%
  dplyr::filter(chromosome %in% 1:20) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(len = max(end), .groups = "drop")

# Genes displaced more than FAR_SNP_KB from their SNP (stage 1) do not stretch
# the interval; a locus left with no gene is the SNP +- SNP_WIN_BP
loci <- readr::read_csv(f_pos, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::mutate(displaced = !is.na(dist_snp_kb) & dist_snp_kb > FAR_SNP_KB,
                start = dplyr::if_else(displaced, snp_pos - SNP_WIN_BP, start),
                end = dplyr::if_else(displaced, snp_pos + SNP_WIN_BP, end)) %>%
  dplyr::group_by(locus_id, chromosome) %>%
  dplyr::summarise(start = min(c(start, snp_pos), na.rm = TRUE), end = max(c(end, snp_pos), na.rm = TRUE),
                   n_displaced = sum(displaced),
                   genes = paste(sort(unique(gene_id)), collapse = "; "), .groups = "drop") %>%
  dplyr::left_join(readr::read_csv(f_cls, show_col_types = FALSE) %>% dplyr::select(locus_id, locus_status),
                   by = "locus_id") %>%
  dplyr::mutate(span_kb = (end - start) / 1000)
nrow(loci)
# expected 21 loci; a locus whose genes have no a4 position cannot be placed

# ---------------------------------------------------------------------------
# 2) Observed overlaps, QTL as published and extended
# ---------------------------------------------------------------------------

# TRUE for each locus [s, e] on chromosome ch that touches at least one QTL
# of the set q after extending the QTL by flank bp on both sides
hits_any <- function(ch, s, e, q, flank) {
  vapply(seq_along(ch), function(i) any(q$chromosome == ch[i] & q$start - flank <= e[i] & q$end + flank >= s[i]), logical(1))
}

# Every locus x QTL pair that overlaps, for the report
pairs_at <- function(flank) {
  loci %>%
    dplyr::inner_join(qtl, by = "chromosome", suffix = c("_locus", "_qtl"), relationship = "many-to-many") %>%
    dplyr::filter(start_qtl - flank <= end_locus, end_qtl + flank >= start_locus) %>%
    dplyr::mutate(flank_mb = flank / 1e6, dist_kb = dist_kb((start_locus + end_locus) / 2, start_qtl, end_qtl)) %>%
    dplyr::select(flank_mb, locus_id, locus_status, chromosome, start_locus, end_locus, qtl_id, source, trait_group,
                  start_qtl, end_qtl, span_mb, dist_kb)
}
(pairs <- dplyr::bind_rows(lapply(QTL_FLANK_BP, pairs_at)))
readr::write_csv(pairs, file.path(dirs$res, "02_locus_qtl_overlaps.csv"))

# ---------------------------------------------------------------------------
# 3) Permutation test, by trait group and for all QTL together
# ---------------------------------------------------------------------------

# Random start for each locus on its own chromosome (or a random chromosome,
# weighted by length, when genome_wide = TRUE), keeping its span
shuffle_loci <- function(genome_wide = FALSE) {
  ch <- if (genome_wide) sample(chr_len$chromosome, nrow(loci), replace = TRUE, prob = chr_len$len) else loci$chromosome
  room <- pmax(chr_len$len[match(ch, chr_len$chromosome)] - (loci$end - loci$start), 1)
  s <- floor(runif(nrow(loci)) * room)
  list(ch = ch, s = s, e = s + (loci$end - loci$start))
}

# Observed number of overlapping loci and its null distribution for one QTL set
perm_test <- function(q, set_name, flank, genome_wide) {
  obs <- sum(hits_any(loci$chromosome, loci$start, loci$end, q, flank))
  null <- replicate(N_PERM_QTL, {
    r <- shuffle_loci(genome_wide)
    sum(hits_any(r$ch, r$s, r$e, q, flank))
  })
  tibble::tibble(qtl_set = set_name, n_qtl = nrow(q), flank_mb = flank / 1e6,
                 null_model = if (genome_wide) "genome_wide" else "same_chromosome", n_loci = nrow(loci),
                 observed = obs, expected = mean(null), p_value = (1 + sum(null >= obs)) / (1 + N_PERM_QTL))
}

sets <- c(list(all = qtl), split(qtl, qtl$trait_group))
grid <- tidyr::expand_grid(set = names(sets), flank = QTL_FLANK_BP, genome_wide = c(FALSE, TRUE))
set.seed(SEED_QTL)
(perm <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) perm_test(sets[[grid$set[i]]], grid$set[i], grid$flank[i], grid$genome_wide[i]))))
readr::write_csv(perm, file.path(dirs$res, "03_overlap_permutation_test.csv"))
# Reading: observed well above expected with a small P, in the same-chromosome
# model, means the loci sit in regions already linked to storability more
# often than chance. Observed close to expected means overlaps are what QTL
# of these sizes give by chance, so a single overlap (including chr04) is not
# independent support. Wide QTL from biparental maps cover several Mb, which
# raises the expected value; compare flank 0 with the extended version.

# ---------------------------------------------------------------------------
# 4) Locus table and figure
# ---------------------------------------------------------------------------
(locus_tab <- loci %>%
    dplyr::mutate(qtl_as_published = hits_any(chromosome, start, end, qtl, 0),
                  qtl_extended = hits_any(chromosome, start, end, qtl, max(QTL_FLANK_BP))) %>%
    dplyr::left_join(pairs %>%
                       dplyr::filter(flank_mb == max(flank_mb)) %>%
                       dplyr::group_by(locus_id) %>%
                       dplyr::summarise(qtl_ids = paste(unique(qtl_id), collapse = "; "),
                                        qtl_sources = paste(unique(source), collapse = "; "), .groups = "drop"),
                     by = "locus_id") %>%
    dplyr::arrange(dplyr::desc(qtl_as_published), dplyr::desc(qtl_extended), chromosome))
readr::write_csv(locus_tab, file.path(dirs$res, "04_locus_qtl_summary.csv"))

p_map <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = chr_len, ggplot2::aes(x = chromosome, xend = chromosome, y = 0, yend = len / 1e6),
                        linewidth = 3, colour = "grey85", lineend = "round") +
  ggplot2::geom_segment(data = qtl, ggplot2::aes(x = chromosome - 0.25, xend = chromosome - 0.25, y = start / 1e6,
                                                 yend = end / 1e6, colour = trait_group), linewidth = 2) +
  ggplot2::geom_point(data = loci, ggplot2::aes(x = chromosome + 0.25, y = (start + end) / 2e6,
                                                shape = locus_status == "convergent_locus"), size = 2.4) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 17), labels = c("other loci", "convergent locus"), name = NULL) +
  ggplot2::scale_x_continuous(breaks = 1:20) +
  ggplot2::scale_y_reverse() +
  ggplot2::labs(x = "Chromosome", y = "Position in Wm82.a4 (Mb)", colour = "QTL trait group") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), legend.position = "bottom")
ggplot2::ggsave(file.path(dirs$fig, "01_loci_vs_published_qtl.png"), p_map, width = 10, height = 6, dpi = 300)

close_stage("07_external_validation_qtl", dirs)
