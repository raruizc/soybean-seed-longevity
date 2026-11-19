##############################################################################
# Soybean seed longevity
# Stage 7a: candidates and modules against the supplementary lists of
# Pereira Lima et al. (2017)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: use the lists published by the authors of GSE98199 as a second reading
# of the candidates and of the network modules.
#   S3 Table: 16,255 genes DE along maturation. The list is biased towards
#             genes that FALL (87.7% fall; header: "only genes with a variation
#             > 1"), so in_S3 is descriptive only, never a validation criterion.
#   S4 Table: 549 transcription factors of the maturation network, with their
#             correlation with longevity (corr P50).
#   S5 Table: DE genes in the LONGEVITY contrast stage 9 vs dried 7.2 (742 up,
#             1,525 down). It compares dried seeds with high vs low longevity,
#             so it separates longevity from drying; stage 2 did not test it.
# Limitation: the three tables come from the SAME RNA-seq used in the temporal
# layer, so this is a complementary reading, not independent validation; the
# independent check is script 02 (Medicago, another species and study).
# IDs: the tables use Wm82.a2; candidates use the stage 2 bridge; for the
# network, a2 IDs are matched directly to a4 (93.9% coincide).
#
# Run from the repository root, after stage 6:
#   source("scripts/07_external_validation/01_pereira_lima_longevity_lists.R")
#
# Packages: readr, readxl, dplyr, tidyr, stringr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "readxl", "dplyr", "tidyr", "stringr", "ggplot2"))

stage <- "07_external_validation"
dirs <- setup_stage(file.path(stage, "01_pereira_lima"))
f_s3 <- file.path(RAW$pl_dir, "S3Table.xlsx")
f_s4 <- file.path(RAW$pl_dir, "S4Table.xlsx")
f_s5 <- file.path(RAW$pl_dir, "S5Table.xlsx")
f_cand <- from_stage("02_maturation_timecourse", "09_candidate_temporal_evidence.csv")
f_de <- from_stage("02_maturation_timecourse", "07_deseq2_all_genes.csv")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_cls <- from_stage("06_evidence_integration", "06_final_classes/01_final_evidence_table_39_genes.csv")
check_inputs(c(f_s3, f_s4, f_s5, f_cand, f_de, f_mod, f_cls))
P50_CUT <- 0.90                                # |corr P50| > 0.9, the authors' criterion (bold in S4)

# ---------------------------------------------------------------------------
# 1) Supplementary tables (headers on different rows: S3 row 5, S4 row 3, S5 row 7)
# ---------------------------------------------------------------------------
s3 <- readxl::read_excel(f_s3, skip = 4) %>%
  dplyr::transmute(gene_a2 = `id Gmax v2`, name_s3 = `Name  (Soybase)`) %>%
  dplyr::filter(stringr::str_detect(gene_a2, "^Glyma\\."))
nrow(s3)                                       # expected: ~16,255
s4 <- readxl::read_excel(f_s4, skip = 2) %>%
  dplyr::transmute(gene_a2 = `id Gmax v2`, tf_name = `Name  (Soybase)`, arabi = `Arabidopsis`,
                   corr_P50 = as.numeric(`corr P50`), max_stage = `max dev stage`) %>%
  dplyr::filter(stringr::str_detect(gene_a2, "^Glyma\\.")) %>%
  dplyr::mutate(P50_class = dplyr::case_when(corr_P50 > P50_CUT ~ "positive_P50", corr_P50 < -P50_CUT ~ "negative_P50",
                                             TRUE ~ "weak"))
table(s4$P50_class)                            # expected: 26 positive

# S5: two sheets; column 1 = a2 ID, column 10 = log2(9/7.2D), column 11 = padj
read_s5 <- function(sheet, dir) {
  x <- readxl::read_excel(f_s5, sheet = sheet, skip = 6, .name_repair = "unique_quiet")
  tibble::tibble(gene_a2 = x[[1]], log2_9_vs_72D = as.numeric(x[[10]]), padj_9_vs_72D = as.numeric(x[[11]]),
                 direction_S5 = dir) %>%
    dplyr::filter(stringr::str_detect(gene_a2, "^Glyma\\."))
}
s5 <- dplyr::bind_rows(read_s5("DEG UP", "up_in_high_longevity"), read_s5("DEG DOWN", "down_in_high_longevity"))
dplyr::count(s5, direction_S5)                 # expected: 742 / 1,525
stopifnot(!anyDuplicated(s5$gene_a2))

# ---------------------------------------------------------------------------
# 2) What the authors' lists say about each candidate
# ---------------------------------------------------------------------------
cand_pl <- readr::read_csv(f_cand, show_col_types = FALSE) %>%
  dplyr::transmute(gene_id = Wm82_a4_id, gene_a2 = strip_version(GSE98199_gene_id)) %>%
  dplyr::left_join(dplyr::select(readr::read_csv(f_cls, show_col_types = FALSE), gene_id, evidence_class, locus_id,
                                 module, maturation_program, maturation_direction), by = "gene_id") %>%
  dplyr::mutate(in_S3_maturation_DE = gene_a2 %in% s3$gene_a2) %>%
  dplyr::left_join(s5, by = "gene_a2") %>%
  dplyr::left_join(dplyr::select(s4, gene_a2, tf_name, corr_P50, P50_class), by = "gene_a2") %>%
  dplyr::mutate(direction_S5 = dplyr::coalesce(direction_S5, "not_DE_9_vs_72D"), is_network_TF = !is.na(tf_name)) %>%
  dplyr::arrange(evidence_class, gene_id)
readr::write_csv(cand_pl, file.path(dirs$res, "01_candidates_PereiraLima.csv"))

# ---------------------------------------------------------------------------
# 3) Are the candidates enriched in longevity genes (S5)?
# Background = genes tested in GSE98199 (non-NA padj in our DESeq2). A
# candidate fraction close to the background means the candidates do not stand
# out (as in stage 6); a much larger one would be a longevity signal (9 vs 7.2D)
# that our contrast did not capture.
# ---------------------------------------------------------------------------
tested <- readr::read_csv(f_de, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::distinct(gene_a2 = strip_version(GSE98199_gene_id))
bg_rate <- mean(tested$gene_a2 %in% s5$gene_a2)
cand_tested <- dplyr::filter(cand_pl, gene_a2 %in% tested$gene_a2)
k <- sum(cand_tested$direction_S5 != "not_DE_9_vs_72D")
(cand_enrich <- tibble::tibble(n = nrow(cand_tested), n_in_S5 = k, frac_cand = k / nrow(cand_tested),
                               frac_background = bg_rate, p_binom = binom.test(k, nrow(cand_tested), bg_rate)$p.value))
readr::write_csv(cand_enrich, file.path(dirs$res, "02_candidate_enrichment_S5.csv"))
# expected: 3/31 (9.7% vs 7.4%, P = 0.50)

# ---------------------------------------------------------------------------
# 4) Network modules x longevity genes (S5) and P50-linked TFs (S4)
# Fisher per module and direction; background = network genes tested in GSE98199.
# A module enriched in "up_in_high_longevity" is active once the seed already
# has high longevity (beyond drying).
# ---------------------------------------------------------------------------
net <- readr::read_csv(f_mod, show_col_types = FALSE) %>%
  dplyr::filter(module != "grey", gene_id %in% tested$gene_a2) %>%
  dplyr::left_join(dplyr::select(s5, gene_id = gene_a2, direction_S5), by = "gene_id") %>%
  dplyr::mutate(up = direction_S5 %in% "up_in_high_longevity", down = direction_S5 %in% "down_in_high_longevity")

# Fisher test of one module for one direction ("up" or "down")
fisher_s5 <- function(mod, dir) {
  in_mod <- net$module == mod
  hit <- net[[dir]]
  ft <- fisher.test(table(factor(in_mod, c(TRUE, FALSE)), factor(hit, c(TRUE, FALSE))))
  tibble::tibble(module = mod, direction = dir, n_module = sum(in_mod), n_hit = sum(in_mod & hit),
                 frac_module = mean(hit[in_mod]), frac_background = mean(hit), odds_ratio = unname(ft$estimate),
                 p = ft$p.value)
}
fisher_s5("black", "down")                     # check with one module
grid <- tidyr::expand_grid(module = sort(unique(net$module)), direction = c("up", "down"))
mod_s5 <- dplyr::bind_rows(Map(fisher_s5, grid$module, grid$direction)) %>%
  dplyr::mutate(padj = p.adjust(p, "BH")) %>%
  dplyr::arrange(padj)
readr::write_csv(mod_s5, file.path(dirs$res, "03_module_enrichment_S5.csv"))
# expected: purple and pink enriched in S5-down

tf_mod <- s4 %>%
  dplyr::inner_join(readr::read_csv(f_mod, show_col_types = FALSE), by = c("gene_a2" = "gene_id")) %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(n_TF = dplyr::n(), n_TF_P50_pos = sum(P50_class == "positive_P50"),
                   n_TF_P50_neg = sum(P50_class == "negative_P50"),
                   TF_P50_pos = paste(tf_name[P50_class == "positive_P50"], collapse = "; ")) %>%
  dplyr::arrange(dplyr::desc(n_TF_P50_pos))
readr::write_csv(tf_mod, file.path(dirs$res, "04_P50_TFs_by_module.csv"))

mods_cand <- setdiff(na.omit(unique(cand_pl$module)), "grey")
p <- mod_s5 %>%
  dplyr::filter(module %in% mods_cand) %>%
  dplyr::mutate(log2_OR = log2(pmax(odds_ratio, 0.05)), sig = ifelse(padj < ALPHA, "BH < 0.05", "n.s."),
                direction = factor(direction, c("up", "down"), c("Up with longevity (9 > 7.2D)", "Down with longevity (9 < 7.2D)"))) %>%
  ggplot2::ggplot(ggplot2::aes(log2_OR, reorder(module, log2_OR), colour = sig)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::geom_point(ggplot2::aes(size = frac_module)) +
  ggplot2::facet_wrap(~ direction) +
  ggplot2::scale_colour_manual(values = c("BH < 0.05" = "#1b7837", "n.s." = "grey60")) +
  ggplot2::labs(x = "log2 odds ratio (module vs network)", y = NULL, colour = NULL, size = "Fraction in module") +
  ggplot2::theme_bw(base_size = 14)
ggplot2::ggsave(file.path(dirs$fig, "01_candidate_modules_vs_S5.png"), p, width = 11, height = 6, dpi = 300)

close_stage("07_external_validation_pereira_lima", dirs)
