##############################################################################
# Soybean seed longevity
# Stage 6h: do the Wm82.a4 genes that really sit at the displaced SNPs have more
# expression support than the genes listed by the authors?
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: for the displaced Thenveettil et al. (2025) SNPs (stage 1), compare the
# published genes with the a4 genes within +-10 kb of the SNP (e.g. the
# lipoxygenase cluster on chr13), using the layers that cover the whole genome:
#   Temporal: GSE98199 DESeq2. To compare genes outside the list, T is rebuilt on
#             the GENOME scale: T = sqrt(rank(padj) x rank(|log2FC|)) among the
#             30,710 tested genes, contrast 9 vs 7.2F; top_pct_T = position in %
#             (1 = best 1% of the transcriptome).
#   Network:  module and its maturation class (script 04).
#   Pereira Lima et al. (2017) S5: direction in the longevity contrast 9 vs 7.2D,
#             if the file is in data/raw/external/pereira_lima_2017/.
#   Atlas:    NOT covered (stage 3 extracted only the 39 genes); 03_ lists the
#             genes to extract.
# This is an EXTENSION: it does not change the main analysis of the 39 genes.
#
# Run from the repository root, after script 06 of this stage:
#   source("scripts/06_evidence_integration/08_a4_genes_at_displaced_snps.R")
#
# Packages: readr, readxl, dplyr, tidyr, stringr, ggplot2, ggrepel
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "readxl", "dplyr", "tidyr", "stringr", "ggplot2"))
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install first: ggrepel")

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "08_displaced_snp_genes"))
f_win <- from_stage("01_candidates", "05_a4_genes_in_snp_window.csv")
f_verd <- from_stage("01_candidates", "06_snp_verdict.csv")
f_pos <- from_stage("01_candidates", "03_candidate_positions_a4.csv")
f_de <- from_stage("02_maturation_timecourse", "07_deseq2_all_genes.csv")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_anchor <- from_stage(stage, "04_maturation_anchor/04_module_maturation_anchor_classes.csv")
f_cls <- from_stage(stage, "06_final_classes/01_final_evidence_table_39_genes.csv")
f_s5 <- file.path(RAW$pl_dir, "S5Table.xlsx")                       # optional
check_inputs(c(f_win, f_verd, f_pos, f_de, f_mod, f_anchor, f_cls, RAW$info_annot))

# ---------------------------------------------------------------------------
# 1) The two gene sets of each displaced SNP (a gene in both stays as "authors")
# ---------------------------------------------------------------------------
displaced <- readr::read_csv(f_verd, show_col_types = FALSE) %>%
  dplyr::filter(verdict == "displaced") %>%
  dplyr::select(locus_id, snp_pos)
genes <- dplyr::bind_rows(
  readr::read_csv(f_pos, show_col_types = FALSE) %>%
    dplyr::semi_join(displaced, by = c("locus_id", "snp_pos")) %>%
    dplyr::transmute(locus_id, snp_pos, gene_id, set = "authors_list", dist_snp_kb),
  readr::read_csv(f_win, show_col_types = FALSE) %>%
    dplyr::semi_join(displaced, by = c("locus_id", "snp_pos")) %>%
    dplyr::transmute(locus_id, snp_pos, gene_id, set = "a4_at_SNP", dist_snp_kb = dist_kb)) %>%
  dplyr::distinct(locus_id, snp_pos, gene_id, .keep_all = TRUE)
print(dplyr::count(genes, locus_id, set))

# ---------------------------------------------------------------------------
# 2) Temporal evidence on the genome scale
# ---------------------------------------------------------------------------
de_dir <- function(p, lfc) dplyr::case_when(p < ALPHA & lfc >= MIN_LFC ~ "up", p < ALPHA & lfc <= -MIN_LFC ~ "down",
                                            TRUE ~ "not_DE")
de <- readr::read_csv(f_de, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(gene_id = strip_version(GSE98199_gene_id),
                ctr = dplyr::if_else(contrast == "stage9_vs_7.2_fresh", "mat", "dry")) %>%
  dplyr::select(gene_id, ctr, baseMean, log2FoldChange, padj) %>%
  tidyr::pivot_wider(names_from = ctr, values_from = c(log2FoldChange, padj), names_sep = "_") %>%
  dplyr::mutate(T_genome = sqrt(rank_frac(padj_mat) * rank_frac(-abs(log2FoldChange_mat))),
                top_pct_T = round(100 * rank_frac(T_genome), 1),
                DE_mat = de_dir(padj_mat, log2FoldChange_mat), DE_dry = de_dir(padj_dry, log2FoldChange_dry))
nrow(de)                                       # expected: 30,710 tested genes

# ---------------------------------------------------------------------------
# 3) Network, Pereira Lima S5 and annotation
# Direct network members come from the module table; for the listed genes the
# module used in the integration (with projection) is used when absent.
# ---------------------------------------------------------------------------
annot <- read_info_annot() %>%
  dplyr::select(gene_id, arabi_hit = `Best-hit-arabi-name`, arabi_defline = `arabi-defline`)
s5 <- if (file.exists(f_s5)) {
  read_s5 <- function(sheet, dir) {
    x <- readxl::read_excel(f_s5, sheet = sheet, skip = 6, .name_repair = "unique_quiet")
    dplyr::filter(tibble::tibble(gene_id = x[[1]], log2_9_vs_72D = as.numeric(x[[10]]), direction_S5 = dir),
                  stringr::str_detect(gene_id, "^Glyma\\."))
  }
  dplyr::bind_rows(read_s5("DEG UP", "up_with_longevity"), read_s5("DEG DOWN", "down_with_longevity"))
} else tibble::tibble(gene_id = character(), log2_9_vs_72D = numeric(), direction_S5 = character())

tab <- genes %>%
  dplyr::left_join(annot, by = "gene_id") %>%
  dplyr::left_join(dplyr::select(de, gene_id, baseMean, log2FoldChange_mat, padj_mat, DE_mat, log2FoldChange_dry,
                                 padj_dry, DE_dry, top_pct_T), by = "gene_id") %>%
  dplyr::left_join(readr::read_csv(f_mod, show_col_types = FALSE), by = "gene_id") %>%
  dplyr::left_join(dplyr::select(readr::read_csv(f_cls, show_col_types = FALSE), gene_id, module_main = module,
                                 class_main = evidence_class), by = "gene_id") %>%
  dplyr::mutate(module = dplyr::coalesce(module, module_main)) %>%
  dplyr::left_join(dplyr::select(readr::read_csv(f_anchor, show_col_types = FALSE), module, trend_class, anchor_class),
                   by = "module") %>%
  dplyr::left_join(s5, by = "gene_id") %>%
  dplyr::mutate(temporal_status = ifelse(is.na(top_pct_T), "not_tested_low_count", "tested"),
                module = dplyr::coalesce(module, "outside_network"),
                direction_S5 = dplyr::coalesce(direction_S5, "not_DE_or_absent")) %>%
  dplyr::arrange(locus_id, snp_pos, set, top_pct_T)
readr::write_csv(tab, file.path(dirs$res, "01_authors_vs_a4_genes.csv"))

# ---------------------------------------------------------------------------
# 4) Comparison per SNP
# A smaller best_top_pct_T = gene more prominent in the maturation
# transcriptome. If the a4 genes beat the listed ones, the locus may have been
# underestimated in the main analysis; if similar, swapping genes would not
# change the reading of the locus.
# ---------------------------------------------------------------------------
(by_snp <- tab %>%
    dplyr::group_by(locus_id, snp_pos, set) %>%
    dplyr::summarise(n_genes = dplyr::n(), n_tested = sum(temporal_status == "tested"),
                     best_top_pct_T = if (any(!is.na(top_pct_T))) min(top_pct_T, na.rm = TRUE) else NA_real_,
                     n_DE_mat = sum(DE_mat %in% c("up", "down")), n_S5 = sum(direction_S5 != "not_DE_or_absent"),
                     modules = paste(unique(module), collapse = "; "), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = set, values_from = c(n_genes, n_tested, best_top_pct_T, n_DE_mat, n_S5, modules)))
readr::write_csv(by_snp, file.path(dirs$res, "02_comparison_by_snp.csv"))
readr::write_csv(dplyr::distinct(dplyr::filter(tab, set == "a4_at_SNP"), gene_id, locus_id, arabi_defline),
                 file.path(dirs$res, "03_a4_genes_to_extract_from_atlas.csv"))

# ---------------------------------------------------------------------------
# 5) Figure: position in the maturation transcriptome, authors vs a4, per SNP
# ---------------------------------------------------------------------------
fig <- tab %>%
  dplyr::filter(!is.na(top_pct_T)) %>%
  dplyr::mutate(snp = paste0(sub("Th_", "", locus_id), " | ", formatC(snp_pos, big.mark = ",", format = "d")),
                set = factor(set, c("authors_list", "a4_at_SNP"), c("Authors' list", "a4 gene at SNP (+-10 kb)")),
                label = sub("Glyma\\.", "", gene_id))
lv <- sort(unique(fig$snp))
fig$y <- match(fig$snp, lv) + ifelse(fig$set == "Authors' list", -0.18, 0.18)
p <- ggplot2::ggplot(fig, ggplot2::aes(top_pct_T, y, colour = set, shape = DE_mat)) +
  ggplot2::geom_vline(xintercept = c(10, 50), linetype = "dashed", colour = "grey70") +
  ggplot2::geom_point(size = 3.2, stroke = 1) +
  ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3.2, show.legend = FALSE, min.segment.length = 0,
                           box.padding = 0.35, max.overlaps = Inf, seed = 1) +
  ggplot2::scale_y_continuous(breaks = seq_along(lv), labels = lv) +
  ggplot2::scale_x_continuous(limits = c(0, 100)) +
  ggplot2::scale_colour_manual(values = c("Authors' list" = "grey45", "a4 gene at SNP (+-10 kb)" = "#1b7837")) +
  ggplot2::scale_shape_manual(values = c(down = 25, up = 24, not_DE = 21)) +
  ggplot2::labs(x = "Position in the maturation transcriptome (% best; GSE98199, 9 vs 7.2F)", y = NULL,
                colour = NULL, shape = "Maturation DE") +
  ggplot2::theme_bw(base_size = 13)
ggplot2::ggsave(file.path(dirs$fig, "01_authors_vs_a4_genes_temporal.png"), p, width = 11, height = 6, dpi = 300)

close_stage("06_evidence_integration_displaced_snps", dirs)
