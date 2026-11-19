##############################################################################
# Soybean seed longevity
# Stage 6d: anchoring the network modules in the maturation program
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: kME says how typical a gene is of its module, not whether the MODULE has
# anything to do with the acquisition of longevity. Each network module (built
# on the atlas) is taken to GSE98199 and tested on two questions:
#   1) do its genes change in a coordinated way along maturation (fresh seeds,
#      stages 7.1 to 9)? COMPETITIVE test: the module is compared with random
#      gene sets of the same size, because about 40% of the transcriptome
#      changes during maturation;
#   2) is it enriched in genes going up or down in the replicated contrasts
#      (stage 9 vs 7.2F, and drying 7.2D vs 7.2F)? Fisher test, BH per contrast.
# Each candidate then inherits the class of its module.
# Limitation: GSE98199 is also the temporal layer, so the anchor is NOT
# independent evidence. It says whether the candidate's coexpression
# neighbourhood follows the maturation program, giving biological meaning to
# the network layer, and it is not another vote in the ranking.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/06_evidence_integration/04_network_maturation_anchor.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2, DESeq2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install first: DESeq2")

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "04_maturation_anchor"))
f_de <- from_stage("02_maturation_timecourse", "07_deseq2_all_genes.csv")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_int <- from_stage(stage, "01_integrated_matrix/02_integrated_evidence_matrix_39.csv")
check_inputs(c(RAW$gse_counts, RAW$gse_meta, f_de, f_mod, f_int))

N_RAND <- 5000                                 # random sets per module
STAGE_ORD <- c("7.1", "7.2", "7.3", "8.1", "8.2", "8.3", "9")
# With hundreds of genes the competitive test is always significant, so the
# class also requires an effect size:
MIN_ABS_RHO <- 0.30                            # |mean rho of the module|
MIN_FRAC_DIR <- 0.70                           # fraction of genes going the same way

# ---------------------------------------------------------------------------
# 1) Counts, metadata and VST (same filter as stage 2: count >= 10 in >= 2 libraries)
# vst(blind = TRUE): the aim here is descriptive (trajectory), not a test.
# ---------------------------------------------------------------------------
meta <- readr::read_csv(RAW$gse_meta, show_col_types = FALSE) %>%
  dplyr::mutate(stage = as.character(stage), stage_index = match(stage, STAGE_ORD))
counts <- read.delim(gzfile(RAW$gse_counts), header = TRUE, row.names = 1, check.names = FALSE)
count_mat <- as.matrix(counts[rownames(counts) != "*", meta$description, drop = FALSE])
storage.mode(count_mat) <- "integer"
dds <- DESeq2::DESeqDataSetFromMatrix(count_mat[rowSums(count_mat >= 10) >= 2, ],
                                      as.data.frame(meta, row.names = meta$description), design = ~ 1)
vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
rownames(vsd) <- strip_version(rownames(vsd))
stopifnot(!anyDuplicated(rownames(vsd)))

# ---------------------------------------------------------------------------
# 2) Network modules measured in GSE98199
# The match rate is computed AFTER the expression filter, so it means "ID found
# and expressed in R7-R9": a low rate (e.g. midnightblue ~20%) indicates a
# module barely expressed at late maturation, not only ID mismatch.
# ---------------------------------------------------------------------------
modules <- readr::read_csv(f_mod, show_col_types = FALSE)
(id_match <- modules %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(n_network = dplyr::n(), n_in_GSE = sum(gene_id %in% rownames(vsd)),
                     match_rate = n_in_GSE / n_network) %>%
    dplyr::arrange(dplyr::desc(n_network)))
readr::write_csv(id_match, file.path(dirs$res, "01_ID_match_by_module.csv"))
bg <- dplyr::filter(modules, gene_id %in% rownames(vsd))     # background = network genes measured
mods_test <- setdiff(sort(unique(bg$module)), "grey")

# ---------------------------------------------------------------------------
# 3) Maturation trajectory: Spearman of each gene with the stage (fresh only)
# rho > 0 = rises along maturation; rho < 0 = falls. Sensitivity: same rho
# without stage 7.1 (single library), since several curves dropped mainly
# between 7.1 and 7.2.
# ---------------------------------------------------------------------------
fresh <- meta$description[meta$treatment == "fresh"]
stage_fresh <- meta$stage_index[match(fresh, meta$description)]
X <- vsd[bg$gene_id, fresh, drop = FALSE]
X <- X[apply(X, 1, sd) > 0, , drop = FALSE]
no71 <- stage_fresh != 1
X71 <- X[, no71, drop = FALSE]
ok71 <- apply(X71, 1, sd) > 0
gene_rho <- tibble::tibble(gene_id = rownames(X), rho_stage = as.vector(cor(t(X), stage_fresh, method = "spearman")),
                           rho_no71 = NA_real_) %>%
  dplyr::left_join(bg, by = "gene_id")
gene_rho$rho_no71[ok71] <- as.vector(cor(t(X71[ok71, , drop = FALSE]), stage_fresh[no71], method = "spearman"))

# Competitive test of one module: mean rho vs means of random same-size sets
# drawn from the background (all modules, grey included); two-sided P by the
# distance to the background mean
mu_bg <- mean(gene_rho$rho_stage)
module_test <- function(mod) {
  r <- gene_rho$rho_stage[gene_rho$module == mod]
  r71 <- gene_rho$rho_no71[gene_rho$module == mod]
  null <- replicate(N_RAND, mean(sample(gene_rho$rho_stage, length(r))))
  tibble::tibble(module = mod, n_genes = length(r), mean_rho = mean(r), frac_rho_pos = mean(r > 0),
                 frac_rho_neg = mean(r < 0), mean_rho_no71 = mean(r71, na.rm = TRUE),
                 p_competitive = (1 + sum(abs(null - mu_bg) >= abs(mean(r) - mu_bg))) / (1 + N_RAND))
}
set.seed(SEED_COMP)
module_trend <- dplyr::bind_rows(lapply(mods_test, module_test)) %>%
  dplyr::mutate(padj_competitive = p.adjust(p_competitive, "BH"),
                trend_class = dplyr::case_when(
                  padj_competitive < ALPHA & mean_rho >= MIN_ABS_RHO & frac_rho_pos >= MIN_FRAC_DIR ~ "rises_in_maturation",
                  padj_competitive < ALPHA & mean_rho <= -MIN_ABS_RHO & frac_rho_neg >= MIN_FRAC_DIR ~ "falls_in_maturation",
                  padj_competitive < ALPHA & mean_rho > mu_bg ~ "weak_rise",
                  padj_competitive < ALPHA & mean_rho < mu_bg ~ "weak_fall",
                  TRUE ~ "no_coordinated_trend"),
                robust_without_7.1 = trend_class %in% c("rises_in_maturation", "falls_in_maturation") &
                  sign(mean_rho_no71) == sign(mean_rho) & abs(mean_rho_no71) >= MIN_ABS_RHO) %>%
  dplyr::arrange(dplyr::desc(abs(mean_rho)))
readr::write_csv(module_trend, file.path(dirs$res, "02_module_maturation_trend.csv"))

# ---------------------------------------------------------------------------
# 4) Module enrichment in DE genes of the replicated contrasts (Fisher)
# ---------------------------------------------------------------------------
de <- readr::read_csv(f_de, show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = strip_version(GSE98199_gene_id), up = significant & log2FoldChange > 0,
                down = significant & log2FoldChange < 0) %>%
  dplyr::inner_join(bg, by = "gene_id")

# Fisher test for one module, contrast and direction ("up" or "down")
fisher_module <- function(mod, ctr, dir) {
  d <- de[de$contrast == ctr, ]
  in_mod <- d$module == mod
  hit <- d[[dir]]
  ft <- fisher.test(table(factor(in_mod, c(TRUE, FALSE)), factor(hit, c(TRUE, FALSE))))
  tibble::tibble(module = mod, contrast = ctr, direction = dir, n_module = sum(in_mod), n_hit_module = sum(in_mod & hit),
                 frac_hit_module = mean(hit[in_mod]), frac_hit_background = mean(hit),
                 odds_ratio = unname(ft$estimate), p = ft$p.value)
}
fisher_module("purple", "stage9_vs_7.2_fresh", "up")   # check with one module
grid <- tidyr::expand_grid(module = mods_test, contrast = unique(de$contrast), direction = c("up", "down"))
module_DE <- dplyr::bind_rows(Map(fisher_module, grid$module, grid$contrast, grid$direction)) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(padj = p.adjust(p, "BH")) %>%
  dplyr::ungroup()
readr::write_csv(module_DE, file.path(dirs$res, "03_module_DE_enrichment.csv"))

mat_DE <- module_DE %>%
  dplyr::filter(contrast == "stage9_vs_7.2_fresh") %>%
  dplyr::select(module, direction, odds_ratio, padj) %>%
  tidyr::pivot_wider(names_from = direction, values_from = c(odds_ratio, padj)) %>%
  dplyr::mutate(DE_class = dplyr::case_when(
    padj_up < ALPHA & odds_ratio_up > 1 & !(padj_down < ALPHA & odds_ratio_down > 1) ~ "enriched_up_in_maturation",
    padj_down < ALPHA & odds_ratio_down > 1 & !(padj_up < ALPHA & odds_ratio_up > 1) ~ "enriched_down_in_maturation",
    padj_up < ALPHA & padj_down < ALPHA ~ "enriched_both_directions",
    TRUE ~ "no_DE_enrichment"))

# ---------------------------------------------------------------------------
# 5) Anchor class of each module and of each candidate
# A module that rises (or falls) coherently AND is enriched in DE genes in the
# same direction belongs to the maturation program. Pereira Lima et al. (2017)
# associated with longevity both genes that rise (chaperones, HSPs, RFO) and
# genes that fall (chloroplast/photosynthesis), so both directions matter.
# A module without trend means high kME is only typicality in a module that
# does not follow maturation.
# ---------------------------------------------------------------------------
module_class <- module_trend %>%
  dplyr::select(module, n_genes, mean_rho, frac_rho_pos, mean_rho_no71, padj_competitive, trend_class, robust_without_7.1) %>%
  dplyr::left_join(mat_DE, by = "module") %>%
  dplyr::mutate(anchor_class = dplyr::case_when(
    trend_class == "rises_in_maturation" & DE_class == "enriched_up_in_maturation" ~ "maturation_up_module",
    trend_class == "falls_in_maturation" & DE_class == "enriched_down_in_maturation" ~ "maturation_down_module",
    trend_class %in% c("weak_rise", "weak_fall") ~ "weak_trend_module",
    trend_class != "no_coordinated_trend" | DE_class != "no_DE_enrichment" ~ "partial_maturation_signal",
    TRUE ~ "no_maturation_signal"))
readr::write_csv(module_class, file.path(dirs$res, "04_module_maturation_anchor_classes.csv"))

cand_anchor <- readr::read_csv(f_int, show_col_types = FALSE) %>%
  dplyr::select(gene_id, module_for_integration, network_support_label, module_coupling,
                network_rank_available_all, temporal_support_label, evidence_coverage_class) %>%
  dplyr::left_join(module_class, by = c("module_for_integration" = "module")) %>%
  dplyr::mutate(anchor_class = dplyr::case_when(is.na(module_for_integration) ~ "not_in_network",
                                                module_for_integration == "grey" ~ "grey_unassigned",
                                                TRUE ~ anchor_class)) %>%
  dplyr::arrange(anchor_class, network_rank_available_all)
readr::write_csv(cand_anchor, file.path(dirs$res, "05_candidate_module_anchor_39.csv"))
print(dplyr::count(cand_anchor, anchor_class))
# expected: 11 candidates in falling modules, 2 in rising, 19 in weak-trend, 7 outside the network

# ---------------------------------------------------------------------------
# 6) Figure: mean z-score trajectory of the candidate modules
# ---------------------------------------------------------------------------
Z <- t(scale(t(vsd[bg$gene_id, meta$description, drop = FALSE])))
Z <- Z[complete.cases(Z), , drop = FALSE]
mods_cand <- setdiff(na.omit(unique(cand_anchor$module_for_integration)), "grey")
prof <- dplyr::bind_rows(lapply(mods_cand, function(m) {
  g <- intersect(bg$gene_id[bg$module == m], rownames(Z))
  tibble::tibble(module = m, description = colnames(Z), mean_z = colMeans(Z[g, , drop = FALSE]))
})) %>%
  dplyr::left_join(meta, by = "description") %>%
  dplyr::left_join(dplyr::select(module_class, module, anchor_class), by = "module")
p_prof <- ggplot2::ggplot(dplyr::filter(prof, treatment == "fresh"), ggplot2::aes(stage_index, mean_z)) +
  ggplot2::stat_summary(fun = mean, geom = "line", colour = "grey40") +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::geom_point(data = dplyr::filter(prof, treatment == "dried"), shape = 2, size = 2) +
  ggplot2::facet_wrap(~ module + anchor_class, ncol = 3) +
  ggplot2::scale_x_continuous(breaks = seq_along(STAGE_ORD), labels = STAGE_ORD) +
  ggplot2::labs(x = "Developmental stage (GSE98199)", y = "Mean gene z-score (VST)",
                subtitle = "Circles: fresh seeds; triangles: 7.2 dried") +
  ggplot2::theme_bw() +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 7))
ggplot2::ggsave(file.path(dirs$fig, "01_candidate_module_trajectories.png"), p_prof, width = 10, height = 8, dpi = 300)
saveRDS(list(module_trend = module_trend, module_class = module_class, profiles = prof),
        file.path("data", "processed", "maturation_anchor_objects.rds"))

close_stage("06_evidence_integration_anchor", dirs)
