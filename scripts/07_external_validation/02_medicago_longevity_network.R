##############################################################################
# Soybean seed longevity
# Stage 7b: external validation with the Medicago truncatula maturation network
# (MatNet) and longevity module of Righetti et al. (2015, Plant Cell)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Question: do the candidates and the soybean network modules have
# counterparts in the Medicago longevity module? This is the first INDEPENDENT
# evidence of the study: another species, experiment and platform.
# Bridge between species: best Arabidopsis hit (AGI, without ".1"), from the
# Wm82.a4 annotation for soybean and BestHit.TAIR10 in Righetti's tables. Two
# genes with the same AGI are treated as counterparts; this is an orthology
# proxy that does not separate paralogs.
# Righetti Supplemental Data Sets: DS1 MatNet (2,912 genes, clusters 1-5) =
# testable universe; DS2 MatNet TFs; DS3 longevity module (499 genes); DS4
# desiccation tolerance genes (200); DS5 ABI3 regulation (Mt-abi3 vs WT).
# Analyses: (1) gene by gene (39 + a4 genes at displaced SNPs); (2) candidate
# enrichment in the longevity module; (3) soybean modules enriched in the
# Medicago longevity and desiccation sets (Fisher, BH), with universe = soybean
# network genes whose AGI is in MatNet, by gene and by unique AGI.
#
# Run from the repository root, after stage 6:
#   source("scripts/07_external_validation/02_medicago_longevity_network.R")
#
# Packages: readr, readxl, dplyr, tidyr, stringr, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "readxl", "dplyr", "tidyr", "stringr", "ggplot2"))

stage <- "07_external_validation"
dirs <- setup_stage(file.path(stage, "02_medicago"))
ds <- function(x) file.path(RAW$righetti_dir, paste0("Righetti2015_", x, ".xls"))
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_cls <- from_stage("06_evidence_integration", "06_final_classes/01_final_evidence_table_39_genes.csv")
f_anchor <- from_stage("06_evidence_integration", "04_maturation_anchor/04_module_maturation_anchor_classes.csv")
f_6h <- from_stage("06_evidence_integration", "08_displaced_snp_genes/01_authors_vs_a4_genes.csv")   # optional
check_inputs(c(ds("DS1_MatNet"), ds("DS2_TF_MatNet"), ds("DS3_longevity_module"), ds("DS4_desiccation_tolerance"),
               ds("DS5_ABI3_regulation"), RAW$info_annot, f_mod, f_cls, f_anchor))
MIN_GENES_MOD <- 30                            # smaller modules are left out of Fisher (little power)

# AT1G55020.1 -> AT1G55020 (upper case, no version)
agi <- function(x) stringr::str_remove(toupper(stringr::str_trim(x)), "\\.\\d+$")

# ---------------------------------------------------------------------------
# 1) Righetti tables (header on a different row in each file), summarized per AGI
# ---------------------------------------------------------------------------
read_righetti <- function(f, skip) {
  readxl::read_excel(f, skip = skip, .name_repair = "unique_quiet") %>%
    dplyr::mutate(AGI = agi(BestHit.TAIR10)) %>%
    dplyr::filter(stringr::str_detect(AGI, "^AT[1-5CM]G\\d{5}$"))
}
matnet <- read_righetti(ds("DS1_MatNet"), 5) %>% dplyr::select(AGI, matnet_cluster = `Cluster in MATNET`)
tf_mt <- read_righetti(ds("DS2_TF_MatNet"), 4) %>% dplyr::distinct(AGI)
longmod <- read_righetti(ds("DS3_longevity_module"), 6) %>% dplyr::select(AGI, longmod_degree = Degree)
dtmod <- read_righetti(ds("DS4_desiccation_tolerance"), 5) %>% dplyr::distinct(AGI)
# DS5: the P of each contrast is the column right after its log2FC (repeated "pvalue" names)
abi3_raw <- read_righetti(ds("DS5_ABI3_regulation"), 5)
col_40 <- which(names(abi3_raw) == "Mt-abi3/WT 40 DAP")
abi3 <- tibble::tibble(AGI = abi3_raw$AGI, lfc = as.numeric(abi3_raw[[col_40]]), p = as.numeric(abi3_raw[[col_40 + 1]])) %>%
  dplyr::group_by(AGI) %>%
  dplyr::summarise(abi3_lfc_40DAP = lfc[which.min(p)], abi3_p_40DAP = min(p, na.rm = TRUE))

mt <- matnet %>%
  dplyr::group_by(AGI) %>%
  dplyr::summarise(matnet_cluster = paste(sort(unique(matnet_cluster)), collapse = "; ")) %>%
  dplyr::full_join(longmod %>% dplyr::group_by(AGI) %>%
                     dplyr::summarise(in_longmod = TRUE, longmod_max_degree = suppressWarnings(max(as.numeric(longmod_degree), na.rm = TRUE))),
                   by = "AGI") %>%
  dplyr::full_join(dplyr::mutate(dtmod, in_DT = TRUE), by = "AGI") %>%
  dplyr::full_join(dplyr::mutate(tf_mt, is_TF_MatNet = TRUE), by = "AGI") %>%
  dplyr::full_join(abi3, by = "AGI") %>%
  dplyr::mutate(longmod_max_degree = ifelse(is.finite(longmod_max_degree), longmod_max_degree, NA),
                dplyr::across(c(in_longmod, in_DT, is_TF_MatNet), ~ dplyr::coalesce(.x, FALSE)),
                in_matnet = !is.na(matnet_cluster) | in_longmod)
c(AGI_MatNet = sum(mt$in_matnet), AGI_longmod = sum(mt$in_longmod), AGI_DT = sum(mt$in_DT))

# ---------------------------------------------------------------------------
# 2) Soybean: AGI of each gene, module and maturation program
# ---------------------------------------------------------------------------
annot <- read_info_annot() %>%
  dplyr::transmute(gene_id, AGI = agi(`Best-hit-arabi-name`), arabi_defline = `arabi-defline`)
anchor <- dplyr::select(readr::read_csv(f_anchor, show_col_types = FALSE), module, trend_class, anchor_class)
network <- readr::read_csv(f_mod, show_col_types = FALSE) %>%
  dplyr::left_join(annot, by = "gene_id") %>%
  dplyr::left_join(mt, by = "AGI") %>%
  dplyr::mutate(dplyr::across(c(in_matnet, in_longmod, in_DT), ~ dplyr::coalesce(.x, FALSE)))

# ---------------------------------------------------------------------------
# 3) Candidate by candidate
# ---------------------------------------------------------------------------
cand <- readr::read_csv(f_cls, show_col_types = FALSE) %>%
  dplyr::transmute(gene_id, evidence_class, locus_id, module, set = "39_published")
if (file.exists(f_6h)) {
  extra <- readr::read_csv(f_6h, show_col_types = FALSE) %>%
    dplyr::filter(set == "a4_at_SNP") %>%
    dplyr::transmute(gene_id, evidence_class = "a4_gene_at_displaced_SNP", locus_id, module, set = "a4_at_SNP")
  cand <- dplyr::distinct(dplyr::bind_rows(cand, extra), gene_id, .keep_all = TRUE)
}
cand_mt <- cand %>%
  dplyr::left_join(annot, by = "gene_id") %>%
  dplyr::left_join(mt, by = "AGI") %>%
  dplyr::mutate(dplyr::across(c(in_matnet, in_longmod, in_DT, is_TF_MatNet), ~ dplyr::coalesce(.x, FALSE)),
                medicago_evidence = dplyr::case_when(in_longmod ~ "longevity_module", in_DT ~ "desiccation_tolerance",
                                                     in_matnet ~ "MatNet_other_cluster", !is.na(AGI) ~ "AGI_outside_MatNet",
                                                     TRUE ~ "no_Arabidopsis_hit")) %>%
  dplyr::arrange(set, evidence_class, gene_id)
readr::write_csv(cand_mt, file.path(dirs$res, "01_candidates_Medicago.csv"))
print(dplyr::count(cand_mt, set, medicago_evidence))   # expected: only 3/39 map into MatNet

# ---------------------------------------------------------------------------
# 4) Are the 39 candidates enriched in the longevity module?
# With few testable candidates the test has little power; it is descriptive.
# ---------------------------------------------------------------------------
univ <- dplyr::filter(network, in_matnet)
c39 <- dplyr::filter(cand_mt, set == "39_published", in_matnet)
(cand_enrich <- tibble::tibble(n = nrow(c39), n_longmod = sum(c39$in_longmod), frac_background = mean(univ$in_longmod),
                               p_binom = if (nrow(c39)) binom.test(sum(c39$in_longmod), nrow(c39), mean(univ$in_longmod))$p.value else NA))
readr::write_csv(cand_enrich, file.path(dirs$res, "02_candidate_enrichment_longevity_module.csv"))

# ---------------------------------------------------------------------------
# 5) Soybean modules x Medicago longevity module and desiccation set
# Enriched in longevity (not desiccation) = conserved program tied to longevity,
# independent support. Enriched only in desiccation = drying program. Grey and
# royalblue (organellar artifact, stage 4) are left out.
# By unique AGI as well: soybean is paleopolyploid, so two paralogs with the
# same AGI would count twice in the per-gene Fisher test.
# ---------------------------------------------------------------------------

# Fisher test, per gene, of one module against one Medicago set
fisher_gene <- function(mod, target) {
  in_mod <- univ$module == mod
  hit <- univ[[target]]
  ft <- fisher.test(table(factor(in_mod, c(TRUE, FALSE)), factor(hit, c(TRUE, FALSE))))
  tibble::tibble(module = mod, target = target, n_module_in_MatNet = sum(in_mod), n_hit = sum(in_mod & hit),
                 frac_module = mean(hit[in_mod]), frac_background = mean(hit), odds_ratio = unname(ft$estimate), p = ft$p.value)
}

# Same test counting unique AGIs
agi_mod <- dplyr::distinct(univ, module, AGI)
agi_bg <- dplyr::distinct(univ, AGI, in_longmod, in_DT)
fisher_agi <- function(mod, target) {
  a <- unique(agi_mod$AGI[agi_mod$module == mod])
  hit <- agi_bg$AGI[agi_bg[[target]]]
  x <- sum(a %in% hit)
  n <- length(a)
  N <- nrow(agi_bg)
  K <- length(hit)
  ft <- fisher.test(matrix(c(x, n - x, K - x, N - n - K + x), 2))
  tibble::tibble(module = mod, target = target, n_AGI_module = n, n_AGI_hit = x, odds_ratio_AGI = unname(ft$estimate),
                 p_AGI = ft$p.value)
}
mods_test <- univ %>%
  dplyr::count(module) %>%
  dplyr::filter(n >= MIN_GENES_MOD, !module %in% c("grey", "royalblue")) %>%
  dplyr::pull(module)
grid <- tidyr::expand_grid(module = mods_test, target = c("in_longmod", "in_DT"))
mod_enrich <- dplyr::bind_rows(Map(fisher_gene, grid$module, grid$target)) %>%
  dplyr::mutate(padj = p.adjust(p, "BH")) %>%
  dplyr::left_join(anchor, by = "module") %>%
  dplyr::arrange(target, padj)
readr::write_csv(mod_enrich, file.path(dirs$res, "03_module_enrichment_longevity_DT.csv"))
mod_enrich_agi <- dplyr::bind_rows(Map(fisher_agi, grid$module, grid$target)) %>%
  dplyr::mutate(padj_AGI = p.adjust(p_AGI, "BH")) %>%
  dplyr::left_join(dplyr::select(mod_enrich, module, target, odds_ratio_gene = odds_ratio, padj_gene = padj),
                   by = c("module", "target")) %>%
  dplyr::arrange(target, padj_AGI)
readr::write_csv(mod_enrich_agi, file.path(dirs$res, "04_module_enrichment_by_unique_AGI.csv"))
print(dplyr::filter(mod_enrich_agi, padj_AGI < ALPHA | padj_gene < ALPHA), n = 40)
# expected: turquoise in longevity (OR 4.2 per gene, 3.15 per AGI), salmon in
# longevity (~4.8), pink in desiccation (~10.6); black not enriched

fig <- mod_enrich %>%
  dplyr::mutate(log2_OR = log2(pmax(odds_ratio, 0.05)), sig = ifelse(padj < ALPHA, "BH < 0.05", "n.s."),
                target = factor(target, c("in_longmod", "in_DT"), c("Longevity module (Medicago)", "Desiccation tolerance (Medicago)")))
ord <- fig %>%
  dplyr::filter(target == "Longevity module (Medicago)") %>%
  dplyr::arrange(log2_OR) %>%
  dplyr::pull(module)
p <- ggplot2::ggplot(fig, ggplot2::aes(log2_OR, factor(module, ord), colour = sig)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::geom_point(ggplot2::aes(size = n_hit)) +
  ggplot2::facet_wrap(~ target) +
  ggplot2::scale_colour_manual(values = c("BH < 0.05" = "#1b7837", "n.s." = "grey60")) +
  ggplot2::labs(x = "log2 odds ratio (module vs network genes with AGI in MatNet)", y = NULL, colour = NULL,
                size = "Genes in target") +
  ggplot2::theme_bw(base_size = 13)
ggplot2::ggsave(file.path(dirs$fig, "01_soy_modules_vs_Medicago.png"), p, width = 11, height = 7, dpi = 300)

close_stage("07_external_validation_medicago", dirs)
