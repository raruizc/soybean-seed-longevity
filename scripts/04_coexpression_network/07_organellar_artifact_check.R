##############################################################################
# Soybean seed longevity
# Stage 4g: organellar-copy check of the modules (royalblue)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: the royalblue module (~44 genes, photosynthesis terms, up after drying)
# could be a nuclear program or an artifact. Nuclear genes whose best
# Arabidopsis hit is ATCG (chloroplast) or ATMG (mitochondrion) are nuclear
# copies of organellar DNA (NUPTs/NUMTs). Reads from plastid transcripts map
# to them, so a module made of such genes tracks the fraction of organellar
# RNA in the library, not a nuclear gene program. For each module we test
# whether the organellar fraction exceeds the network background (one-sided
# binomial, BH across modules).
#
# Run from the repository root, after script 05 of this stage and stage 2:
#   source("scripts/04_coexpression_network/07_organellar_artifact_check.R")
#
# Packages: readr, dplyr, tidyr, stringr
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "stringr"))

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_mod <- file.path(dirs$res, "21_module_assignments.csv")
f_de <- from_stage("02_maturation_timecourse", "07_deseq2_all_genes.csv")
check_inputs(c(f_mod, f_de, RAW$info_annot))

annot <- read_info_annot() %>%
  dplyr::transmute(gene_id, arabi_hit = `Best-hit-arabi-name`, arabi_defline = `arabi-defline`,
                   organelle = dplyr::case_when(stringr::str_detect(arabi_hit, "^ATCG") ~ "plastid_like",
                                                stringr::str_detect(arabi_hit, "^ATMG") ~ "mito_like",
                                                TRUE ~ "nuclear"))
mods <- readr::read_csv(f_mod, show_col_types = FALSE) %>%
  dplyr::left_join(annot, by = "gene_id") %>%
  dplyr::mutate(organelle = dplyr::coalesce(organelle, "nuclear"))

# ---------------------------------------------------------------------------
# 1) Organellar fraction per module against the network background
# ---------------------------------------------------------------------------
bg <- mean(mods$organelle != "nuclear")
(frac <- mods %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(n_genes = dplyr::n(), n_plastid = sum(organelle == "plastid_like"),
                     n_mito = sum(organelle == "mito_like"), .groups = "drop") %>%
    dplyr::mutate(frac_organelle = (n_plastid + n_mito) / n_genes,
                  p_binom = mapply(function(k, n) binom.test(k, n, bg, alternative = "greater")$p.value,
                                   n_plastid + n_mito, n_genes),
                  padj = p.adjust(p_binom, "BH")) %>%
    dplyr::arrange(dplyr::desc(frac_organelle)))
readr::write_csv(frac, file.path(dirs$res, "34_organelle_fraction_by_module.csv"))

# ---------------------------------------------------------------------------
# 2) Royalblue gene by gene, with its GSE98199 responses
# ---------------------------------------------------------------------------
de <- readr::read_csv(f_de, show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = strip_version(GSE98199_gene_id)) %>%
  dplyr::select(gene_id, contrast, log2FoldChange, padj) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = c(log2FoldChange, padj))
royal <- mods %>%
  dplyr::filter(module == "royalblue") %>%
  dplyr::left_join(de, by = "gene_id") %>%
  dplyr::arrange(organelle, gene_id)
readr::write_csv(royal, file.path(dirs$res, "35_royalblue_genes.csv"))
print(dplyr::count(royal, organelle))          # expected: 34 plastid_like, 3-4 mito_like of 44

# Reading: if most of royalblue is plastid_like/mito_like, the module is an
# organellar-RNA artifact (the share of chloroplast transcripts in the library
# changes with drying). It is reported as such and NOT as a photosynthesis
# program. Any other module far above the background (padj < 0.05) needs the
# same caveat.

close_stage("04_coexpression_network_organelle", dirs)
