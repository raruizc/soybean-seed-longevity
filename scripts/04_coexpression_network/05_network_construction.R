##############################################################################
# Soybean seed longevity
# Stage 4e: coexpression network (BioNERO/WGCNA, blockwise) and modules
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: build the main network on the frozen input (231 samples, 23,648 genes,
# beta 16, signed hybrid, Pearson) with BioNERO::exp2gcn_blockwise(), which
# calls WGCNA::blockwiseModules() in blocks of 5,000 genes and so never holds a
# 23k x 23k TOM in memory. Modules with eigengene correlation >= 0.80 are
# merged; minimum module size 30. Then: module sizes, eigengenes, where the 39
# candidates fall (in the network, eligible but removed by the variance
# filter, or not eligible) and module x body part association.
# Grey = genes not assigned to any module; it is kept in the tables but never
# interpreted as a module. Hubs, neighbours and post hoc projection are in
# script 06. Runtime: about 10 to 30 minutes.
#
# Run from the repository root, after script 04 of this stage:
#   source("scripts/04_coexpression_network/05_network_construction.R")
#
# Packages: readr, dplyr, tibble, ggplot2, WGCNA, BioNERO
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tibble", "ggplot2"))
for (p in c("WGCNA", "BioNERO")) if (!requireNamespace(p, quietly = TRUE)) stop("Install first: ", p)

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_log <- file.path("data", "processed", "network_prevalence_log2TPM.rds")
f_genes <- file.path(dirs$res, "20_network_genes.txt")
f_grid <- file.path(dirs$res, "19_candidate_retention_grid.csv")
f_dec <- file.path(dirs$res, "03_sample_curation_decisions.csv")
check_inputs(c(f_log, f_genes, f_grid, f_dec))

log_prev <- readRDS(f_log)
mat_net <- log_prev[readr::read_lines(f_genes), , drop = FALSE]
dim(mat_net)                                   # expected: 23,648 genes x 231 samples
stopifnot(WGCNA::goodSamplesGenes(t(mat_net), verbose = 0)$allOK)
saveRDS(mat_net, file.path("data", "processed", "network_input_log2TPM.rds"))
meta <- readr::read_csv(f_dec, show_col_types = FALSE)
meta <- meta[match(colnames(mat_net), meta$sample_id), ]
stopifnot(identical(meta$sample_id, colnames(mat_net)))

# ---------------------------------------------------------------------------
# 1) Network construction (the long step)
# ---------------------------------------------------------------------------
t0 <- Sys.time()
net <- BioNERO::exp2gcn_blockwise(exp = as.data.frame(mat_net), net_type = NET$type,
                                  module_merging_threshold = NET$merge_cor, SFTpower = NET$power,
                                  cor_method = NET$cor, max_block_size = NET$block,
                                  min_module_size = NET$min_module, verbose = 3)
elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(net, file.path("data", "processed", "network_BioNERO_blockwise.rds"))

# ---------------------------------------------------------------------------
# 2) Modules and eigengenes
# ---------------------------------------------------------------------------
modules <- tibble::tibble(gene_id = as.character(net$genes_and_modules$Genes),
                          module = as.character(net$genes_and_modules$Modules))
readr::write_csv(modules, file.path(dirs$res, "21_module_assignments.csv"))
(module_sizes <- modules %>%
    dplyr::count(module, name = "n_genes") %>%
    dplyr::mutate(is_grey = module == "grey") %>%
    dplyr::arrange(is_grey, dplyr::desc(n_genes)))
readr::write_csv(module_sizes, file.path(dirs$res, "22_module_sizes.csv"))
cat("Modules (without grey):", sum(!module_sizes$is_grey), "| grey genes:",
    sum(module_sizes$n_genes[module_sizes$is_grey]), "| minutes:", round(elapsed_min, 1), "\n")
# expected: 21 modules, 6,614 grey genes

MEs <- as.data.frame(net$MEs)
readr::write_csv(tibble::as_tibble(MEs, rownames = "sample_id"), file.path(dirs$res, "23_module_eigengenes.csv"))

# ---------------------------------------------------------------------------
# 3) Where the 39 candidates are
# ---------------------------------------------------------------------------
cand_net <- readr::read_csv(f_grid, show_col_types = FALSE) %>%
  dplyr::select(gene_id, pass_prevalence, variance_log2TPM, variance_percentile,
                retained_network = dplyr::all_of(paste0("retained_top", 100 * NET$var_keep))) %>%
  dplyr::left_join(modules, by = "gene_id") %>%
  dplyr::left_join(dplyr::select(module_sizes, module, module_size = n_genes), by = "module") %>%
  dplyr::mutate(network_status = dplyr::case_when(
    !pass_prevalence ~ "not_eligible_prevalence",
    !retained_network ~ "eligible_but_removed_by_variance",
    !is.na(module) ~ "direct_member_of_network",
    TRUE ~ "review")) %>%
  dplyr::arrange(network_status, module, gene_id)
readr::write_csv(cand_net, file.path(dirs$res, "24_candidate_module_assignment.csv"))
print(dplyr::count(cand_net, network_status))  # expected: 18 direct, 14 removed by variance, 7 not eligible

# ---------------------------------------------------------------------------
# 4) Module x body part (body part as categorical trait, Spearman)
# BioProject is not used as a trait because it is partly confounded with
# anatomy (script 03). "Strongest" body part is descriptive, not exclusive.
# ---------------------------------------------------------------------------
trait <- data.frame(body_part = meta$body_part, row.names = meta$sample_id)
me_trait <- BioNERO::module_trait_cor(exp = as.data.frame(mat_net), metadata = trait, MEs = MEs,
                                      cor_method = "spearman") %>%
  as.data.frame() %>%
  dplyr::mutate(padj_BH_global = p.adjust(pvalue, "BH")) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(padj_BH_within_trait = p.adjust(pvalue, "BH")) %>%
  dplyr::ungroup()
readr::write_csv(me_trait, file.path(dirs$res, "25_module_bodypart_correlations.csv"))

module_context <- me_trait %>%
  dplyr::group_by(ME) %>%
  dplyr::slice_max(abs(cor), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(ME, strongest_bodypart = trait, strongest_cor = cor, strongest_pvalue = pvalue,
                   strongest_padj_BH_global = padj_BH_global)
readr::write_csv(module_context, file.path(dirs$res, "26_module_strongest_bodypart.csv"))

p_trait <- me_trait %>%
  dplyr::mutate(sig = dplyr::case_when(padj_BH_global < 0.001 ~ "***", padj_BH_global < 0.01 ~ "**",
                                       padj_BH_global < 0.05 ~ "*", TRUE ~ "")) %>%
  ggplot2::ggplot(ggplot2::aes(trait, ME, fill = cor)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sig), size = 3) +
  ggplot2::scale_fill_gradient2(midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(x = "Body part", y = "Module eigengene", fill = "Spearman",
                subtitle = "* BH P < 0.05; ** < 0.01; *** < 0.001") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(dirs$fig, "06_module_bodypart_heatmap.png"), p_trait, width = 9,
                height = max(7, 0.25 * ncol(MEs)), dpi = 300)
tryCatch(ggplot2::ggsave(file.path(dirs$fig, "07_dendrogram_modules.png"), BioNERO::plot_dendro_and_colors(net),
                         bg = "white", width = 10, height = 6, dpi = 300),
         error = function(e) message("Dendrogram not saved: ", conditionMessage(e)))

readr::write_csv(tibble::tibble(parameter = c("n_genes", "n_samples", names(NET)[-(1:3)], "elapsed_minutes"),
                                value = c(nrow(mat_net), ncol(mat_net), unlist(NET[-(1:3)]), round(elapsed_min, 1))),
                 file.path(dirs$res, "27_network_parameters.csv"))

close_stage("04_coexpression_network_build", dirs)
