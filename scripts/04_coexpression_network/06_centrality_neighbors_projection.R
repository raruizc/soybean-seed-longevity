##############################################################################
# Soybean seed longevity
# Stage 4f: module membership, intramodular connectivity, hubs, candidate
# neighbours and post hoc projection of candidates outside the network
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: describe the position of each candidate in the frozen network.
#   kME     = Pearson correlation between a gene and its module eigengene;
#   kWithin = sum of signed hybrid adjacencies (max(r, 0)^beta) to the other
#             genes of the same module, computed in chunks (no 23k x 23k matrix);
#   hub     = kME > 0.8 AND top 10% kWithin in the module (BioNERO definition).
# Neighbours: the 30 genes of the same module with the highest adjacency to
# each candidate (used for GO in stage 5).
# Candidates removed by the variance filter, or left in grey, are NOT put back
# into the network: they are only correlated post hoc with the existing
# eigengenes (strongest positive r, since edges are positive in signed hybrid).
# A hub is a structural position, not evidence of causality.
#
# Run from the repository root, after script 05 of this stage:
#   source("scripts/04_coexpression_network/06_centrality_neighbors_projection.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))

stage <- "04_coexpression_network"
dirs <- setup_stage(stage)
f_net <- file.path("data", "processed", "network_BioNERO_blockwise.rds")
f_in <- file.path("data", "processed", "network_input_log2TPM.rds")
f_prev <- file.path("data", "processed", "network_prevalence_log2TPM.rds")
f_cand <- file.path(dirs$res, "24_candidate_module_assignment.csv")
f_ctx <- file.path(dirs$res, "26_module_strongest_bodypart.csv")
check_inputs(c(f_net, f_in, f_prev, f_cand, f_ctx))

net <- readRDS(f_net)
mat_net <- readRDS(f_in)
log_prev <- readRDS(f_prev)
cand <- readr::read_csv(f_cand, show_col_types = FALSE)
module_context <- readr::read_csv(f_ctx, show_col_types = FALSE)

modules <- tibble::tibble(gene_id = as.character(net$genes_and_modules$Genes),
                          module = as.character(net$genes_and_modules$Modules))
modules <- modules[match(rownames(mat_net), modules$gene_id), ]
MEs <- as.data.frame(net$MEs)[colnames(mat_net), , drop = FALSE]
MEs_bio <- MEs[, setdiff(colnames(MEs), "MEgrey"), drop = FALSE]
dat <- t(mat_net)                              # samples x genes, WGCNA orientation

# ---------------------------------------------------------------------------
# 1) Module membership (kME) of every gene in its own module
# ---------------------------------------------------------------------------

# Two-sided P of a Pearson correlation with n observations
cor_p <- function(r, n) {
  r <- pmin(pmax(r, -0.999999999), 0.999999999)
  2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), df = n - 2)
}
kme_all <- cor(dat, as.matrix(MEs_bio), method = NET$cor)
own_me <- paste0("ME", modules$module)
idx <- cbind(seq_len(nrow(modules)), match(own_me, colnames(kme_all)))   # NA column for grey
membership <- modules %>%
  dplyr::mutate(ME = ifelse(module == "grey", NA, own_me), kME = kme_all[idx],
                kME_pvalue = cor_p(kME, nrow(dat)), kME_BH = p.adjust(kME_pvalue, "BH"))

# ---------------------------------------------------------------------------
# 2) kWithin per module, in chunks of 250 genes
# ---------------------------------------------------------------------------

# Sum of signed hybrid adjacencies of each gene to the rest of its module
kwithin_module <- function(genes, chunk = 250) {
  x_all <- dat[, genes, drop = FALSE]
  starts <- seq(1, length(genes), by = chunk)
  unlist(lapply(starts, function(a) {
    ids <- genes[a:min(a + chunk - 1, length(genes))]
    adj <- pmax(cor(dat[, ids, drop = FALSE], x_all, method = NET$cor), 0)^NET$power
    adj[cbind(seq_along(ids), match(ids, genes))] <- 0          # no self-connection
    setNames(rowSums(adj), ids)
  }))
}
bio_modules <- setdiff(sort(unique(modules$module)), "grey")
kwithin <- unlist(lapply(bio_modules, function(m) kwithin_module(modules$gene_id[modules$module == m])))

centrality <- membership %>%
  dplyr::mutate(kWithin = unname(kwithin[gene_id])) %>%
  dplyr::group_by(module) %>%
  dplyr::mutate(kWithin_percentile = ifelse(module == "grey", NA, rank(kWithin, na.last = "keep") / sum(!is.na(kWithin))),
                top10_kWithin = module != "grey" & kWithin_percentile >= 1 - NET$hub_top,
                hub = module != "grey" & !is.na(kME) & kME > NET$hub_kme & top10_kWithin) %>%
  dplyr::ungroup()
readr::write_csv(centrality, file.path(dirs$res, "28_gene_centrality_kME_kWithin.csv"))
(hub_summary <- centrality %>%
    dplyr::filter(module != "grey") %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(module_size = dplyr::n(), n_kME_gt_0.8 = sum(kME > NET$hub_kme), n_hubs = sum(hub)) %>%
    dplyr::arrange(dplyr::desc(n_hubs)))
readr::write_csv(hub_summary, file.path(dirs$res, "29_hub_summary_by_module.csv"))

# ---------------------------------------------------------------------------
# 3) Candidates inside informative modules and their top-30 neighbours
# ---------------------------------------------------------------------------
cand_cent <- centrality %>%
  dplyr::filter(gene_id %in% cand$gene_id) %>%
  dplyr::mutate(in_informative_module = module != "grey") %>%
  dplyr::arrange(dplyr::desc(hub), dplyr::desc(kME))
readr::write_csv(cand_cent, file.path(dirs$res, "30_candidate_centrality.csv"))

# Top neighbours of one candidate by adjacency within its module
top_neighbors <- function(gene, module) {
  genes_mod <- setdiff(modules$gene_id[modules$module == module], gene)
  r <- as.numeric(cor(dat[, gene], dat[, genes_mod], method = NET$cor))
  tibble::tibble(candidate_gene = gene, module = module, neighbor_gene = genes_mod, correlation = r,
                 adjacency = pmax(r, 0)^NET$power) %>%
    dplyr::slice_max(adjacency, n = NET$n_neighbors, with_ties = FALSE) %>%
    dplyr::mutate(neighbor_rank = dplyr::row_number())
}
direct <- dplyr::filter(cand_cent, in_informative_module)
neighbors <- dplyr::bind_rows(Map(top_neighbors, direct$gene_id, direct$module)) %>%
  dplyr::left_join(dplyr::select(centrality, neighbor_gene = gene_id, neighbor_kME = kME,
                                 neighbor_kWithin_percentile = kWithin_percentile, neighbor_is_hub = hub),
                   by = "neighbor_gene")
readr::write_csv(neighbors, file.path(dirs$res, "31_candidate_top30_neighbors.csv"))

# ---------------------------------------------------------------------------
# 4) Post hoc projection: candidates removed by variance or left in grey
# The 7 candidates that fail the prevalence filter are not projected.
# ---------------------------------------------------------------------------
to_project <- cand$gene_id[cand$network_status == "eligible_but_removed_by_variance" |
                             (cand$network_status == "direct_member_of_network" & cand$module %in% "grey")]
to_project <- intersect(to_project, rownames(log_prev))

# Correlation (and BH within gene) of one gene with every biological eigengene
project_gene <- function(gene) {
  x <- log_prev[gene, colnames(mat_net)]
  ct <- lapply(MEs_bio, function(y) suppressWarnings(cor.test(x, y, method = NET$cor)))
  tibble::tibble(gene_id = gene, module = sub("^ME", "", names(ct)),
                 correlation = vapply(ct, function(z) unname(z$estimate), numeric(1)),
                 pvalue = vapply(ct, function(z) z$p.value, numeric(1))) %>%
    dplyr::mutate(padj_BH_within_gene = p.adjust(pvalue, "BH"))
}
projection <- dplyr::bind_rows(lapply(to_project, project_gene))
readr::write_csv(projection, file.path(dirs$res, "32_posthoc_projection_long.csv"))

projection_best <- projection %>%
  dplyr::group_by(gene_id) %>%
  dplyr::slice_max(correlation, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(gene_id, strongest_positive_module = module, strongest_positive_r = correlation,
                   strongest_positive_padj = padj_BH_within_gene, positive_r_ge_0_8 = correlation >= 0.8,
                   positive_r_ge_0_6 = correlation >= 0.6) %>%
  dplyr::left_join(module_context %>% dplyr::transmute(strongest_positive_module = sub("^ME", "", ME),
                                                       projected_module_bodypart = strongest_bodypart,
                                                       projected_module_bodypart_cor = strongest_cor,
                                                       projected_module_bodypart_padj = strongest_padj_BH_global),
                   by = "strongest_positive_module")

# ---------------------------------------------------------------------------
# 5) Network evidence of the 39 candidates (input of stages 5 and 6)
# ---------------------------------------------------------------------------
cand_network <- cand %>%
  dplyr::select(gene_id, network_status, direct_module = module, module_size) %>%
  dplyr::left_join(dplyr::select(cand_cent, gene_id, direct_kME = kME, direct_kME_BH = kME_BH, direct_kWithin = kWithin,
                                 direct_kWithin_percentile = kWithin_percentile, direct_hub = hub), by = "gene_id") %>%
  dplyr::left_join(projection_best, by = "gene_id") %>%
  dplyr::mutate(network_evidence_type = dplyr::case_when(
    network_status == "direct_member_of_network" & direct_module != "grey" ~ "direct_module_member",
    network_status == "direct_member_of_network" ~ "grey_with_posthoc_projection",
    network_status == "eligible_but_removed_by_variance" ~ "posthoc_projection_only",
    network_status == "not_eligible_prevalence" ~ "not_evaluated_in_global_network",
    TRUE ~ "review"))
readr::write_csv(cand_network, file.path(dirs$res, "33_candidate_network_evidence_39.csv"))
print(dplyr::count(cand_network, network_evidence_type))

p_cent <- ggplot2::ggplot(direct, ggplot2::aes(kME, kWithin_percentile, label = gene_id, shape = hub)) +
  ggplot2::geom_vline(xintercept = NET$hub_kme, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 1 - NET$hub_top, linetype = "dashed") +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(check_overlap = TRUE, nudge_y = 0.025, size = 3) +
  ggplot2::labs(x = "Module membership (kME)", y = "Intramodular connectivity percentile", shape = "Hub") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "08_candidate_kME_kWithin.png"), p_cent, width = 9, height = 7, dpi = 300)

close_stage("04_coexpression_network_centrality", dirs)
