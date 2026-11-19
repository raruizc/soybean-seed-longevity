##############################################################################
# Soybean seed longevity
# Stage 6b: ranking sensitivity (RP, Borda, RRA, leave-one-layer-out, top-k)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: check how much the order of the candidates depends on the way layers are
# aggregated. Everything uses the SAME complete-case universe (temporal + atlas
# + network), with each layer re-ranked inside that universe.
# Caution: with only 2 lists, RP and RRA reward a gene that dominates ONE
# layer; this script is therefore a sensitivity analysis, and convergence is
# tested formally in script 03.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/06_evidence_integration/02_ranking_sensitivity.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2, RobustRankAggreg
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))
if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) stop("Install first: RobustRankAggreg")

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "02_ranking_sensitivity"))
f_int <- from_stage(stage, "01_integrated_matrix/02_integrated_evidence_matrix_39.csv")
check_inputs(f_int)

all_genes <- readr::read_csv(f_int, show_col_types = FALSE)
rank_pos <- function(x) rank(x, ties.method = "average", na.last = "keep")
geom2 <- function(a, b) sqrt(a * b)
geom3 <- function(a, b, c) (a * b * c)^(1/3)
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# ---------------------------------------------------------------------------
# 1) Complete-case universe and re-ranked layers
# Expected N = 30. Stage 6a scores used different denominators (31 temporal,
# 39 atlas, 32 network); here they are re-ranked within the 30.
# ---------------------------------------------------------------------------
cc <- all_genes %>%
  dplyr::filter(!is.na(temporal_domain_score), !is.na(atlas_context_score), !is.na(network_rank_coupling)) %>%
  dplyr::mutate(rankfrac_T = rank_frac(temporal_domain_score), rankfrac_A = rank_frac(atlas_context_score),
                rankfrac_N = rank_frac(network_rank_coupling), SEA_score_CC = geom2(rankfrac_A, rankfrac_N),
                rankfrac_SEA = rank_frac(SEA_score_CC))   # SEA is ranked again before combining
(N_CC <- nrow(cc))

# ---------------------------------------------------------------------------
# 2) Aggregation strategies
# Primary (2 domains): RP and Borda of temporal + SEA. Sensitivity (3 layers):
# RP3 and Borda3 with atlas and network apart. Leave-one-layer-out: pairs.
# These strategies are deterministic functions of the same 3 vectors, so the
# top-k frequency among them is NOT statistical stability.
# ---------------------------------------------------------------------------
cc <- cc %>%
  dplyr::mutate(rank_primary_RP = rank_pos(geom2(rankfrac_T, rankfrac_SEA)),
                rank_primary_Borda = rank_pos((rankfrac_T + rankfrac_SEA) / 2),
                rank_expanded_RP3 = rank_pos(geom3(rankfrac_T, rankfrac_A, rankfrac_N)),
                rank_expanded_Borda3 = rank_pos((rankfrac_T + rankfrac_A + rankfrac_N) / 3),
                rank_LOO_no_temporal = rank_pos(geom2(rankfrac_A, rankfrac_N)),
                rank_LOO_no_atlas = rank_pos(geom2(rankfrac_T, rankfrac_N)),
                rank_LOO_no_network = rank_pos(geom2(rankfrac_T, rankfrac_A)))

# ---------------------------------------------------------------------------
# 3) Robust rank aggregation (Kolde et al., 2012)
# Each list = genes ordered best to worst; N = N_CC. exact = FALSE applies the
# Bonferroni correction to the rho score, so Score saturates at 1 with ties.
# If no gene has Score < 0.05, RRA is not presented as a significance test.
# ---------------------------------------------------------------------------
ordered_ids <- function(score) cc$gene_id[order(score)]
# RRA of a list of ordered gene vectors; returns score and rank per gene
rra <- function(glist, label) {
  x <- RobustRankAggreg::aggregateRanks(glist = glist, N = N_CC, method = "RRA", full = TRUE, exact = FALSE)
  out <- tibble::tibble(gene_id = x$Name, score = x$Score, rank = rank_pos(x$Score))
  names(out)[2:3] <- c(paste0(label, "_score"), paste0("rank_", label))
  out
}
rra_2 <- rra(list(Temporal = ordered_ids(cc$rankfrac_T), SEA = ordered_ids(cc$rankfrac_SEA)), "RRA2")
rra_3 <- rra(list(Temporal = ordered_ids(cc$rankfrac_T), Atlas = ordered_ids(cc$rankfrac_A),
                  Network = ordered_ids(cc$rankfrac_N)), "RRA3")
cc <- cc %>%
  dplyr::left_join(rra_2, by = "gene_id") %>%
  dplyr::left_join(rra_3, by = "gene_id")
min(cc$RRA2_score)                             # previous run: 0.08 (no gene < 0.05)

strategy_ranks <- dplyr::select(cc, gene_id, rank_primary_RP, rank_primary_Borda, rank_expanded_RP3,
                                rank_expanded_Borda3, rank_RRA2, rank_RRA3, rank_LOO_no_temporal,
                                rank_LOO_no_atlas, rank_LOO_no_network)
readr::write_csv(cc, file.path(dirs$res, "01_complete_case_layer_scores_30.csv"))
readr::write_csv(strategy_ranks, file.path(dirs$res, "02_strategy_rankings.csv"))

# ---------------------------------------------------------------------------
# 4) Agreement among strategies, top-k frequency and leave-one-layer-out shifts
# High correlation among RP, Borda and RRA is partly mechanical (same inputs).
# Low correlation between leave-one-out pairs means the layers order the genes
# differently, i.e. little convergence.
# ---------------------------------------------------------------------------
spearman <- cor(as.matrix(strategy_ranks[, -1]), method = "spearman")
readr::write_csv(tibble::as_tibble(spearman, rownames = "strategy"), file.path(dirs$res, "03_rank_spearman_matrix.csv"))

rank_long <- tidyr::pivot_longer(strategy_ranks, -gene_id, names_to = "strategy", values_to = "rank")
topk <- rank_long %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(strategies_n = dplyr::n(), top5_fraction = mean(rank <= 5), top10_fraction = mean(rank <= 10),
                   median_rank = median(rank), best_rank = min(rank), worst_rank = max(rank)) %>%
  dplyr::arrange(dplyr::desc(top10_fraction), median_rank)
readr::write_csv(topk, file.path(dirs$res, "04_topk_stability.csv"))

loo <- cc %>%
  dplyr::transmute(gene_id, baseline_rank = rank_expanded_RP3,
                   shift_no_temporal = rank_LOO_no_temporal - baseline_rank,
                   shift_no_atlas = rank_LOO_no_atlas - baseline_rank,
                   shift_no_network = rank_LOO_no_network - baseline_rank,
                   max_abs_LOO_shift = pmax(abs(shift_no_temporal), abs(shift_no_atlas), abs(shift_no_network))) %>%
  dplyr::arrange(dplyr::desc(max_abs_LOO_shift))
readr::write_csv(loo, file.path(dirs$res, "05_leave_one_layer_out_shifts.csv"))

(audit <- tibble::tibble(
  metric = c("complete_case_candidates", "strategies", "top10_in_all_strategies", "top5_in_all_strategies",
             "median_pairwise_spearman", "min_RRA2_score"),
  value = c(N_CC, ncol(strategy_ranks) - 1, sum(topk$top10_fraction == 1), sum(topk$top5_fraction == 1),
            median(spearman[upper.tri(spearman)]), min(cc$RRA2_score))))
readr::write_csv(audit, file.path(dirs$res, "06_audit.csv"))
# expected: 30, 9, 1, 0, 0.738. Few genes in the top 10 of every strategy and
# low leave-one-out correlations mean the fine order is not robust: the text
# should talk about evidence classes, not exact positions.

p_cor <- tibble::as_tibble(spearman, rownames = "s1") %>%
  tidyr::pivot_longer(-s1, names_to = "s2", values_to = "rho") %>%
  ggplot2::ggplot(ggplot2::aes(s2, s1, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.8) +
  ggplot2::scale_fill_gradient(low = "grey95", high = "grey20", limits = c(0, 1)) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Spearman rho") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(dirs$fig, "01_rank_correlation_heatmap.png"), p_cor, width = 11, height = 9, dpi = 300)

close_stage("06_evidence_integration_ranking_sensitivity", dirs)
