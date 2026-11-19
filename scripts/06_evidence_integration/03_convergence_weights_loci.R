##############################################################################
# Soybean seed longevity
# Stage 6c: genomic background, convergence test, weakest link, weight
# sensitivity and within-locus prioritization
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: answer the central question of the integration: do the GWAS candidates
# that rank high in one independent domain (temporal, GSE98199) also rank high
# in the other (SEA = atlas + network), more often than by chance?
#   1) DE rate of the candidates vs the whole tested transcriptome;
#   2) Spearman between layers in the 30 complete cases;
#   3) permutation test of convergence in the top K = 10;
#   4) weakest link: worst_rank = max(rank T, rank SEA), so a gene ranks well
#      only if it is good in BOTH domains (the opposite of RP/RRA with 2 lists);
#   5) Dirichlet weights for (T, A, N) and for (T, SEA);
#   6) prioritization inside each GWAS locus (locus map from stage 1).
#
# Run from the repository root, after scripts 01 and 02 of this stage:
#   source("scripts/06_evidence_integration/03_convergence_weights_loci.R")
#
# Packages: readr, dplyr, tidyr, tibble, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2"))

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "03_convergence_loci"))
f_int <- from_stage(stage, "01_integrated_matrix/02_integrated_evidence_matrix_39.csv")
f_cc <- from_stage(stage, "02_ranking_sensitivity/01_complete_case_layer_scores_30.csv")
f_de <- from_stage("02_maturation_timecourse", "07_deseq2_all_genes.csv")
f_locus <- from_stage("01_candidates", "02_locus_map.csv")
check_inputs(c(f_int, f_cc, f_de, f_locus))

K_TOP <- 10                                    # "top" for convergence (1/3 of the 30)
N_PERM <- 20000
N_WEIGHTS <- 10000
all_genes <- readr::read_csv(f_int, show_col_types = FALSE)
cc <- readr::read_csv(f_cc, show_col_types = FALSE)
(N_CC <- nrow(cc))                             # expected: 30

# ---------------------------------------------------------------------------
# 1) Differential expression background in GSE98199
# Exact two-sided binomial test of the candidate DE fraction against the
# genome rate. A small P with a higher candidate rate would be enrichment to
# report; a large P means "10 responded to maturation" does not set the
# candidates apart from the rest of the transcriptome.
# ---------------------------------------------------------------------------
genome_rate <- readr::read_csv(f_de, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(genome_tested = dplyr::n(), genome_DE = sum(significant)) %>%
  dplyr::mutate(genome_rate = genome_DE / genome_tested)
tested <- dplyr::filter(all_genes, tested_DESeq2 %in% TRUE)
(background_DE <- tibble::tibble(contrast = c("stage9_vs_7.2_fresh", "7.2_dried_vs_7.2_fresh"),
                                 cand_tested = nrow(tested),
                                 cand_DE = c(sum(tested$significant_maturation %in% TRUE),
                                             sum(tested$significant_dry %in% TRUE))) %>%
    dplyr::left_join(genome_rate, by = "contrast") %>%
    dplyr::mutate(cand_rate = cand_DE / cand_tested,
                  binom_p = mapply(function(x, n, p) binom.test(x, n, p)$p.value, cand_DE, cand_tested, genome_rate)))
readr::write_csv(background_DE, file.path(dirs$res, "01_background_DE_rates.csv"))
# expected: maturation 39.4% genome vs 10/31 (P ~ 0.47); drying 30.2% vs 6/31 (P ~ 0.24)

# ---------------------------------------------------------------------------
# 2) Agreement between layers (Spearman, 30 complete cases)
# Values near zero or negative = layers order the genes almost independently.
# ---------------------------------------------------------------------------
layer_mat <- as.matrix(dplyr::select(cc, T = rankfrac_T, A = rankfrac_A, N = rankfrac_N, SEA = rankfrac_SEA))
(layer_cor <- round(cor(layer_mat, method = "spearman"), 3))
readr::write_csv(tibble::as_tibble(layer_cor, rownames = "layer"), file.path(dirs$res, "02_layer_spearman.csv"))

# ---------------------------------------------------------------------------
# 3) Convergence permutation test
# Statistic: genes in the top K of temporal AND of SEA. Null: shuffle the
# temporal vector (breaks only the temporal-SEA link and keeps the atlas-network
# dependence inside SEA). 3-layer version: shuffle T and A.
# P = (1 + permutations >= observed) / (1 + N_PERM).
# High P = convergence no greater than chance; P < 0.05 would mean more genes
# than expected at the top of both domains at once.
# ---------------------------------------------------------------------------
rT <- cc$rankfrac_T * N_CC                     # back to positions 1..30
rA <- cc$rankfrac_A * N_CC
rN <- cc$rankfrac_N * N_CC
rS <- cc$rankfrac_SEA * N_CC
conv2 <- function(t, s, k) sum(t <= k & s <= k)
conv3 <- function(t, a, n, k) sum(t <= k & a <= k & n <= k)
set.seed(SEED_PERM)
null2 <- vapply(seq_len(N_PERM), function(i) conv2(sample(rT), rS, K_TOP), numeric(1))
null3 <- vapply(seq_len(N_PERM), function(i) conv3(sample(rT), sample(rA), rN, K_TOP), numeric(1))
obs2 <- conv2(rT, rS, K_TOP)
obs3 <- conv3(rT, rA, rN, K_TOP)
(convergence <- tibble::tibble(test = c("Temporal_and_SEA_topK", "Temporal_Atlas_Network_topK"), K = K_TOP,
                               observed = c(obs2, obs3), expected_null = c(mean(null2), mean(null3)),
                               p_perm = c((1 + sum(null2 >= obs2)), (1 + sum(null3 >= obs3))) / (1 + N_PERM),
                               n_perm = N_PERM))
readr::write_csv(convergence, file.path(dirs$res, "03_convergence_permutation_test.csv"))
# expected: observed 2, expected 3.3, P ~ 0.94

# ---------------------------------------------------------------------------
# 4) Weakest link (temporal + SEA)
# p_indep = P(max(R1, R2) <= w) = (floor(w) / N)^2 if domains were
# independent; BH over the 30 genes. A gene that is not significant after BH
# is reported as a descriptive priority, not as a test result.
# ---------------------------------------------------------------------------
weakest <- cc %>%
  dplyr::transmute(gene_id, rank_T = rT, rank_SEA = rS, rank_A = rA, rank_N = rN,
                   worst_rank_2dom = pmax(rank_T, rank_SEA), sum_rank_2dom = rank_T + rank_SEA,
                   n_domains_topK = (rank_T <= K_TOP) + (rank_SEA <= K_TOP),
                   n_layers_topK = (rank_T <= K_TOP) + (rank_A <= K_TOP) + (rank_N <= K_TOP),
                   p_indep = (floor(worst_rank_2dom) / N_CC)^2) %>%
  dplyr::mutate(p_indep_BH = p.adjust(p_indep, "BH"),
                rank_weakest_link = rank(worst_rank_2dom + sum_rank_2dom / 1000)) %>%   # sum breaks ties
  dplyr::arrange(rank_weakest_link)
readr::write_csv(weakest, file.path(dirs$res, "04_weakest_link_2domains.csv"))
# expected top: 04G192100 (worst rank 6) and 04G192900 (7), both from Wu rs46324094

# ---------------------------------------------------------------------------
# 5) Weight sensitivity
# Score = weighted geometric mean of rank fractions, exp(sum w_i log R_i).
# The primary ranking equals w = (0.50, 0.25, 0.25). Weights are drawn from a
# Dirichlet (normalized gammas); we count how often each gene is in the top 5/10.
# Dirichlet(1,1,1) over T, A, N gives atlas + network 2/3 of the weight on
# average, treating as two pieces of evidence what was defined as one domain;
# hence the second version over (T, SEA) with w_T ~ U(0, 1).
# ---------------------------------------------------------------------------
rdirichlet_cols <- function(n, alpha) {
  g <- matrix(rgamma(n * length(alpha), shape = alpha), nrow = length(alpha))
  sweep(g, 2, colSums(g), "/")
}

# Frequency in top 5/10 and rank quantiles over random weights
weight_freq <- function(L, alpha) {
  set.seed(SEED_PERM)
  r <- apply(log(L) %*% rdirichlet_cols(N_WEIGHTS, alpha), 2, rank)
  tibble::tibble(gene_id = cc$gene_id, freq_top5 = rowMeans(r <= 5), freq_top10 = rowMeans(r <= 10),
                 median_rank = apply(r, 1, median), q05_rank = apply(r, 1, quantile, 0.05),
                 q95_rank = apply(r, 1, quantile, 0.95))
}
weights3 <- weight_freq(cbind(cc$rankfrac_T, cc$rankfrac_A, cc$rankfrac_N), c(1, 1, 1)) %>%
  dplyr::mutate(rank_ref_50_25_25 = rank(as.vector(log(cbind(cc$rankfrac_T, cc$rankfrac_A, cc$rankfrac_N)) %*%
                                                     c(0.50, 0.25, 0.25)))) %>%
  dplyr::arrange(dplyr::desc(freq_top10), median_rank)
weights2 <- weight_freq(cbind(cc$rankfrac_T, cc$rankfrac_SEA), c(1, 1)) %>%
  dplyr::rename_with(~ paste0(.x, "_2dom"), -gene_id) %>%
  dplyr::left_join(dplyr::select(weights3, gene_id, freq_top10_3layer = freq_top10), by = "gene_id") %>%
  dplyr::arrange(dplyr::desc(freq_top10_2dom), median_rank_2dom)
readr::write_csv(weights3, file.path(dirs$res, "05_weight_sensitivity_3layers.csv"))
readr::write_csv(weights2, file.path(dirs$res, "06_weight_sensitivity_2domains.csv"))
# A gene stable in both versions has a position that hardly depends on the
# weights; a gene that drops in the 2-domain version relied on the extra weight
# the 3-layer version gives to SEA.

# ---------------------------------------------------------------------------
# 6) Prioritization inside each GWAS locus
# Complete cases are ordered by weakest link; incomplete genes enter with their
# coverage class and available ranks, without a position in the locus.
# A locus where each gene is good in one layer only has no preferred candidate,
# and the text says so.
# ---------------------------------------------------------------------------
locus_map <- readr::read_csv(f_locus, show_col_types = FALSE)
within_locus <- locus_map %>%
  dplyr::select(gene_id, study, chromosome, locus_id, n_genes_locus) %>%
  dplyr::left_join(dplyr::select(all_genes, gene_id, evidence_coverage_class, temporal_support_label,
                                 project_seed_context, network_support_label, functional_support_class,
                                 temporal_rank_available_all, atlas_rank_all_39, network_rank_available_all,
                                 sea_rank_available_all), by = "gene_id") %>%
  dplyr::left_join(dplyr::select(weakest, gene_id, rank_T, rank_SEA, worst_rank_2dom, n_domains_topK,
                                 rank_weakest_link), by = "gene_id") %>%
  dplyr::group_by(locus_id) %>%
  dplyr::mutate(n_complete_locus = sum(!is.na(worst_rank_2dom)),
                rank_in_locus = rank(rank_weakest_link, na.last = "keep", ties.method = "min"),
                best_in_locus = rank_in_locus %in% 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(chromosome, locus_id, rank_in_locus)
readr::write_csv(within_locus, file.path(dirs$res, "07_within_locus_prioritization.csv"))
(locus_summary <- within_locus %>%
    dplyr::group_by(locus_id, study, chromosome) %>%
    dplyr::summarise(n_genes = dplyr::n(), n_complete = sum(!is.na(worst_rank_2dom)),
                     genes = paste(gene_id, collapse = "; "),
                     best_gene = if (any(best_in_locus)) gene_id[best_in_locus][1] else NA_character_,
                     best_worst_rank = if (any(!is.na(worst_rank_2dom))) min(worst_rank_2dom, na.rm = TRUE) else NA_real_,
                     any_gene_both_domains_topK = any(n_domains_topK %in% 2), .groups = "drop") %>%
    dplyr::arrange(best_worst_rank))
readr::write_csv(locus_summary, file.path(dirs$res, "08_locus_summary.csv"))
# expected: Wu chr04 rs46324094 is the only locus with a gene in the top 10 of both domains

# ---------------------------------------------------------------------------
# 7) Figures
# ---------------------------------------------------------------------------
lab <- dplyr::filter(weakest, worst_rank_2dom <= 15 | rank_T <= 3 | rank_SEA <= 3)
p_wl <- ggplot2::ggplot(weakest, ggplot2::aes(rank_T, rank_SEA)) +
  ggplot2::annotate("rect", xmin = 0.5, xmax = K_TOP + 0.5, ymin = 0.5, ymax = K_TOP + 0.5, fill = "grey85", alpha = 0.6) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(data = lab, ggplot2::aes(label = gene_id), size = 2.8, nudge_y = 0.7, check_overlap = TRUE) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.05))) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = "Temporal rank (1 = strongest)", y = "SEA rank (atlas + network)",
                subtitle = paste0("Shaded: top ", K_TOP, " in both domains; complete cases (n = ", N_CC, ")")) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "01_weakest_link_temporal_vs_SEA.png"), p_wl, width = 8, height = 8, dpi = 300)

p_null <- dplyr::bind_rows(tibble::tibble(test = "Temporal and SEA", value = null2),
                           tibble::tibble(test = "Temporal, atlas and network", value = null3)) %>%
  ggplot2::ggplot(ggplot2::aes(value)) +
  ggplot2::geom_bar(fill = "grey60") +
  ggplot2::geom_vline(data = tibble::tibble(test = c("Temporal and SEA", "Temporal, atlas and network"),
                                            observed = c(obs2, obs3)),
                      ggplot2::aes(xintercept = observed), linetype = 2) +
  ggplot2::facet_wrap(~ test, scales = "free_y") +
  ggplot2::labs(x = paste0("Genes in the top ", K_TOP, " of all compared rankings"), y = "Permutations",
                subtitle = "Dashed line: observed value") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "02_convergence_permutation_null.png"), p_null, width = 10, height = 4.5, dpi = 300)

saveRDS(list(null2 = null2, null3 = null3, convergence = convergence, weakest = weakest),
        file.path("data", "processed", "integration_convergence_objects.rds"))
close_stage("06_evidence_integration_convergence", dirs)
