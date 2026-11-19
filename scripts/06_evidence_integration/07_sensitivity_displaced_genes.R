##############################################################################
# Soybean seed longevity
# Stage 6g: sensitivity to the assembly mix-up of the Thenveettil candidates
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: the Thenveettil et al. (2025) SNPs are on Wm82.a4, but the genes in
# their +-10 kb window were looked up on a2/a1/ZH13 (stage 1, section 4). For
# 8 SNPs our genes lie 0.24 to 1.95 Mb from the SNP on a4. The integration is
# redone without those genes:
#   S0_all     : the 39 genes (reproduces scripts 03 and 06; check)
#   S1_200kb   : drops genes > 200 kb from their own SNP
#   S2_50kb    : also drops genes at 50-200 kb
# Wu et al. (2024) genes and genes without SNP position always stay.
# In each scenario: re-rank the layers within the remaining complete cases,
# Spearman T x SEA, convergence test (20,000 permutations) and classes A-E with
# K fixed (10/15) or proportional (n/3 and n/2, the rule that gave 10/15 for 30).
#
# Run from the repository root, after script 06 of this stage:
#   source("scripts/06_evidence_integration/07_sensitivity_displaced_genes.R")
#
# Packages: readr, dplyr, tidyr, tibble
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "tibble"))

stage <- "06_evidence_integration"
dirs <- setup_stage(file.path(stage, "07_displaced_gene_sensitivity"))
f_int <- from_stage(stage, "01_integrated_matrix/02_integrated_evidence_matrix_39.csv")
f_cls <- from_stage(stage, "06_final_classes/01_final_evidence_table_39_genes.csv")
f_pos <- from_stage("01_candidates", "03_candidate_positions_a4.csv")
check_inputs(c(f_int, f_cls, f_pos))
N_PERM <- 20000

dist_snp <- dplyr::select(readr::read_csv(f_pos, show_col_types = FALSE), gene_id, study, snp_pos, dist_snp_kb)
SCENARIOS <- c(S0_all = Inf, S1_200kb = FAR_SNP_KB, S2_50kb = NEAR_SNP_KB)
dropped <- function(lim) dist_snp$gene_id[!is.na(dist_snp$dist_snp_kb) & dist_snp$dist_snp_kb > lim]
lengths(lapply(SCENARIOS, dropped))            # expected: 0, 13 and 16 genes
all_genes <- readr::read_csv(f_int, show_col_types = FALSE)
class_main <- dplyr::select(readr::read_csv(f_cls, show_col_types = FALSE), gene_id, class_main = evidence_class)

# ---------------------------------------------------------------------------
# 1) One scenario: re-rank, convergence, weakest link and classes
# ---------------------------------------------------------------------------
run_scenario <- function(name, lim) {
  cc <- all_genes %>%
    dplyr::filter(!gene_id %in% dropped(lim), !is.na(temporal_domain_score), !is.na(atlas_context_score),
                  !is.na(network_rank_coupling)) %>%
    dplyr::mutate(rankfrac_T = rank_frac(temporal_domain_score), rankfrac_A = rank_frac(atlas_context_score),
                  rankfrac_N = rank_frac(network_rank_coupling), rankfrac_SEA = rank_frac(sqrt(rankfrac_A * rankfrac_N)))
  n <- nrow(cc)
  rT <- cc$rankfrac_T * n
  rS <- cc$rankfrac_SEA * n
  test_k <- function(k, rule) {
    set.seed(SEED_PERM)
    null <- vapply(seq_len(N_PERM), function(i) sum(sample(rT) <= k & rS <= k), numeric(1))
    obs <- sum(rT <= k & rS <= k)
    tibble::tibble(scenario = name, n_complete_cases = n, K_rule = rule, K = k, observed = obs,
                   expected_null = mean(null), p_perm = (1 + sum(null >= obs)) / (1 + N_PERM))
  }
  k_prop <- round(n / 3)
  conv <- dplyr::bind_rows(test_k(10, "fixed_10"), if (k_prop != 10) test_k(k_prop, "proportional_n/3")) %>%
    dplyr::mutate(spearman_T_SEA = cor(rT, rS, method = "spearman"))
  cls <- function(kt, km) dplyr::case_when(rT <= kt & rS <= kt ~ "A_convergent", rT <= km & rS <= km ~ "B_balanced_moderate",
                                           rT <= kt & rS > km ~ "C_temporal_dominant", rS <= kt & rT > km ~ "D_SEA_dominant",
                                           TRUE ~ "E_low_support")
  genes <- tibble::tibble(scenario = name, gene_id = cc$gene_id, rank_T = rT, rank_SEA = rS, worst_rank = pmax(rT, rS),
                          class_fixed_10_15 = cls(10, 15), class_prop = cls(k_prop, round(n / 2))) %>%
    dplyr::arrange(worst_rank, rank_T + rank_SEA)
  list(conv = conv, genes = genes)
}
res <- Map(run_scenario, names(SCENARIOS), SCENARIOS)

# ---------------------------------------------------------------------------
# 2) Results
# If P stays high in S1/S2, the lack of convergence did not depend on the
# displaced genes. Genes that change class appear in 03.
# ---------------------------------------------------------------------------
(conv_all <- dplyr::bind_rows(lapply(res, `[[`, "conv")))
readr::write_csv(conv_all, file.path(dirs$res, "01_convergence_by_scenario.csv"))
# expected (proportional K): S0 P = 0.94, S1 P = 0.79, S2 P = 0.35
genes_all <- dplyr::bind_rows(lapply(res, `[[`, "genes")) %>%
  dplyr::left_join(class_main, by = "gene_id") %>%
  dplyr::left_join(dplyr::select(dist_snp, gene_id, dist_snp_kb), by = "gene_id")
readr::write_csv(genes_all, file.path(dirs$res, "02_genes_ranks_classes_by_scenario.csv"))
(changes <- genes_all %>%
    dplyr::filter(scenario != "S0_all", class_fixed_10_15 != class_main | class_prop != class_main) %>%
    dplyr::select(scenario, gene_id, class_main, class_fixed_10_15, class_prop, rank_T, rank_SEA))
readr::write_csv(changes, file.path(dirs$res, "03_class_changes.csv"))
# expected: class A (04G192100, 04G192900) stable; 06G153700 rises to A in S2

close_stage("06_evidence_integration_displaced", dirs)
