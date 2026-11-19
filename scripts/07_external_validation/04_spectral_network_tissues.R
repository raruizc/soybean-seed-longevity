##############################################################################
# Soybean seed longevity
# Stage 7d: spectral comparison of coexpression graphs between seed tissues
# and random graph model selection (statGraph)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Question 1: is the coexpression program that carries the convergent
# chromosome 4 locus (module of its best gene, black) wired the same way in
# seed coat, embryo and whole seed? Each tissue is split at random into
# disjoint subsets of SG$subset_n samples; each subset gives one graph of the
# module genes (Pearson on log2 TPM, the SG$density fraction of gene pairs
# with the highest positive correlation become edges, as in a signed network).
# Fixing the edge density makes the comparison about how edges are arranged,
# not about how strong correlations are in each tissue. The graphs are
# compared by their spectral density (eigenvalues of the adjacency matrix)
# with ANOGVA (Fujita et al., 2017) and, two tissues at a time, with the
# Takahashi et al. (2012) test, both from statGraph (Guzman et al., 2026).
# The ANOGVA of statGraph 1.0.6 does not divide the permutation statistic by
# the number of groups, as it does for the observed one, which pushes P
# towards 1 (checked by simulation, see docs/statgraph_anogva_note.md). The
# function anogva_fixed() below uses the same statistic on both sides.
# Specificity: the same statistic for random gene sets of the same size,
# drawn from the network genes, on the same sample subsets.
# Question 2: which random graph model (Erdos-Renyi, geometric, Watts-Strogatz,
# Barabasi-Albert) is closest to the module graph built with all samples,
# chosen by the Graph Information Criterion (GIC, Takahashi et al., 2012).
# Limit: tissue is confounded with BioProject (stage 4c); subsets of one
# tissue share BioProjects, so a tissue effect can still carry a study effect.
#
# Run from the repository root, after stages 4 and 6:
#   source("scripts/07_external_validation/04_spectral_network_tissues.R")
#
# Packages: readr, dplyr, tidyr, ggplot2, igraph, statGraph (CRAN or
# github.com/fujita-team/statGraph)
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "tidyr", "ggplot2"))
for (p in c("igraph", "statGraph")) if (!requireNamespace(p, quietly = TRUE)) stop("Install first: ", p)

stage <- "07_external_validation"
dirs <- setup_stage(file.path(stage, "04_statgraph"))
f_expr <- file.path("data", "processed", "network_input_log2TPM.rds")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_dec <- from_stage("04_coexpression_network", "03_sample_curation_decisions.csv")
f_locus <- from_stage("06_evidence_integration", "06_final_classes/03_final_locus_table.csv")
f_anchor <- from_stage("06_evidence_integration", "04_maturation_anchor/05_candidate_module_anchor_39.csv")
check_inputs(c(f_expr, f_mod, f_dec, f_locus, f_anchor))

expr <- readRDS(f_expr)
mods <- readr::read_csv(f_mod, show_col_types = FALSE)
meta <- readr::read_csv(f_dec, show_col_types = FALSE) %>%
  dplyr::filter(sample_id %in% colnames(expr))

# Module of the best gene of the convergent locus (direct or post-hoc projection)
best_gene <- readr::read_csv(f_locus, show_col_types = FALSE) %>%
  dplyr::filter(locus_status == "convergent_locus") %>%
  dplyr::pull(best_gene)
(target_module <- readr::read_csv(f_anchor, show_col_types = FALSE) %>%
    dplyr::filter(gene_id %in% best_gene) %>%
    dplyr::pull(module_for_integration))
module_genes <- mods$gene_id[mods$module == target_module]
length(module_genes)
# expected: Glyma.04G192100, projected to black (661 genes)

# ---------------------------------------------------------------------------
# 1) Graph construction and tissue subsets
# ---------------------------------------------------------------------------

# Unweighted graph of a gene set in a set of samples: the top `density`
# fraction of gene pairs by positive Pearson correlation become edges
build_graph <- function(genes, samples, density) {
  r <- suppressWarnings(cor(t(expr[genes, samples, drop = FALSE])))
  r[is.na(r)] <- 0                              # gene without variance in the subset: no edges
  thr <- quantile(r[upper.tri(r)], 1 - density)
  A <- (r >= thr & r > 0) * 1
  diag(A) <- 0
  igraph::graph_from_adjacency_matrix(A, mode = "undirected")
}

# Disjoint random subsets of SG$subset_n samples per tissue; tissues with fewer
# than SG$min_graphs subsets are dropped
split_tissues <- function() {
  sets <- lapply(SG$tissues, function(t) {
    s <- sample(meta$sample_id[meta$body_part == t])
    k <- floor(length(s) / SG$subset_n)
    if (k < SG$min_graphs) return(NULL)
    split(s[seq_len(k * SG$subset_n)], rep(seq_len(k), each = SG$subset_n))
  })
  names(sets) <- SG$tissues
  sets <- sets[!vapply(sets, is.null, logical(1))]
  list(samples = unlist(sets, recursive = FALSE), tissue = rep(names(sets), lengths(sets)))
}

# ANOGVA statistic and permutation P with both sides on the same scale:
# mean distance (KL) between each tissue's mean spectral density and the
# mean of the tissue means; labels permuted among graphs
anogva_fixed <- function(graphs, labels, n_perm, dist = "KL") {
  graphs <- statGraph:::set.list.spectral.density(graphs)
  stat <- function(lab) {
    means <- lapply(split(seq_along(graphs), lab), function(i) statGraph:::get.mean.spectral.density(graphs[i]))
    overall <- means[[1]]
    overall$y <- Reduce("+", lapply(means, function(d) d$y)) / length(means)
    mean(vapply(means, function(d) statGraph:::distance(d, overall, dist = dist), numeric(1)))
  }
  obs <- stat(labels)
  null <- if (n_perm > 0) replicate(n_perm, stat(sample(labels))) else numeric(0)
  c(statistic = obs, p_value = (1 + sum(null >= obs)) / (1 + n_perm))
}

set.seed(SEED_SG)
partitions <- replicate(SG$n_partitions, split_tissues(), simplify = FALSE)
table(partitions[[1]]$tissue)
# expected with 10 samples per subset: embryo 6, seed 7, seed coat 5
# (cotyledon has 2 subsets and is dropped)

# ---------------------------------------------------------------------------
# 2) Tissue comparison of the target module, by partition and edge density
# ---------------------------------------------------------------------------

# ANOGVA and pairwise Takahashi tests for one partition and one density
test_partition <- function(i, density) {
  part <- partitions[[i]]
  graphs <- lapply(part$samples, function(s) build_graph(module_genes, s, density))
  an <- anogva_fixed(graphs, part$tissue, SG$n_perm)
  pairs <- combn(unique(part$tissue), 2, simplify = FALSE)
  tk <- dplyr::bind_rows(lapply(pairs, function(pr) {
    t2 <- statGraph::graph.takahashi.test(graphs[part$tissue == pr[1]], graphs[part$tissue == pr[2]], maxBoot = SG$n_perm)
    tibble::tibble(comparison = paste(pr, collapse = " vs "), JS = unname(t2$statistic), p_value = t2$p.value)
  }))
  dplyr::bind_rows(tibble::tibble(comparison = "ANOGVA (all tissues)", JS = NA_real_, KL_statistic = an[["statistic"]],
                                  p_value = an[["p_value"]]), tk) %>%
    dplyr::mutate(partition = i, density = density, .before = 1)
}

grid <- tidyr::expand_grid(i = seq_len(SG$n_partitions), density = SG$density)
(tissue_tests <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(k) test_partition(grid$i[k], grid$density[k]))))
readr::write_csv(tissue_tests, file.path(dirs$res, "01_tissue_tests_by_partition.csv"))

(tissue_summary <- tissue_tests %>%
    dplyr::group_by(density, comparison) %>%
    dplyr::summarise(n_partitions = dplyr::n(), median_p = median(p_value), frac_p_lt_0.05 = mean(p_value < 0.05),
                     median_KL = median(KL_statistic), median_JS = median(JS), .groups = "drop"))
readr::write_csv(tissue_summary, file.path(dirs$res, "02_tissue_tests_summary.csv"))
# Reading: P < 0.05 in most partitions means the module graph is arranged
# differently between tissues (the pairwise tests say which ones). P spread
# over [0, 1] means no difference beyond what random subsets of 10 samples
# give. Both readings are conditional on BioProject: see the header.

# ---------------------------------------------------------------------------
# 3) Specificity: the same ANOGVA statistic for random gene sets
# ---------------------------------------------------------------------------
net_genes <- setdiff(rownames(expr), module_genes)

# ANOGVA statistic (no permutation) of one random gene set on one partition
random_stat <- function(k) {
  genes <- sample(net_genes, length(module_genes))
  part <- partitions[[(k - 1) %% SG$n_partitions + 1]]
  graphs <- lapply(part$samples, function(s) build_graph(genes, s, SG$density_main))
  anogva_fixed(graphs, part$tissue, 0)[["statistic"]]
}
set.seed(SEED_SG + 1)
null_stats <- vapply(seq_len(SG$n_random_sets), random_stat, numeric(1))
module_stats <- tissue_tests %>%
  dplyr::filter(density == SG$density_main, comparison == "ANOGVA (all tissues)") %>%
  dplyr::pull(KL_statistic)
(specificity <- tibble::tibble(module = target_module, n_genes = length(module_genes), density = SG$density_main,
                               module_median_KL = median(module_stats), random_median_KL = median(null_stats),
                               random_q95_KL = quantile(null_stats, 0.95), n_random_sets = SG$n_random_sets,
                               p_more_rewired = (1 + sum(null_stats >= median(module_stats))) / (1 + SG$n_random_sets),
                               p_more_conserved = (1 + sum(null_stats <= median(module_stats))) / (1 + SG$n_random_sets)))
readr::write_csv(specificity, file.path(dirs$res, "03_module_vs_random_gene_sets.csv"))
readr::write_csv(tibble::tibble(random_set = seq_along(null_stats), KL_statistic = null_stats),
                 file.path(dirs$res, "04_random_gene_set_statistics.csv"))
# Reading: if the module statistic sits inside the random distribution, the
# tissue difference (if any) is a property of the tissues, shared by any gene
# set, and not of this program. A small p_more_rewired means the module is
# rewired between tissues more than genes in general; a small
# p_more_conserved means its wiring is kept across tissues more than for
# random genes, i.e. the program is shared by seed coat, embryo and seed.
# Caution: module genes were grouped because they correlate across these same
# samples, so some stability is built in. The fairer reference is the other
# modules, below.

# Same statistic for every other module at least as large as the target,
# with SG$module_draws random draws of the same number of genes, over the
# first SG$module_partitions partitions (to keep the run time reasonable)
other_mods <- mods %>%
  dplyr::filter(module != "grey", module != target_module) %>%
  dplyr::count(module) %>%
  dplyr::filter(n >= length(module_genes)) %>%
  dplyr::pull(module)

# ANOGVA statistic of a draw of one module over the partitions, median
module_stat <- function(mod) {
  genes <- sample(mods$gene_id[mods$module == mod], length(module_genes))
  median(vapply(seq_len(SG$module_partitions), function(i) {
    part <- partitions[[i]]
    graphs <- lapply(part$samples, function(s) build_graph(genes, s, SG$density_main))
    anogva_fixed(graphs, part$tissue, 0)[["statistic"]]
  }, numeric(1)))
}
set.seed(SEED_SG + 3)
(module_ref <- tidyr::expand_grid(module = other_mods, draw = seq_len(SG$module_draws)) %>%
    dplyr::mutate(median_KL = vapply(module, module_stat, numeric(1))) %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(median_KL = median(median_KL), .groups = "drop") %>%
    dplyr::bind_rows(tibble::tibble(module = target_module, median_KL = median(module_stats[seq_len(SG$module_partitions)]))) %>%
    dplyr::mutate(is_target = module == target_module, rank_low_to_high = rank(median_KL)) %>%
    dplyr::arrange(median_KL))
readr::write_csv(module_ref, file.path(dirs$res, "04b_module_vs_other_modules.csv"))
# Reading: the target ranked first (lowest statistic) means its wiring is
# the most similar across tissues among modules of comparable size; a middle
# rank means its stability is what any coexpression module shows.

# ---------------------------------------------------------------------------
# 4) Random graph model closest to the module graph (all samples)
# ---------------------------------------------------------------------------
g_all <- build_graph(module_genes, meta$sample_id, SG$density_main)
set.seed(SEED_SG + 2)
ms <- statGraph::graph.model.selection(g_all, models = SG$models, search = "ternary", eps = 0.01, ngraphs = SG$ngraphs_gic)
(model_sel <- tibble::tibble(module = target_module, n_genes = igraph::vcount(g_all), n_edges = igraph::ecount(g_all),
                             model = rownames(ms$estimates), parameter = ms$estimates[, "param"],
                             GIC = ms$estimates[, "GIC"]) %>%
    dplyr::arrange(GIC))
readr::write_csv(model_sel, file.path(dirs$res, "05_graph_model_selection.csv"))
# Reading: the lowest GIC is the model whose spectrum is closest. BA (scale
# free) agrees with the premise behind the soft threshold of the network;
# GRG (geometric) is common for correlation graphs, where genes behave as
# points in a latent expression space; ER would mean no structure beyond the
# edge count.

# ---------------------------------------------------------------------------
# 5) Figure: mean spectral density per tissue (first partition, main density)
# ---------------------------------------------------------------------------
part <- partitions[[1]]
graphs <- statGraph:::set.list.spectral.density(lapply(part$samples, function(s) build_graph(module_genes, s, SG$density_main)))
dens <- dplyr::bind_rows(lapply(seq_along(graphs), function(i) {
  tibble::tibble(tissue = part$tissue[i], graph = i, eigenvalue = graphs[[i]]$density$x, density = graphs[[i]]$density$y)
}))
p_sd <- ggplot2::ggplot(dens, ggplot2::aes(eigenvalue, density, group = graph, colour = tissue)) +
  ggplot2::geom_line(alpha = 0.6) +
  ggplot2::scale_y_sqrt() +                     # the spike at 0 (isolated genes) hides the rest otherwise
  ggplot2::coord_cartesian(xlim = c(-3, 6)) +
  ggplot2::labs(x = "Eigenvalue of the adjacency matrix", y = "Spectral density", colour = "Tissue",
                title = paste0("Module ", target_module, ", ", SG$subset_n, " samples per graph, edge density ",
                               SG$density_main)) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "01_spectral_density_by_tissue.png"), p_sd, width = 8, height = 5, dpi = 300)
readr::write_csv(dens, file.path(dirs$res, "06_spectral_density_partition1.csv"))   # used by the figure of script 05

close_stage("07_external_validation_statgraph", dirs)
