##############################################################################
# Soybean seed longevity
# Stage 7e: figure of the chromosome 4 locus programs as coexpression graphs
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: draw, for each gene of the convergent locus, the part of its module
# that surrounds it: the gene, its NET$n_neighbors most correlated genes of
# the module (the neighbourhood used in stage 5) and the top FIG$n_hubs hubs
# of the module by kWithin. Edges come from the signed hybrid adjacency of
# the network (cor^beta for positive correlations, 231 samples); only the
# FIG$edge_density strongest pairs of each subgraph are drawn, so that the
# picture stays readable. Node size and shade = degree in the drawn graph.
# Panel C repeats the statGraph result (script 04): mean spectral density of
# the module graph in each tissue, with the range over the subsets.
# The drawing is a summary, not a test; the tests are in stages 5 and 7d.
#
# Run from the repository root, after scripts 04 of this stage:
#   source("scripts/07_external_validation/05_locus_program_graphs.R")
#
# Packages: readr, dplyr, ggplot2, ggrepel, patchwork, igraph
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("readr", "dplyr", "ggplot2", "ggrepel", "patchwork"))
if (!requireNamespace("igraph", quietly = TRUE)) stop("Install first: igraph")

stage <- "07_external_validation"
dirs <- setup_stage(file.path(stage, "05_locus_graphs"))
f_expr <- file.path("data", "processed", "network_prevalence_log2TPM.rds")
f_mod <- from_stage("04_coexpression_network", "21_module_assignments.csv")
f_cent <- from_stage("04_coexpression_network", "28_gene_centrality_kME_kWithin.csv")
f_pos <- from_stage("01_candidates", "03_candidate_positions_a4.csv")
f_locus <- from_stage("06_evidence_integration", "06_final_classes/03_final_locus_table.csv")
f_class <- from_stage("06_evidence_integration", "06_final_classes/01_final_evidence_table_39_genes.csv")
f_anchor <- from_stage("06_evidence_integration", "04_maturation_anchor/05_candidate_module_anchor_39.csv")
f_dens <- from_stage(file.path(stage, "04_statgraph"), "06_spectral_density_partition1.csv")
check_inputs(c(f_expr, f_mod, f_cent, f_pos, f_locus, f_class, f_anchor, f_dens, RAW$info_annot))

# Prevalence-filtered matrix (all 231 samples), so that genes removed by the
# variance filter and projected onto a module can be drawn too
expr <- readRDS(f_expr)
mods <- readr::read_csv(f_mod, show_col_types = FALSE)
cent <- readr::read_csv(f_cent, show_col_types = FALSE)
annot <- read_info_annot() %>%
  dplyr::transmute(gene_id, label = dplyr::if_else(is.na(`arabi-symbol`) | `arabi-symbol` == "",
                                                   sub("\\.[0-9]+$", "", `Best-hit-arabi-name`),
                                                   sub(",.*", "", `arabi-symbol`)))

locus <- readr::read_csv(f_locus, show_col_types = FALSE) %>%
  dplyr::filter(locus_status == "convergent_locus") %>%
  dplyr::pull(locus_id)
# Genes of the locus in class A (the third gene, Glyma.04G192600, is class E)
(locus_genes <- readr::read_csv(f_pos, show_col_types = FALSE) %>%
    dplyr::filter(locus_id %in% locus) %>%
    dplyr::select(gene_id) %>%
    dplyr::inner_join(readr::read_csv(f_class, show_col_types = FALSE) %>% dplyr::select(gene_id, evidence_class),
                      by = "gene_id") %>%
    dplyr::filter(evidence_class == "A_convergent") %>%
    dplyr::inner_join(readr::read_csv(f_anchor, show_col_types = FALSE) %>% dplyr::select(gene_id, module_for_integration),
                      by = "gene_id"))
# expected: Glyma.04G192100 (black, projected) and Glyma.04G192900 (turquoise)

# ---------------------------------------------------------------------------
# 1) Subgraph around one candidate
# ---------------------------------------------------------------------------

# Nodes: candidate + its most correlated module genes + module hubs.
# Edges: signed hybrid adjacency, strongest FIG$edge_density of the pairs
locus_graph <- function(gene, module) {
  mg <- setdiff(mods$gene_id[mods$module == module], gene)
  r_g <- cor(expr[gene, ], t(expr[mg, ]))[1, ]
  neigh <- names(sort(r_g, decreasing = TRUE))[seq_len(NET$n_neighbors)]
  hubs <- cent %>%
    dplyr::filter(module == !!module) %>%
    dplyr::slice_max(kWithin, n = FIG$n_hubs) %>%
    dplyr::pull(gene_id)
  nodes <- unique(c(gene, neigh, hubs))
  r <- cor(t(expr[nodes, ]))
  adj <- ifelse(r > 0, r^NET$power, 0)
  diag(adj) <- 0
  thr <- quantile(adj[upper.tri(adj)], 1 - FIG$edge_density)
  g <- igraph::graph_from_adjacency_matrix((adj >= thr) * adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  igraph::V(g)$role <- dplyr::case_when(nodes == gene ~ "candidate", nodes %in% hubs ~ "module hub", TRUE ~ "neighbour")
  igraph::V(g)$degree <- igraph::degree(g)
  g
}

# ggplot of a subgraph in the style of BioNERO::plot_gcn (Fruchterman-Reingold
# layout, grey edges, nodes shaded by degree, labels on the candidate and on
# the FIG$n_labels genes of highest degree)
draw_graph <- function(g, title, max_degree) {
  set.seed(SEED_FIG)
  xy <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  nd <- tibble::tibble(gene_id = igraph::V(g)$name, x = xy[, 1], y = xy[, 2], degree = igraph::V(g)$degree,
                       role = igraph::V(g)$role) %>%
    dplyr::left_join(annot, by = "gene_id") %>%
    dplyr::mutate(text = paste0(sub("Glyma\\.", "", gene_id), dplyr::if_else(is.na(label), "", paste0("\n", label))),
                  show = role == "candidate" | rank(-degree, ties.method = "first") <= FIG$n_labels)
  el <- igraph::as_edgelist(g, names = FALSE)
  ed <- tibble::tibble(x = xy[el[, 1], 1], y = xy[el[, 1], 2], xend = xy[el[, 2], 1], yend = xy[el[, 2], 2])
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed, ggplot2::aes(x, y, xend = xend, yend = yend), colour = "grey80", alpha = 0.6,
                          linewidth = 0.3) +
    ggplot2::geom_point(data = dplyr::filter(nd, role != "candidate"),
                        ggplot2::aes(x, y, size = degree, alpha = degree, shape = role), colour = "midnightblue") +
    ggplot2::geom_point(data = dplyr::filter(nd, role == "candidate"), ggplot2::aes(x, y), colour = "firebrick",
                        size = 5, shape = 18) +
    ggrepel::geom_label_repel(data = dplyr::filter(nd, show), ggplot2::aes(x, y, label = text), size = 2.6,
                              colour = dplyr::if_else(dplyr::filter(nd, show)$role == "candidate", "firebrick", "grey30"),
                              lineheight = 0.9, min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, seed = SEED_FIG) +
    ggplot2::scale_shape_manual(values = c(`module hub` = 17, neighbour = 16), name = NULL) +
    ggplot2::scale_size_continuous(range = c(1, 6), limits = c(0, max_degree), name = "Degree") +
    ggplot2::scale_alpha_continuous(range = c(0.35, 0.95), limits = c(0, max_degree), name = "Degree") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10))
}

graphs <- lapply(seq_len(nrow(locus_genes)), function(i) locus_graph(locus_genes$gene_id[i], locus_genes$module_for_integration[i]))
names(graphs) <- locus_genes$gene_id
(graph_summary <- tibble::tibble(gene_id = locus_genes$gene_id, module = locus_genes$module_for_integration,
                                 n_nodes = vapply(graphs, igraph::vcount, numeric(1)),
                                 n_edges = vapply(graphs, igraph::ecount, numeric(1)),
                                 candidate_degree = vapply(seq_along(graphs), function(i) igraph::V(graphs[[i]])$degree[1], numeric(1)),
                                 max_degree = vapply(graphs, function(g) max(igraph::V(g)$degree), numeric(1))))
readr::write_csv(graph_summary, file.path(dirs$res, "01_locus_subgraph_summary.csv"))
edges <- dplyr::bind_rows(lapply(names(graphs), function(n) {
  igraph::as_data_frame(graphs[[n]], what = "edges") %>% dplyr::mutate(candidate = n, .before = 1)
}))
readr::write_csv(edges, file.path(dirs$res, "02_locus_subgraph_edges.csv"))
# Reading: a candidate with degree near the maximum sits in the core of its
# program; a low degree means it is attached to the program by few genes,
# which is expected for a gene projected after the variance filter. Expected
# here: Glyma.04G192100 degree 1 and Glyma.04G192900 degree 4, against a
# maximum of 21; both are at the edge of their programs, not hubs.

# ---------------------------------------------------------------------------
# 2) Spectral density by tissue (from script 04) and final figure
# ---------------------------------------------------------------------------
dens <- readr::read_csv(f_dens, show_col_types = FALSE) %>%
  dplyr::group_by(tissue, eigenvalue) %>%
  dplyr::summarise(mean = mean(density), lo = min(density), hi = max(density), .groups = "drop")
p_c <- ggplot2::ggplot(dens, ggplot2::aes(eigenvalue, mean, colour = tissue, fill = tissue)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::scale_y_sqrt() +
  ggplot2::coord_cartesian(xlim = c(-1, 4)) +
  ggplot2::labs(x = "Eigenvalue of the adjacency matrix", y = "Spectral density (square root scale)", colour = "Tissue",
                fill = "Tissue",
                title = paste0("C  ", locus_genes$module_for_integration[1], " module graph by tissue (spectral density)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "right")

titles <- paste0(c("A  ", "B  "), locus_genes$gene_id, " in ", locus_genes$module_for_integration)
top_degree <- max(graph_summary$max_degree)          # one degree scale for both panels
p_ab <- draw_graph(graphs[[1]], titles[1], top_degree) + draw_graph(graphs[[2]], titles[2], top_degree) +
  patchwork::plot_layout(guides = "collect")
fig <- p_ab / p_c + patchwork::plot_layout(heights = c(1.6, 1))
ggplot2::ggsave(file.path(dirs$fig, "01_locus_program_graphs.pdf"), fig, width = 10, height = 9)
ggplot2::ggsave(file.path(dirs$fig, "01_locus_program_graphs.png"), fig, width = 10, height = 9, dpi = 300)

close_stage("07_external_validation_locus_graphs", dirs)
