##############################################################################
# Soybean seed longevity
# Stage 2a: import, quality control and descriptive profiles (GSE98199)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: audit the GSE98199 count matrix (Pereira Lima et al., 2017; cultivar
# BRS 284, stages R7.1 to R9, plus dried R7.2), check the quality of the 11
# libraries (PCA, distances, replicate correlation) and build the normalized
# profiles and descriptive metrics of the 39 candidates.
# Only 7.2F, 7.2D and 9 have replicates: the full trajectory is DESCRIPTIVE.
# Inference is done in script 02 (replicated contrasts).
#
# Run from the repository root:
#   source("scripts/02_maturation_timecourse/01_qc_profiles_GSE98199.R")
#
# Packages: DESeq2, readr, dplyr, tibble, ggplot2
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("DESeq2", "readr", "dplyr", "tibble", "ggplot2"))

stage <- "02_maturation_timecourse"
dirs <- setup_stage(stage)
check_inputs(c(RAW$gse_counts, RAW$gse_meta, RAW$aliases))

SAMPLES <- c("Gm_7.1f", "Gm_7.2fa", "Gm_7.2fb", "Gm_7.2da", "Gm_7.2db",
             "Gm_7.3", "Gm_8.1", "Gm_8.2", "Gm_8.3", "Gm_9a", "Gm_9b")

# ---------------------------------------------------------------------------
# 1) Count matrix: reading and audit
# Checks: 11 libraries in the expected order, no NA, integers only (raw counts,
# required by DESeq2), technical row "*" removed.
# ---------------------------------------------------------------------------
counts <- read.delim(gzfile(RAW$gse_counts), header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[rownames(counts) != "*", SAMPLES, drop = FALSE]
stopifnot(!anyNA(counts), all(vapply(counts, function(x) all(x == round(x)), logical(1))))
count_mat <- as.matrix(counts)
storage.mode(count_mat) <- "integer"
dim(count_mat)                                 # expected: 56,044 genes x 11
(audit <- tibble::tibble(sample = SAMPLES, total_reads = colSums(count_mat),
                         genes_detected = colSums(count_mat > 0)))
readr::write_csv(audit, file.path(dirs$res, "01_library_audit.csv"))

meta <- readr::read_csv(RAW$gse_meta, show_col_types = FALSE) %>%
  dplyr::rename(sample = description) %>%
  dplyr::mutate(stage = factor(as.character(stage), c("7.1", "7.2", "7.3", "8.1", "8.2", "8.3", "9")),
                treatment = factor(treatment, c("fresh", "dried")))
stopifnot(identical(meta$sample, SAMPLES))

# ---------------------------------------------------------------------------
# 2) Global QC: VST, PCA, library distances and replicate correlation
# Reading: replicates should cluster in the PCA and correlate highly (> 0.95).
# If 7.2D separates from 7.2F on PC1/PC2, drying changes the transcriptome
# broadly, which justifies the drying contrast in script 02.
# ---------------------------------------------------------------------------
dds_all <- DESeq2::DESeqDataSetFromMatrix(count_mat, as.data.frame(meta, row.names = meta$sample), design = ~ 1)
dds_all <- DESeq2::estimateSizeFactors(dds_all)
vsd <- DESeq2::vst(dds_all, blind = TRUE)

pca <- DESeq2::plotPCA(vsd, intgroup = c("stage", "treatment"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, label = name, shape = treatment)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(nudge_y = 1.2, size = 3) +
  ggplot2::labs(x = paste0("PC1: ", pct[1], "%"), y = paste0("PC2: ", pct[2], "%"), shape = "Condition") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(dirs$fig, "01_PCA_libraries.png"), p_pca, width = 8, height = 6, dpi = 300)

dist_mat <- as.matrix(dist(t(SummarizedExperiment::assay(vsd))))
png(file.path(dirs$fig, "02_library_distance_heatmap.png"), width = 1700, height = 1500, res = 180)
heatmap(dist_mat, symm = TRUE)
dev.off()
readr::write_csv(tibble::as_tibble(dist_mat, rownames = "sample"), file.path(dirs$res, "02_library_distances.csv"))

log_norm <- log2(DESeq2::counts(dds_all, normalized = TRUE) + 1)
(rep_cor <- tibble::tibble(
  group = c("7.2 fresh", "7.2 dried", "9 fresh"),
  pearson_log2 = c(cor(log_norm[, "Gm_7.2fa"], log_norm[, "Gm_7.2fb"]),
                   cor(log_norm[, "Gm_7.2da"], log_norm[, "Gm_7.2db"]),
                   cor(log_norm[, "Gm_9a"], log_norm[, "Gm_9b"]))))
readr::write_csv(rep_cor, file.path(dirs$res, "03_replicate_correlation.csv"))

# ---------------------------------------------------------------------------
# 3) Candidates in GSE98199 (Wm82.a2 IDs) through the pangene aliases
# ---------------------------------------------------------------------------
aliases <- readr::read_csv(RAW$aliases, show_col_types = FALSE) %>%
  dplyr::mutate(alias_gene_id = sub("\\.m[0-9]+$", "", alias_gene_id)) %>%
  dplyr::distinct(Wm82_a4_id, alias_gene_id)
cand_map <- tibble::tibble(GSE98199_gene_id = rownames(count_mat), alias_gene_id = strip_version(rownames(count_mat))) %>%
  dplyr::inner_join(aliases, by = "alias_gene_id") %>%
  dplyr::distinct(Wm82_a4_id, .keep_all = TRUE) %>%
  dplyr::select(Wm82_a4_id, GSE98199_gene_id) %>%
  dplyr::arrange(Wm82_a4_id)
stopifnot(dplyr::n_distinct(cand_map$Wm82_a4_id) == 39)
readr::write_csv(cand_map, file.path(dirs$res, "04_candidate_map_GSE98199.csv"))

# ---------------------------------------------------------------------------
# 4) Normalized profiles and descriptive metrics of the candidates
# Counts normalized by size factors; replicates averaged.
# rho_stage = Spearman between expression and the order of the 7 fresh stages
# (n = 7, descriptive, no P). Values near +1 or -1 = monotonic trend.
# ---------------------------------------------------------------------------
norm <- DESeq2::counts(dds_all, normalized = TRUE)[cand_map$GSE98199_gene_id, , drop = FALSE]
profiles <- tibble::as_tibble(norm) %>%
  dplyr::transmute(Wm82_a4_id = cand_map$Wm82_a4_id,
                   `7.1` = Gm_7.1f, `7.2F` = (Gm_7.2fa + Gm_7.2fb) / 2, `7.2D` = (Gm_7.2da + Gm_7.2db) / 2,
                   `7.3` = Gm_7.3, `8.1` = Gm_8.1, `8.2` = Gm_8.2, `8.3` = Gm_8.3, `9` = (Gm_9a + Gm_9b) / 2)
readr::write_csv(profiles, file.path(dirs$res, "05_candidate_profiles_normalized.csv"))

fresh_stages <- c("7.1", "7.2F", "7.3", "8.1", "8.2", "8.3", "9")
m_fresh <- as.matrix(profiles[, fresh_stages])
metrics <- profiles %>%
  dplyr::mutate(mean_norm_all = rowMeans(as.matrix(profiles[, -1])),
                log2FC_9_vs_7.2F = log2((`9` + 1) / (`7.2F` + 1)),
                log2FC_7.2D_vs_7.2F = log2((`7.2D` + 1) / (`7.2F` + 1)),
                rho_stage = suppressWarnings(apply(m_fresh, 1, function(v) cor(1:7, v, method = "spearman"))))
readr::write_csv(metrics, file.path(dirs$res, "06_candidate_descriptive_metrics.csv"))

close_stage("02_maturation_timecourse_qc", dirs)
