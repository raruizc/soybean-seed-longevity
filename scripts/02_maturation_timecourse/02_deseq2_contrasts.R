##############################################################################
# Soybean seed longevity
# Stage 2b: differential expression in the replicated contrasts (GSE98199)
# Author: Ricardo Antonio Ruiz Cardozo
#
# Aim: test the only two contrasts with biological replicates on both sides,
# drying (7.2D vs 7.2F) and late maturation (9 vs 7.2F), with DESeq2 (Love et
# al., 2014), and build the temporal evidence table of the 39 candidates used
# in the integration (stage 6).
# DE criterion: padj (BH) < 0.05 and |log2FC| >= 1 in the original test. The
# shrunken log2FC (lfcShrink, type = "normal") is kept as effect size; the
# original padj is not replaced.
#
# Run from the repository root, after script 01 of this stage:
#   source("scripts/02_maturation_timecourse/02_deseq2_contrasts.R")
#
# Packages: DESeq2, readr, dplyr, tibble
##############################################################################

source("R/config.R")
source("R/utils.R")
load_packages(c("DESeq2", "readr", "dplyr", "tibble"))

stage <- "02_maturation_timecourse"
dirs <- setup_stage(stage)
f_map <- file.path(dirs$res, "04_candidate_map_GSE98199.csv")
check_inputs(c(RAW$gse_counts, f_map))

# ---------------------------------------------------------------------------
# 1) Counts of the replicated libraries and minimum filter
# Tested gene = count >= 10 in at least 2 of the 6 libraries. A candidate that
# does not pass becomes "not tested, low count" (information, not NA).
# ---------------------------------------------------------------------------
REPLICATED <- c("Gm_7.2fa", "Gm_7.2fb", "Gm_7.2da", "Gm_7.2db", "Gm_9a", "Gm_9b")
counts <- read.delim(gzfile(RAW$gse_counts), header = TRUE, row.names = 1, check.names = FALSE)
mat <- as.matrix(counts[rownames(counts) != "*", REPLICATED])
storage.mode(mat) <- "integer"

coldata <- data.frame(row.names = REPLICATED,
                      group = factor(rep(c("stage7.2_fresh", "stage7.2_dried", "stage9_fresh"), each = 2),
                                     c("stage7.2_fresh", "stage7.2_dried", "stage9_fresh")))
dds <- DESeq2::DESeqDataSetFromMatrix(mat, coldata, design = ~ group)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
nrow(dds)                                      # expected: 30,710 tested genes
dds_fit <- DESeq2::DESeq(dds)

# ---------------------------------------------------------------------------
# 2) Contrasts, shrinkage and MA plots
# ---------------------------------------------------------------------------

# Wald test + shrunken log2FC for one contrast against fresh 7.2
contrast_vs_fresh <- function(level, label) {
  res <- DESeq2::results(dds_fit, contrast = c("group", level, "stage7.2_fresh"), alpha = ALPHA)
  coef <- paste0("group_", level, "_vs_stage7.2_fresh")
  shr <- DESeq2::lfcShrink(dds_fit, coef = coef, res = res, type = "normal")
  png(file.path(dirs$fig, paste0("03_MA_", label, ".png")), width = 1500, height = 1200, res = 180)
  DESeq2::plotMA(shr, alpha = ALPHA, main = label)
  dev.off()
  tibble::as_tibble(as.data.frame(res), rownames = "GSE98199_gene_id") %>%
    dplyr::mutate(shrunk_log2FoldChange = shr$log2FoldChange, contrast = label,
                  significant = !is.na(padj) & padj < ALPHA & abs(log2FoldChange) >= MIN_LFC)
}
res_dry <- contrast_vs_fresh("stage7.2_dried", "7.2_dried_vs_7.2_fresh")
res_mat <- contrast_vs_fresh("stage9_fresh", "stage9_vs_7.2_fresh")

all_genes <- dplyr::bind_rows(res_dry, res_mat)
readr::write_csv(all_genes, file.path(dirs$res, "07_deseq2_all_genes.csv"))
(genome_summary <- all_genes %>%
    dplyr::group_by(contrast) %>%
    dplyr::summarise(tested = sum(!is.na(padj)), DE = sum(significant),
                     up = sum(significant & log2FoldChange > 0), down = sum(significant & log2FoldChange < 0)))
readr::write_csv(genome_summary, file.path(dirs$res, "08_DE_genome_summary.csv"))

# ---------------------------------------------------------------------------
# 3) Temporal evidence table of the 39 candidates (input of stage 6)
# ---------------------------------------------------------------------------

# Add the contrast suffix to the columns (padj -> padj_dry etc.)
wide <- function(res, suffix) {
  x <- dplyr::select(res, GSE98199_gene_id, baseMean, log2FoldChange, shrunk_log2FoldChange,
                     lfcSE, stat, pvalue, padj, significant)
  names(x)[-1] <- paste0(names(x)[-1], "_", suffix)
  x
}

cand_map <- readr::read_csv(f_map, show_col_types = FALSE)
cand <- cand_map %>%
  dplyr::mutate(tested_DESeq2 = GSE98199_gene_id %in% rownames(dds_fit), filtered_low_count = !tested_DESeq2) %>%
  dplyr::left_join(wide(res_dry, "dry"), by = "GSE98199_gene_id") %>%
  dplyr::left_join(wide(res_mat, "maturation"), by = "GSE98199_gene_id") %>%
  dplyr::mutate(evidence_pattern = dplyr::case_when(
    filtered_low_count ~ "low_count_not_tested",
    significant_dry & significant_maturation ~ "both",
    significant_maturation ~ "maturation_only",
    significant_dry ~ "drying_only",
    TRUE ~ "neither"))
readr::write_csv(cand, file.path(dirs$res, "09_candidate_temporal_evidence.csv"))
print(dplyr::count(cand, evidence_pattern))

# Reading: the number of DE candidates only makes sense against the genome
# (08_). If the DE fraction among candidates is close to the genome fraction,
# the GWAS list does not stand out in the maturation program; the formal test
# is done in the integration (stage 6, genomic background).

close_stage("02_maturation_timecourse_deseq2", dirs)
