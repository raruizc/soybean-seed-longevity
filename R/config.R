##############################################################################
# Soybean seed longevity: GWAS candidates confronted with expression evidence
# Shared configuration: input paths and fixed study parameters
# Author: Ricardo Antonio Ruiz Cardozo
#
# Every script runs from the REPOSITORY ROOT (opening soybean-seed-longevity.Rproj
# in RStudio sets the working directory) and starts with:
#   source("R/config.R")
#   source("R/utils.R")
##############################################################################

if (!file.exists("R/config.R")) stop("Run from the repository root (open the .Rproj file)")

# ---------------------------------------------------------------------------
# Raw inputs (data/raw/, not versioned; see data/raw/README.md)
# ---------------------------------------------------------------------------
RAW <- list(
  gwas_xlsx    = file.path("data", "raw", "gwas", "gwas_candidates_curated.xlsx"),
  aliases      = file.path("data", "raw", "gwas", "candidate_aliases.csv"),
  info_annot   = file.path("data", "raw", "annotation", "glyma.Wm82.gnm4.ann1.T8TQ.info_annot.txt.gz"),
  gff          = file.path("data", "raw", "annotation", "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz"),
  gse_counts   = file.path("data", "raw", "GSE98199", "GSE98199_processed.txt.gz"),
  gse_meta     = file.path("data", "raw", "GSE98199", "metadata_GSE98199.csv"),
  atlas_tpm    = file.path("data", "raw", "atlas", "SEA_39genes_TPM.tsv"),
  atlas_meta   = file.path("data", "raw", "atlas", "SEA_39genes_metadata.tsv"),
  network_dir  = file.path("data", "raw", "network"),
  pl_dir       = file.path("data", "raw", "external", "pereira_lima_2017"),
  righetti_dir = file.path("data", "raw", "external", "righetti_2015"),
  qtl_table    = file.path("data", "raw", "external", "qtl", "storability_qtl_curated.csv")
)

# ---------------------------------------------------------------------------
# Fixed parameters (set before looking at the results)
# ---------------------------------------------------------------------------
ALPHA        <- 0.05       # significance level after BH
MIN_LFC      <- 1          # minimum |log2FC| for differential expression
LOCUS_WIN_BP <- 250000     # SNPs on the same chromosome closer than this = same locus
SNP_WIN_BP   <- 10000      # gene search window around each SNP (Thenveettil et al., 2025)
NEAR_SNP_KB  <- 50         # listed gene within 50 kb of its SNP in Wm82.a4 = consistent position
FAR_SNP_KB   <- 200        # beyond this = displaced gene (sensitivity scenario S1)
SEED_PERM    <- 20261101   # permutations and weight draws (integration and sensitivity)
SEED_COMP    <- 20261102   # competitive test of modules along maturation
SEED_QTL     <- 20261116   # QTL colocalization permutations (stage 7c)
N_PERM_QTL   <- 10000     # random placements of the loci (stage 7c)
QTL_FLANK_BP <- c(0, 1e6)  # QTL as published and extended by 1 Mb on each side
SEED_SG      <- 20261117   # statGraph tissue comparison (stage 7d)

# Expression atlas: body parts treated as seed-related (radicle excluded)
SEED_PARTS <- c("seed", "embryo", "cotyledon", "endosperm", "seed coat", "suspensor")

# Coexpression network (chosen in 04_coexpression_network/04, frozen afterwards)
NET <- list(
  min_tpm = 1, min_prop = 0.25,       # prevalence filter: TPM >= 1 in >= 25% of samples
  zk_cut = -2,                        # sample connectivity outliers (BioNERO default)
  var_keep = 0.75,                    # keep the 75% most variable genes
  power = 16, type = "signed hybrid", cor = "pearson",
  merge_cor = 0.80, min_module = 30, block = 5000,
  hub_kme = 0.80, hub_top = 0.10, n_neighbors = 30
)

# Spectral comparison of module graphs between tissues (stage 7d, statGraph)
SG <- list(
  tissues = c("seed coat", "embryo", "seed", "cotyledon"),
  subset_n = 10, min_graphs = 3,      # samples per graph; tissues with < 3 graphs are dropped
  density = c(0.02, 0.05, 0.10),      # fraction of gene pairs kept as edges
  density_main = 0.05,
  n_partitions = 20, n_perm = 1000,   # random splits of the samples; permutations per test
  n_random_sets = 100,                # random gene sets for the specificity check
  module_draws = 3, module_partitions = 5,  # draws of each larger module at the target size
  models = c("ER", "GRG", "WS", "BA"), ngraphs_gic = 10
)

# Figure of the chromosome 4 locus programs as graphs (stage 7e)
SEED_FIG <- 20261119
FIG <- list(n_hubs = 5, edge_density = 0.25, n_labels = 3)
