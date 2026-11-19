# soybean-seed-longevity

**Prioritizing GWAS candidates for soybean seed longevity with independent expression evidence**

Ricardo Antonio Ruiz Cardozo, Graduate Program in Statistics and Agricultural
Experimentation (PPGEEA), ESALQ/University of São Paulo, Brazil.

## What this study does

Two genome-wide association studies (GWAS) of soybean seed longevity proposed
39 candidate genes in 21 loci (Wu et al., 2024; Thenveettil et al., 2025). A
GWAS locus points to a region, not to a gene, and the candidate lists were
built with different methods and genome assemblies. This repository checks
those candidates against public expression data that the GWAS did not use:

| Layer | Data | Question |
|---|---|---|
| Temporal | GSE98199, seed maturation R7 to R9 and drying (Pereira Lima et al., 2017) | Does the gene respond to late maturation, when longevity is acquired? |
| Atlas | Soybean Expression Atlas v2 (Almeida-Silva et al., 2023) | Is the gene expressed in seed tissues, consistently across studies? |
| Network | Coexpression network of 231 curated Williams 82 seed samples (BioNERO/WGCNA) | Which program (module) does the gene belong to, and does that module follow maturation? |
| Function | GO enrichment of modules and candidate neighbourhoods | What do those programs do? |
| External | Pereira Lima et al. (2017) longevity lists; *Medicago truncatula* longevity network (Righetti et al., 2015) | Is the program conserved in another legume and another study? |

The integration asks whether candidates that rank high in one independent
domain (temporal) also rank high in the other (atlas + network, "SEA") more
often than by chance, and classifies each gene and locus by where its support
comes from.

**Main findings (to be read with the manuscript):**
- Convergence between domains is not greater than chance (permutation P = 0.94);
  the result holds after removing genes displaced by the assembly mix-up (P = 0.79 and 0.35).
- One locus is convergent: Wu et al. (2024) QTN rs46324094 on chromosome 4
  (Glyma.04G192100 and Glyma.04G192900), stable in every sensitivity scenario.
- The Thenveettil et al. (2025) SNPs are on assembly Wm82.a4, but the genes in
  their windows were taken from Wm82.a2/a1/ZH13: 8 of 18 SNPs have the listed
  gene more than 200 kb away. On chromosome 13 the a4 window falls in a
  lipoxygenase cluster.
- The module of Glyma.04G192900 (turquoise; RNA processing and translation) is
  enriched in the *Medicago* longevity module; a seed-coat module (black,
  Glyma.04G192100) has no counterpart there.

## How to run

1. Open `soybean-seed-longevity.Rproj` in RStudio (the working directory becomes
   the repository root). R >= 4.3; packages listed below.
2. Download the raw data into `data/raw/` as described in [`data/raw/README.md`](data/raw/README.md).
3. Run everything with `source("run_all.R")`, or one script at a time in the
   numbered order.

Each script reads only from `data/raw/`, `data/processed/` and the `results/`
of earlier stages, writes to `results/<stage>/` and `figures/<stage>/`, and saves
its `sessionInfo()` in `logs/`. Paths and fixed parameters (random seeds,
thresholds, network settings) live in `R/config.R`; helpers used by more than
one stage live in `R/utils.R`. Each script states in its comments the values
expected at key steps, so a re-run can be checked line by line.

Runtime on a laptop: stages 1 to 3 take minutes; stage 4 can take about an hour
(the soft-threshold grid and the blockwise network are the long steps; saved
soft-threshold diagnostics are reused on re-runs); stages 5 to 7 take minutes.

## Structure

```
R/
  config.R        input paths and fixed parameters
  utils.R         shared helpers (stage folders, rank fractions, annotation, GO ORA)
scripts/
  01_candidates/              candidate table, GWAS loci, Wm82.a4 positions, genes at each SNP
  02_maturation_timecourse/   GSE98199: QC, profiles, DESeq2 (drying and late maturation)
  03_expression_atlas/        Soybean Expression Atlas: body part summary, Tau, seed context, study support
  04_coexpression_network/    sample curation, Z.K QC, confounding, soft threshold, network, hubs, organellar check
  05_functional_enrichment/   GO of modules and neighbourhoods, BH-family sensitivity, redundancy reduction
  06_evidence_integration/    ranks, sensitivity, convergence test, weakest link, weights, loci,
                              maturation anchor, tissue within BioProject, final classes,
                              assembly sensitivity, a4 genes at displaced SNPs
  07_external_validation/     Pereira Lima et al. (2017) lists; Medicago MatNet (Righetti et al., 2015); colocalization with published QTL (optional, curated table); spectral comparison of module graphs between tissues (statGraph); figure of the chromosome 4 locus programs as coexpression graphs
data/raw/         public inputs (not versioned; see data/raw/README.md)
data/processed/   intermediate objects (not versioned; regenerated)
results/, figures/, logs/   outputs (not versioned; regenerated)
docs/             methods draft (Portuguese, LaTeX) and discussion references
run_all.R         runs all stages in order
```

## Reproducibility notes

- Random seeds: 20261101 (convergence permutations and weight draws) and
  20261102 (competitive test of modules), set in `R/config.R`.
- The network was built with BioNERO 1.18.0 and WGCNA 1.74 on R 4.5.3. Other
  versions can change module colours and the borders of small modules.
- GO results depend on the GO.db release (GO.db 3.22.0 was used). With an older
  release the ancestor map changes and so do term counts.
- Soybean IDs: candidates and atlas use Wm82.a4.v1; GSE98199 uses Wm82.a2,
  bridged through the pangene table of the Legume Information System.

## Software

R 4.5.3. Bioconductor: DESeq2, BioNERO, WGCNA, GO.db, AnnotationDbi,
SummarizedExperiment. CRAN: RobustRankAggreg, dplyr, tidyr, readr, readxl,
tibble, stringr, purrr, matrixStats, ggplot2, ggrepel.

## How to cite

Manuscript in preparation. Until then, please cite this repository and the data
sources listed in `data/raw/README.md`.

## License

Code: MIT License (see `LICENSE`). The data belong to their original authors
and are not redistributed here.
