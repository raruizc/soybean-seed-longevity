# Note on `statGraph::anogva()` (version 1.0.6)

## What the function does
In `anogva()` (file `R/anogva.R` of github.com/fujita-team/statGraph, master branch, and the CRAN 1.0.6 release), the observed statistic is the mean distance between each group's mean spectral density and the mean of the group means:

```r
distOrig <- distOrig/length(ulabels)
```

The statistic of each permutation is the **sum** of the same distances, without that division:

```r
distBoot[boot] <- distBoot[boot] + distance(sample_group_mean_den[[label]], all_mean_den, dist = dist)
```

The permutation statistic is therefore k times larger than the observed one, for k groups, which pushes the P value towards 1.

## Check by simulation
Erdős–Rényi graphs with 50 vertices, 20 graphs per group, 3 groups; script below.

| Scenario | Original `anogva()` | Same statistic on both sides |
|---|---|---|
| Null (p = 0.50 in all groups), mean P over 30 runs | 0.88 (minimum 0.47) | 0.52 (3% below 0.05) |
| Small effect (0.50 / 0.50 / 0.51), power at 5% over 30 runs | 0 | 0.57 |
| Large effect (0.50 / 0.50 / 0.60), P | 0 | 0 |

The original test is still valid in the sense that it rejects strong effects. It is conservative, though, and has almost no power for moderate ones.

## What stage 7d does
`anogva_fixed()` in `scripts/07_external_validation/04_spectral_network_tissues.R` uses the same statistic for the observed and the permuted labels. Within each permutation it also recomputes the mean of the group means and reports P = (1 + #{permuted ≥ observed}) / (1 + B). The spectral densities come from `statGraph` itself (internal functions `set.list.spectral.density`, `get.mean.spectral.density` and `distance`).

## Reproduce
```r
library(statGraph)
set.seed(1)
mk <- function(p) lapply(1:20, function(i) igraph::sample_gnp(50, p))
lab <- rep(1:3, each = 20)
fix <- statGraph::anogva
body(fix) <- parse(text = sub("pvalue <- length", "distBoot <- distBoot/length(ulabels); pvalue <- length",
                              paste(deparse(body(fix)), collapse = "\n")))[[1]]
environment(fix) <- asNamespace("statGraph")
null_orig <- replicate(30, anogva(c(mk(0.5), mk(0.5), mk(0.5)), lab, maxBoot = 100)$p.value)
null_fix <- replicate(30, fix(c(mk(0.5), mk(0.5), mk(0.5)), lab, maxBoot = 100)$p.value)
c(mean(null_orig), mean(null_fix))
```

It may be worth reporting this to the maintainer (issue on GitHub).
