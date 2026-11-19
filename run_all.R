##############################################################################
# Soybean seed longevity
# Runs the whole study in order, from the repository root
# Author: Ricardo Antonio Ruiz Cardozo
#
# Use: source("run_all.R")
# To run a subset, define STAGES first, e.g. STAGES <- "06_evidence_integration"
# Every script also runs on its own; this file only fixes the order.
##############################################################################

STAGES <- if (exists("STAGES")) STAGES else c("01_candidates", "02_maturation_timecourse", "03_expression_atlas",
                                              "04_coexpression_network", "05_functional_enrichment",
                                              "06_evidence_integration", "07_external_validation")

# Scripts that need an input curated by hand are skipped while it is missing
OPTIONAL <- c("03_qtl_colocalization.R" = file.path("data", "raw", "external", "qtl", "storability_qtl_curated.csv"))
# Scripts that need a package used only by them are skipped while it is not installed
OPTIONAL_PKG <- c("04_spectral_network_tissues.R" = "statGraph")

for (st in STAGES) {
  for (s in sort(list.files(file.path("scripts", st), pattern = "\\.R$", full.names = TRUE))) {
    if (basename(s) %in% names(OPTIONAL) && !file.exists(OPTIONAL[[basename(s)]])) {
      cat("\n>>> skipped (input not found: ", OPTIONAL[[basename(s)]], "): ", s, "\n", sep = "")
      next
    }
    if (basename(s) %in% names(OPTIONAL_PKG) && !requireNamespace(OPTIONAL_PKG[[basename(s)]], quietly = TRUE)) {
      cat("\n>>> skipped (package not installed: ", OPTIONAL_PKG[[basename(s)]], "): ", s, "\n", sep = "")
      next
    }
    cat("\n>>> ", s, "\n", sep = "")
    source(s, echo = FALSE)
  }
}
