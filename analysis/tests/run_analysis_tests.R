#!/usr/bin/env Rscript
# Run analysis/repository integration tests.
#
# These tests cover `scripts/r/` helpers and `analysis/` QMDs.
# They are intentionally kept separate from the package test suite
# (`tests/testthat/`) to avoid inflating `R CMD check` run times and to
# permit sourcing scripts not bundled in the installed package.
#
# The package is loaded from the current checkout via devtools::load_all()
# before running tests, so the analysis suite always tests the current source
# rather than any previously installed version.
#
# Usage (from the repository root):
#   Rscript analysis/tests/run_analysis_tests.R
#
# Or from an R session:
#   source("analysis/tests/run_analysis_tests.R")
#
# Or via testthat directly (after devtools::load_all()):
#   devtools::load_all()
#   testthat::test_dir("analysis/tests/testthat")

for (pkg in c("testthat", "devtools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("The '", pkg, "' package is required to run analysis tests.")
  }
}

cat("Loading stimgate from current checkout...\n")
devtools::load_all(quiet = TRUE)

test_dir <- file.path("analysis", "tests", "testthat")

cat("Running analysis integration tests from:", test_dir, "\n")
testthat::test_dir(test_dir, reporter = testthat::default_reporter())
