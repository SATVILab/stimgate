#!/usr/bin/env Rscript
# Run analysis/repository integration tests.
#
# These tests cover `scripts/r/` helpers and `analysis/` QMDs.
# They are intentionally kept separate from the package test suite
# (`tests/testthat/`) to avoid inflating `R CMD check` run times and to
# permit sourcing scripts not bundled in the installed package.
#
# Usage (from the repository root):
#   Rscript analysis/tests/run_analysis_tests.R
#
# Or from an R session:
#   source("analysis/tests/run_analysis_tests.R")
#
# Or via testthat directly:
#   testthat::test_dir("analysis/tests/testthat")

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The 'testthat' package is required to run analysis tests.")
}

test_dir <- file.path("analysis", "tests", "testthat")

cat("Running analysis integration tests from:", test_dir, "\n")
testthat::test_dir(test_dir, reporter = testthat::default_reporter())
