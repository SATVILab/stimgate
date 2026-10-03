#!/usr/bin/env Rscript
# Run from the repository root. Examples:
# Rscript analysis/tests/run_qmd_tests.R --list
# Rscript analysis/tests/run_qmd_tests.R 1,3 9
# Rscript analysis/tests/run_qmd_tests.R all
# These targets test each analysis's scientific/API contracts, not full renders.
source(file.path("analysis", "tests", "qmd-test-targets.R"))
targets <- .qmd_test_targets()
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  for (document in names(targets)) {
    cat(document, " -> ", paste(targets[[document]], collapse = ", "), "\n", sep = "")
  }
} else {
  selected <- .select_qmd_test_targets(args, targets)
  for (pkg in c("testthat", "devtools")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The '", pkg, "' package is required to run QMD tests.")
    }
  }
  devtools::load_all(quiet = TRUE)
  testthat::set_max_fails(Inf)
  for (document in names(selected)) {
    cat("Running QMD target:", document, "\n")
    for (test_file in selected[[document]]) {
      testthat::test_file(
        file.path("analysis", "tests", "testthat", test_file),
        env = testthat::test_env(),
        reporter = "summary",
        stop_on_failure = TRUE
      )
    }
  }
}
