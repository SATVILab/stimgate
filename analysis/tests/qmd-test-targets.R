# Scientific/API test targets for the top-level analysis documents.
.qmd_test_targets <- function(root_dir = ".") {
  targets <- list(
    "1-sim-trans.qmd" = "test-sim-trans.R",
    "2a-sim-bw-freq_bs-global.qmd" = c(
      "test-sim-bw-freq_bs-global-simcyto.R", "test-sim-bandwidth-analysis-run.R",
      "test-bandwidth-bias-plots-by-cells.R", "test-qmd-2-cache-guidance.R",
      "test-threshold-density-plots.R", "test-analysis-mcse.R", "test-ratio-companion-plots.R",
      "test-sim-debug-loc.R", "test-sim-debug-compare.R"
    ),
    "2b-sim-bias_uns-freq_bs.qmd" = c(
      "test-sim-bias-uns-freq.R", "test-sim-bandwidth-analysis-run.R",
      "test-bandwidth-bias-plots-by-cells.R", "test-qmd-2-cache-guidance.R",
      "test-analysis-mcse.R", "test-ratio-companion-plots.R"
    ),
    "3-sim-bw-est-base.qmd" = c(
      "test-sim-bw-est-base-simcyto.R", "test-sim-bw-est-base-run.R",
      "test-analysis-mcse.R"
    ),
    "4-sim-bw-est-norm.qmd" = c(
      "test-sim-bw-est-norm-simcyto.R", "test-analysis-mcse.R"
    ),
    "5-sim-bw-est-adaptive.qmd" = "test-sim-bw-est-adaptive-simcyto.R",
    "6-sim-bw-freq_bs-adaptive.qmd" = "test-sim-bw-freq_bs-adaptive-simcyto.R",
    "7-sim-compare-freq_bs.qmd" = c(
      "test-sim-compare-freq_bs-simcyto.R", "test-analysis-7-transactional-multichunk.R",
      "test-analysis-mcse.R", "test-sim-compare-dataset-differences.R", "test-ratio-companion-plots.R",
      "test-sim-performance-estimands.R", "test-sim-performance-figures.R",
      "test-analysis-mcse-modes.R", "test-sim-performance-failures.R",
      "test-sim-compare-recorded-errors.R",
      "test-analysis-sim-size.R"
    ),
    "8-sim-compare-freq_bs-batch.qmd" = c(
      "test-sim-compare-freq_bs-batch.R", "test-sim-compare-classification.R",
      "test-analysis-mcse.R", "test-sim-compare-dataset-differences.R", "test-ratio-companion-plots.R",
      "test-sim-performance-estimands.R", "test-sim-performance-figures.R",
      "test-analysis-mcse-modes.R", "test-sim-performance-failures.R",
      "test-sim-compare-recorded-errors.R",
      "test-analysis-sim-size.R"
    ),
    "9-real-compare-acs-cytof.qmd" = c(
      "test-acs-cytof-gate.R", "test-acs-cytof-methods.R", "test-ratio-companion-plots.R",
      "test-acs-cytof-paths.R"
    ),
    "10-real-compare-acs-cytof-validation.qmd" = c(
      "test-acs-cytof-validation-plots.R", "test-acs-cytof-paths.R"
    ),
    "2c-sim-test.qmd" = c(
      "test-sim-test-qmd.R", "test-sim-debug-loc.R", "test-sim-debug-compare.R"
    )
  )
  documents <- list.files(file.path(root_dir, "analysis"), pattern = "[.]qmd$")
  if (!setequal(documents, names(targets))) {
    stop("QMD test registry must match every top-level analysis QMD. Missing: ",
      paste(setdiff(names(targets), documents), collapse = ", "),
      "; unregistered: ", paste(setdiff(documents, names(targets)), collapse = ", ")
    )
  }
  tests <- unique(unlist(targets, use.names = FALSE))
  missing <- tests[!file.exists(file.path(root_dir, "analysis", "tests", "testthat", tests))]
  if (length(missing)) {
    stop("Missing QMD test files: ", paste(missing, collapse = ", "))
  }
  targets
}

.select_qmd_test_targets <- function(args, targets) {
  if (!length(args)) {
    return(targets)
  }
  if (any(!nzchar(trimws(args)))) {
    stop("QMD selection must not be empty.")
  }
  selection <- unlist(strsplit(trimws(args), "[,[:space:]]+"), use.names = FALSE)
  selection <- selection[nzchar(selection)]
  if (!length(selection)) {
    stop("QMD selection must not be empty.")
  }
  if (identical(selection, "all")) {
    return(targets)
  }
  filenames <- names(targets)
  stems <- sub("[.]qmd$", "", filenames)
  numbers <- sub("-.*", "", stems)
  selection <- sub("^analysis/", "", gsub("\\\\", "/", selection))
  indices <- vapply(selection, function(value) {
    match(value, filenames, nomatch = match(value, stems, nomatch = match(value, numbers)))
  }, integer(1))
  if (anyNA(indices)) {
    stop("Unknown QMD selection: ", paste(selection[is.na(indices)], collapse = ", "),
      ". Use --list to see available targets."
    )
  }
  targets[unique(indices)]
}
