.performance_failure_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.performance_failure_fixture <- function() {
  tidyr::expand_grid(iter = 1:20, sample = as.character(1:20)) |>
    dplyr::mutate(scenario = "a", sim_seed = 12345L,
      approach = "stimgate", method = "stimgate_error", error = "gating crashed",
      gateReturnPoint = "stimgate_error", thresholdOrigin = "error",
      thresholdFallbackUsed = NA, threshold = NA_real_,
      propRespTruth = 0.01, propRespEst = NA_real_,
      propStim = NA_real_, propUns = NA_real_,
      nTruePos = NA_integer_, nFalsePos = NA_integer_,
      nFalseNeg = NA_integer_, nTrueNeg = NA_integer_)
}

test_that("all StimGate failures remain explicit in pooled performance coverage", {
  env <- .performance_failure_env()
  raw <- .performance_failure_fixture()
  keys <- c("scenario", "method")
  frequency <- env$.simCompareSummariseFreqBs(raw, keys)
  expect_identical(frequency$method, "stimgate")
  expect_equal(frequency$n, 400L)
  expect_equal(frequency$n_run_error, 400L)
  expect_equal(frequency$n_est, 0L)
  expect_equal(frequency$n_dataset_total, 20L)
  expect_equal(frequency$n_dataset_med_abs_rel_error, 0L)
  expect_true(is.na(frequency$med_abs_rel_error))
  unsigned <- env$.simCompareUnsignedErrorSummary(raw, keys)
  signed <- env$.simCompareSignedErrorSummary(raw, keys)
  expect_equal(unsigned$n_dataset_total, 20L)
  expect_equal(unsigned$n_dataset_median, 0L)
  expect_true(all(is.na(signed$median)))
  expect_true(all(is.na(signed$prop)))
  expect_equal(signed$n_dataset_total, c(20L, 20L))
  classification <- env$.simCompareClassificationSummary(raw, keys)
  expect_equal(classification$n_failed, 400L)
  expect_equal(classification$n_run_error, 400L)
  expect_equal(classification$n_valid, 0L)
  expect_true(is.na(classification$fdp_median))
  # Summaries must not relabel the saved diagnostic/provenance rows in place.
  expect_true(all(raw$method == "stimgate_error"))
  expect_true(all(raw$gateReturnPoint == "stimgate_error"))
})

test_that("all-failure StimGate maximum cohort survives alongside valid competitors", {
  env <- .performance_failure_env()
  failure <- .performance_failure_fixture()
  competitor <- failure |>
    dplyr::mutate(method = "fbeta", approach = "fbeta", error = NA_character_,
      propRespEst = 0.011)
  raw <- dplyr::bind_rows(failure, competitor)
  result <- env$.simCompareDatasetMaxSummary(raw, c("scenario", "method"),
    expected_samples = 20L, expected_datasets = 1:20, mcse = TRUE)
  failed <- result[result$method == "stimgate", ]
  expect_equal(nrow(failed), 2L)
  expect_equal(failed$n_eligible, c(0L, 0L))
  expect_equal(failed$n_incomplete, c(20L, 20L))
  expect_equal(failed$n_affected, c(0L, 0L))
  expect_true(all(is.na(failed$occurrence)))
  expect_true(all(is.na(failed$severity)))
  expect_equal(failed$bootstrap_coverage_occurrence, c(0, 0))
  expect_false(any(failed$interval_available_occurrence))
  expect_true(all(failed$cohort_matches_stimgate))
  expect_false(any(result$cohort_matches_stimgate[result$method == "fbeta"]))
  mixed <- failure
  mixed$method[mixed$iter > 5] <- "stimgate"
  mixed$error[mixed$iter > 5] <- NA_character_
  mixed$propRespEst[mixed$iter > 5] <- 0.01
  cohort <- env$.simCompareDatasetMaxSummary(mixed, c("scenario", "method"),
    expected_samples = 20L, expected_datasets = 1:20)
  expect_equal(cohort$n_eligible, c(15L, 15L))
  expect_equal(cohort$n_incomplete, c(5L, 5L))
  expect_equal(cohort$occurrence, c(0, 0))
})
