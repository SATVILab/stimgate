.recorded_error_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R")) {
    source(file.path(root, "scripts/r", file), local = env)
  }
  env
}

.recorded_error_fixture <- function(n_sample = 20L, n_iter = 20L) {
  tidyr::expand_grid(iter = seq_len(n_iter), sample = as.character(seq_len(n_sample)),
    method = c("stimgate", "fbeta", "tailgate")) |>
    dplyr::mutate(sim_id = 1L, sim_seed = 100L, propRespTruth = 0.1,
      propRespEst = 0.1, nCellStim = 100L, nPosStim = 12L,
      nTruePos = 10L, nFalsePos = 2L, nFalseNeg = 1L, nTrueNeg = 87L,
      threshold = 1, thresholdMetric = 0.5, thresholdOrigin = "calculated",
      gateReturnPoint = paste0(method, "_calculated"), thresholdFallbackUsed = FALSE,
      propStim = 0.12, propUns = 0.02, nPosUns = 2L,
      unsExprSum = 1.5, error = NA_character_)
}

.recorded_error_mark <- function(raw, method = "tailgate", iter = 2L, sample = "10") {
  selected <- raw$method == method & raw$iter == iter & raw$sample == sample
  raw$error[selected] <- "Tailgate bandwidth must be a finite positive scalar."
  raw$thresholdOrigin[selected] <- paste0("error: ", raw$error[selected])
  raw$gateReturnPoint[selected] <- paste0(method, "_error")
  for (col in c("threshold", "thresholdMetric", "propRespEst", "propStim", "propUns",
    "nPosStim", "nPosUns", "nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg")) {
    raw[[col]][selected] <- NA
  }
  raw
}

test_that("complete simulations retain explicitly recorded comparator failures", {
  env <- .recorded_error_env()
  raw <- .recorded_error_mark(.recorded_error_fixture())
  expect_true(env$.simComparePrimaryOutputComplete(raw, 20L, 20L))
  status <- env$.simCompareGridOutputStatus(raw, data.frame(sim_id = 1L), 20L, 20L)
  expect_true(status$validation_ok)
  expect_equal(status$completed_ids, 1L)
  expect_length(status$failed_ids, 0L)
  for (retry in c(FALSE, TRUE)) {
    expect_true(env$.simCompareValidateScenarioCache(raw, data.frame(sim_id = 1L),
      nSample = 20L, nIter = 20L, retryErrors = retry))
  }
  # The failure remains unscored, rather than becoming a zero-response estimate.
  summary <- env$.simCompareSummariseFreqBs(raw, "method")
  failed <- summary[summary$method == "tailgate", ]
  expect_equal(failed$n_run_error, 1L)
  expect_equal(failed$n_est, 399L)
  expect_true(all(is.na(raw$propRespEst[env$.simCompareRecordedComparatorErrors(raw)])))
})

test_that("recorded errors cannot conceal absent, unlabelled or malformed outcomes", {
  env <- .recorded_error_env()
  raw <- .recorded_error_mark(.recorded_error_fixture())
  failed <- which(env$.simCompareRecordedComparatorErrors(raw))
  cases <- list(missing_row = raw[-1L, ], duplicate_row = dplyr::bind_rows(raw, raw[1L, ]))
  for (field in c("error", "gateReturnPoint", "thresholdOrigin")) {
    bad <- raw
    bad[[field]][failed] <- NA_character_
    cases[[field]] <- bad
  }
  for (field in c("threshold", "propRespEst", "nPosStim", "nTruePos")) {
    bad <- raw
    bad[[field]][failed] <- 0
    cases[[field]] <- bad
  }
  for (field in c("propRespTruth", "unsExprSum", "nCellStim")) {
    bad <- raw
    bad[[field]][failed] <- NA
    cases[[field]] <- bad
  }
  for (field in c("propRespTruth", "unsExprSum", "nCellStim")) {
    bad <- raw
    bad[[field]][failed] <- bad[[field]][failed] + 0.1
    cases[[paste0(field, "_differs")]] <- bad
  }
  cases$missing_method_column <- dplyr::select(raw, -"method")
  bad <- raw
  bad$thresholdFallbackUsed[failed] <- TRUE
  cases$fallback_disguised_as_error <- bad
  bad <- raw
  bad$sample[bad$sample == "20"] <- "21"
  cases$missing_intended_sample <- bad
  bad <- raw
  bad$iter[bad$iter == 20L] <- 21L
  cases$missing_intended_iteration <- bad
  bad <- raw
  bad$propRespEst[1L] <- NA_real_
  cases$unlabelled_estimate <- bad
  bad <- raw
  bad$nPosStim[1L] <- 99L
  cases$bad_gate_counts <- bad
  bad <- raw
  other <- which(bad$method == "fbeta")[1L]
  bad$nTruePos[other] <- 9L
  bad$nFalsePos[other] <- 3L
  cases$mismatched_biological_counts <- bad
  bad <- raw
  bad$error[1L] <- "StimGate failed"
  cases$stimgate_failure <- bad
  bad <- raw[1L, ]
  bad$method <- NA_character_
  bad$iter <- NA_integer_
  bad$sample <- NA_character_
  bad$error <- "disk write failed"
  cases$infrastructure_failure <- bad
  for (name in names(cases)) {
    expect_false(env$.simComparePrimaryOutputComplete(cases[[name]], 20L, 20L), info = name)
    status <- env$.simCompareGridOutputStatus(cases[[name]], data.frame(sim_id = 1L), 20L, 20L)
    expect_false(status$validation_ok, info = name)
    expect_true(nzchar(status$failure_reasons[["1"]]), info = name)
    expect_false(env$.simCompareValidateScenarioCache(cases[[name]], data.frame(sim_id = 1L),
      nSample = 20L, nIter = 20L, retryErrors = TRUE), info = name)
  }
})

test_that("mismatch checks preserve pairing and finite zero-mismatch coverage", {
  env <- .recorded_error_env()
  raw <- .recorded_error_fixture(20L, 2L) |>
    dplyr::mutate(base_scenario_id = 1L, mismatch_type = "mean_shift_all", mismatch_val = 0)
  shifted <- dplyr::mutate(raw, sim_id = 2L, mismatch_type = "mean_shift_negative")
  shifted <- .recorded_error_mark(shifted)
  combined <- dplyr::bind_rows(raw, shifted)
  out <- env$.simCompareValidateMismatch(combined)
  expect_true(all(out$pairing_check$paired))
  failed <- out$zero_agreement[out$zero_agreement$method == "tailgate", ]
  expect_equal(failed$n_pairs, 40L)
  expect_equal(failed$n_failed_pairs, 1L)
  expect_equal(failed$n_compared, 39L)
  expect_equal(failed$n_same_threshold, 39L)
  expect_equal(failed$n_same_counts, 39L)
  mismatch <- combined
  index <- which(mismatch$sim_id == 2L & mismatch$method == "fbeta")[1L]
  mismatch$threshold[index] <- 2
  expect_error(env$.simCompareValidateMismatch(mismatch), "did not reproduce")
  unpaired <- combined
  unpaired$unsExprSum[index] <- 2
  expect_error(env$.simCompareValidateMismatch(unpaired), "not paired")
})
