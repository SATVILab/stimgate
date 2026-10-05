root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

.classification_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

test_that("confusion counts follow the package's strict gate convention", {
  env <- .classification_env()
  x <- c(0.1, 0.5, 1, 2, 3)
  labels <- c("gn", "gn", "gp", "gp", "gp")

  # Perfect gate between the components.
  perfect <- env$.simCompareConfusionCounts(x, labels, 0.75)
  expect_equal(
    unlist(perfect),
    c(nTruePos = 3, nFalsePos = 0, nFalseNeg = 0, nTrueNeg = 2)
  )

  # A cell exactly at the gate is negative (`x > gate`, as in R/pos_ind.R).
  at_value <- env$.simCompareConfusionCounts(x, labels, 1)
  expect_equal(
    unlist(at_value),
    c(nTruePos = 2, nFalsePos = 0, nFalseNeg = 1, nTrueNeg = 2)
  )
  at_negative <- env$.simCompareConfusionCounts(x, labels, 0.1)
  expect_equal(
    unlist(at_negative),
    c(nTruePos = 3, nFalsePos = 1, nFalseNeg = 0, nTrueNeg = 1)
  )

  # Empty gate above every cell.
  empty <- env$.simCompareConfusionCounts(x, labels, 10)
  expect_equal(
    unlist(empty),
    c(nTruePos = 0, nFalsePos = 0, nFalseNeg = 3, nTrueNeg = 2)
  )

  # No genuine positives.
  none <- env$.simCompareConfusionCounts(x, rep("gn", 5), 0.75)
  expect_equal(
    unlist(none),
    c(nTruePos = 0, nFalsePos = 3, nFalseNeg = 0, nTrueNeg = 2)
  )

  # No finite gate: counts are unknown, not zero.
  expect_true(all(is.na(unlist(env$.simCompareConfusionCounts(x, labels, NA)))))
  expect_error(env$.simCompareConfusionCounts(x, labels[-1], 1), "same length")
})

test_that("counts reproduce the gated stimulated count of the estimate helper", {
  env <- .classification_env()
  x_stim <- c(0.1, 0.5, 1, 2, 3)
  est <- env$.simCompareEstimateFromThreshold(
    xStim = x_stim, xUns = c(0, 0.2, 1.5), threshold = 1,
    labelsStim = c("gn", "gn", "gp", "gp", "gp")
  )
  expect_equal(est$nTruePos + est$nFalsePos, est$nPosStim)
  expect_equal(
    est$nTruePos + est$nFalsePos + est$nFalseNeg + est$nTrueNeg,
    est$nCellStim
  )
  # Without labels the estimate is unchanged and has no counts.
  plain <- env$.simCompareEstimateFromThreshold(
    xStim = x_stim, xUns = c(0, 0.2, 1.5), threshold = 1
  )
  expect_false("nTruePos" %in% names(plain))
  expect_equal(plain$propRespEst, est$propRespEst)
})

test_that("classification metrics handle empty, perfect and failed gates", {
  env <- .classification_env()
  tbl <- tibble::tibble(
    case = c("perfect", "empty", "no_pos", "mixed", "fallback_empty", "failed"),
    nTruePos = c(3L, 0L, 0L, 2L, 0L, NA),
    nFalsePos = c(0L, 0L, 3L, 1L, 0L, NA),
    nFalseNeg = c(0L, 3L, 0L, 1L, 3L, NA),
    nTrueNeg = c(2L, 2L, 2L, 1L, 2L, NA),
    thresholdFallbackUsed = c(FALSE, FALSE, FALSE, FALSE, TRUE, NA),
    error = c(NA, NA, NA, NA, NA, "failed")
  )
  out <- env$.simCompareClassificationMetrics(tbl)

  perfect <- out[out$case == "perfect", ]
  expect_equal(perfect$fdp, 0)
  expect_equal(perfect$sensitivity, 1)
  expect_equal(perfect$false_positive_rate, 0)
  expect_equal(perfect$selected_fraction, 3 / 5)

  # Empty gate: FDP undefined, zero sensitivity when positives exist.
  empty <- out[out$case == "empty", ]
  expect_true(is.na(empty$fdp))
  expect_equal(empty$sensitivity, 0)
  expect_true(empty$gate_empty)
  expect_equal(empty$gate_status, "calculated_empty")

  # No genuine positives: sensitivity undefined.
  no_pos <- out[out$case == "no_pos", ]
  expect_true(is.na(no_pos$sensitivity))
  expect_equal(no_pos$fdp, 1)

  mixed <- out[out$case == "mixed", ]
  expect_equal(mixed$fdp, 1 / 3)
  expect_equal(mixed$sensitivity, 2 / 3)
  expect_equal(mixed$false_positive_rate, 1 / 2)

  expect_equal(
    out$gate_status[out$case == "fallback_empty"], "fallback_empty"
  )
  expect_equal(out$gate_status[out$case == "failed"], "failed")
  expect_true(is.na(out$fdp[out$case == "failed"]))

  props <- unlist(out[, c(
    "fdp", "sensitivity", "false_positive_rate", "selected_fraction"
  )])
  props <- props[!is.na(props)]
  expect_true(all(props >= 0 & props <= 1))
})

test_that("classification summary uses replicate values and excludes undefined FDP", {
  env <- .classification_env()
  raw <- tibble::tibble(
    scenario = "a",
    method = "stimgate",
    nTruePos = c(10L, 5L, 0L),
    nFalsePos = c(0L, 5L, 0L),
    nFalseNeg = c(0L, 5L, 10L),
    nTrueNeg = c(90L, 85L, 90L),
    thresholdFallbackUsed = c(FALSE, FALSE, TRUE),
    error = NA_character_
  )
  out <- env$.simCompareClassificationSummary(raw, scenarioCols = c("scenario", "method"))
  expect_equal(out$n, 3L)
  expect_equal(out$n_valid, 3L)
  expect_equal(out$n_fdp_defined, 2L)
  expect_equal(out$n_empty, 1L)
  expect_equal(out$n_fallback, 1L)
  expect_equal(out$n_fallback_empty, 1L)
  # FDP over the two defined replicates only (0 and 0.5).
  expect_equal(out$fdp_median, 0.25)
  expect_equal(out$fdp_q90, stats::quantile(c(0, 0.5), 0.9, names = FALSE))
  # Sensitivity includes the empty gate's zero.
  expect_equal(out$sensitivity_median, 0.5)
  expect_equal(out$sensitivity_q10, stats::quantile(c(1, 0.5, 0), 0.1, names = FALSE))
})

test_that("pairing and zero-mismatch checks flag unpaired data", {
  env <- .classification_env()
  base <- tidyr::expand_grid(
    base_scenario_id = 1L,
    mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = c(0, 0.1),
    iter = 1:2,
    sample = "1",
    method = "stimgate"
  ) |>
    dplyr::mutate(
      sim_id = as.integer(factor(paste(mismatch_type, mismatch_val))),
      unsExprSum = 10,
      nCellStim = 100,
      threshold = 1,
      nTruePos = 10L, nFalsePos = 0L, nFalseNeg = 0L, nTrueNeg = 90L
    )
  expect_true(all(env$.simComparePairingCheck(base)$paired))
  agreement <- env$.simCompareZeroMismatchAgreement(base)
  expect_equal(agreement$n_compared, 2L)
  expect_equal(agreement$n_same_counts, 2L)

  unpaired <- base
  unpaired$unsExprSum[unpaired$iter == 2L & unpaired$mismatch_val == 0.1] <- 11
  check <- env$.simComparePairingCheck(unpaired)
  expect_equal(check$paired, c(TRUE, FALSE))

  differ <- base
  differ$threshold[differ$mismatch_type == "mean_shift_negative" &
    differ$mismatch_val == 0 & differ$iter == 1L] <- 1.1
  agreement <- env$.simCompareZeroMismatchAgreement(differ)
  expect_equal(agreement$n_same_threshold, 1L)
})

test_that("mismatch settings are paired across replicates and agree at zero", {
  skip_if_not_installed("simcyto")
  withr::local_preserve_seed()
  env <- .classification_env()
  # The package is already loaded from this checkout; do not reload it.
  env$.simCompareEnsureCurrentCheckout <- function() invisible(NULL)
  path_fbeta <- file.path(root_dir, "scripts", "python", "fbeta.py")
  base <- tibble::tibble(
    base_scenario_id = 1L, transformation = "gaussian", mean_pos = 5,
    prob_response = 0.1, n_cell = 400, bias_uns = 0, bw = 0.1,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.1,
    n_cell_uns_relative_to_stim = 1, stim_sd_multiplier = 1, sim_seed = 101L
  )
  rows <- dplyr::bind_rows(
    dplyr::mutate(base, sim_id = 1L, mismatch_type = "mean_shift_all",
      mismatch_val = 0, stim_mean_shift = 0, stim_mean_shift_clusters = NA_character_),
    dplyr::mutate(base, sim_id = 2L, mismatch_type = "mean_shift_negative",
      mismatch_val = 0, stim_mean_shift = 0, stim_mean_shift_clusters = "gn"),
    dplyr::mutate(base, sim_id = 3L, mismatch_type = "mean_shift_all",
      mismatch_val = 1, stim_mean_shift = 1, stim_mean_shift_clusters = NA_character_)
  )
  runs <- lapply(seq_len(nrow(rows)), function(i) {
    env$.simCompareRunScenario(
      row = rows[i, , drop = FALSE], nSample = 1L, nIter = 2L,
      resume = FALSE, tailgateAutoTol = TRUE, pathFbeta = path_fbeta,
      keepCells = TRUE
    )
  })
  cells <- lapply(runs, attr, "cells")
  out <- dplyr::bind_rows(lapply(runs, function(r) {
    attr(r, "cells") <- NULL
    r
  }))
  expect_true(all(is.na(out$error)), info = paste(unique(out$error), collapse = "; "))

  # The second replicate's unstimulated tube and labels are identical in every
  # setting, so the settings are paired beyond the first replicate.
  for (it in 1:2) {
    uns <- lapply(cells, function(d) d$expr[d$iter == it & d$condition == "unstim"])
    lab <- lapply(cells, function(d) d$label[d$iter == it & d$condition == "stim"])
    expect_identical(uns[[1]], uns[[3]])
    expect_identical(lab[[1]], lab[[3]])
    stim0 <- cells[[1]]$expr[cells[[1]]$iter == it & cells[[1]]$condition == "stim"]
    stim1 <- cells[[3]]$expr[cells[[3]]$iter == it & cells[[3]]$condition == "stim"]
    expect_equal(stim1 - stim0, rep(1, length(stim0)), tolerance = 1e-12)
  }
  expect_true(all(env$.simComparePairingCheck(out)$paired))

  primary <- out[out$method %in% c("stimgate", "fbeta", "tailgate"), ]
  expect_true(env$.simCompareCountsConsistent(primary))

  # Zero-shift variants give identical data and identical results.
  expect_identical(cells[[1]], cells[[2]])
  agreement <- env$.simCompareZeroMismatchAgreement(out)
  expect_equal(agreement$n_same_threshold, agreement$n_compared)
  expect_equal(agreement$n_same_counts, agreement$n_compared)
})

test_that("classification plot fixes 0-100% scales and frees the x range", {
  env <- .classification_env()
  tbl <- tidyr::expand_grid(
    base_scenario_id = 1:2,
    method = c("stimgate", "fbeta"),
    mismatch_val = c(0, 0.05, 0.1)
  ) |>
    dplyr::mutate(
      transformation = ifelse(base_scenario_id == 1L, "gaussian", "gamma"),
      scenario_desc = ifelse(base_scenario_id == 1L, "Clean", "Very clean"),
      mismatch_val = ifelse(base_scenario_id == 1L, mismatch_val * 10, mismatch_val),
      fdp_median = 0.1, fdp_q90 = 0.3,
      sensitivity_median = 0.9, sensitivity_q10 = 0.5,
      fpr_median = 0.001, fpr_q90 = 0.002
    )
  p <- env$.simComparePlotClassification(tbl, x_label = "Shift")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$layout$panel_scales_y[[1]]$limits, c(0, 1))
  expect_setequal(
    levels(built$layout$layout$outcome),
    c("False discovery proportion", "Sensitivity")
  )
  # Gaussian first, then gamma, each column with its own range.
  expect_equal(
    levels(built$layout$layout$scenario),
    c("Gaussian: Clean", "Gamma: Very clean")
  )
  x_ranges <- lapply(built$layout$panel_scales_x, function(s) s$dimension())
  expect_false(isTRUE(all.equal(x_ranges[[1]], x_ranges[[2]])))

  fpr <- env$.simComparePlotClassification(tbl, outcomes = "fpr", unit_scale = FALSE)
  expect_s3_class(fpr, "ggplot")
  expect_no_error(ggplot2::ggplot_build(fpr))
})

test_that("gate diagnostic plot draws every method's gate per setting", {
  env <- .classification_env()
  cells <- tidyr::expand_grid(
    mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = c(0, 0.1),
    condition = c("unstim", "stim"),
    i = 1:20
  ) |>
    dplyr::mutate(
      expr = i / 10 + mismatch_val,
      label = ifelse(i > 15, "gp", "gn")
    )
  gates <- tidyr::expand_grid(
    mismatch_type = c("mean_shift_all", "mean_shift_negative"),
    mismatch_val = c(0, 0.1),
    method = c("stimgate", "fbeta", "tailgate")
  ) |>
    dplyr::mutate(threshold = 1.5)
  p <- env$.simComparePlotGateDiagnostic(cells, gates)
  built <- ggplot2::ggplot_build(p)
  vline <- built$data[[3]]
  # Coincident gates stay separate lines.
  expect_equal(nrow(vline), nrow(gates))
  expect_equal(nrow(built$layout$layout), 4L)
})

test_that("analysis 8 leads with FDP and sensitivity and keeps error results", {
  content <- paste(readLines(
    file.path(root_dir, "analysis", "8-sim-compare-freq_bs-batch.qmd"),
    warn = FALSE
  ), collapse = "\n")
  expect_true(grepl("## Gate purity and detection", content, fixed = TRUE))
  expect_true(grepl(
    "## Secondary results: background-subtracted frequency error",
    content, fixed = TRUE
  ))
  expect_lt(
    regexpr("## Gate purity and detection", content, fixed = TRUE),
    regexpr("## Secondary results", content, fixed = TRUE)
  )
  expect_true(grepl("validate_full = .simCompareValidateMismatch", content, fixed = TRUE))
  expect_true(grepl("mismatch_validation <- .simCompareValidateMismatch(compare_raw)", content, fixed = TRUE))
  expect_true(grepl("gate_diagnostic_spec = gate_diagnostic_spec", content, fixed = TRUE))
  expect_true(grepl('c("collated", "gate_diagnostic.rds")', content, fixed = TRUE))
  # The diagnostic chunk reads saved results and never simulates at plot time.
  diag_chunk <- regmatches(
    content,
    regexpr("(?s)#\\| label: gate-diagnostic-data.*?```", content, perl = TRUE)
  )
  expect_length(diag_chunk, 1L)
  expect_false(grepl(".simCompareGateDiagnosticRun", diag_chunk, fixed = TRUE))
})
