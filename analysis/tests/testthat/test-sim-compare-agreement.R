.sim_compare_agreement_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R", "acs_cytof-plot_cyt.R", "sim-compare-freq_bs.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("precision and specificity retain undefined outcomes and tube coverage", {
  env <- .sim_compare_agreement_env()
  raw <- tibble::tibble(
    scenario = "a", method = "stimgate", iter = 1:6,
    nTruePos = c(6L, 0L, 0L, 3L, 0L, 6L),
    nFalsePos = c(2L, 0L, 2L, 0L, 0L, 2L),
    nFalseNeg = c(4L, 5L, 0L, 1L, 0L, 4L),
    nTrueNeg = c(8L, 10L, 8L, 0L, 0L, 8L),
    error = c(rep(NA_character_, 5), "runtime failure")
  )
  metrics <- env$.simCompareClassificationMetrics(raw)
  expect_equal(metrics$precision, c(0.75, NA, 0, 1, NA, NA))
  expect_equal(metrics$specificity, c(0.8, 1, 0.8, NA, NA, NA))
  expect_equal(metrics$precision, 1 - metrics$fdp)
  expect_equal(metrics$specificity, 1 - metrics$false_positive_rate)
  expect_equal(metrics$sensitivity, c(0.6, 0, NA, 0.75, NA, NA))
  expect_equal(metrics$f1, c(12 / 18, 0, 0, 6 / 7, NA, NA))
  expect_equal(metrics$n_selected, c(8, 0, 2, 3, 0, NA))
  expect_equal(metrics$n_genuine_neg, c(10, 10, 10, 0, 0, NA))

  summary <- env$.simCompareClassificationSummary(raw, c("scenario", "method"))
  expect_equal(summary$n, 6L)
  expect_equal(summary$n_valid, 5L)
  expect_equal(summary$n_failed, 1L)
  expect_equal(summary$n_fdp_defined, 3L)
  expect_equal(summary$n_precision_defined, 3L)
  expect_equal(summary$n_specificity_defined, 3L)
  expect_equal(summary$n_sensitivity_defined, 3L)
  expect_equal(summary$n_f1_defined, 4L)
  for (metric in c("precision", "specificity", "sensitivity", "f1")) {
    finite <- metrics[[metric]][is.finite(metrics[[metric]])]
    expect_equal(summary[[paste0(metric, "_median")]], stats::median(finite))
    expect_equal(summary[[paste0(metric, "_q10")]],
      unname(stats::quantile(finite, 0.1)))
  }
  expect_equal(env$.simCompareClassificationOutcomes$precision[["tail"]], "precision_q10")
  expect_equal(env$.simCompareClassificationOutcomes$specificity[["tail"]], "specificity_q10")
})

test_that("frequency correlations reuse population-moment CCC and retain finite pairs", {
  env <- .sim_compare_agreement_env()
  raw <- tibble::tibble(
    scenario = "a", mismatch_val = 0, method = "stimgate",
    propRespTruth = c(1, 2, 3, 4, NA, 5),
    propRespEst = c(3, 5, 7, Inf, 9, 11),
    # The final finite estimate is a runtime error, not a valid pair.
    error = c(rep(NA_character_, 5), "runtime failure")
  )
  identity <- tibble::tibble(
    scenario = "a", mismatch_val = 1, method = "fbeta",
    propRespTruth = 1:3, propRespEst = 1:3, error = NA_character_
  )
  out <- env$.simCompareCorrelationTable(dplyr::bind_rows(raw, identity),
    c("scenario", "mismatch_val", "method"))
  shifted <- out[out$method == "stimgate", ]
  expect_equal(shifted$n, 6L)
  expect_equal(shifted$n_valid, 4L)
  expect_equal(shifted$n_pairs, 3L)
  expect_equal(shifted$pearson, 1)
  # x = 2y + 1: population covariance 4/3, variances 8/3 and 2/3,
  # squared mean difference 9. The sample-moment convention gives 1/3.
  expect_equal(shifted$ccc, 8 / 37)
  expect_equal(shifted$ccc, env$.acsCytofValidationCcc(c(3, 5, 7), 1:3))
  expect_equal(out$ccc[out$method == "fbeta"], 1)
  expect_equal(out$n_pairs[out$method == "fbeta"], 3L)
})

test_that("frequency correlations are undefined for too few pairs or either zero variance", {
  env <- .sim_compare_agreement_env()
  raw <- dplyr::bind_rows(
    tibble::tibble(scenario = "short", propRespEst = c(1, 2, NA), propRespTruth = 1:3),
    tibble::tibble(scenario = "constant_estimate", propRespEst = 2, propRespTruth = 1:3),
    tibble::tibble(scenario = "constant_truth", propRespEst = 1:3, propRespTruth = 2),
    tibble::tibble(scenario = "missing", propRespEst = NA_real_, propRespTruth = 1:3)
  ) |>
    dplyr::mutate(method = "stimgate")
  out <- env$.simCompareCorrelationTable(raw, "scenario")
  expect_true(all(is.na(out$pearson)))
  expect_true(all(is.na(out$ccc)))
  expect_equal(out$n_pairs[match(c("short", "constant_estimate", "constant_truth", "missing"), out$scenario)],
    c(2L, 3L, 3L, 0L))
  # Historical failure labels count in the primary method's cohort.
  raw$method <- "stimgate_error"
  raw$error <- "runtime failure"
  failed <- env$.simCompareCorrelationTable(raw, "scenario")
  expect_true(all(failed$method == "stimgate"))
  expect_true(all(failed$n_pairs == 0L))
  expect_true(all(is.na(failed$ccc)))
})

test_that("both QMDs export classification and agreement tables only with plot results", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  for (qmd in c("7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd")) {
    env <- .sim_compare_agreement_env()
    source(file.path(root, "scripts", "r", "sim-compare-performance-plot.R"), local = env)
    lines <- readLines(file.path(root, "analysis", qmd), warn = FALSE)
    chunk <- function(label) {
      start <- which(lines == paste0("#| label: ", label))
      end <- which(lines == "```" & seq_along(lines) > start)[1L]
      parse(text = lines[seq.int(start + 1L, end - 1L)])
    }
    env$sim_grid <- tibble::tibble(
      base_scenario_id = 1L, scenario_desc = "synthetic", transformation = "gaussian",
      mean_pos_setting = "high", prob_response = 0.1, n_cell = 100L,
      mismatch_type = "mean_shift_all", mismatch_val = 0
    )
    env$compare_tbl <- tidyr::expand_grid(
      env$sim_grid, approach = "a", method = c("stimgate", "fbeta", "tailgate"), iter = 1:3
    ) |>
      dplyr::mutate(nTruePos = 8L, nFalsePos = 2L, nFalseNeg = 2L, nTrueNeg = 88L,
        propRespEst = .data$iter / 10, propRespTruth = .data$iter / 10)
    env$compare_raw <- env$compare_tbl
    env$scenario_cols <- c(names(env$sim_grid), "approach", "method")
    env$show_mcse <- FALSE
    env$mcse_mode <- "off"
    env$fig_key <- "synthetic"
    env$root_dir <- root
    env$run_plots <- TRUE
    env$results_available <- TRUE
    tables <- list()
    paths <- list()
    env$.analysis_report_table <- function(tbl, path_parts, ...) {
      tables[[length(tables) + 1L]] <<- tbl
      paths[[length(paths) + 1L]] <<- path_parts
    }
    eval(chunk("classification-summary"), env)
    eval(chunk("frequency-correlation"), env)
    expect_length(tables, 3L)
    expect_equal(paths[[1]], c("synthetic", "mcse_off", "classification-summary.csv"))
    expect_equal(paths[[3]], c("synthetic", "frequency-correlation.csv"))
    expect_true(all(c("n", "n_valid", "n_precision_defined", "n_specificity_defined",
      "n_sensitivity_defined", "n_f1_defined", "precision_median", "precision_q10",
      "specificity_median", "specificity_q10", "f1_median", "f1_q10") %in% names(tables[[1]])))
    expect_true(all(tables[[1]]$n_precision_defined == 3L))
    expect_true(all(tables[[3]]$n_pairs == 3L))
    expect_true(all(tables[[3]]$ccc == 1))

    tables <- list()
    env$run_plots <- FALSE
    eval(chunk("classification-summary"), env)
    eval(chunk("frequency-correlation"), env)
    expect_length(tables, 0L)
    env$run_plots <- TRUE
    env$results_available <- FALSE
    env$compare_tbl <- tibble::tibble()
    eval(chunk("classification-summary"), env)
    eval(chunk("frequency-correlation"), env)
    expect_length(tables, 0L)
  }
})

test_that("classification mismatch panels include F1 beside FDP and sensitivity without titles", {
  env <- .sim_compare_agreement_env()
  tbl <- tibble::tibble(
    base_scenario_id = 1L, scenario_desc = "synthetic", transformation = "gaussian",
    method = "stimgate", mismatch_val = c(0, 1),
    fdp_median = c(0.1, 0.2), fdp_q90 = c(0.2, 0.3),
    sensitivity_median = c(0.9, 0.8), sensitivity_q10 = c(0.8, 0.7),
    f1_median = c(0.9, 0.8), f1_q10 = c(0.8, 0.7)
  )
  plot <- env$.simComparePlotClassification(tbl, outcomes = c("fdp", "sensitivity", "f1"))
  expect_setequal(as.character(plot$data$outcome),
    c("False discovery proportion", "Sensitivity", "F1 score"))
  expect_null(plot$labels$title)
  expect_null(plot$labels$subtitle)
  layout <- ggplot2::ggplot_build(plot)$layout$layout
  expect_equal(length(unique(layout$ROW)), 1L)
  expect_equal(length(unique(layout$COL)), 3L)

  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  lines <- readLines(file.path(root, "analysis", "8-sim-compare-freq_bs-batch.qmd"), warn = FALSE)
  for (label in c("dataset-paired-differences", "fdp-sensitivity")) {
    start <- which(lines == paste0("#| label: ", label))
    end <- which(lines == "```" & seq_along(lines) > start)[1L]
    expect_true(any(grepl('outcomes = c("fdp", "sensitivity", "f1")',
      lines[seq.int(start, end)], fixed = TRUE)))
  }
})
