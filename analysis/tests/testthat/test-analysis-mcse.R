root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

.mcse_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  env
}

.has_errorbar <- function(p) {
  any(vapply(p$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1)))
}

.errorbar_data <- function(p) {
  idx <- which(vapply(
    p$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1)
  ))
  ggplot2::layer_data(p, idx[[1]])
}

test_that("mean MCSE is sd / sqrt(n), with the n - 1 standard deviation", {
  env <- .mcse_env()
  x <- c(1, 4, 2, 8, 5, 7, NA, Inf)
  res <- env$.analysis_mcse_mean(x)
  xf <- c(1, 4, 2, 8, 5, 7)
  expect_equal(res$estimate, mean(xf))
  expect_equal(res$mcse, stats::sd(xf) / sqrt(6))
  expect_equal(res$n, 6L)
  expect_equal(res$lower, mean(xf) - stats::qnorm(0.975) * res$mcse)
  expect_equal(res$upper, mean(xf) + stats::qnorm(0.975) * res$mcse)
  # Shifting and scaling the values scales the MCSE only.
  shifted <- env$.analysis_mcse_mean(10 + 3 * xf)
  expect_equal(shifted$mcse, 3 * res$mcse)
  # Too few values for a standard error.
  expect_true(is.na(env$.analysis_mcse_mean(2)$mcse))

  prop <- env$.analysis_mcse_proportion(c(TRUE, FALSE, TRUE, TRUE, NA))
  expect_equal(prop$estimate, 0.75)
  expect_equal(prop$mcse, sqrt(0.75 * 0.25 / 4))
  expect_lte(prop$upper, 1)
})

test_that("percentile intervals use binomial order statistics", {
  env <- .mcse_env()
  idx <- env$.analysis_mcse_quantile_index(100, 0.5)
  expect_equal(
    unname(idx),
    c(stats::qbinom(0.025, 100, 0.5), stats::qbinom(0.975, 100, 0.5) + 1)
  )
  expect_equal(unname(idx), c(40L, 61L))

  # With the values 1..100 (in any order), x_(i) = i.
  x <- c(37:100, 1:36)
  res <- env$.analysis_mcse_quantile(x, 0.5)
  expect_equal(res$estimate, stats::median(x))
  expect_equal(c(res$lower, res$upper), c(40, 61))
  expect_equal(res$mcse, (61 - 40) / (2 * stats::qnorm(0.975)))
  # No randomness: the same interval every time.
  expect_identical(res, env$.analysis_mcse_quantile(x, 0.5))

  # Coverage is at least 95% whenever an interval is returned.
  for (n in c(10, 25, 50, 100, 200)) {
    for (p in c(0.1, 0.5, 0.9, 0.95)) {
      i <- env$.analysis_mcse_quantile_index(n, p)
      if (anyNA(i)) next
      coverage <- stats::pbinom(i[["upper"]] - 1, n, p) -
        stats::pbinom(i[["lower"]] - 1, n, p)
      expect_gte(coverage, 0.95)
    }
  }
})

test_that("tiny samples, out-of-sample bounds and maxima give no interval", {
  env <- .mcse_env()
  expect_true(all(is.na(env$.analysis_mcse_quantile_index(4, 0.5))))
  tiny <- env$.analysis_mcse_quantile(c(1, 2, 3, 4), 0.5)
  expect_true(is.na(tiny$lower) && is.na(tiny$upper) && is.na(tiny$mcse))
  expect_equal(tiny$estimate, 2.5)
  # The 90th percentile of 25 values would need x_(26).
  expect_true(all(is.na(env$.analysis_mcse_quantile_index(25, 0.9))))
  mx <- env$.analysis_mcse_max(c(1, 5, 3))
  expect_equal(mx$estimate, 5)
  expect_true(is.na(mx$mcse) && is.na(mx$lower) && is.na(mx$upper))
  empty <- env$.analysis_mcse_quantile(numeric(), 0.5)
  expect_true(is.na(empty$estimate))
})

test_that("scenario averages combine MCSEs as independent scenarios", {
  env <- .mcse_env()
  res <- env$.analysis_mcse_average(c(1, 2, 3), c(0.1, 0.2, 0.2))
  expect_equal(res$estimate, 2)
  expect_equal(res$mcse, sqrt(0.01 + 0.04 + 0.04) / 3)
  expect_equal(res$lower, 2 - stats::qnorm(0.975) * res$mcse)
  # Missing estimates are left out; a missing MCSE gives no interval.
  left_out <- env$.analysis_mcse_average(c(1, NA, 3), c(0.1, NA, 0.1))
  expect_equal(left_out$estimate, 2)
  expect_equal(left_out$mcse, sqrt(0.02) / 2)
  expect_true(is.na(env$.analysis_mcse_average(c(1, 3), c(0.1, NA))$mcse))

  tbl <- tibble::tibble(
    g = c("a", "a", "b", "b"),
    stat = c(1, 3, 2, 4),
    stat_mcse = c(0.3, 0.4, 0.1, NA)
  )
  avg <- env$.analysis_mcse_average_cols(tbl, "g", "stat")
  expect_equal(avg$stat, c(2, 3))
  expect_equal(avg$stat_mcse, c(0.25, NA))
  expect_equal(avg$stat_upper[[1]], 2 + stats::qnorm(0.975) * 0.25)
})

test_that("dataset-level MCSE uses the spread of per-dataset statistics", {
  env <- .mcse_env()
  set.seed(1)
  x <- stats::rnorm(60)
  unit <- rep(1:6, each = 10)
  per_unit <- vapply(split(x, unit), stats::median, numeric(1))
  expect_equal(
    env$.analysis_mcse_between_units(x, unit, stats::median),
    stats::sd(per_unit) / sqrt(6)
  )
  # Fewer than five datasets: no interval.
  keep <- unit <= 4
  expect_true(is.na(
    env$.analysis_mcse_between_units(x[keep], unit[keep], stats::median)
  ))
  # Rows without a dataset are ignored.
  expect_true(is.na(env$.analysis_mcse_between_units(
    x, c(unit[1:40], rep(NA, 20)), stats::median
  )))
})

test_that("signed summaries add per-direction intervals on their own side", {
  env <- .mcse_env()
  rel_error <- c(seq(0.01, 1, length.out = 60), -seq(0.01, 0.9, length.out = 40))
  plain <- env$.simBandwidthSignedErrorSides(rel_error)
  expect_named(plain, c("direction", "prop", "median", "q90", "q95", "max"))

  sides <- env$.simBandwidthSignedErrorSides(rel_error, mcse = TRUE)
  expect_equal(sides[, names(plain)], plain)
  over <- sort(rel_error[rel_error > 0])
  i <- env$.analysis_mcse_quantile_index(length(over), 0.5)
  expect_equal(sides$median_lower[[1]], over[[i[["lower"]]]])
  expect_equal(sides$median_upper[[1]], over[[i[["upper"]]]])
  under <- sides[sides$direction == "under", ]
  expect_lte(under$median_lower, under$median)
  expect_lte(under$median, under$median_upper)
  expect_lte(under$median_upper, 0)
  expect_true(all(is.na(sides$max_lower)))

  # Per dataset: the MCSE is the spread of the datasets' own medians.
  unit <- rep(1:5, length.out = length(rel_error))
  by_unit <- env$.simBandwidthSignedErrorSides(rel_error, mcse = TRUE, unit = unit)
  expect_equal(by_unit$median, plain$median)
  pos <- rel_error > 0
  expect_equal(
    by_unit$median_mcse[[1]],
    env$.analysis_mcse_between_units(rel_error[pos], unit[pos], stats::median)
  )
  expect_gte(by_unit$median_lower[[1]], 0)

  summ <- env$.simBandwidthSignedErrorSummary(
    tibble::tibble(g = rep(c("a", "b"), each = 100), rel_error = c(rel_error, rel_error)),
    "g",
    mcse = TRUE
  )
  avg <- env$.simBandwidthSignedErrorAverage(summ, character())
  expect_equal(nrow(avg), 2L)
  expect_equal(
    avg$median_mcse[avg$direction == "over"],
    sqrt(2 * sides$median_mcse[[1]]^2) / 2
  )
  expect_lte(avg$median_upper[avg$direction == "under"], 0)
})

test_that("plots draw interval layers only when the toggle is on", {
  env <- .mcse_env()
  set.seed(2)
  raw <- tidyr::expand_grid(
    transformation = "gaussian", bw = c(0.1, 0.2), i = 1:200
  ) |>
    dplyr::mutate(rel_error = stats::rnorm(dplyr::n(), 0.2, 1))
  sides <- env$.simBandwidthSignedErrorSummary(
    raw, c("transformation", "bw"),
    mcse = TRUE
  )
  off <- env$.simBandwidthGlobalSignedErrorPlot(sides)
  on <- env$.simBandwidthGlobalSignedErrorPlot(sides, mcse = TRUE)
  expect_false(.has_errorbar(off))
  expect_true(.has_errorbar(on))
  expect_no_error(ggplot2::ggplotGrob(on))

  # Bars never pass the +1500% (16x) cap, drawn at log2(16) on the scale.
  capped <- sides |>
    dplyr::mutate(median_upper = ifelse(direction == "over", 40, median_upper))
  bars <- .errorbar_data(
    env$.simBandwidthGlobalSignedErrorPlot(capped, mcse = TRUE)
  )
  expect_lte(max(bars$ymax), log2(16) + 1e-8)

  bias_tbl <- tidyr::expand_grid(
    bias_uns_multiplier = c(0.5, 1), bw = 0.1, mismatch_label = "none"
  ) |>
    dplyr::mutate(
      bias_uns_basis = "bandwidth",
      median_abs_rel_error = 0.2, median_abs_rel_error_lower = 0.1,
      median_abs_rel_error_upper = 0.3,
      q90_abs_rel_error = 0.5, max_abs_rel_error = 1
    )
  expect_false(.has_errorbar(env$.simBandwidthBiasRelativeErrorPlot(bias_tbl)))
  bias_on <- env$.simBandwidthBiasRelativeErrorPlot(bias_tbl, mcse = TRUE)
  expect_true(.has_errorbar(bias_on))
  expect_no_error(ggplot2::ggplotGrob(bias_on))

  cell_tbl <- tidyr::expand_grid(
    n_cell = c(1e3, 1e4), method = c("stimgate", "fbeta"),
    transformation = "gaussian", prob_response = 0.01
  ) |>
    dplyr::mutate(
      mean_abs_error = 0.01, mean_abs_error_lower = 0.005,
      mean_abs_error_upper = 0.015
    )
  expect_false(.has_errorbar(env$.simComparePlotByCell(cell_tbl, "mean_abs_error", "MAE")))
  cell_on <- env$.simComparePlotByCell(cell_tbl, "mean_abs_error", "MAE", mcse = TRUE)
  expect_true(.has_errorbar(cell_on))
  expect_no_error(ggplot2::ggplotGrob(cell_on))

  class_tbl <- tidyr::expand_grid(
    base_scenario_id = 1L, method = c("stimgate", "fbeta"),
    mismatch_val = c(0, 0.1)
  ) |>
    dplyr::mutate(
      transformation = "gaussian", scenario_desc = "Clean",
      fdp_median = 0.1, fdp_median_lower = 0, fdp_median_upper = 0.2,
      fdp_q90 = 0.3, sensitivity_median = 0.9,
      sensitivity_median_lower = 0.8, sensitivity_median_upper = 1,
      sensitivity_q10 = 0.5
    )
  expect_false(.has_errorbar(env$.simComparePlotClassification(class_tbl)))
  class_on <- env$.simComparePlotClassification(class_tbl, mcse = TRUE)
  expect_true(.has_errorbar(class_on))
  expect_no_error(ggplot2::ggplotGrob(class_on))
  bars <- .errorbar_data(class_on)
  expect_true(all(bars$ymin >= 0 & bars$ymax <= 1))
})

test_that("comparison summaries take MCSEs from the spread between datasets", {
  env <- .mcse_env()
  set.seed(3)
  make_raw <- function(n_iter) {
    tidyr::expand_grid(
      method = c("stimgate", "fbeta"), iter = seq_len(n_iter), sample = 1:4
    ) |>
      dplyr::mutate(
        n_cell = 1000,
        propRespTruth = 0.01,
        propRespEst = 0.01 * (1 + stats::rnorm(dplyr::n(), 0, 0.3)),
        thresholdFallbackUsed = stats::runif(dplyr::n()) < 0.2,
        propStim = 0.02, propUns = stats::runif(dplyr::n(), 0, 0.01),
        threshold = 1,
        nTruePos = 8L, nFalsePos = sample(0:3, dplyr::n(), TRUE),
        nFalseNeg = sample(0:4, dplyr::n(), TRUE), nTrueNeg = 980L
      )
  }
  raw <- make_raw(6)
  plain <- env$.simCompareSummariseFreqBs(raw, c("n_cell", "method"))
  summ <- env$.simCompareSummariseFreqBs(raw, c("n_cell", "method"), mcse = TRUE)
  expect_equal(summ$med_abs_rel_error, plain$med_abs_rel_error)
  stim <- raw[raw$method == "stimgate", ]
  rel <- abs(stim$propRespEst - stim$propRespTruth) / stim$propRespTruth
  expect_equal(
    summ$med_abs_rel_error_mcse[summ$method == "stimgate"],
    env$.analysis_mcse_between_units(rel, stim$iter, stats::median)
  )
  expect_equal(
    summ$med_abs_rel_error_upper - summ$med_abs_rel_error,
    stats::qnorm(0.975) * summ$med_abs_rel_error_mcse
  )
  expect_true(all(summ$fallback_rate_lower >= 0 & summ$fallback_rate_upper <= 1))

  few <- env$.simCompareSummariseFreqBs(make_raw(4), c("n_cell", "method"), mcse = TRUE)
  expect_true(all(is.na(few$med_abs_rel_error_mcse)))
  expect_true(all(is.na(few$med_abs_rel_error_lower)))

  class_summ <- env$.simCompareClassificationSummary(
    raw, c("n_cell", "method"),
    mcse = TRUE
  )
  expect_true(all(class_summ$fdp_median_lower >= 0))
  expect_true(all(class_summ$sensitivity_median_upper <= 1))
  expect_false(anyNA(class_summ$sensitivity_median_mcse))
  expect_error(
    env$.simCompareClassificationSummary(
      dplyr::select(raw, -"iter"), c("n_cell", "method"),
      mcse = TRUE
    ),
    "dataset column"
  )

  err <- env$.simCompareUnsignedErrorSummary(raw, c("n_cell", "method"), mcse = TRUE)
  expect_true(all(is.na(err$max_lower)))
  expect_false(anyNA(err$median_mcse))
  plain_err <- env$.simCompareUnsignedErrorSummary(raw, c("n_cell", "method"))
  expect_named(plain_err, c("n_cell", "method", "median", "q95", "max"))
  avg <- env$.simCompareErrorAverage(err, "method")
  expect_equal(avg$median_mcse, err$median_mcse[match(avg$method, err$method)])
})

test_that("summary-plot QMDs declare and read the show_mcse toggle", {
  for (file in c(
    "2a-sim-bw-freq_bs-global.qmd", "2b-sim-bias_uns-freq_bs.qmd",
    "3-sim-bw-est-base.qmd", "4-sim-bw-est-norm.qmd",
    "7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd"
  )) {
    lines <- readLines(file.path(root_dir, "analysis", file), warn = FALSE)
    expect_true("  show_mcse: true" %in% lines, info = file)
    expect_true(any(grepl(
      'show_mcse <- .as_flag(.get_qmd_param_env("show_mcse", "SHOW_MCSE", TRUE))',
      lines,
      fixed = TRUE
    )), info = file)
    style <- grep("analysis-plot-style.R", lines, fixed = TRUE)
    mcse <- grep("analysis-mcse.R", lines, fixed = TRUE)
    expect_length(mcse, 1L)
    expect_gt(mcse, style[[1]])
    expect_true(any(grepl("show_mcse", lines[-seq_len(mcse)], fixed = TRUE)), info = file)
  }
})
