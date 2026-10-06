.signed_percentile_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R", "acs_cytof-manual.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("unconditional percentiles retain zeros and match raw type-7 quantiles", {
  env <- .signed_percentile_env()
  x <- c(-2, -1, 0, 0, 0.5, 1, 7, NA, Inf)
  probs <- env$.simBandwidthSignedErrorProbs
  out <- env$.simBandwidthSignedErrorPercentiles(x, mcse = TRUE)
  expect_equal(unname(unlist(out[1, names(probs)], use.names = FALSE)),
    unname(stats::quantile(x[is.finite(x)], probs, names = FALSE)))
  expect_equal(out$median, 0)
  expect_equal(out$n_finite, 7L)
  expect_false("direction" %in% names(out))
  expect_false("prop" %in% names(out))
  expect_true(all(is.na(env$.simBandwidthSignedErrorPercentiles(c(NA, Inf))[
    names(probs)])))
})

test_that("independent percentile bounds and scenario averages use existing MC rules", {
  env <- .signed_percentile_env()
  x <- seq(-2, 4, length.out = 300)
  out <- env$.simBandwidthSignedErrorPercentiles(x, mcse = TRUE)
  for (name in names(env$.simBandwidthSignedErrorProbs)) {
    q <- env$.analysis_mcse_quantile(x, env$.simBandwidthSignedErrorProbs[[name]])
    expect_equal(out[[paste0(name, "_lower")]], q$lower)
    expect_equal(out[[paste0(name, "_upper")]], q$upper)
    expect_equal(out[[paste0(name, "_mcse")]], q$mcse)
  }
  # Bounds may cross zero; these are not conditional side summaries.
  around_zero <- env$.simBandwidthSignedErrorPercentiles(seq(-1, 1, length.out = 300))
  expect_lt(around_zero$median_lower, 0)
  expect_gt(around_zero$median_upper, 0)
  tbl <- dplyr::bind_rows(out, out)
  avg <- env$.simBandwidthSignedPercentileAverage(tbl, character())
  expect_equal(avg$median, out$median)
  expect_equal(avg$median_mcse, out$median_mcse / sqrt(2))
  expect_equal(avg$n_scenario_median, 2L)
})

test_that("outer-pair fallback is a whole-figure decision independent of display", {
  env <- .signed_percentile_env()
  ample <- env$.simBandwidthSignedErrorPercentiles(seq(-1, 2, length.out = 300))
  small <- env$.simBandwidthSignedErrorPercentiles(seq(-1, 2, length.out = 100))
  expect_identical(env$.simBandwidthSignedPercentileOuter(ample), c("q025", "q975"))
  expect_identical(env$.simBandwidthSignedPercentileOuter(small), c("q05", "q95"))
  expect_true(is.finite(small$q05_lower))
  suppressed <- ample
  suppressed$q975_upper <- NA_real_
  expect_identical(env$.simBandwidthSignedPercentileOuter(
    dplyr::bind_rows(ample, suppressed)), c("q05", "q95"))
  tbl <- dplyr::bind_rows(ample, small)
  tbl$bw <- c(1, 2)
  tbl$transformation <- "gaussian"
  off <- env$.simBandwidthSignedPercentilePlot(tbl)
  on <- env$.simBandwidthSignedPercentilePlot(tbl, mcse = TRUE)
  expect_equal(off$data$value, on$data$value)
  expect_null(on$labels$title)
  expect_null(on$labels$subtitle)
  expect_match(on$labels$caption, "5th/95th", fixed = TRUE)
  expect_setequal(as.character(on$data$percentile),
    c("5th", "10th", "25th", "50th (median)", "75th", "90th", "95th"))
  unavailable <- env$.simBandwidthSignedErrorPercentiles(c(-1, 0, 1))
  unavailable <- dplyr::bind_rows(unavailable, unavailable)
  unavailable$bw <- c(1, 2)
  unavailable$transformation <- "gaussian"
  expect_no_error(ggplot2::ggplot_build(env$.simBandwidthSignedPercentilePlot(
    unavailable, mcse = TRUE)))
})

test_that("comparison percentiles resample complete iter blocks and retain missing units", {
  env <- .signed_percentile_env()
  raw <- tidyr::expand_grid(iter = c(1L, 2L, 4L, 5L, 6L), sample = 1:3,
    method = c("stimgate", "fbeta")) |>
    dplyr::mutate(scenario = "a", sim_seed = 23L, propRespTruth = 1,
      propRespEst = 1 + iter^2 + sample, error = NA_character_)
  out <- env$.simCompareSignedPercentileSummary(raw, c("scenario", "method"), mcse = TRUE)
  rows <- raw[raw$method == "stimgate", ]
  x <- c(rows$propRespEst - 1, NA_real_)
  unit <- c(as.character(rows$iter), "3")
  direct <- env$.analysis_mcse_pooled_cols(x, unit,
    function(v) env$.analysis_mcse_quantile_finite(v, 0.5), "median", "sim_seed:23")
  expect_equal(out$median, rep(stats::median(rows$propRespEst - 1), 2))
  expect_equal(out$.boot_median[[1]], direct$.boot_median[[1]])
  expect_equal(out$.boot_median[[1]], out$.boot_median[[2]])
  expect_equal(out$median_mcse, rep(direct$median_mcse, 2))
  # Same biological draws in two scenarios: averaging must preserve covariance.
  paired <- dplyr::bind_rows(raw, dplyr::mutate(raw, scenario = "b"))
  avg <- env$.simCompareSignedPercentileAverage(paired, c("scenario", "method"), "method", mcse = TRUE)
  expect_equal(avg$median_mcse, out$median_mcse)
  expect_equal(avg$.boot_median, out$.boot_median)
  expect_true(all(avg$n_scenario_median == 2L))
  failed <- dplyr::mutate(raw, error = "comparator failed")
  expect_true(all(is.na(env$.simCompareSignedPercentileSummary(
    failed, c("scenario", "method"))$median)))
})

test_that("percentile plots build with single-method and method-colour encodings", {
  env <- .signed_percentile_env()
  q <- env$.simBandwidthSignedErrorPercentiles(seq(-2, 20, length.out = 300))
  tbl <- dplyr::bind_rows(q, q)
  tbl$bw <- c(1, 2)
  tbl$transformation <- "gaussian"
  single <- env$.simBandwidthSignedPercentilePlot(tbl, mcse = TRUE)
  expect_no_error(ggplot2::ggplot_build(single))
  expect_equal(length(unique(single$data$percentile)), 7L)
  # Brown below and teal above the median, as in the over/under views.
  colours <- unname(single$scales$get_scales("colour")$palette(7))
  expect_equal(colours, c(rep(env$.simBandwidthSignedErrorColours[["under_q95"]], 3), "#333333",
    rep(env$.simBandwidthSignedErrorColours[["over_q95"]], 3)))
  # Percentile aesthetics share one legend title, so ggplot2 merges them.
  titles <- vapply(c("colour", "alpha", "linewidth", "linetype"),
    function(a) single$scales$get_scales(a)$name, character(1))
  expect_true(all(titles == "Percentile"))
  multi_tbl <- dplyr::bind_rows(dplyr::mutate(tbl, method = "stimgate"),
    dplyr::mutate(tbl, method = "fbeta"))
  multi <- env$.simBandwidthSignedPercentilePlot(multi_tbl,
    alphas = c(outer = 0.2, tail = 0.4, quartile = 0.6, median = 1),
    linewidths = c(outer = 0.2, tail = 0.4, quartile = 0.6, median = 1))
  built <- ggplot2::ggplot_build(multi)
  expect_setequal(built$data[[2]]$linetype, c("dotted", "solid", "dashed"))
  expect_setequal(built$data[[2]]$alpha, c(0.2, 0.4, 0.6, 1))
  bias <- dplyr::mutate(tbl, bias_uns_multiplier = c(0.25, 0.5),
    bias_uns_basis = "bandwidth", mismatch_label = "No mismatch", n_cell = 1000, bw = 1)
  expect_no_error(ggplot2::ggplot_build(env$.simBandwidthSignedPercentilePlot(
    bias, x = "bias_uns_multiplier")))
})

test_that("ACS percentile intervals use common donor draws across methods and strata", {
  env <- .signed_percentile_env()
  tbl <- tidyr::expand_grid(SampleID = letters[1:6], cyt = c("a", "b"),
    pop = "CD4", method = c("stimgate", "fbeta")) |>
    dplyr::mutate(rel_error = match(SampleID, letters) - 3, thresholdFailed = FALSE)
  out <- env$.acsCytofManualSignedPercentiles(tbl)
  expect_true(all(vapply(out$.boot_median, function(v) identical(v, out$.boot_median[[1]]), logical(1))))
  expect_true(all(out$n_dataset_median == 6L))
  plot <- env$.simBandwidthSignedPercentilePlot(out, x = "cyt", mcse = TRUE,
    facet = ggplot2::facet_wrap(~pop))
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("all requested percentile QMD chunks parse and obey plot controls", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  labels <- list(
    "2a-sim-bw-freq_bs-global.qmd" = c("signed-percentiles-averaged-n-cell",
      "signed-percentiles-by-n-cell", "signed-percentiles-by-n-cell-prob"),
    "2b-sim-bias_uns-freq_bs.qmd" = c("signed-percentiles",
      "signed-percentiles-averaged-n-cell", "signed-percentiles-by-n-cell"),
    "7-sim-compare-freq_bs.qmd" = "signed-percentiles-cell-count",
    "8-sim-compare-freq_bs-batch.qmd" = "signed-percentiles-averaged",
    "9-real-compare-acs-cytof.qmd" = "manual-signed-percentiles"
  )
  for (doc in names(labels)) {
    lines <- readLines(file.path(root, "analysis", doc), warn = FALSE)
    if (startsWith(doc, "9-")) {
      expect_true(any(grepl('source(file.path(scripts_r_dir, "analysis-mcse.R"))', lines, fixed = TRUE)))
    }
    for (label in labels[[doc]]) {
      start <- which(lines == paste0("#| label: ", label))
      expect_length(start, 1L)
      end <- start + which(lines[seq.int(start + 1L, length(lines))] == "```")[[1]]
      code <- paste(lines[seq.int(start + 1L, end - 1L)], collapse = "\n")
      expect_no_error(parsed <- parse(text = code))
      env <- new.env(parent = baseenv())
      env$run_plots <- FALSE
      env$results_available <- FALSE
      expect_no_error(eval(parsed, env))
      expect_match(code, "signed_percentiles", fixed = TRUE)
      expect_false(grepl("RatioTwin|ratio_twins = TRUE", code))
    }
  }
})
