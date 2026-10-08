.signed_percentile_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R", "acs_cytof-manual.R",
    "acs_cytof-plot_cyt.R")) {
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
  # The fallback pair is named in the band legend rather than the caption.
  expect_identical(levels(on$layers[[1]]$data$band),
    c("5th\u201395th", "10th\u201390th", "25th\u201375th"))
  expect_false(grepl("Outer pair|Unavailable intervals", paste(on$labels$caption, "")))
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
  # Nested bands share one hue, darker towards the median; the median is a line.
  built_single <- ggplot2::ggplot_build(single)
  ribbon <- built_single$data[[1]]
  expect_setequal(ribbon$fill, c("#C7EAE5", "#80CDC1", "#35978F"))
  band_q <- list(c("q025", "q975"), c("q10", "q90"), c("q25", "q75"))
  for (i in seq_along(band_q)) {
    rows <- ribbon[ribbon$fill == c("#C7EAE5", "#80CDC1", "#35978F")[[i]], ]
    expect_equal(sort(rows$ymin), sort(env$.simBandwidthSignedErrorTrans()$transform(pmin(tbl[[band_q[[i]][1]]], env$.simBandwidthSignedErrorCap))))
    expect_equal(sort(rows$ymax), sort(env$.simBandwidthSignedErrorTrans()$transform(pmin(tbl[[band_q[[i]][2]]], env$.simBandwidthSignedErrorCap))))
  }
  expect_identical(single$scales$get_scales("fill")$name, "Percentile range")
  expect_identical(single$labels$y, "Percentage deviation from true response")
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
    "8-sim-compare-freq_bs-batch.qmd" = c("signed-percentiles-averaged",
      "signed-percentiles-by-n-cell"),
    "9-real-compare-acs-cytof.qmd" = c("manual-signed-percentiles",
      "manual-signed-percentiles-by-stimulus"),
    "10-real-compare-acs-cytof-validation.qmd" = "validation-signed-percentiles-by-stimulus"
  )
  for (doc in names(labels)) {
    lines <- readLines(file.path(root, "analysis", doc), warn = FALSE)
    if (startsWith(doc, "9-") || startsWith(doc, "10-")) {
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


test_that("stimulus percentile summaries retain zeros, failures and common donor draws", {
  env <- .signed_percentile_env()
  tbl <- tidyr::expand_grid(SampleID = letters[1:6], cyt = c("a", "b"),
    pop = "CD4", method = c("stimgate", "fbeta"), stim = c("s1", "s2")) |>
    dplyr::mutate(rel_error = match(SampleID, letters) - 3, thresholdFailed = FALSE)
  groups <- c("method", "pop", "cyt", "stim")
  out <- env$.acsCytofManualSignedPercentiles(tbl, group_cols = groups)
  expect_equal(nrow(out), 8L)
  expect_true(all(out$n_finite == 6L))
  expect_true(all(vapply(out$.boot_median,
    function(v) identical(v, out$.boot_median[[1]]), logical(1))))
  failed <- tbl
  failed$thresholdFailed[failed$stim == "s2"] <- TRUE
  failed_out <- env$.acsCytofManualSignedPercentiles(failed, group_cols = groups)
  expect_true(all(is.na(failed_out$median[failed_out$stim == "s2"])))
  expect_true(all(failed_out$n_finite[failed_out$stim == "s2"] == 0L))
  # A donor absent from one stimulus remains an empty block in its bootstrap.
  sparse <- tbl[!(tbl$stim == "s2" & tbl$SampleID == "f"), ]
  sparse_out <- env$.acsCytofManualSignedPercentiles(sparse, group_cols = groups)
  direct <- env$.simBandwidthSignedErrorPercentiles(c(-2, -1, 0, 1, 2, NA),
    mcse = TRUE, unit = letters[1:6], bootstrap_family = "acs-donors")
  expect_equal(sparse_out$.boot_median[[which(sparse_out$stim == "s2")[[1]]]],
    direct$.boot_median[[1]])
})

test_that("ACS stimulus percentile chunks print every requested method and keep strata separate", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  docs <- c("9-real-compare-acs-cytof.qmd", "10-real-compare-acs-cytof-validation.qmd")
  labels <- c("manual-signed-percentiles-by-stimulus", "validation-signed-percentiles-by-stimulus")
  for (i in seq_along(docs)) {
    env <- .signed_percentile_env()
    env$manual_comparison_tbl <- tidyr::expand_grid(SampleID = letters[1:6],
      cyt = c("a", "b"), pop = c("CD4", "CD8"),
      method = c("stimgate", "fbeta", "tailgate"), stim = c("s1", "s2")) |>
      dplyr::mutate(rel_error = match(SampleID, letters) - 3, thresholdFailed = FALSE)
    env$validation_methods <- env$.acsCytofValidationMethods()
    env$fig_key <- "test"
    env$root_dir <- root
    # Exercise real print orchestration, intercepting filesystem saves only.
    env$.analysis_fig_dir <- function(parts, ...) file.path("fig", paste(parts, collapse = "/"))
    env$.simBandwidthDisplayNote <- function(...) invisible(NULL)
    plots <- list()
    paths <- character()
    env$print <- function(x, ...) {
      if (inherits(x, "ggplot")) plots[[length(plots) + 1L]] <<- x
      invisible(x)
    }
    env$.analysis_save_fig <- function(plot, path, ...) {
      paths <<- c(paths, path)
      invisible(NULL)
    }
    tables <- list()
    env$.analysis_report_table <- function(tbl, path_parts, ...) {
      tables[[length(tables) + 1L]] <<- path_parts
      invisible(tbl)
    }
    lines <- readLines(file.path(root, "analysis", docs[[i]]), warn = FALSE)
    start <- which(lines == paste0("#| label: ", labels[[i]]))
    end <- start + which(lines[seq.int(start + 1L, length(lines))] == "```")[[1]]
    code <- parse(text = lines[seq.int(start + 1L, end - 1L)])
    env$run_plots <- TRUE
    capture.output(eval(code, env))
    expect_length(plots, if (i == 1L) 4L else 6L)
    expect_length(unique(paths), length(plots))
    # One bootstrap-validity CSV per figure, each at its own path
    expect_length(unique(tables), length(plots))
    for (plot in plots) {
      expect_equal(dplyr::n_distinct(plot$data$stim), 1L)
      expect_equal(dplyr::n_distinct(plot$data$percentile), 7L)
      expect_setequal(plot$data$pop, c("CD4", "CD8"))
      expect_no_error(ggplot2::ggplot_build(plot))
    }
    expect_setequal(unlist(lapply(plots, function(p) as.character(p$data$method))),
      c("stimgate", "fbeta", "tailgate"))
    if (i == 1L) {
      expect_setequal(plots[[1]]$data$method, c("stimgate", "fbeta", "tailgate"))
      expect_setequal(plots[[2]]$data$method, c("stimgate", "fbeta"))
    } else {
      expect_true(all(vapply(plots, function(p) dplyr::n_distinct(p$data$method) == 1L,
        logical(1))))
    }
    plots <- list()
    paths <- character()
    env$run_plots <- FALSE
    capture.output(eval(code, env))
    expect_length(plots, 0L)
    expect_length(paths, 0L)
  }
})

test_that("main signed percentile views precede conditional severity diagnostics", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (doc in c("7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd",
    "9-real-compare-acs-cytof.qmd")) {
    lines <- readLines(file.path(root, "analysis", doc), warn = FALSE)
    main <- which(grepl("^#\\| label: .*signed-percentiles", lines))
    conditional <- which(grepl("^#\\| label: (signed-error|manual-signed-error|dataset-max-severity)", lines))
    expect_true(max(main) < min(conditional), info = doc)
  }
})


test_that("Analysis 8 per-cell percentile chunk keeps cell counts and mismatch types separate", {
  env <- .signed_percentile_env()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  lines <- readLines(file.path(root, "analysis", "8-sim-compare-freq_bs-batch.qmd"),
    warn = FALSE)
  chunk <- function(label) {
    start <- which(lines == paste0("#| label: ", label))
    end <- start + which(lines[seq.int(start + 1L, length(lines))] == "```")[[1]]
    parse(text = lines[seq.int(start + 1L, end - 1L)])
  }
  env$run_plots <- TRUE
  env$results_available <- TRUE
  env$show_mcse <- FALSE
  env$mcse_mode <- "off"
  eval(chunk("figure-helpers"), env)
  env$.compare_fresh_fig_dir <- function(name) name
  q <- env$.simBandwidthSignedErrorPercentiles(seq(-1, 1, length.out = 300))
  keys <- tidyr::expand_grid(method = c("stimgate", "fbeta", "tailgate"),
    n_cell = c(1000, 10000), mismatch_type = c("mean_shift_negative", "sd_inflation"),
    mismatch_val = c(0, 1)) |>
    dplyr::mutate(transformation = "gaussian", mean_pos_setting = "high")
  env$percentile_by_cell <- dplyr::bind_cols(keys, q[rep(1L, nrow(keys)), ])
  plots <- list()
  paths <- character()
  env$.analysis_print_save_fig <- function(plot, path, ..., mcse_mode) {
    expect_identical(mcse_mode, "off")
    plots[[length(plots) + 1L]] <<- plot
    paths <<- c(paths, path)
    invisible(plot)
  }
  capture.output(eval(chunk("signed-percentiles-by-n-cell"), env))
  expect_length(plots, 8L)
  expect_length(unique(paths), 8L)
  for (plot in plots) {
    expect_equal(dplyr::n_distinct(plot$data$n_cell), 1L)
    expect_equal(dplyr::n_distinct(plot$data$mismatch_type), 1L)
    expect_equal(dplyr::n_distinct(plot$data$percentile), 7L)
    expect_no_error(ggplot2::ggplot_build(plot))
  }
  expect_setequal(plots[[1]]$data$method, c("stimgate", "fbeta", "tailgate"))
  expect_setequal(plots[[5]]$data$method, c("stimgate", "fbeta"))
  plots <- list()
  env$results_available <- FALSE
  capture.output(eval(chunk("signed-percentiles-by-n-cell"), env))
  expect_length(plots, 0L)
})
