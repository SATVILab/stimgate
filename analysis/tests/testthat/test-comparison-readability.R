.readability_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R",
    "sim-compare-performance-plot.R", "sim-compare-qmd7-presentation.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("comparison panels have independent ranges and shared method encodings", {
  env <- .readability_env()
  tbl <- tidyr::expand_grid(transformation = c("gaussian", "gamma"),
    prob_response = c(0.01, 0.2), method = c("stimgate", "tailgate", "fbeta"),
    n_cell = c(1000, 10000)) |>
    dplyr::mutate(value = ifelse(.data$transformation == "gamma", 1000, 1) *
      .data$n_cell / 1000)
  p <- env$.simComparePlotByCell(tbl, "value", "Error", free_y = TRUE)
  built <- ggplot2::ggplot_build(p)
  expect_length(built$layout$panel_scales_y, 4L)
  expect_equal(p$data$value, tbl$value)
  expect_length(unique(built$data[[1]]$linetype), 3L)
  expect_length(unique(built$data[[2]]$shape), 3L)
  for (aes in c("colour", "shape", "linetype")) {
    expect_identical(p$scales$get_scales(aes)$name, "Method")
    expect_identical(p$scales$get_scales(aes)$labels, env$.analysis_method_labels)
  }
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)

  facets <- ggplot2::ggplot(tbl, ggplot2::aes(n_cell, value)) +
    ggplot2::geom_point() + env$.simCompareMismatchFacet(TRUE, ncol = 4)
  facets$data$statistic <- "Median"
  expect_length(ggplot2::ggplot_build(facets)$layout$panel_scales_y, 4L)
})

test_that("classification tails name the displayed outcomes and use readable wraps", {
  env <- .readability_env()
  tbl <- tidyr::expand_grid(base_scenario_id = 1:3,
    method = c("stimgate", "tailgate", "fbeta"), mismatch_val = c(0, 0.1)) |>
    dplyr::mutate(transformation = "gaussian", scenario_desc = paste("Scenario", .data$base_scenario_id),
      fdp_median = 0.1, fdp_q90 = 0.2, sensitivity_median = 0.9,
      sensitivity_q10 = 0.8, fpr_median = .data$base_scenario_id / 1000,
      fpr_q90 = .data$base_scenario_id / 100)
  paired <- env$.simComparePlotClassification(tbl)
  fpr <- env$.simComparePlotClassification(tbl, outcomes = "fpr", unit_scale = FALSE)
  expect_identical(unname(paired$scales$get_scales("linetype")$labels[["tail"]]),
    "90th FDP / 10th sensitivity")
  expect_identical(unname(fpr$scales$get_scales("linetype")$labels[["tail"]]), "90th percentile")
  expect_match(fpr$labels$x, "square-root spacing", fixed = TRUE)
  expect_equal(paired$facet$params$ncol, 2)
  pb <- ggplot2::ggplot_build(paired)
  expect_equal(pb$layout$panel_scales_y[[1]]$limits, c(0, 1))
  expect_length(ggplot2::ggplot_build(fpr)$layout$panel_scales_y, 3L)
  expect_length(unique(pb$data[[2]]$shape), 3L)
})

test_that("upper-tail percentages and placement methods retain the original values", {
  env <- .readability_env()
  tbl <- tidyr::expand_grid(scenario_desc = c("Clean", "Very clean"),
    method = c("stimgate", "tailgate", "fbeta"), mismatch_val = c(0, 0.1)) |>
    dplyr::mutate(mismatch_type = "mean_shift_negative", q90_abs_rel_error = 25,
      propUns_mean = 0.001)
  p <- env$.simComparePlotUpperTail(tbl)
  expect_identical(p$data$q90_abs_rel_error, tbl$q90_abs_rel_error)
  expect_match(p$scales$get_scales("y")$labels(25), "%", fixed = TRUE)
  expect_equal(p$facet$params$ncol, 2)
  expect_length(ggplot2::ggplot_build(p)$layout$panel_scales_y, 2L)
  placement <- ggplot2::ggplot_build(env$.simComparePlotPlacement(tbl))
  expect_length(unique(placement$data[[1]]$linetype), 3L)
  expect_length(unique(placement$data[[2]]$shape), 3L)
})

test_that("floor tables count zero negative and small positive estimates separately", {
  env <- .readability_env()
  tbl <- tidyr::expand_grid(method = c("stimgate", "tailgate"),
    prob_response = c(0.01, 0.2), sample = 1:6) |>
    dplyr::mutate(propRespEst = c(-0.01, 0, 0.00001, 0.0001, 0.01, 0.2)[.data$sample])
  out <- env$.simCompareQmd7FloorTable(tbl, 0.0001)
  expect_equal(nrow(out), 4L)
  expect_equal(out$n_method_sample_estimates, rep(6L, 4))
  expect_equal(out$zero, rep("1 / 6", 4))
  expect_equal(out$negative, rep("1 / 6", 4))
  expect_equal(out$positive_at_or_below_floor, rep("2 / 6", 4))
})

test_that("gate patterns distinguish methods without moving coincident thresholds", {
  env <- .readability_env()
  cells <- tidyr::expand_grid(condition = c("stim", "unstim"), i = 1:20) |>
    dplyr::mutate(expr = .data$i / 10, label = ifelse(.data$i > 15, "gp", "gn"),
      mismatch_type = "mean_shift_negative", mismatch_val = 0)
  gates <- tibble::tibble(method = c("stimgate", "tailgate", "fbeta"), threshold = 1.5,
    mismatch_type = "mean_shift_negative", mismatch_val = 0)
  p <- env$.simComparePlotGateDiagnostic(cells, gates)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[3]]$xintercept, rep(1.5, 3))
  expect_length(unique(built$data[[3]]$linetype), 3L)
  expect_identical(p$scales$get_scales("alpha")$name, "Reference tube")
})

test_that("dataset maxima keep occurrence fixed and free each severity panel", {
  env <- .readability_env()
  tbl <- tidyr::expand_grid(transformation = c("gaussian", "gamma"),
    direction = c("over", "under"), n_cell = c(1000, 10000)) |>
    dplyr::mutate(method = "stimgate", prob_response = 0.01,
      occurrence = 0.5, occurrence_lower = 0.2, occurrence_upper = 0.8,
      severity = ifelse(.data$transformation == "gamma", 10, 0.2),
      severity_lower = .data$severity / 2, severity_upper = .data$severity * 1.2,
      n_eligible = 20, n_incomplete = 0, n_affected = 10)
  occurrence <- env$.simCompareDatasetMaxPlot(tbl, "occurrence", by_prob = FALSE)
  severity <- env$.simCompareDatasetMaxPlot(tbl, "severity", by_prob = FALSE)
  expect_equal(ggplot2::ggplot_build(occurrence)$layout$panel_scales_y[[1]]$limits, c(0, 1))
  expect_length(ggplot2::ggplot_build(severity)$layout$panel_scales_y, 4L)
  expect_equal(severity$data$value, ifelse(tbl$direction == "over", 1, -1) * tbl$severity)
})

test_that("figure loop places floor tables immediately after each figure", {
  env <- .readability_env()
  events <- character()
  env$.analysis_print_save_fig <- function(...) events <<- c(events, "figure")
  tbl <- tibble::tibble(method = c("stimgate", "tailgate", "fbeta"),
    mean_pos_setting = "high", x = 1, y = 1)
  invisible(utils::capture.output(env$.simCompareFigureLoop(tbl,
    make_plot = function(d) ggplot2::ggplot(d, ggplot2::aes(x, y)),
    dir = tempdir(), file_fn = function(pos, extra) "unused.png", height = 9, level = 4L,
    after_plot = function(d, p) {
      expect_equal(nrow(p$data), nrow(d))
      events <<- c(events, "table")
    })))
  expect_gt(length(events), 0L)
  expect_identical(events, rep(c("figure", "table"), length(events) / 2))
})

test_that("QMD split headings and filenames retain response and mismatch settings", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  q7 <- paste(readLines(file.path(root, "analysis", "7-sim-compare-freq_bs.qmd")), collapse = "\n")
  q8 <- paste(readLines(file.path(root, "analysis", "8-sim-compare-freq_bs-batch.qmd")), collapse = "\n")
  for (quantity in c("occurrence", "severity")) {
    expect_match(q7, paste0('.simCompareDatasetMaxPlot(d, "', quantity,
      '", by_prob = FALSE'), fixed = TRUE)
    expect_match(q7, paste0('"dataset_max_', quantity, '_", pos,\n      "_prob_response_", extra'), fixed = TRUE)
  }
  expect_match(q7, "after_plot = function(d, p)", fixed = TRUE)
  expect_match(q7, "denominators count method/sample estimates", fixed = TRUE)
  for (label in c("threshold-placement-unstim", "upper-tail-and-fallbacks")) {
    block <- strsplit(q8, paste0("#| label: ", label), fixed = TRUE)[[1]][2]
    block <- strsplit(block, "```", fixed = TRUE)[[1]][1]
    expect_match(block, 'extra_col = "mismatch_type"', fixed = TRUE)
    expect_match(block, "file_fn = .compare_mismatch_file", fixed = TRUE)
    expect_match(block, "extra_heading = .compare_mismatch_heading", fixed = TRUE)
  }
})
