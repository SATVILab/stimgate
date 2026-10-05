.performance_figure_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R",
    "sim-compare-performance-plot.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.performance_figure_fixture <- function() {
  tibble::tibble(
    method = "stimgate", transformation = "gaussian", prob_response = 0.01,
    n_cell = c(1000, 10000, 1000, 10000), direction = c("over", "over", "under", "under"),
    occurrence = c(0.5, 0.8, 0.3, 0.6), occurrence_lower = c(0.3, 0.6, 0.1, 0.4),
    occurrence_upper = c(0.7, 1, 0.5, 0.8),
    severity = c(0.4, 0.6, 0.2, 0.5), severity_lower = c(0.2, 0.3, 0.1, 0.3),
    severity_upper = c(0.7, 0.9, 0.4, 0.8),
    n_dataset_total = 20L, n_eligible = 20L, n_incomplete = 0L,
    n_affected = c(10L, 16L, 6L, 12L),
    eligible_dataset_ids = paste(1:20, collapse = ","), cohort_matches_stimgate = TRUE
  )
}

test_that("dedicated maximum plots preserve row-specific intervals and ratio geometry", {
  env <- .performance_figure_env()
  data <- .performance_figure_fixture()
  occurrence <- env$.simCompareDatasetMaxPlot(data, "occurrence", mcse = TRUE)
  expect_equal(occurrence$data$lower, data$occurrence_lower)
  expect_equal(occurrence$data$upper, data$occurrence_upper)
  expect_equal(occurrence$data$value_shown, data$occurrence)
  expect_no_error(ggplot2::ggplotGrob(occurrence))
  severity <- env$.simCompareDatasetMaxPlot(data, "severity", mcse = TRUE)
  expect_equal(severity$data$value_shown, c(0.4, 0.6, -0.2, -0.5))
  expect_equal(severity$data$lower, c(0.2, 0.3, -0.4, -0.8))
  expect_equal(severity$data$upper, c(0.7, 0.9, -0.1, -0.3))
  plain <- env$.simCompareDatasetMaxPlot(data, "severity")
  expect_equal(severity$data$value_shown, plain$data$value_shown)
  original <- ggplot2::ggplot_build(severity)
  ratio <- ggplot2::ggplot_build(env$.simBandwidthRatioPlot(severity))
  expect_length(ratio$data, length(original$data))
  for (i in seq_along(original$data)) {
    coordinates <- intersect(c("x", "y", "ymin", "ymax", "xend", "yend"), names(original$data[[i]]))
    expect_equal(ratio$data[[i]][coordinates], original$data[[i]][coordinates])
  }
  expect_match(severity$labels$caption, "affected")
  expect_match(severity$labels$caption, "95%", fixed = TRUE)
})

test_that("main percentile and separate maximum figure contracts remain distinct", {
  env <- .performance_figure_env()
  data <- tidyr::expand_grid(transformation = "gaussian", method = "stimgate",
    mismatch_val = c(0, 0.1)) |>
    dplyr::mutate(median = 0.1, q95 = 0.3, max = 20)
  main <- env$.simComparePlotMismatchError(data)
  expect_setequal(as.character(main$data$statistic), c("Median", "95th percentile"))
  expect_true(max(main$data$value) < 20)
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (qmd in c("7-sim-compare-freq_bs.qmd", "8-sim-compare-freq_bs-batch.qmd")) {
    text <- paste(readLines(file.path(root, "analysis", qmd), warn = FALSE), collapse = "\n")
    expect_match(text, "label: dataset-max-occurrence", fixed = TRUE)
    expect_match(text, "label: dataset-max-severity", fixed = TRUE)
    expect_match(text, "label: dataset-max-coverage", fixed = TRUE)
    expect_match(text, "dataset_max_occurrence_by_n_cell", fixed = TRUE)
    expect_match(text, "dataset_max_signed_error_severity_by_n_cell", fixed = TRUE)
    expect_match(text, "ratio_twins = TRUE", fixed = TRUE)
    expect_false(grepl("label: max-rel-error-cell-count", text, fixed = TRUE))
  }
})
