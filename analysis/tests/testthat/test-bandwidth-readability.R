.bandwidth_readability_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
                 "sim-bandwidth-analysis-plot.R")) {
    source(file.path(testthat::test_path(), "../../../scripts/r", file), local = env)
  }
  env
}

test_that("signed display transform is continuous, monotone and invertible", {
  env <- .bandwidth_readability_env()
  trans <- env$.simBandwidthSignedErrorTrans()
  x <- c(-1e6, -100, -8, -2, -1 - 1e-8, -1, -1 + 1e-8, -0.5, 0, 1e-8, 1, 15, 100)
  expect_true(all(diff(trans$transform(x)) > 0))
  expect_equal(trans$inverse(trans$transform(x)), x, tolerance = 1e-8)
  expect_equal(trans$transform(c(-1 - 1e-8, -1 + 1e-8)), c(-1, -1), tolerance = 1e-7)
  expect_equal(trans$transform(c(-8, -1, -0.5, 0, 3)), c(-4, -1, -0.5, 0, 2))
  expect_true(all(c(-8, -4, -2, -1) %in% trans$breaks(c(-8, 15))))
})

test_that("absolute error scale has nonnegative percent breaks", {
  env <- .bandwidth_readability_env()
  plot <- ggplot2::ggplot() + env$.simBandwidthAbsErrorLayers()
  scale <- plot$scales$get_scales("y")
  expect_true(all(scale$breaks(c(0, 15)) >= 0))
  expect_identical(scale$labels(0.25), "25%")
  tbl <- tibble::tibble(bw = 0.1, bias_uns_multiplier = 1, bias_uns_basis = "bandwidth",
    mismatch_label = "None", median_abs_rel_error = 0.25, q90_abs_rel_error = 0.5, max_abs_rel_error = 1,
    median_abs_rel_error_lower = 0.1, median_abs_rel_error_upper = 0.4)
  bias <- env$.simBandwidthBiasRelativeErrorPlot(tbl, mcse = TRUE)
  expect_identical(bias$scales$get_scales("y")$labels(0.25), "25%")
  expect_s3_class(bias$layers[[1]]$geom, "GeomErrorbar")
  expect_equal(bias$layers[[1]]$aes_params$linewidth, 0.3)
  expect_equal(bias$layers[[1]]$aes_params$alpha, 0.6)
  built <- ggplot2::ggplot_build(bias)
  expect_equal(built$data[[1]]$ymin, 0.1)
  expect_equal(built$data[[1]]$ymax, 0.4)
  expect_gt(built$data[[1]]$xmax - built$data[[1]]$xmin, 0)
  expect_s3_class(bias$facet, "FacetWrap")
  expect_no_error(ggplot2::ggplotGrob(bias))
})

test_that("bandwidth ranks preserve unequal transformation grids and merge legends", {
  env <- .bandwidth_readability_env()
  tbl <- tibble::tibble(transformation = c("gaussian", "gaussian", "skew", "skew", "gamma", "gamma", "gamma"),
    bw = c(0.05, 0.1, 0.05, 0.1, 0.001, 0.0025, 0.005))
  ranked <- env$.simBandwidthRankData(tbl)
  expect_equal(as.integer(ranked$bw_rank), c(1, 2, 1, 2, 1, 2, 3))
  expect_identical(ranked$bw, tbl$bw)
  missing_middle <- ranked[as.character(ranked$bw_rank) != "2", ]
  sparse_scales <- env$.simBandwidthRankScales(missing_middle)
  expect_identical(sparse_scales[[1]]$breaks, c("1", "3"))
  scales <- env$.simBandwidthRankScales(tbl)
  expect_identical(scales[[1]]$name, scales[[2]]$name)
  expect_identical(scales[[1]]$labels, scales[[2]]$labels)
  expect_match(scales[[1]]$labels[1], "0.001 (Gamma)", fixed = TRUE)
  expect_match(scales[[1]]$labels[1], "0.05 (Gaussian/Skew)", fixed = TRUE)
  palette <- scales[[1]]$palette(3)
  expect_equal(length(unique(palette)), 3L)
  expect_lt(max(grDevices::col2rgb(palette)), 220)
})

test_that("threshold figures retain values with horizontal IQR and vertical jitter", {
  env <- .bandwidth_readability_env()
  tbl <- tidyr::expand_grid(transformation = "gaussian", prob_response = 0.01,
    n_cell = 1000, bw = c(0.05, 0.1), threshold = 1:5)
  med <- env$.simBandwidthThresholdPlot(tbl)
  expect_equal(med$data$threshold_median, c(3, 3))
  expect_equal(med$data$threshold_iqr_lower, c(2, 2))
  expect_equal(med$data$threshold_iqr_upper, c(4, 4))
  expect_s3_class(med$layers[[1]]$geom, "GeomSegment")
  invalid <- dplyr::mutate(tbl, valid_estimate = threshold < 5)
  expect_equal(env$.simBandwidthThresholdPlot(invalid)$data$threshold_median, c(2.5, 2.5))
  expect_equal(env$.simBandwidthThresholdPlot(invalid, samples = TRUE)$data$threshold_median, c(3, 3))
  sample <- env$.simBandwidthThresholdPlot(tbl, samples = TRUE)
  built <- ggplot2::ggplot_build(sample)
  expect_equal(built$data[[1]]$x, tbl$threshold)
  expect_equal(nrow(built$data[[1]]), nrow(tbl))
  expect_no_error(ggplot2::ggplotGrob(sample))
})

test_that("cap triangles and negative-estimate notes identify displayed figures", {
  env <- .bandwidth_readability_env()
  tbl <- tibble::tibble(transformation = "gaussian", bw = 0.1, direction = "under",
    prop = 1, median = -2, q95 = -4, max = 40)
  plot <- env$.simBandwidthGlobalSignedErrorPlot(tbl)
  built <- ggplot2::ggplot_build(plot)
  triangle <- which(vapply(plot$layers, function(l) identical(l$aes_params$shape, 17), logical(1)))
  expect_length(triangle, 1L)
  expect_equal(built$data[[triangle]]$y, log2(16))
  expect_equal(nrow(built$data[[triangle]]), 1L)
  text <- NULL
  expect_message(text <- paste(utils::capture.output(env$.simBandwidthDisplayNote(plot, "example.png")), collapse = "\n"),
    "example.png:.*below -100%")
  expect_match(text, "triangle: above display cap", fixed = TRUE)
  expect_match(text, "compressed log scale", fixed = TRUE)
  ordinary <- env$.simBandwidthGlobalSignedErrorPlot(dplyr::mutate(tbl, median = -0.2, q95 = -0.5, max = -1))
  expect_message(env$.simBandwidthDisplayNote(ordinary, "ordinary.png"), NA)
  expect_no_error(ggplot2::ggplotGrob(plot))
})

test_that("readability QMD contracts retain primary views without density overlays or cropping", {
  root <- file.path(testthat::test_path(), "../../..")
  read <- function(file) paste(readLines(file.path(root, "analysis", file), warn = FALSE), collapse = "\n")
  bias <- read("2b-sim-bias_uns-freq_bs.qmd")
  expect_false(grepl('#| label: relative-error\n', bias, fixed = TRUE))
  expect_false(grepl('#| label: signed-error\n', bias, fixed = TRUE))
  for (label in c("relative-error-averaged-n-cell", "relative-error-by-n-cell", "signed-error-averaged-n-cell", "signed-error-by-n-cell")) {
    expect_match(bias, paste0("#| label: ", label), fixed = TRUE)
  }
  trans <- read("1-sim-trans.qmd")
  expect_false(grepl("coord_cartesian(ylim = c(0, 0.5))", trans, fixed = TRUE))
  global <- read("2a-sim-bw-freq_bs-global.qmd")
  expect_match(global, ".simBandwidthThresholdPlot(threshold_tbl_curr, samples = TRUE)", fixed = TRUE)
  expect_match(global, ".simBandwidthThresholdDensityPlot(NULL", fixed = TRUE)
  base <- read("3-sim-bw-est-base.qmd")
  expect_match(base, "shape = .data$bw_mtd", fixed = TRUE)
})


test_that("owned readability QMD R chunks parse", {
  root <- file.path(testthat::test_path(), "../../..", "analysis")
  for (file in c("1-sim-trans.qmd", "2a-sim-bw-freq_bs-global.qmd",
    "2b-sim-bias_uns-freq_bs.qmd", "3-sim-bw-est-base.qmd")) {
    lines <- readLines(file.path(root, file), warn = FALSE)
    starts <- which(grepl("^```\\{r", lines))
    for (start in starts) {
      end <- which(lines == "```" & seq_along(lines) > start)[1L]
      expect_no_error(parse(text = lines[seq.int(start + 1L, end - 1L)]))
    }
  }
})

test_that("coverage notes stay compact and export every plotted setting", {
  env <- .bandwidth_readability_env()
  panels <- tibble::tibble(mean_pos_setting = "low", n_cell = 100,
    transformation = env$.analysis_trans_factor("gaussian"), bw = seq(0.1, 1, by = 0.1))
  plot <- ggplot2::ggplot(panels, ggplot2::aes(bw, n_cell)) + ggplot2::geom_point()
  summary <- tibble::tibble(mean_pos_setting = "low", n_cell = 100,
    transformation = "gaussian", bw = seq(0.1, 1, by = 0.1), n_sample = 25L,
    n_valid = c(24L, rep(25L, 9)), n_failed = c(1L, rep(0L, 9)),
    n_provenance = 25L, n_fallback = c(0L, 2L, rep(0L, 8)))
  root <- .local_projr_root()
  out <- paste(utils::capture.output(tbl <- env$.simBandwidthPrintCoverage(plot, summary,
    table_parts = c("analysis-name", "coverage.csv"), path_root = root)),
    collapse = "\n")
  expect_match(out, "10 plotted settings; 25 samples each (250 in total)", fixed = TRUE)
  expect_match(out, "Failed estimates: 1 (0.4%)", fixed = TRUE)
  expect_match(out, "Threshold fallbacks: 2 of 250 (0.8%)", fixed = TRUE)
  # No scenario rows are printed; the CSV and returned table retain all settings.
  expect_equal(nrow(tbl), 10L)
  expect_false(grepl("\n|", out, fixed = TRUE))
  saved <- readr::read_csv(.projr_output_path(root, "table", "analysis-name", "coverage.csv"),
    show_col_types = FALSE)
  expect_equal(nrow(saved), nrow(tbl))
  expect_equal(saved$n_failed, tbl$n_failed)
  expect_equal(saved$n_fallback, tbl$n_fallback)
  expect_match(out, "output/table/analysis-name/coverage.csv", fixed = TRUE)
  clean <- summary; clean$n_failed <- 0L; clean$n_fallback <- 0L
  out <- paste(utils::capture.output(env$.simBandwidthPrintCoverage(plot, clean)), collapse = "\n")
  expect_false(grepl("Settings with failures", out, fixed = TRUE))
})
