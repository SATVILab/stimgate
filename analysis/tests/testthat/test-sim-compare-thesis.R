root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.thesisCompareEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c("analysis-runtime.R", "analysis-plot-style.R", "sim-bandwidth-analysis-plot.R",
    "sim-compare-thesis.R")) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

.thesisGrid <- function() {
  g <- expand.grid(method = c("stimgate", "fbeta", "tailgate", "tailgate_default"),
    transformation = c("gaussian", "skew", "gamma"), prob_response = c(0.2, 0.002),
    n_cell = c(1e3, 1e5), mean_pos_setting = c("low", "high"), stringsAsFactors = FALSE)
  g$q10 <- -0.5; g$median <- 0.1; g$q90 <- 3
  g$f1_median <- seq(0, 1, length.out = nrow(g))
  g$f1_q10 <- rev(g$f1_median)
  g$f1_q90 <- pmax(g$f1_median, g$f1_q10)
  g
}

test_that("thesis comparison figures build for one separation and leave out default Tailgate", {
  env <- .thesisCompareEnv()
  g <- .thesisGrid()
  f1 <- env$.simCompareThesisF1Plot(g, "low")
  expect_s3_class(f1, "ggplot")
  expect_false("tailgate_default" %in% as.character(f1$data$method))
  expect_setequal(unique(f1$data$mean_pos_setting), "low")
  expect_no_error(ggplot2::ggplot_build(f1))
  # the 10th-percentile version fills from f1_q10, not the median
  f1_q10 <- env$.simCompareThesisF1Plot(g, "low", stat = "q10")
  expect_equal(f1_q10$data$value, f1_q10$data$f1_q10)
  expect_equal(f1$data$value, f1$data$f1_median)
  expect_no_error(ggplot2::ggplot_build(f1_q10))
  expect_error(env$.simCompareThesisF1Plot(g, "low", stat = "q90"))
  freq <- env$.simCompareThesisFreqPlot(g, "high")
  expect_s3_class(freq, "ggplot")
  dir <- withr::local_tempdir()
  paths <- env$.analysis_save_thesis_fig(freq, file.path(dir, "freq"), height = 215)
  expect_true(all(file.exists(paths)))
})

test_that("frequency axis labels are actual frequencies on the signed scale", {
  env <- .thesisCompareEnv()
  expect_identical(env$.simCompareThesisPercent(c(0.2, 0.002, 1e-4, 0)), c("20%", "0.2%", "0.01%", "0%"))
  expect_identical(env$.simCompareThesisProbLab("0.002"), "Response 0.2%")
})

test_that("Analysis 7 saves its thesis figures to output-thesis results/qmd7", {
  qmd <- paste(readLines(file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")), collapse = "\n")
  expect_match(qmd, 'thesis_parts <- c("results", "qmd7", fig_key[-1L])', fixed = TRUE)
  for (kind in c('"f1")', '"f1-q10")', 'paste0("f1-percentiles", suffix)', 'paste0("frequency", suffix)',
    '"thresholds")', '"thresholds-by-cells")',
    'paste0("thresholds-", tr)', 'paste0("thresholds-", .analysis_label_cells(nc))')) {
    expect_match(qmd, kind, fixed = TRUE)
  }
  expect_match(qmd, "rug_n = 100L", fixed = TRUE)
  expect_false(grepl("no_tailgate", qmd, fixed = TRUE))
})

test_that("threshold figures pool, split by cell count and change columns", {
  env <- .thesisCompareEnv()
  withr::local_seed(1)
  grid <- expand.grid(approach = c("stimgate", "fbeta", "tailgate"),
    transformation = c("gaussian", "skew"), prob_response = c(0.01, 1e-4),
    n_cell = c(5e3, 1e5), rep = 1:10, stringsAsFactors = FALSE)
  # the lowest response only exists at the largest cell count, as in the grid
  grid <- grid[!(grid$prob_response == 1e-4 & grid$n_cell == 5e3), ]
  grid$threshold <- stats::rnorm(nrow(grid), 4 + (grid$n_cell == 5e3))
  grid$mean_pos_setting <- "low"
  grid$mean_pos <- 4.5
  dens <- expand.grid(transformation = c("gaussian", "skew"), prob_response = c(0.01, 1e-4),
    mean_pos = 4.5, condition = c("stimulated", "positive_rug"), expression = seq(0, 8, 0.5),
    stringsAsFactors = FALSE)
  dens$density <- stats::dnorm(dens$expression, 2)

  pooled <- env$.simCompareThesisThresholdPlot(grid, dens, "low")
  expect_no_error(ggplot2::ggplot_build(pooled))
  steps <- pooled$layers[[2]]$data
  expect_true(all(is.na(steps$cells)))
  expect_setequal(levels(steps$column), c("Gaussian", "Skew"))

  split <- env$.simCompareThesisThresholdPlot(grid, dens, "low", split_cells = TRUE)
  expect_no_error(ggplot2::ggplot_build(split))
  s <- split$layers[[2]]$data
  # one outline per method and cell count where both exist
  expect_setequal(as.character(unique(s$cells)), c("5k", "100k"))
  expect_equal(nrow(dplyr::distinct(s[s$prob_response == "0.01" & s$transformation == "Gaussian", ],
    approach, cells)), 6L)

  one <- grid[grid$transformation == "gaussian", ]
  by_cells <- env$.simCompareThesisThresholdPlot(one, dens, "low", cols = "n_cell")
  b <- ggplot2::ggplot_build(by_cells)
  expect_setequal(levels(by_cells$layers[[2]]$data$column), c("5000", "1e+05"))
  # 2 response probabilities x 2 cell counts, with the missing one left blank
  expect_equal(nrow(b$layout$layout), 4L)

  per_cell <- env$.simCompareThesisThresholdPlot(grid[grid$n_cell == 1e5, ], dens, "low")
  # 100k cells: both response probabilities for both transformations
  expect_equal(nrow(ggplot2::ggplot_build(per_cell)$layout$layout), 4L)
  per_small <- env$.simCompareThesisThresholdPlot(grid[grid$n_cell == 5e3, ], dens, "low")
  expect_setequal(as.character(unique(per_small$layers[[2]]$data$prob_response)), "0.01")
  expect_equal(env$.simCompareThesisThresholdHeight(1), 79)
  expect_equal(env$.simCompareThesisThresholdHeight(10), 215)
})

test_that("percentile figures show background-subtracted frequency and F1 with error bars", {
  env <- .thesisCompareEnv()
  g <- .thesisGrid()
  # one response probability present at one cell count only, as 0.01% is
  g <- rbind(g, transform(g[g$prob_response == 0.2 & g$n_cell == 1e5, ], prob_response = 1e-4))
  freq <- env$.simCompareThesisPctPlot(g, "low", "frequency", "prob")
  b <- ggplot2::ggplot_build(freq)
  expect_identical(freq$labels$y, "Background-subtracted frequency (%)")
  # percent of cells: the median 0.1 relative error at 20% response is 22%
  expect_equal(unique(freq$data$mid[freq$data$prob_response == 0.2]), 22)
  bars <- b$data[[2]]
  # the single-cell-count row still gets error bars
  expect_true(all(c("ymin", "ymax") %in% names(bars)))
  expect_true(any(abs(bars$ymin - 100 * 1e-4 * 0.5) < 1e-9))
  cells <- env$.simCompareThesisPctPlot(g, "low", "frequency", "cells")
  expect_equal(unique(cells$data$mid), 1.1)
  expect_no_error(ggplot2::ggplot_build(cells))
  f1 <- env$.simCompareThesisPctPlot(g, "low", "f1", "prob")
  expect_equal(f1$data$hi, pmax(f1$data$mid, f1$data$lo))
  expect_no_error(ggplot2::ggplot_build(f1))
  expect_false("tailgate_default" %in% as.character(f1$data$method))
})


test_that("frequency by cells has even probability spacing and a shared linear ratio scale", {
  env <- .thesisCompareEnv()
  g <- .thesisGrid()
  g <- rbind(g, transform(g[g$prob_response == 0.2, ], prob_response = 0.01))
  g$q10 <- -2; g$q90 <- 5
  p <- env$.simCompareThesisPctPlot(g, "low", "frequency", "cells")
  expect_identical(levels(p$data$row), c("1,000 cells", "100,000 cells"))
  expect_equal(sort(unique(p$data$x)), 1:3)
  expect_equal(p$data$x, match(p$data$prob_response, c(0.002, 0.01, 0.2)))
  b <- ggplot2::ggplot_build(p)
  for (panel in b$layout$panel_params) {
    expect_equal(panel$y.range, c(0, 4))
    expect_equal(panel$y$get_breaks(), 0:4)
    expect_equal(panel$y$get_labels(), paste0(0:4, "x"))
  }
  expect_equal(unique(p$data$lo), -1)
  expect_equal(unique(p$data$hi), 6)
  bands <- env$.simCompareThesisPctPlot(g, "low", "frequency", "cells", intervals = "bands")
  expect_equal(bands$data, p$data)
  expect_true(any(vapply(bands$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))))
  expect_false(any(vapply(bands$layers, function(l) inherits(l$geom, "GeomErrorbar"), logical(1))))
  expect_no_error(ggplot2::ggplot_build(bands))
})

test_that("percentile figures default to thin bands with room above the truth", {
  env <- .thesisCompareEnv()
  g <- .thesisGrid()
  g <- rbind(g, transform(g[g$prob_response == 0.2 & g$n_cell == 1e5, ], prob_response = 1e-4))
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  bars <- env$.simCompareThesisPctPlot(g, "low", "frequency", "prob", intervals = "bars")
  expect_true("GeomErrorbar" %in% geoms(bars))
  for (measure in c("frequency", "f1")) for (rows in c("prob", "cells")) {
    p <- env$.simCompareThesisPctPlot(g, "low", measure, rows)
    expect_true("GeomRibbon" %in% geoms(p))
    expect_false("GeomErrorbar" %in% geoms(p))
    # median plus the two thin percentile lines
    expect_equal(sum(geoms(p) == "GeomLine"), 3L)
    lines <- p$layers[geoms(p) == "GeomLine"]
    expect_equal(lines[[1]]$aes_params$alpha, env$.analysis_band_line_styles$thin$alpha)
    expect_no_error(ggplot2::ggplot_build(p))
  }
  p <- env$.simCompareThesisPctPlot(g, "low", "frequency", "prob")
  # same values as the error-bar figure: actual frequencies in percent
  expect_equal(p$data, bars$data)
  # the response simulated at one cell count gets bars, not an empty ribbon
  expect_true("GeomLinerange" %in% geoms(p))
  # each row runs from 0 to at least 1.5 times its truth
  b <- ggplot2::ggplot_build(p)
  truth <- c(20, 0.2, 0.01)
  for (i in seq_len(nrow(b$layout$layout))) {
    row <- b$layout$layout$ROW[[i]]
    yr <- b$layout$panel_params[[i]]$y.range
    expect_lte(yr[1], 0)
    expect_gte(yr[2], 1.5 * truth[[row]])
  }
  qmd <- paste(readLines(file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")), collapse = "\n")
  expect_false(grepl("bands-dotted|frequency-by-cells-bands|band_lines", qmd))
})

test_that("Analysis 8 band figures put the shift on a free x axis per transformation", {
  env <- .thesisCompareEnv()
  g <- expand.grid(method = c("stimgate", "fbeta", "tailgate", "tailgate_default"),
    transformation = c("gaussian", "skew", "gamma"), prob_response = c(0.2, 0.002),
    n_cell = c(5e3, 1e5), mean_pos_setting = c("low", "high"),
    mismatch_type = "mean_shift_all", mismatch_val = c(0, 0.1, 0.35), stringsAsFactors = FALSE)
  g$mismatch_val[g$transformation == "gamma"] <- g$mismatch_val[g$transformation == "gamma"] / 2
  g$q10 <- -0.5; g$median <- 0; g$q90 <- 1
  g$f1_q10 <- 0.2; g$f1_median <- 0.5; g$f1_q90 <- 0.9
  freq <- env$.simCompareThesisMismatchPlot(g, "low", "mean_shift_all", 5e3, "frequency")
  expect_false("tailgate_default" %in% as.character(freq$data$method))
  expect_setequal(unique(freq$data$x[freq$data$transformation == "Gamma"]), c(0, 0.05, 0.175))
  # actual frequency in percent, truth from the response probability
  expect_equal(unique(freq$data$mid[freq$data$row == "Response 20%"]), 20)
  expect_equal(unique(freq$data$hi[freq$data$row == "Response 0.2%"]), 0.4)
  expect_identical(levels(freq$data$row), c("Response 20%", "Response 0.2%"))
  expect_s3_class(freq$facet, "FacetWrap")
  expect_identical(freq$facet$params$free, list(x = TRUE, y = TRUE))
  # every panel runs from 0 to at least 1.5 times its own truth
  bf <- ggplot2::ggplot_build(freq)
  rows <- as.character(bf$layout$layout$.row_label)
  for (i in seq_along(rows)) {
    yr <- bf$layout$panel_params[[i]]$y.range
    expect_gte(yr[2], 1.5 * if (rows[[i]] == "Response 20%") 20 else 0.2)
    expect_lte(yr[1], 0)
  }
  b <- ggplot2::ggplot_build(freq)
  ranges <- vapply(b$layout$panel_params, function(p) p$x.range[2], numeric(1))
  expect_true(length(unique(round(ranges, 6))) > 1L)
  f1 <- env$.simCompareThesisMismatchPlot(g, "high", "mean_shift_all", 1e5, "f1")
  expect_s3_class(f1$facet, "FacetGrid")
  expect_identical(f1$facet$params$free, list(x = TRUE, y = FALSE))
  expect_equal(unique(f1$data$mid), 0.5)
  expect_no_error(ggplot2::ggplot_build(f1))
  expect_equal(env$.simCompareThesisBandHeight(4), 208)
  qmd <- paste(readLines(file.path(root_dir, "analysis", "8-sim-compare-freq_bs-batch.qmd")), collapse = "\n")
  expect_match(qmd, 'thesis_dir(thesis_types[[type]], "bands", measure)', fixed = TRUE)
  expect_match(qmd, "nrd0 bandwidth", fixed = TRUE)
  expect_false(grepl("band_lines|bands-dotted", qmd))
})
