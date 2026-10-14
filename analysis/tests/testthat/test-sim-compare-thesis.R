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

test_that("Analysis 7 saves four thesis figures per separation to output-thesis results", {
  qmd <- paste(readLines(file.path(root_dir, "analysis", "7-sim-compare-freq_bs.qmd")), collapse = "\n")
  expect_match(qmd, '"output-thesis", c("results", fig_key[-1L])', fixed = TRUE)
  for (kind in c('stem("f1")', 'stem("f1-q10")', 'stem("frequency")', 'stem("thresholds")')) {
    expect_match(qmd, kind, fixed = TRUE)
  }
  expect_match(qmd, "rug_n = 100L", fixed = TRUE)
  expect_false(grepl("no_tailgate", qmd, fixed = TRUE))
})
