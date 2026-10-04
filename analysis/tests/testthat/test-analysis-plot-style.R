.plot_style_env <- function() {
  env <- new.env(parent = globalenv())
  source(
    file.path(testthat::test_path(), "../../../scripts/r/analysis-plot-style.R"),
    local = env
  )
  env
}

test_that("transformations use one label set in the standard order", {
  env <- .plot_style_env()
  f <- env$.analysis_trans_factor(c("gamma", "gaussian", "skew", "other"))
  expect_identical(levels(f), c("Gaussian", "Skew", "Gamma", "other"))
  expect_identical(as.character(f), c("Gamma", "Gaussian", "Skew", "other"))
  expect_identical(
    env$.analysis_trans_order(c("gamma", "skew", "gaussian")),
    c("gaussian", "skew", "gamma")
  )
})

test_that("numbers and proportions avoid scientific notation", {
  env <- .plot_style_env()
  expect_identical(
    env$.analysis_label_number(c(0.0025, 0.05, 1.5, 1e5, 1000)),
    c("0.0025", "0.05", "1.5", "100,000", "1,000")
  )
  expect_identical(
    env$.analysis_label_percent(c(1e-4, 0.002, 0.05, 0.2)),
    c("0.01%", "0.2%", "5%", "20%")
  )
  expect_identical(env$.analysis_labeller_percent("p = ")(1e-4), "p = 0.01%")
})

test_that("figures save at the standard width with page-capped height", {
  skip_if_not_installed("png")
  env <- .plot_style_env()
  path <- file.path(withr::local_tempdir(), "fig", "p.png")
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    env$.analysis_theme()
  env$.analysis_save_fig(p, path, height = 40)
  expect_true(file.exists(path))
  dims <- dim(png::readPNG(path))
  expect_equal(dims[[1]] / dims[[2]], 22 / 16, tolerance = 0.01)
})

test_that("strips have a black border and no fill", {
  env <- .plot_style_env()
  strip <- env$.analysis_theme()$strip.background
  expect_true(is.na(strip$fill))
  expect_identical(strip$colour, "black")
})

test_that("loop headings are markdown headings at the given level", {
  env <- .plot_style_env()
  expect_output(env$.analysis_heading("Mean position: low", 4), "#### Mean position: low")
})

test_that("method colours cover every method with readable labels", {
  env <- .plot_style_env()
  expect_named(env$.analysis_method_colours, c("stimgate", "tailgate", "fbeta"))
  expect_identical(unname(env$.analysis_method_labels), c("StimGate", "Tailgate", "F-beta"))
  expect_s3_class(env$.analysis_scale_method(), "ScaleDiscrete")
})
