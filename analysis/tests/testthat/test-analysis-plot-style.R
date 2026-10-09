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
  expect_identical(env$.analysis_method_colours,
    c(stimgate = "#E69F00", tailgate = "#009E73", fbeta = "#0072B2",
      tailgate_default = "#CC79A7"))
  expect_identical(unname(env$.analysis_method_labels),
    c("StimGate", "Tailgate", "F-beta", "Tailgate (default settings)"))
  for (encoding in list(env$.analysis_method_labels, env$.analysis_method_shapes,
    env$.analysis_method_linetypes)) {
    expect_setequal(names(encoding), names(env$.analysis_method_colours))
  }
  expect_s3_class(env$.analysis_scale_method(), "ScaleDiscrete")
})


test_that("axis floors reach every free panel and retain larger data and intervals", {
  env <- .plot_style_env()
  data <- data.frame(panel = c("tiny", "large"), x = c(1, 10),
    fdp = c(0.002, 0.3), lower = c(0.001, -0.05), upper = c(0.003, 0.6))
  p <- ggplot2::ggplot(data, ggplot2::aes(x, fdp, colour = panel)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0) +
    ggplot2::facet_wrap(~panel, scales = "free")
  before <- ggplot2::ggplot_build(p)
  built <- ggplot2::ggplot_build(p + env$.analysis_y_floor())
  for (panel in built$layout$panel_params) {
    expect_lte(panel$y.range[[1]], 0)
    expect_gte(panel$y.range[[2]], 0.1)
  }
  large <- which(built$layout$layout$panel == "large")
  expect_lte(built$layout$panel_params[[large]]$y.range[[1]], -0.05)
  expect_gte(built$layout$panel_params[[large]]$y.range[[2]], 0.6)
  expect_equal(built$data[1:2], before$data[1:2])
  expect_equal(lapply(built$layout$panel_scales_x, function(s) s$dimension()),
    lapply(before$layout$panel_scales_x, function(s) s$dimension()))
})

test_that("signed floors train in data space and preserve the signed transform", {
  env <- .plot_style_env()
  source(file.path(testthat::test_path(),
    "../../../scripts/r/sim-bandwidth-analysis-plot.R"), local = env)
  trans <- env$.simBandwidthSignedErrorTrans()
  data <- data.frame(panel = c("tiny", "negative"), x = 1, y = c(0.002, -2))
  p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::scale_y_continuous(transform = trans) +
    env$.analysis_y_floor(c(-0.1, 0.1))
  built <- ggplot2::ggplot_build(p)
  for (panel in built$layout$panel_params) {
    expect_lte(panel$y.range[[1]], trans$transform(-0.1))
    expect_gte(panel$y.range[[2]], trans$transform(0.1))
  }
  expect_equal(built$data[[1]]$y, trans$transform(data$y))
  expect_equal(trans$transform(c(-2, -1, 0, 1)), c(-2, -1, 0, 1))
})

test_that("fitted figure height follows facet rows for print and save", {
  env <- .plot_style_env()
  data <- data.frame(panel = letters[1:6], x = 1, y = 0.002)
  p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_equal(env$.analysis_facet_height(p + ggplot2::facet_wrap(~panel, ncol = 3)), 15)
  p <- p + ggplot2::facet_wrap(~panel, ncol = 2)
  expect_equal(env$.analysis_facet_height(p), 20)
  env$.analysis_print_fig <- function(plot, mcse_mode, width, height) {
    env$printed <- c(width = width, height = height)
  }
  env$.analysis_save_fig <- function(plot, path, mcse_mode, width, height, allow_tall) {
    env$saved <- c(width = width, height = height)
    env$allow_tall <- allow_tall
  }
  env$.analysis_print_save_fig(p, "unused.pdf", height = 6, fit_panels = TRUE)
  expect_equal(env$printed, c(width = 16, height = 20))
  expect_identical(env$saved, env$printed)
  expect_true(env$allow_tall)
})

test_that("fitted HTML figures use their own device dimensions", {
  skip_if_not_installed("png")
  skip_if_not_installed("base64enc")
  env <- .plot_style_env()
  env$plot <- ggplot2::ggplot(data.frame(x = 1, y = 0.002), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  rendered <- knitr::knit(text = c(
    "```{r fitted, echo=FALSE, results='asis', fig.width=7, fig.height=3}",
    ".analysis_print_fig(plot, width = 16, height = 20)",
    "```"
  ), envir = env, quiet = TRUE)
  # An embedded Markdown image, so Quarto cannot drop or unlink its file.
  uri <- regmatches(rendered, regexpr("data:image/png;base64,[A-Za-z0-9+/=]+", rendered))
  expect_length(uri, 1L)
  path <- withr::local_tempfile(fileext = ".png")
  writeBin(base64enc::base64decode(sub("^data:image/png;base64,", "", uri)), path)
  dims <- dim(png::readPNG(path))
  expect_equal(dims[[1]] / dims[[2]], 20 / 16, tolerance = 0.01)
})
