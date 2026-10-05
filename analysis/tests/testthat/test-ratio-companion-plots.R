test_that("ratio twins preserve signed geometry, intervals and original scales", {
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-plot-style.R", "sim-bandwidth-analysis-plot.R")) {
    source(file.path(testthat::test_path(), "../../../scripts/r", file), local = env)
  }
  data <- tibble::tibble(x = 1:6, y = c(-2, -1, 0, 1, 15, 15),
    lower = c(-2.1, -1.1, -0.1, 0.9, 14, 14),
    upper = c(-1.9, -0.9, 0.1, 1.1, 15, 15))
  original <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper)) +
    env$.simBandwidthSignedErrorLayers(capped = TRUE)
  before <- ggplot2::ggplot_build(original)
  twin <- env$.simBandwidthRatioPlot(original)
  after <- ggplot2::ggplot_build(twin)
  expect_equal(after$data, before$data)
  expect_equal(after$layout$panel_params[[1]]$y.range,
    before$layout$panel_params[[1]]$y.range)
  expect_equal(after$layout$panel_params[[1]]$y$get_breaks(),
    before$layout$panel_params[[1]]$y$get_breaks())
  expect_identical(twin$layers, original$layers)
  expect_identical(twin$data, original$data)
  labels <- twin$scales$get_scales("y")$labels(c(-2, -1, 0, 1, 15))
  expect_identical(labels, c("-1x", "0x", "1x (exact)", "2x", "\u2265 16x"))
  expect_identical(original$labels$y, "Relative error")
  expect_identical(original$scales$get_scales("y")$labels(c(-1, 0, 1, 15)),
    c("-100% (0x)", "0%", "+100% (2x)", "\u2265 +1500% (16x)"))
  expect_equal(ggplot2::ggplot_build(original)$data, before$data)
  uncapped <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_point() + env$.simBandwidthSignedErrorLayers()
  expect_identical(env$.simBandwidthRatioPlot(uncapped)$scales$get_scales("y")$labels(15),
    "16x")
})

test_that("ratio companions save beside originals and print their own heading", {
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-plot-style.R", "sim-bandwidth-analysis-plot.R")) {
    source(file.path(testthat::test_path(), "../../../scripts/r", file), local = env)
  }
  env$.analysis_save_fig <- function(plot, path, height, allow_tall, mcse_mode = NULL) {
    env$saved <- list(plot = plot, path = path, height = height, allow_tall = allow_tall)
  }
  env$.analysis_print_fig <- function(plot, mcse_mode = NULL) env$printed <- plot
  p <- ggplot2::ggplot() + env$.simBandwidthSignedErrorLayers()
  md <- utils::capture.output(env$.simBandwidthPrintRatioTwin(p,
    "output/fig/draft/signed_error_by_n_cell/sd_inflation/a.pdf", 22))
  expect_identical(env$saved$path,
    "output/fig/draft/ratio_by_n_cell/sd_inflation/a.pdf")
  expect_identical(env$saved$height, 22)
  expect_identical(env$saved$allow_tall, FALSE)
  expect_identical(env$printed, env$saved$plot)
  expect_true(any(grepl("Estimate / reference ratio", md, fixed = TRUE)))
})
