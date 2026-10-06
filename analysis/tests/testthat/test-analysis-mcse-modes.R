.mcse_modes_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.mcse_modes_plot <- function(env) {
  data <- tibble::tibble(x = 1:2, y = c(0.1, 0.2),
    lower = c(0.05, 0.1), upper = c(0.2, 0.4))
  ggplot2::ggplot(data, ggplot2::aes(x, y)) + ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.1) +
    env$.analysis_mcse_errorbar(dplyr::mutate(data, lower = -2, upper = 3), width = 2)
}

test_that("MCSE modes default to on and retain Boolean compatibility", {
  env <- .mcse_modes_env()
  expect_identical(env$.analysis_mcse_mode(), "on")
  for (value in list(TRUE, "true", "1", "on")) {
    expect_identical(env$.analysis_mcse_mode(value), "on")
  }
  for (value in list(FALSE, "false", "0", "off")) {
    expect_identical(env$.analysis_mcse_mode(value), "off")
  }
  for (value in list(NA, character(), c("on", "off"), "nonsense", " BOTH ")) {
    expect_error(env$.analysis_mcse_mode(value), "off or on", fixed = TRUE)
  }
})

test_that("off variants remove only MC intervals without mutating the source plot", {
  env <- .mcse_modes_env()
  plot <- .mcse_modes_plot(env)
  variants <- c(env$.analysis_mcse_plot_variants(plot, "off"),
    env$.analysis_mcse_plot_variants(plot, "on"))
  expect_named(variants, c("off", "on"))
  expect_length(plot$layers, 4L)
  expect_length(variants$on$layers, 4L)
  expect_length(variants$off$layers, 4L)
  expect_true(inherits(variants$off$layers[[4]]$geom, "GeomBlank"))
  expect_true(inherits(variants$off$layers[[3]]$geom, "GeomErrorbar"))
  expect_false(isTRUE(attr(variants$off$layers[[3]], "analysis_mcse")))
  before <- ggplot2::ggplot_build(plot)
  off <- ggplot2::ggplot_build(variants$off)
  for (i in 1:4) expect_equal(off$data[[i]], before$data[[i]])
  expect_equal(off$layout$panel_params[[1]]$y.range,
    before$layout$panel_params[[1]]$y.range)
  expect_equal(off$layout$panel_params[[1]]$x.range,
    before$layout$panel_params[[1]]$x.range)
  expect_no_error(ggplot2::ggplotGrob(variants$off))
  expect_length(env$.analysis_mcse_plot_variants(plot, FALSE), 1L)
  expect_length(env$.analysis_mcse_plot_variants(plot, TRUE), 1L)
})

test_that("each render saves and prints one mode and preserves its sibling", {
  env <- .mcse_modes_env()
  plot <- .mcse_modes_plot(env)
  dir <- withr::local_tempdir()
  captured <- list()
  env$print <- function(x) captured[[length(captured) + 1L]] <<- x
  for (mode in c("off", "on")) {
    paths <- env$.analysis_save_fig(plot, file.path(dir, "fixture.png"),
      height = 3, width = 4, mcse_mode = mode)
    expect_length(paths, 1L)
    expect_identical(normalizePath(paths, winslash = "/"),
      normalizePath(file.path(dir, paste0("mcse_", mode), "fixture.png"), winslash = "/"))
    output <- utils::capture.output(env$.analysis_print_fig(plot, mode))
    expect_true(any(grepl(paste0("Monte Carlo error bars: ", mode), output, fixed = TRUE)))
  }
  expect_true(all(file.exists(file.path(dir, c("mcse_off", "mcse_on"), "fixture.png"))))
  expect_length(captured, 2L)
  expect_true(inherits(captured[[1]]$layers[[4]]$geom, "GeomBlank"))
  expect_true(inherits(captured[[2]]$layers[[4]]$geom, "GeomErrorbar"))
  expect_length(env$.analysis_mcse_plot_variants(plot), 1L)
})

test_that("performance QMDs default on and explicitly forward plot modes", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  for (stem in c("2a-", "2b-", "3-", "4-", "7-", "8-")) {
    file <- list.files(file.path(root, "analysis"), pattern = paste0("^", stem), full.names = TRUE)
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(text, '  show_mcse: "on"', fixed = TRUE)
    expect_match(text, 'show_mcse <- mcse_mode != "off"', fixed = TRUE)
    expect_match(text, "mcse_mode = mcse_mode", fixed = TRUE)
  }
})

test_that("method outcome validation accepts multiple biological bootstrap families", {
  env <- .mcse_modes_env()
  raw <- tidyr::expand_grid(method = c("stimgate", "fbeta"), sim_seed = c(11L, 22L)) |>
    dplyr::mutate(iter = 1L, sample = "1", nTruePos = 1L, nFalsePos = 0L,
      nFalseNeg = 0L, nTrueNeg = 9L, thresholdFallbackUsed = FALSE,
      thresholdOrigin = "calculated", error = NA_character_)
  raw$error[raw$method == "fbeta" & raw$sim_seed == 22L] <- "failed"
  counts <- env$.simCompareMethodOutcomeCounts(raw)
  expect_equal(counts$n, c(2L, 2L))
  expect_equal(counts$n_valid[counts$method == "stimgate"], 2L)
  expect_equal(counts$n_run_error[counts$method == "fbeta"], 1L)
  expect_equal(counts$n_valid[counts$method == "fbeta"], 1L)
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  text <- paste(readLines(file.path(root, "analysis/8-sim-compare-freq_bs-batch.qmd")), collapse = "\n")
  expect_match(text, ".simCompareMethodOutcomeCounts(compare_raw)", fixed = TRUE)
})

test_that("comparison loops print before saving and save ratios without printing", {
  env <- .mcse_modes_env()
  constructed <- 0L
  calls <- character()
  printed <- list()
  saved <- list()
  env$print <- function(x) {
    calls <<- c(calls, "print")
    printed[[length(printed) + 1L]] <<- x
  }
  env$.analysis_save_fig <- function(plot, path, height, allow_tall, mcse_mode) {
    calls <<- c(calls, "save")
    saved[[length(saved) + 1L]] <<- list(path = path, mode = mcse_mode)
  }
  data <- tibble::tibble(method = c("stimgate", "tailgate", "fbeta"),
    mean_pos_setting = "high")
  utils::capture.output(env$.simCompareFigureLoop(data,
    make_plot = function(d) {
      constructed <<- constructed + 1L
      .mcse_modes_plot(env) + env$.simBandwidthSignedErrorLayers()
    }, dir = "signed_error_fixture", file_fn = function(pos, extra) "fixture.png",
    height = 3, level = 3L, ratio_twins = TRUE, mcse_mode = "on"))
  # Each method subset prints its original, then saves the original and ratio.
  expect_identical(calls, rep(c("print", "save", "save"), 2L))
  expect_equal(constructed, 2L)
  expect_length(printed, 2L)
  expect_length(saved, 4L)
  expect_true(all(vapply(saved, function(x) identical(x$mode, "on"), logical(1))))
  expect_true(any(vapply(saved, function(x) grepl("ratio_fixture", x$path), logical(1))))
})

test_that("shared asis helper prints and separates the figure before saving", {
  env <- .mcse_modes_env()
  calls <- character()
  env$print <- function(x) calls <<- c(calls, "print")
  env$cat <- function(...) calls <<- c(calls, "separator")
  env$.analysis_save_fig <- function(plot, path, ..., mcse_mode = NULL) {
    calls <<- c(calls, "save")
  }
  env$.analysis_print_save_fig(ggplot2::ggplot(), "fixture.pdf", height = 3)
  expect_identical(calls, c("print", "separator", "save"))
})

test_that("MC error bars have visible caps and strokes without changing estimates", {
  env <- .mcse_modes_env()
  data <- tibble::tibble(x = 1, y = 0.2, lower = 0.1, upper = 0.3)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x, y)) + ggplot2::geom_point() +
    env$.analysis_mcse_errorbar(data)
  expect_equal(plot$layers[[2]]$aes_params$alpha, 0.6)
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1]]$y, data$y)
  expect_true(all(built$data[[2]]$xmax > built$data[[2]]$xmin))
  expect_equal(built$data[[2]]$ymin, data$lower)
  expect_equal(built$data[[2]]$ymax, data$upper)
})

test_that("analysis HTML keeps retina size bounded and preserves sibling MC files", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  files <- list.files(file.path(root, "analysis"), pattern = "\\.qmd$", full.names = TRUE)
  for (file in files) {
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(text, "knitr:\n  opts_chunk:\n    fig.retina: 1", fixed = TRUE, info = basename(file))
    expect_match(text, "    embed-resources: true", fixed = TRUE, info = basename(file))
    expect_false(grepl("unlink(path, recursive = TRUE)", text, fixed = TRUE), info = basename(file))
  }
})
