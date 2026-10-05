test_that("threshold reference densities preserve RNG and reveal stimulated response mass", {
  withr::local_preserve_seed()
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R", "sim-bandwidth-analysis-plot.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  panels <- tibble::tibble(
    transformation = "gaussian", mean_pos = 8, prob_response = 0.3,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell = 1000, threshold = 4
  )
  settings <- list(probExact = TRUE, covEvMin = 1.5, covEvMax = 1.5)
  set.seed(92)
  before <- .Random.seed
  densities <- env$.simBandwidthThresholdDensities(
    dplyr::bind_rows(panels, panels), settings, n_cell = 500, density_n = 256
  )
  expect_identical(.Random.seed, before)
  expect_s3_class(densities, "tbl_df")
  expect_setequal(unique(densities$condition), c("unstimulated", "stimulated"))
  expect_equal(nrow(densities), 512L)
  expect_true(all(is.finite(densities$density)))
  expect_true(all(densities$density >= 0))
  upper_mass <- function(condition) {
    d <- densities[densities$condition == condition, ]
    sum(d$density[d$expression > 4]) * diff(d$expression)[1]
  }
  expect_gt(upper_mass("stimulated"), 0.2)
  expect_lt(upper_mass("unstimulated"), 0.03)
  expect_identical(
    env$.simBandwidthThresholdDensities(panels, settings, n_cell = 500, density_n = 256),
    densities
  )
  expect_identical(.Random.seed, before)

  panels$transformation <- env$.analysis_trans_factor(panels$transformation)
  # Reuse the same curves in different cell-count panels without simulating again.
  panels <- dplyr::bind_rows(panels, dplyr::mutate(panels, n_cell = 2000))
  original <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      data = panels, ggplot2::aes(xintercept = threshold),
      linetype = "dashed", colour = "black", alpha = 0.75
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(prob_response, n_cell),
      cols = ggplot2::vars(transformation), scales = "free_x",
      labeller = ggplot2::labeller(prob_response = env$.analysis_labeller_percent())
    ) + env$.analysis_theme()
  plot <- env$.simBandwidthThresholdDensityPlot(original, panels, densities)
  expect_length(original$layers, 1L)
  expect_identical(plot$layers[[3]], original$layers[[1]])
  expect_s3_class(plot$layers[[1]]$geom, "GeomArea")
  expect_s3_class(plot$layers[[2]]$geom, "GeomLine")
  built <- ggplot2::ggplot_build(plot)
  expect_equal(length(unique(built$data[[1]]$PANEL)), 2L)
  expect_equal(length(unique(built$data[[1]]$group)), 2L)
  expect_equal(length(unique(built$data[[2]]$group)), 2L)
  expect_equal(built$data[[3]]$xintercept, c(4, 4))
  expect_no_error(ggplot2::ggplotGrob(plot))
})
