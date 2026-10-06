.qmd7_presentation_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R",
    "sim-compare-qmd7-presentation.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("QMD 7 appends 1 percent without changing legacy scenario IDs or seeds", {
  env <- .qmd7_presentation_env()
  keys <- c("transformation", "prob_response", "n_cell", "mean_pos_setting",
    "mean_pos", "sample_perturbation_sd", "condition_perturbation_sd",
    "cluster_perturbation_sd", "background_relative_to_response",
    "n_cell_uns_relative_to_stim")
  old <- tidyr::expand_grid(
    transformation = c("skew", "gaussian", "gamma"), mean_pos_setting = c("low", "high"),
    prob_response = c(0.000005, 0.00002, 0.0001, 0.0005, 0.002, 0.2),
    n_cell = c(1000, 5000, 20000, 100000), condition_perturbation_sd = c(0, 0.5)
  ) |>
    dplyr::mutate(mean_pos = 5, sample_perturbation_sd = 0,
      cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
      n_cell_uns_relative_to_stim = 1) |>
    dplyr::filter(!(.data$n_cell * .data$prob_response < 5 & .data$prob_response < 0.04))
  legacy <- old |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::mutate(base_scenario_id = dplyr::cur_group_id(),
      sim_seed = as.integer(12345L + .data$base_scenario_id - 1L)) |>
    dplyr::ungroup() |>
    dplyr::mutate(sim_id = dplyr::row_number())
  added <- old |>
    dplyr::select(-"prob_response") |>
    dplyr::distinct() |>
    dplyr::mutate(prob_response = 0.01)
  # Put new rows first to ensure identity preservation does not depend on
  # their position in expand_grid()'s response vector.
  actual <- env$.simCompareQmd7SeedGrid(dplyr::bind_rows(added, old), 12345L)
  expect_equal(dplyr::select(dplyr::filter(actual, .data$prob_response != 0.01),
    dplyr::all_of(names(legacy))), legacy)
  new <- dplyr::filter(actual, .data$prob_response == 0.01)
  expect_true(all(new$sim_id > max(legacy$sim_id)))
  expect_true(all(new$sim_seed > max(legacy$sim_seed)))
  expect_equal(anyDuplicated(actual$sim_id), 0L)
  expect_equal(anyDuplicated(actual$sim_seed), 0L)
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  qmd <- paste(readLines(file.path(root, "analysis", "7-sim-compare-freq_bs.qmd")), collapse = "\n")
  expect_match(qmd, ")[c(2, 4, 6, 8, 10, 12, 16)]", fixed = TRUE)
  expect_match(qmd, "analysis_grid_spec = analysis_grid_spec", fixed = TRUE)
})

test_that("frequency plots separate methods and preserve nonpositive estimates", {
  skip_if_not_installed("ggforce")
  withr::local_preserve_seed()
  env <- .qmd7_presentation_env()
  data <- tidyr::expand_grid(prob_response = c(0.0001, 0.01, 0.2),
    method = c("stimgate", "tailgate", "fbeta"), sample = 1:6) |>
    dplyr::mutate(transformation = "gaussian", mean_pos_setting = "high",
      condition_perturbation_sd = 0, n_cell = 100000,
      propRespTruth = .data$prob_response,
      propRespEst = dplyr::case_when(.data$sample == 1L ~ 0,
        .data$sample == 2L ~ -0.01, TRUE ~ .data$prob_response * .data$sample / 4))
  plot <- env$.simComparePlotEstVsTruth(data, maxwidth = 0.1, lower_limit = 0.00001)
  expect_equal(plot$data$response_position, match(data$prob_response, c(0.0001, 0.01, 0.2)))
  expect_equal(sort(unique(plot$data$method_position - plot$data$response_position)), c(-0.25, 0, 0.25))
  expect_true(all(plot$data$estimate_shown > 0))
  expect_equal(plot$data$estimate_shown[data$propRespEst <= 0], rep(0.00001, 18))
  expect_match(plot$labels$caption, "zero and negative")
  built <- ggplot2::ggplot_build(plot)
  # Sina points first; truth segments are drawn on top.
  expect_equal(nrow(built$data[[1]]), nrow(data))
  expect_true(all(is.finite(built$data[[1]]$y)))
  expect_equal(built$data[[2]]$y, log10(c(0.0001, 0.01, 0.2)))
  expect_equal(built$data[[2]]$xend - built$data[[2]]$x, rep(0.8, 3))
  without <- env$.simComparePlotEstVsTruth(dplyr::filter(data, .data$method != "tailgate"))
  expect_equal(sort(unique(without$data$method_position - without$data$response_position)), c(-0.25, 0.25))
  expect_error(env$.simComparePlotEstVsTruth(data, maxwidth = 0.3), "method spacing")
  labelled <- env$.simCompareQmd7Subtitle(plot, data)
  expect_match(labelled$labels$subtitle, "All methods")
  expect_match(labelled$labels$subtitle, "Mean position: high", fixed = TRUE)
  labelled <- env$.simCompareQmd7Subtitle(without, without$data)
  expect_match(labelled$labels$subtitle, "Without Tailgate")
  expect_no_error(ggplot2::ggplotGrob(labelled))
})

test_that("provenance table counts errors in the all-sample denominator", {
  env <- .qmd7_presentation_env()
  summary <- tibble::tibble(method = "tailgate", transformation = "gamma",
    mean_pos_setting = "high", condition_perturbation_sd = c(0, 0.5, 0.5),
    n = c(400L, 200L, 200L), n_threshold_fallback = c(0L, 1L, 1L),
    n_run_error = c(0L, 70L, 79L), n_no_cutpoint = c(0L, 1L, 1L))
  table <- env$.simCompareQmd7FallbackTable(summary)
  expect_equal(table$n_samples, c(400L, 400L))
  expect_equal(table$condition_perturbation_sd, c(0, 0.5))
  expect_equal(table$n_run_error[[2]], paste0("149 / 400 (", env$.analysis_label_percent(149 / 400), ")"))
  expect_equal(table$n_threshold_fallback[[2]], paste0("2 / 400 (", env$.analysis_label_percent(2 / 400), ")"))
  expect_equal(table$n_no_cutpoint[[2]], table$n_threshold_fallback[[2]])
})

test_that("threshold panels retain linear histograms and square root only references", {
  env <- .qmd7_presentation_env()
  data <- tidyr::expand_grid(transformation = c("gaussian", "gamma"),
    prob_response = c(0.01, 0.2), approach = c("stimgate", "fbeta"), sample = 1:5) |>
    dplyr::mutate(threshold = ifelse(.data$transformation == "gamma", .data$sample / 10, .data$sample * 3),
      mean_pos = 5, mean_pos_setting = "high")
  densities <- tidyr::expand_grid(transformation = c("gaussian", "gamma"),
    prob_response = c(0.01, 0.2), condition = c("stimulated", "unstimulated"),
    expression = c(0, 1, 2)) |>
    dplyr::mutate(mean_pos = 5, density = (.data$expression + 1)^2 / 100)
  plot <- env$.simComparePlotThresholdDensity(data, densities, reference_scale = 2)
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1]]$y, sqrt(plot$layers[[1]]$data$density) * 2)
  expect_equal(built$data[[2]]$y, built$data[[2]]$density)
  expect_equal(length(unique(built$layout$layout$SCALE_X)), 4L)
  expect_equal(length(unique(built$layout$layout$SCALE_Y)), 4L)
  expect_equal(built$data[[1]]$colour, rep("gray25", nrow(densities)))
  expect_match(plot$labels$y, "Threshold density", fixed = TRUE)
  expect_no_error(ggplot2::ggplotGrob(plot))
  expect_error(env$.simComparePlotThresholdDensity(data, densities, reference_scale = 0), "positive")
})

test_that("QMD 7 reference tubes retain unstimulated positive cells and preserve RNG", {
  skip_if_not_installed("simcyto")
  withr::local_preserve_seed()
  env <- .qmd7_presentation_env()
  panels <- tibble::tibble(transformation = "gaussian", mean_pos = 8,
    prob_response = 0.3, sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2)
  settings <- list(probExact = TRUE, covEvMin = 1.5, covEvMax = 1.5)
  set.seed(13)
  before <- .Random.seed
  densities <- env$.simBandwidthThresholdDensities(dplyr::bind_rows(panels, panels),
    settings, n_cell = 1000, density_n = 256, unstimulated_negative_only = FALSE)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(densities), 512L)
  uns <- dplyr::filter(densities, .data$condition == "unstimulated")
  mass <- sum(uns$density[uns$expression > 4]) * diff(uns$expression)[[1]]
  expect_gt(mass, 0.04)
  expect_lt(mass, 0.09)
  expect_identical(env$.simBandwidthThresholdDensities(panels,
    settings, n_cell = 1000, density_n = 256, unstimulated_negative_only = FALSE), densities)
  expect_identical(.Random.seed, before)
})
