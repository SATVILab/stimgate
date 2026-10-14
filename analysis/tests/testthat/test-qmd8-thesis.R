root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

test_that("shifted reference densities carry each panel's stimulated-tube shift", {
  skip_if_not_installed("simcyto")
  withr::local_preserve_seed()
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c("analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R")) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  base <- tibble::tibble(
    transformation = "gaussian", mean_pos = 8, prob_response = 0.3,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    stim_mean_shift = 0, stim_mean_shift_clusters = NA_character_,
    stim_sd_multiplier = 1, stim_sd_multiplier_clusters = NA_character_
  )
  panels <- dplyr::bind_rows(
    base,
    dplyr::mutate(base, stim_mean_shift = 2),
    dplyr::mutate(base, stim_mean_shift = 2, stim_mean_shift_clusters = "gn")
  )
  settings <- list(probExact = TRUE, covEvMin = 1.5, covEvMax = 1.5)
  set.seed(5)
  before <- .Random.seed
  dens <- env$.simBandwidthThresholdDensities(panels, settings, n_cell = 2000,
    density_n = 256, unstimulated_negative_only = FALSE, rug_n = 50L, mismatch = TRUE)
  expect_identical(.Random.seed, before)
  stim_mean <- function(shift, clusters) {
    d <- dens[dens$condition == "stimulated" & dens$stim_mean_shift == shift &
      dplyr::coalesce(dens$stim_mean_shift_clusters, "") == clusters, ]
    sum(d$expression * d$density) / sum(d$density)
  }
  rug_mean <- function(shift, clusters) {
    d <- dens[dens$condition == "positive_rug" & dens$stim_mean_shift == shift &
      dplyr::coalesce(dens$stim_mean_shift_clusters, "") == clusters, ]
    mean(d$expression)
  }
  # Shifting every component moves the whole stimulated tube and its
  # responders; shifting the negatives moves the tube by their share only
  # (1 - 30% response - 6% background).
  expect_equal(stim_mean(2, "") - stim_mean(0, ""), 2, tolerance = 0.05)
  expect_equal(stim_mean(2, "gn") - stim_mean(0, ""), 2 * 0.64, tolerance = 0.05)
  expect_equal(rug_mean(2, "") - rug_mean(0, ""), 2, tolerance = 0.05)
  expect_equal(rug_mean(2, "gn"), rug_mean(0, ""))
  # The unstimulated tube never shifts (its grid spans both tubes, so compare
  # means rather than curves).
  uns_mean <- function(shift) {
    d <- dens[dens$condition == "unstimulated" & dens$stim_mean_shift == shift &
      is.na(dens$stim_mean_shift_clusters), ]
    sum(d$expression * d$density) / sum(d$density)
  }
  expect_equal(uns_mean(2), uns_mean(0), tolerance = 0.01)
  # Without mismatch the result is unchanged by the new argument.
  plain <- env$.simBandwidthThresholdDensities(base, settings, n_cell = 2000,
    density_n = 256, unstimulated_negative_only = FALSE)
  expect_false("stim_mean_shift" %in% names(plain))
})

test_that("Analysis 8 saves Analysis 7's thesis figures for every shift level", {
  qmd <- readLines(file.path(root_dir, "analysis", "8-sim-compare-freq_bs-batch.qmd"))
  chunk <- paste(qmd, collapse = "\n")
  expect_match(chunk, '"sim-compare-thesis.R"', fixed = TRUE)
  expect_match(chunk, 'thesis_parts <- c("results", "qmd8", fig_key[-1L])', fixed = TRUE)
  for (fn in c(".simCompareThesisF1Plot", ".simCompareThesisPctPlot",
    ".simCompareThesisThresholdPlot", "mismatch = TRUE")) {
    expect_match(chunk, fn, fixed = TRUE)
  }
})
