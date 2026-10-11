.null_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
    "sim-misc.R", "sim-compare-freq_bs.R", "sim-compare-null.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

.null_capture_background <- function(env, ...) {
  captured <- NULL
  env$.simCompareSimCytExperiment <- function(...) {
    captured <<- list(...)
    stop("captured")
  }
  withr::with_seed(1, expect_error(env$.simCompareFreqBs(
    nSample = 2, nIter = 1, nMarker = 1, nCondition = 2, nCluster = 2, nCellStim = 1000, meanPos = 5, transformation = "gaussian",
    samplePerturbationSd = 0, conditionPerturbationSd = 0, clusterPerturbationSd = 0,
    backgroundRelativeToResponse = 0.2, ncellUnsRelativeToStim = 1,
    biasUns = NULL, bw = NULL, ...
  ), "captured"))
  captured
}

test_that("an absolute background replaces the response-relative one", {
  env <- .null_env()
  relative <- .null_capture_background(env, probResponse = 0.05)
  expect_equal(relative$probVecUns, c(0.99, 0.01))
  expect_equal(relative$probResponseVecByStimCondition, list(c(-0.05, 0.05)))
  for (bg in c(0, 1e-4, 1e-2)) {
    null <- .null_capture_background(env, probResponse = 0, probBackground = bg)
    expect_equal(null$probVecUns, c(1 - bg, bg))
    expect_equal(null$probResponseVecByStimCondition, list(c(0, 0)))
  }
  # NA (no prob_background column value) keeps the relative background.
  na <- .null_capture_background(env, probResponse = 0.05, probBackground = NA_real_)
  expect_equal(na$probVecUns, c(0.99, 0.01))
})

test_that("null seeds depend only on the biological scenario", {
  env <- .null_env()
  grid <- tidyr::expand_grid(
    transformation = c("gaussian", "gamma"), mean_pos_setting = c("low", "high"),
    prob_background = env$.simCompareNullBackgrounds, n_cell = c(5000, 20000),
    bw_mtd = c("nrd0", "hpi1")
  ) |>
    dplyr::mutate(mean_pos = 5, sample_perturbation_sd = 0, condition_perturbation_sd = 0,
      cluster_perturbation_sd = 0, n_cell_uns_relative_to_stim = 1)
  seeded <- env$.simCompareNullSeedGrid(grid, 100L)
  expect_equal(seeded$sim_id, seq_len(nrow(grid)))
  # Rows differing only by a method setting share a seed; scenarios do not.
  per_scenario <- seeded |>
    dplyr::group_by(.data$transformation, .data$mean_pos_setting,
      .data$prob_background, .data$n_cell) |>
    dplyr::summarise(n_seed = dplyr::n_distinct(.data$sim_seed), seed = .data$sim_seed[[1]],
      .groups = "drop")
  expect_true(all(per_scenario$n_seed == 1L))
  expect_equal(anyDuplicated(per_scenario$seed), 0L)
  # Filtering the grid afterwards does not change any row's seed.
  kept <- seeded[seeded$n_cell == 20000, ]
  expect_equal(kept$sim_seed, seeded$sim_seed[seeded$n_cell == 20000])
})

test_that("null summaries pool samples and keep missing estimates missing", {
  env <- .null_env()
  raw <- tidyr::expand_grid(iter = 1:5, sample = 1:4, method = c("stimgate", "fbeta")) |>
    dplyr::mutate(
      sim_seed = 7L, transformation = "gaussian", mean_pos_setting = "high", mean_pos = 8,
      prob_background = 0.01, n_cell = 1000, propRespTruth = 0,
      propRespEst = (.data$iter + .data$sample - 5) / 1e4,
      locGenerated = .data$sample != 1, locGeneratedDirect = .data$sample > 2
    )
  raw$propRespEst[raw$method == "fbeta" & raw$iter == 1] <- NA_real_
  errors <- env$.simCompareNullErrors(raw)
  expect_equal(errors$error_pp, 100 * raw$propRespEst)
  cols <- c("transformation", "mean_pos_setting", "prob_background", "n_cell", "method")
  out <- env$.simCompareNullSummary(errors, cols, mcse = FALSE)
  sg <- errors$error_pp[errors$method == "stimgate"]
  fb <- errors$error_pp[errors$method == "fbeta"]
  s <- out[out$method == "stimgate", ]
  f <- out[out$method == "fbeta", ]
  expect_equal(s$error_median, stats::quantile(sg, 0.5, names = FALSE))
  expect_equal(s$error_p90, stats::quantile(sg, 0.9, names = FALSE))
  expect_equal(s$abs_error_p90, stats::quantile(abs(sg), 0.9, names = FALSE))
  expect_equal(s$share_positive, mean(sg > 0))
  expect_equal(f$n_sample, 20L)
  expect_equal(f$n_finite, 16L)
  expect_equal(f$error_p10, stats::quantile(fb, 0.1, na.rm = TRUE, names = FALSE))
  expect_false(any(startsWith(names(out), ".boot_")))

  with_mcse <- env$.simCompareNullSummary(errors, cols, mcse = TRUE)
  expect_equal(with_mcse$error_median, out$error_median)
  # Intervals need five datasets with an estimate: F-beta has four here.
  sg_mcse <- with_mcse[with_mcse$method == "stimgate", ]
  expect_true(sg_mcse$interval_available_error_median)
  expect_true(sg_mcse$error_median_lower <= sg_mcse$error_median)
  expect_false(with_mcse$interval_available_error_median[with_mcse$method == "fbeta"])

  prov <- env$.simCompareNullProvenance(errors, setdiff(cols, "method"))
  expect_equal(c(prov$n_sample, prov$n_direct, prov$n_shared, prov$n_fallback), c(20, 10, 5, 5))

  p <- env$.simCompareNullPlot(with_mcse, "error_median", "Median error (pp)")
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(as.character(ggplot2::ggplot_build(p)$plot$data$method)),
    c("stimgate", "fbeta"))
  expect_s3_class(env$.simCompareNullPlot(out, "abs_error_p90", "x",
    zero_line = FALSE, by_separation = FALSE), "ggplot")
})

test_that("Analysis 7b reuses Analysis 7's method settings", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  settings <- function(file) {
    text <- readLines(file.path(root, "analysis", file))
    start <- grep("^# StimGate bandwidth settings[.]", text)
    end <- grep('^tailgate_x <- "stim"', text)
    text[start:end]
  }
  expect_identical(settings("7b-sim-compare-null.qmd"), settings("7-sim-compare-freq_bs.qmd"))
})
