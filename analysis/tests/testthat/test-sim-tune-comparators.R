root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.tuneEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-compare-freq_bs.R", "sim-compare-tune.R"
  )) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

# Tube rows for two settings in one scenario: `good` gates every responder
# exactly; `loose` lets negatives through.
.tuneToyScores <- function(env, n_iter = 3L) {
  scen <- tibble::tibble(
    sim_id = 1L, transformation = "gaussian", mean_pos_setting = "low",
    mean_pos = 4.5, prob_response = 0.002, n_cell = 5000
  )
  rows <- tidyr::expand_grid(iter = seq_len(n_iter), sample = c("1", "2")) |>
    dplyr::mutate(propRespTruth = 0.002, propStimTruth = 0.0024)
  good <- dplyr::mutate(rows,
    method = "tailgate", tol_type = "relative", tol = 0.01, bias = 0.1, beta = NA_real_,
    tol_value = 1, cut = 2, threshold = 2.1, error = NA_character_,
    thresholdFallbackUsed = FALSE, nTruePos = 12L, nFalsePos = 0L,
    nFalseNeg = 0L, nTrueNeg = 4988L, propRespEst = 0.002
  )
  loose <- dplyr::mutate(good, bias = 0, threshold = 2, nFalsePos = 12L,
    nTrueNeg = 4976L, propRespEst = 0.004)
  fbeta <- dplyr::mutate(good, method = "fbeta", tol_type = NA_character_,
    tol = NA_real_, bias = NA_real_, beta = 0.8)
  dplyr::bind_cols(scen, dplyr::bind_rows(good, loose, fbeta))
}

test_that("the grid covers 0.2% response, every transformation and separation, at 5,000 and 100,000 cells", {
  env <- .tuneEnv()
  grid <- env$.simTuneGrid(60600L)
  expect_equal(nrow(grid), 12L)
  expect_equal(grid$sim_id, 1:12)
  expect_setequal(grid$transformation, c("gaussian", "skew", "gamma"))
  expect_setequal(as.character(grid$mean_pos_setting), c("low", "high"))
  expect_equal(unique(grid$prob_response), 0.002)
  expect_setequal(grid$n_cell, c(5e3, 1e5))
  expect_identical(grid, env$.simTuneGrid(60600L))
  expect_equal(anyDuplicated(grid$sim_seed), 0L)
  # Separations are Analysis 7's.
  pos <- env$.simMiscGetMeanPosTbl()
  joined <- dplyr::inner_join(grid, pos, by = c("transformation", "mean_pos_setting", "mean_pos"))
  expect_equal(nrow(joined), 12L)
})

test_that("the QMD selects the requested scenarios and keeps a rerun chunk", {
  qmd <- readLines(file.path(root_dir, "analysis", "6-sim-tune-comparators.qmd"))
  txt <- paste(qmd, collapse = "\n")
  expect_match(txt, ".simTuneGrid(simulation_seed, prob_response = 0.002, n_cell = c(5e3, 1e5))", fixed = TRUE)
  expect_match(txt, "label: rerun-one-simulation\n#| include: false\n#| eval: false", fixed = TRUE)
  expect_false(grepl("stimgate:::", txt, fixed = TRUE))
  # Trailing empty chunk, used to run every chunk above.
  tail_lines <- utils::tail(qmd[nzchar(trimws(qmd))], 2L)
  expect_equal(tail_lines, c("```{r}", "```"))
})

test_that("the relative tolerance 1% reproduces Tailgate's automatic tolerance", {
  testthat::skip_if_not_installed("cytoUtils")
  testthat::skip_if_not_installed("ks")
  env <- .tuneEnv()
  x <- withr::with_seed(1, c(stats::rnorm(4000), stats::rnorm(40, 4)))
  cuts <- env$.simTuneTailgateCuts(x)
  auto <- env$.simCompareTailgateThreshold(x, autoTol = TRUE, bias = 0)$threshold
  expect_equal(cuts$cut[cuts$tol_type == "relative" & cuts$tol == 0.01], auto)
  # Half steps on a log10 scale, from 1e-6 to 1e-1.
  expect_equal(log10(env$.simTuneSettings()$tailgate$tol_rel), seq(-6, -1, by = 0.5))
  expect_equal(env$.simTuneTolText(c(1e-6, 10^-5.5, 1e-2)), c("1e-6", "3.2e-6", "1e-2"))
  absolute <- env$.simCompareTailgateThreshold(x, autoTol = FALSE, tol = 0.01, bias = 0)$threshold
  expect_equal(cuts$cut[cuts$tol_type == "absolute"], absolute)
  # A smaller tolerance never moves the cutpoint left.
  rel <- cuts[cuts$tol_type == "relative", ]
  rel <- rel[order(rel$tol), ]
  expect_true(all(diff(rel$cut) <= 1e-12))
  gates <- env$.simTuneTailgateGates(x)
  expect_equal(gates$threshold, gates$cut + gates$bias)
  expect_equal(nrow(gates), nrow(cuts) * length(env$.simTuneSettings()$tailgate$bias))
})

test_that("F1 and estimate percentiles keep errors missing", {
  env <- .tuneEnv()
  scores <- .tuneToyScores(env)
  scores$error[1] <- "boom"
  m <- env$.simTuneMetrics(scores)
  expect_true(is.na(m$f1[1]))
  expect_true(is.na(m$propRespEst[1]))
  expect_equal(m$f1[2], 1)
  loose <- m$bias %in% 0
  expect_equal(unique(m$f1[loose]), 24 / 36)
  s <- env$.simTuneSummary(m)
  good <- s[s$method == "tailgate" & s$bias == 0.1, ]
  expect_equal(good$n_error, 1L)
  expect_equal(good$n_f1_finite, 5L)
  expect_equal(good$est_band_error, 0)
  expect_true(good$truth_in_band)
  expect_equal(s$est_band_error[s$method == "tailgate" & s$bias == 0], 1)
})

test_that("selection keeps near-best F1 settings and then prefers accurate estimates", {
  env <- .tuneEnv()
  s <- env$.simTuneSummary(env$.simTuneMetrics(.tuneToyScores(env)))
  r <- env$.simTuneRank(s, f1_tolerance = 0.02)
  sel <- r[r$selected, ]
  expect_setequal(sel$method, c("fbeta", "tailgate"))
  expect_equal(sel$bias[sel$method == "tailgate"], 0.1)
  # With a wide F1 tolerance both Tailgate settings are eligible; the
  # estimate error still picks the accurate one.
  r_wide <- env$.simTuneRank(s, f1_tolerance = 1)
  expect_equal(r_wide$bias[r_wide$selected & r_wide$method == "tailgate"], 0.1)
  # A setting worse on F1 than the tolerance is not chosen for its estimates.
  s2 <- s
  s2$est_band_error[s2$method == "tailgate" & s2$bias == 0] <- 0
  s2$f1_median[s2$method == "tailgate" & s2$bias == 0.1] <- 1
  r2 <- env$.simTuneRank(s2, f1_tolerance = 0.02)
  expect_equal(r2$bias[r2$selected & r2$method == "tailgate"], 0.1)
})

test_that("compared settings include every reference and survive a selected setting equal to one", {
  env <- .tuneEnv()
  m <- env$.simTuneMetrics(.tuneToyScores(env))
  s <- env$.simTuneSummary(m)
  r <- env$.simTuneRank(s)
  cmp <- env$.simTuneCompareSettings(r)
  expect_true(all(env$.simTuneReferenceSettings()$setting_ref %in% cmp$setting_ref))
  expect_true(all(c("fbeta_selected", "tailgate_selected") %in% cmp$setting_ref))
  joined <- env$.simTuneJoinSettings(s, cmp)
  # Selected Tailgate equals Analysis 7's current setting; F-beta selected
  # equals the published one. Both rows are kept.
  expect_setequal(
    as.character(joined$setting_ref),
    c("fbeta_published", "fbeta_selected", "tailgate_default", "tailgate_current", "tailgate_selected")
  )
  d <- env$.simTuneDatasetDifferences(m, cmp)
  expect_equal(unique(d$n_datasets), 3L)
  expect_equal(d$f1_diff_mean[d$setting_ref == "tailgate_default"], 1 - 24 / 36)
  expect_equal(d$f1_diff_mean[d$setting_ref == "tailgate_current"], 0)
})

test_that("the plots build from saved summaries and histograms", {
  env <- .tuneEnv()
  m <- env$.simTuneMetrics(.tuneToyScores(env))
  s <- env$.simTuneSummary(m)
  r <- env$.simTuneRank(s)
  cmp <- env$.simTuneCompareSettings(r)
  hist <- dplyr::bind_cols(
    dplyr::select(m[1, ], dplyr::all_of(env$.simTuneScenarioCols)),
    env$.simTuneHist(c(0, 1, 2.5), c(0, 1, 1.5), c("gn", "gn", "gp"), n_bins = 5L)
  ) |>
    dplyr::mutate(iter = 1L, sample = "1", .after = "n_cell")
  expect_equal(sum(hist$stim_positive), 1L)
  expect_equal(sum(hist$unstim), 3L)
  plots <- list(
    env$.simTunePlotTailgateHeat(s, "f1_median", "Median F1", r[r$selected & r$method == "tailgate", ]),
    env$.simTunePlotTailgateLines(s, "f1"),
    env$.simTunePlotTailgateLines(s, "est"),
    env$.simTunePlotTailgateTol(s, "f1"),
    env$.simTunePlotTailgateTol(s, "est"),
    env$.simTunePlotFbeta(s, "f1"),
    env$.simTunePlotFbeta(s, "est"),
    env$.simTunePlotCompare(s, cmp, "f1"),
    env$.simTunePlotCompare(s, cmp, "est"),
    env$.simTunePlotGateSpread(m, cmp),
    env$.simTunePlotGates(hist, m, cmp)
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    expect_null(p$labels$title)
  }
})

test_that("a scenario run is reproducible, restores the RNG and scores every setting", {
  testthat::skip_if_not_installed("cytoUtils")
  testthat::skip_if_not_installed("simcyto")
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("numpy"))
  env <- .tuneEnv()
  row <- env$.simTuneGrid(60600L)[1, ]
  row$n_cell <- 2000
  fbeta_env <- env$.simCompareFbetaEnvironment(
    pathFbeta = file.path(root_dir, "scripts", "python", "fbeta.py")
  )
  set.seed(5)
  before <- .Random.seed
  a <- env$.simTuneRunScenario(row, n_sample = 2L, n_iter = 1L, fbeta_env = fbeta_env, n_hist_sample = 1L)
  expect_identical(.Random.seed, before)
  b <- env$.simTuneRunScenario(row, n_sample = 2L, n_iter = 1L, fbeta_env = fbeta_env, n_hist_sample = 1L)
  expect_identical(a, b)
  settings <- env$.simTuneSettings()
  n_settings <- (length(settings$tailgate$tol_rel) + 1L) * length(settings$tailgate$bias) +
    length(settings$fbeta$beta)
  expect_equal(nrow(a$scores), 2L * n_settings)
  expect_equal(unique(a$hist$sample), "1")
  expect_true(all(a$scores$propRespTruth > 0))
  ok <- is.na(a$scores$error) & is.finite(a$scores$threshold)
  expect_equal(
    a$scores$nTruePos[ok] + a$scores$nFalsePos[ok],
    a$scores$nPosStim[ok]
  )
})

test_that("parallel workers reproduce the serial run exactly", {
  testthat::skip_if_not_installed("cytoUtils")
  testthat::skip_if_not_installed("simcyto")
  testthat::skip_if_not_installed("furrr")
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("numpy"))
  env <- .tuneEnv()
  grid <- env$.simTuneGrid(60600L)[c(1L, 9L), ]
  grid$n_cell <- c(1500, 1000)
  path_fbeta <- file.path(root_dir, "scripts", "python", "fbeta.py")
  serial <- env$.simTuneRunGrid(grid, n_sample = 2L, n_iter = 2L,
    path_fbeta = path_fbeta, n_hist_sample = 1L)
  parallel <- env$.simTuneRunGrid(grid, n_sample = 2L, n_iter = 2L,
    path_fbeta = path_fbeta, n_hist_sample = 1L, workers = 2L, root_dir = root_dir)
  expect_identical(parallel, serial)
  expect_equal(sort(unique(serial$scores$iter)), 1:2)
  expect_equal(unique(serial$hist$iter), 1L)
  # A scenario run alone gives the same rows as within the grid.
  one <- env$.simTuneRunScenario(grid[2L, ], n_sample = 2L, n_iter = 2L,
    fbeta_env = env$.simCompareFbetaEnvironment(pathFbeta = path_fbeta),
    n_hist_sample = 1L)
  expect_identical(one$scores, dplyr::filter(serial$scores, .data$sim_id == grid$sim_id[[2L]]))
})
