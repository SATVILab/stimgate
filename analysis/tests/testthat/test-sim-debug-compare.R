root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

test_that(".simDebugCompare matches the comparison analyses' gates", {
  testthat::skip_if_not_installed("simcyto")
  testthat::skip_if_not_installed("cytoUtils")
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("numpy"))
  # flowStats warns about replacing an import when it is first loaded.
  suppressWarnings(loadNamespace("flowStats"))
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-run.R", "sim-compare-freq_bs.R",
    "sim-debug-loc.R", "sim-debug-compare.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  row <- tibble::tibble(
    transformation = "skew", prob_response = 0.02, n_cell = 2000,
    mean_pos_setting = "high", mean_pos = 8.5, bw = 0.25,
    bias_uns_setting = "low", bias_uns = 0.05,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1, sim_seed = 123L, sim_id = 1L
  )
  settings <- list(
    nSample = 2L, nMarker = 1, nCondition = 2, nCluster = 2, nIter = 1,
    bwMin = "none", bwMax = "none", probExact = TRUE, covEvMin = 1.5,
    covEvMax = 1.5, clusterGates = FALSE, locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE
  )
  dbg <- NULL
  utils::capture.output(dbg <- suppressMessages(env$.simDebugLoc(
    env$.simBandwidthRunRow(
      row, env$.simBandwidthFreqBsGlobalScenario, settings
    )
  )))
  # The comparison methods get the simulated values that StimGate stored.
  expect_equal(
    sort(dbg$sim$stim), sort(dbg$inputs$exTblStimOrig$F1),
    tolerance = 1e-6
  )

  cmp <- env$.simDebugCompare(dbg)
  expect_null(cmp$fbeta$error)
  expect_null(cmp$tailgate$error)
  expect_equal(
    cmp$fbeta$threshold,
    env$.simCompareFbetaThreshold(dbg$sim$uns, dbg$sim$stim)$threshold
  )
  expect_equal(
    cmp$tailgate$threshold,
    env$.simCompareTailgateThreshold(dbg$sim$stim, autoTol = TRUE)$threshold
  )
  expect_true(cmp$tailgate$reproduced)
  expect_equal(cmp$tailgate$tolUsed / cmp$tailgate$derivMax, 0.01)

  info <- env$.simDebugCompareInfo(cmp, dbg$sim)
  expect_named(
    info,
    c("fbetaSettings", "fbetaResult", "tailgateSettings", "tailgateResult")
  )
  expect_true("response frequency (estimated)" %in% info$fbetaResult$name)

  fig <- env$.simDebugFigure(dbg, row, cmp)
  expect_s3_class(fig, "ggplot")
  expect_gt(attr(fig, "height_cm"), 60)
  stimgate_only <- env$.simDebugFigure(dbg, row)
  expect_lt(attr(stimgate_only, "height_cm"), attr(fig, "height_cm"))

  # A failing method is reported in its section; the figure still renders.
  failed <- env$.simDebugCompare(dbg, tailgate = list(x = "neither"))
  expect_match(failed$tailgate$error, "Unknown Tailgate x")
  expect_s3_class(env$.simDebugFigure(dbg, row, failed), "ggplot")
})
