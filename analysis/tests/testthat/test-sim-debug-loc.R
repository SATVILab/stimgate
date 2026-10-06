root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

test_that(".simDebugLoc records one sample without changing the rerun", {
  testthat::skip_if_not_installed("simcyto")
  env <- new.env(parent = getNamespace("stimgate"))
  for (fn in c(
    "analysis-runtime.R", "analysis-plot-style.R", "sim-misc.R",
    "sim-bandwidth.R", "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-run.R", "sim-debug-loc.R"
  )) {
    source(file.path(root_dir, "scripts", "r", fn), local = env)
  }
  row <- tibble::tibble(
    transformation = "skew", prob_response = 0.02, n_cell = 2000,
    mean_pos_setting = "high", mean_pos = env$.simMiscGetMeanPosTbl()$mean_pos[1],
    bw = 0.25, bias_uns_setting = "low", bias_uns = 0.05,
    sample_perturbation_sd = 0, condition_perturbation_sd = 0,
    cluster_perturbation_sd = 0, background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1, sim_seed = 123L, sim_id = 1L
  )
  settings <- list(
    nSample = 3L, nMarker = 1, nCondition = 2, nCluster = 2, nIter = 1,
    bwMin = "none", bwMax = "none", probExact = TRUE, covEvMin = 1.5,
    covEvMax = 1.5, clusterGates = FALSE, locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE
  )
  rerun <- function() {
    env$.simBandwidthRunRow(
      row, env$.simBandwidthFreqBsGlobalScenario, settings
    )
  }
  quiet <- function(code) {
    out <- NULL
    utils::capture.output(out <- suppressMessages(code))
    out
  }
  ref <- quiet(rerun())

  withr::local_seed(1)
  seed_before <- .Random.seed
  dbg <- quiet(env$.simDebugLoc(rerun(), sample = 2, stopAfter = FALSE))
  expect_identical(.Random.seed, seed_before)
  expect_s3_class(dbg, "simDebugLoc")
  expect_equal(dbg$result, ref)
  expect_identical(dbg$ind, "4")
  stored <- ref$threshold[ref$method == "loc_sample" & ref$ind == "4"]
  expect_equal(dbg$cp$cp, stored)
  expect_false(inherits(
    get(".getCpUnsLocCondition", envir = asNamespace("stimgate")),
    "functionWithTrace"
  ))

  summary <- env$.simDebugLocSummary(dbg)
  expect_equal(
    summary$propRespEst,
    ref$propRespEst[ref$method == "loc_sample" & ref$ind == "4"]
  )
  expect_true(all(c("fdp", "sensitivity", "propRespTruth") %in% names(summary)))
  plots <- env$.simDebugLocPlots(dbg)
  expect_setequal(
    names(plots),
    c("density", "prob", "deriv", "respCells", "taut", "truth")
  )
  expect_true(all(vapply(plots, ggplot2::is.ggplot, logical(1L))))
  xr <- vapply(plots, function(p) {
    paste(p$coordinates$limits$x, collapse = ",")
  }, character(1L))
  expect_length(unique(xr), 1L)
  expect_true(is.data.frame(attr(plots, "lines")))
  zoomed <- env$.simDebugLocPlots(dbg, xlim = c(3, 8))
  expect_equal(zoomed$density$coordinates$limits$x, c(3, 8))
  expect_s3_class(env$.simDebugLocPlotGrid(plots), "ggplot")
  info <- env$.simDebugLocInfo(dbg, row)
  expect_named(info, c("gating", "simulation", "estimate"))
  value <- function(tbl, nm) tbl$value[tbl$name == nm]
  expect_identical(value(info$gating, "bw (grid)"), "0.25")
  expect_identical(value(info$gating, "bandwidth used"), "0.25")
  expect_identical(value(info$gating, "biasUns used"), "0.05")
  expect_identical(value(info$simulation, "sample"), "2")
  expect_identical(value(info$simulation, "transformation"), "skew")
  expect_false("bw" %in% info$simulation$name)
  expect_true(all(c(
    "response frequency (estimated)", "response frequency (true)",
    "relative error"
  ) %in% info$estimate$name))
  expect_s3_class(env$.simDebugLocPlotGrid(plots, info), "ggplot")

  first <- quiet(env$.simDebugLoc(rerun()))
  expect_identical(first$ind, "2")
  expect_null(first$result)

  later <- quiet(env$.simDebugLoc(rerun(), sample = 2, subsequent = TRUE))
  expect_s3_class(later, "simDebugLocList")
  expect_identical(names(later), c("dataset1_ind4", "dataset1_ind6"))
  expect_equal(
    unname(vapply(later, function(x) x$cp$cp, numeric(1L))),
    ref$threshold[ref$method == "loc_sample" & ref$ind %in% c("4", "6")]
  )
})
