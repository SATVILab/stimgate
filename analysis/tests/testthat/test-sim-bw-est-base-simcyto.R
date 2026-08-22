root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")
script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

test_that("base bandwidth simulation helpers source cleanly without legacy functionsForBenchmarking-Cyt.R", {
  for (f in c(script_misc, script_bw, script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  expect_no_error(source(script_misc, local = env))
  expect_no_error(source(script_bw, local = env))
  expect_no_error(source(script_bw_io, local = env))
  expect_no_error(source(script_bw_plot, local = env))

  # Ensure no local unexported simCytExperiment is introduced into the helper environment
  expect_false(exists("simCytExperiment", envir = env, inherits = FALSE))
})

test_that("analysis/3-sim-bw-est-base.qmd does not source functionsForBenchmarking-Cyt.R", {
  qmd_path <- file.path(root_dir, "analysis", "3-sim-bw-est-base.qmd")
  expect_true(file.exists(qmd_path))

  lines <- readLines(qmd_path, warn = FALSE)
  expect_false(
    any(grepl("functionsForBenchmarking-Cyt\\.R", lines)),
    info = "analysis/3-sim-bw-est-base.qmd should not source functionsForBenchmarking-Cyt.R"
  )
})

test_that(".simBandwidthEstBwDirect calls simcyto::simCytExperiment and produces valid Gaussian bandwidth results", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  orig_simcyto_experiment <- simcyto::simCytExperiment
  called_simcyto <- FALSE

  testthat::with_mocked_bindings(
    simCytExperiment = function(...) {
      called_simcyto <<- TRUE
      orig_simcyto_experiment(...)
    },
    .package = "simcyto",
    {
      set.seed(42)
      res_gauss <- env$.simBandwidthEstBwDirect(
        nSample = 2L,
        nMarker = 1L,
        nCondition = 2L,
        nCluster = 2L,
        nIter = 1L,
        biasUns = 0.05,
        bwMtd = "hpi1",
        bwFallback = 0.234,
        bwMin = -Inf,
        bwMax = Inf,
        bwCluster = 0.5,
        nCellStim = 500L,
        probResponse = 0.002,
        probExact = TRUE,
        meanPos = 8,
        transformation = "gaussian",
        backgroundRelativeToResponse = 0.2,
        ncellUnsRelativeToStim = 1,
        covEvMin = 1.5,
        covEvMax = 1.5,
        tolClust = NULL,
        summarise = FALSE
      )

      expect_true(called_simcyto)
      expect_s3_class(res_gauss, "tbl_df")
      expect_equal(nrow(res_gauss), 2L)
      expect_equal(res_gauss$n_cell_stim, c(500L, 500L))
      expect_equal(res_gauss$n_cell_uns, c(500L, 500L))
      expect_true(all(is.finite(res_gauss$bw)))
      expect_true(all(res_gauss$bw > 0))
      expect_true(all(c("bw_uns", "bw_stim", "bw") %in% colnames(res_gauss)))
    }
  )
})

test_that(".simBandwidthEstBwDirect works with gamma and skew transformations from simcyto", {
  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  for (tr in c("gamma", "skew")) {
    set.seed(123)
    res <- env$.simBandwidthEstBwDirect(
      nSample = 2L,
      nMarker = 1L,
      nCondition = 2L,
      nCluster = 2L,
      nIter = 1L,
      biasUns = 0.05,
      bwMtd = "hpi1",
      bwFallback = 0.234,
      bwMin = -Inf,
      bwMax = Inf,
      bwCluster = 0.5,
      nCellStim = 500L,
      probResponse = 0.002,
      probExact = TRUE,
      meanPos = 7,
      transformation = tr,
      backgroundRelativeToResponse = 0.2,
      ncellUnsRelativeToStim = 1,
      covEvMin = 1.5,
      covEvMax = 1.5,
      tolClust = NULL,
      summarise = FALSE
    )

    expect_s3_class(res, "tbl_df")
    expect_equal(nrow(res), 2L)
    expect_equal(res$transformation, c(tr, tr))
    expect_equal(res$n_cell_stim, c(500L, 500L))
    expect_equal(res$n_cell_uns, c(500L, 500L))
    expect_true(all(is.finite(res$bw)))
    expect_true(all(res$bw > 0))
  }
})

test_that("simcyto::simCytExperiment generates expected flowFrames and cluster labels", {
  nCellStim <- 300L
  nCellUns <- 300L
  nCellByCondition <- c(nCellUns, nCellStim)
  transformationFunc <- simcyto::simCytTransformGaussian()
  meanExprMat <- matrix(c(0, 8), byrow = TRUE, ncol = 1)
  clusterLabelVec <- c("gn", "gp")
  probResponse <- 0.05
  probResponseUns <- probResponse * 0.2
  probVecUns <- c(1 - probResponseUns, probResponseUns)
  probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

  set.seed(42)
  sim_res <- simcyto::simCytExperiment(
    nSample = 2L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nCellByCondition = nCellByCondition,
    transformationFunc = transformationFunc,
    mixtureType = "gaussianOnly",
    meanExprMat = meanExprMat,
    clusterLabelVec = clusterLabelVec,
    probVecUns = probVecUns,
    probExact = TRUE,
    probResponseVecByStimCondition = probResponseVecByStimCondition,
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    covEvMin = 1.5,
    covEvMax = 1.5
  )

  expect_named(sim_res, c("flowFrameList", "labelsList"))
  expect_length(sim_res$flowFrameList, 4L)
  expect_length(sim_res$labelsList, 4L)

  # Check cell counts and labels
  for (i in seq_along(sim_res$flowFrameList)) {
    ff <- sim_res$flowFrameList[[i]]
    expect_s4_class(ff, "flowFrame")
    expect_equal(nrow(flowCore::exprs(ff)), 300L)
    expect_equal(ncol(flowCore::exprs(ff)), 1L)

    labels <- sim_res$labelsList[[i]]
    expect_length(labels, 300L)
    expect_true(all(labels %in% c("gn", "gp")))
  }
})
