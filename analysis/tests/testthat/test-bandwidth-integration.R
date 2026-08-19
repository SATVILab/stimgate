root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")

test_that(".simBandwidthBwOne agrees numerically with stimgate:::.bwCalcOne", {
  for (f in c(script_misc, script_bw)) {
    if (!file.exists(f)) stop("Expected analysis helper not found: ", f)
  }

  env <- new.env(parent = getNamespace("stimgate"))
  source(script_misc, local = env)
  source(script_bw, local = env)

  set.seed(42)
  x_sample <- rnorm(200, mean = 5, sd = 1.5)

  bw_params <- list(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normPeakMinRel = 0.75,
    normExtraFrac = 0.25,
    normExtraMax = Inf,
    normExtraJitterFrac = 0.25,
    normLambda = seq(-2, 2, length.out = 81),
    normDensityN = 512L,
    normExcessBwMtd = "hpi3",
    normExcessNcell = 10000L,
    normAdaptiveNcell = 2500L,
    bwAdaptiveCore = NULL,
    bwAdaptiveExtra = NULL,
    bwAdaptiveCrossover = NULL,
    bwAdaptiveTransitionWidth = 0,
    normMtd = "moments",
    adaptive = FALSE
  )

  set.seed(99)
  bw_direct <- do.call(
    env$.simBandwidthBwOne,
    c(bw_params, list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL))
  )

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  expect_equal(as.numeric(bw_direct), as.numeric(bw_pkg))
})

test_that(".simBandwidthBwOne resolves the current stimgate namespace in future workers", {
  skip_if_not_installed("future")

  if (!exists("stimgate", mode = "package")) {
    skip("stimgate package not attached")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2L)

  x_sample <- rnorm(200, mean = 5, sd = 1.5)
  bw_params <- list(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normPeakMinRel = 0.75,
    normExtraFrac = 0.25,
    normExtraMax = Inf,
    normExtraJitterFrac = 0.25,
    normLambda = seq(-2, 2, length.out = 81),
    normDensityN = 512L,
    normExcessBwMtd = "hpi3",
    normExcessNcell = 10000L,
    normAdaptiveNcell = 2500L,
    bwAdaptiveCore = NULL,
    bwAdaptiveExtra = NULL,
    bwAdaptiveCrossover = NULL,
    bwAdaptiveTransitionWidth = 0,
    normMtd = "moments",
    adaptive = FALSE,
    bwMin = 1e-10,
    bwMax = 1e10,
    bwFallback = NULL
  )

  result <- future::future({
    env <- new.env(parent = getNamespace("stimgate"))
    source(script_misc, local = env)
    source(script_bw, local = env)
    set.seed(42)
    do.call(env$.simBandwidthBwOne, bw_params)
  }, seed = TRUE)

  bw_worker <- future::value(result)
  expect_true(is.finite(as.numeric(bw_worker)))
})
