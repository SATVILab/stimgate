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

  bw_wrapper_params <- c(
    bw_params,
    list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL)
  )

  set.seed(99)
  bw_direct <- do.call(env$.simBandwidthBwOne, bw_wrapper_params)

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  expect_equal(as.numeric(bw_direct), as.numeric(bw_pkg))
})

test_that(".simBandwidthBwOne resolves the current stimgate namespace in future workers", {
  skip_if_not_installed("future")
  skip_if_not_installed("pkgload")

  if (!requireNamespace("stimgate", quietly = TRUE)) {
    skip("stimgate package not available")
  }

  root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

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
  bw_wrapper_params <- c(
    bw_params,
    list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL)
  )

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2L)

  result <- future::future(
    {
      env <- new.env(parent = globalenv())
      source(script_misc, local = env)
      source(script_bw, local = env)
      env$.simBandwidthEnsureCurrentCheckout(root_dir)
      ns <- asNamespace("stimgate")

      suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
      set.seed(99)
      bw_worker <- do.call(env$.simBandwidthBwOne, bw_wrapper_params)

      set.seed(99)
      bw_worker_pkg <- do.call(
        get(".bwCalcOne", mode = "function", envir = ns),
        bw_params
      )

      list(
        ns_path = normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = FALSE),
        bw_worker = bw_worker,
        bw_pkg = bw_worker_pkg
      )
    },
    seed = TRUE
  )

  worker_out <- future::value(result)

  expect_equal(
    normalizePath(worker_out$ns_path, winslash = "/", mustWork = FALSE),
    normalizePath(root_dir, winslash = "/", mustWork = FALSE)
  )
  expect_true(is.finite(as.numeric(worker_out$bw_worker)))
  expect_equal(as.numeric(worker_out$bw_worker), as.numeric(worker_out$bw_pkg))
  expect_equal(as.numeric(worker_out$bw_worker), as.numeric(bw_pkg))
})
