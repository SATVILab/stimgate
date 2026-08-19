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

  # Shared params passed explicitly to both calls so they are identical.
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

  # Use the same seed before each call: .bwCalcOne uses internal randomness
  # (e.g., when subsampling), and .simBandwidthBwOne delegates to it without
  # introducing any additional RNG calls before the delegation.
  set.seed(99)
  bw_direct <- do.call(
    env$.simBandwidthBwOne,
    c(bw_params, list(bwMin = 1e-10, bwMax = 1e10, bwFallback = NULL))
  )

  set.seed(99)
  bw_pkg <- do.call(stimgate:::.bwCalcOne, bw_params)

  expect_equal(as.numeric(bw_direct), as.numeric(bw_pkg))
})
