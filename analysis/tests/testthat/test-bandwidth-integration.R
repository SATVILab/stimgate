root_dir <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = FALSE)
if (!file.exists(file.path(root_dir, "scripts", "r", "sim-bandwidth.R"))) {
  root_dir <- normalizePath(
    file.path(testthat::test_path(), "../../.."),
    mustWork = FALSE
  )
}

script_misc <- file.path(root_dir, "scripts", "r", "sim-misc.R")
script_bw <- file.path(root_dir, "scripts", "r", "sim-bandwidth.R")

test_that(".simBandwidthBwOne agrees numerically with stimgate:::.bwCalcOne", {
  skip_if_not_installed("stimgate")
  skip_if_not(file.exists(script_bw), "sim-bandwidth.R not found")

  env <- new.env(parent = baseenv())
  source(script_misc, local = env)
  source(script_bw, local = env)

  set.seed(42)
  x_sample <- rnorm(200, mean = 5, sd = 1.5)

  bw_direct <- env$.simBandwidthBwOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwMin = 1e-10,
    bwMax = 1e10,
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    bwFallback = NULL,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25
  )

  bw_pkg <- stimgate:::.bwCalcOne(
    x = x_sample,
    bwMtd = "hpi1Norm",
    bwAdj = 1,
    bwNcellMin = 10L,
    bwNcellMax = 1000L,
    normPeakFrac = 0.15,
    normExtraFrac = 0.25,
    adaptive = FALSE
  )

  expect_equal(as.numeric(bw_direct), as.numeric(bw_pkg))
})
