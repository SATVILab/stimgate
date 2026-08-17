test_that("native tailgate preserves the relative-derivative shoulder rule", {
  x <- seq(0, 8, length.out = 500)
  y <- stats::dnorm(x, mean = 2.5, sd = 0.8)

  tail <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = x, y = y),
    peakX = 2.5,
    fraction = 1 / 200
  )

  expect_true(is.finite(tail$lowerBoundX))
  expect_gt(tail$lowerBoundX, 2.5)
  expect_lt(tail$lowerBoundX, max(x))
  expect_equal(tail$info$reason, "identified_stimulated_density_shoulder_lower_bound")
})

test_that("native tailgate returns NA when the descending shoulder never reaches the target", {
  x <- seq(0, 1, length.out = 300)
  y <- seq(0.1, 1, length.out = 300)

  tail <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = x, y = y),
    peakX = max(x),
    fraction = 1 / 200
  )

  expect_false(is.finite(tail$lowerBoundX))
  expect_equal(tail$info$reason, "no_negative_derivative_right_of_stimulated_peak")
})

test_that("native tailgate has a deterministic multimodal fixture and differs from legacy tol gating", {
  set.seed(1)
  x <- c(rnorm(200, mean = 0, sd = 0.6), rnorm(200, mean = 4, sd = 0.8))
  density <- stats::density(x, n = 512)

  native <- stimgate:::.getStimGateTailgate(
    density = data.frame(x = density$x, y = density$y),
    fraction = 1 / 200
  )
  wrapper <- stimgate:::.getCpUnsLocMarginalDensityLowerBound(
    density = data.frame(x = density$x, y = density$y),
    peakX = density$x[which.max(density$y)],
    fraction = 1 / 200
  )
  legacy <- .cytokineCutpoint(
    x = x,
    numPeaks = 1,
    refPeak = 1,
    tol = 1e-2,
    side = "right",
    strict = FALSE,
    plot = FALSE
  )

  expect_true(is.finite(native$lowerBoundX))
  expect_true(is.finite(legacy))
  expect_equal(wrapper$lowerBoundX, native$lowerBoundX, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(native$lowerBoundX, legacy, tolerance = 1e-8)))
  expect_true(native$lowerBoundX > min(x))
  expect_true(native$lowerBoundX < max(x))
})
