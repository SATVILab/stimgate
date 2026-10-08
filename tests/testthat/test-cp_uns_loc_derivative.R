.getCpUnsLocDerivThreshold <- get(
  ".getCpUnsLocDerivThreshold",
  envir = asNamespace("stimgate"),
  mode = "function"
)

test_that("invalid psi falls back to the stage default psi", {
  x <- seq_len(20)
  deriv <- c(0, 0, 0, 1, 2, 4, 6, 8, 9, 10, 9, 8, 6, 4, 2, 1, 0, 0, 0, 0)
  prob <- seq(0.5, 1, length.out = 20)

  out <- .getCpUnsLocDerivThreshold(
    x = x,
    prob = prob,
    deriv = deriv,
    alpha = 0.5,
    omega = 0.15,
    psi = NA_real_,
    stage = "marginal"
  )

  expect_identical(out$info$psi, -0.2)
  expect_identical(out$info$thresholdBasis, "right_fall_to_psi_times_peak")
  expect_equal(out$thresholdX, 15)
})

# A small early rise that levels off at 0.2, then the real rise to 0.95: the
# early rise is as steep, so without the check it starts the region.
test_that("a derivative peak counts only if its rise reaches locMinRiseProb", {
  x <- seq(0, 2, length.out = 401)
  prob <- 0.2 * stats::plogis((x - 0.5) / 0.007) +
    0.75 * stats::plogis((x - 1.5) / 0.03)
  deriv <- c(diff(prob) / diff(x), 0)
  peakX <- function(minRiseProb) {
    out <- stimgate:::.getCpUnsLocDerivPeak(
      x, prob, deriv,
      alpha = 0.5, minRiseProb = minRiseProb
    )
    out$data$x[[out$index]]
  }
  expect_equal(peakX(NA_real_), 0.5, tolerance = 0.01)
  expect_equal(peakX(0), 0.5, tolerance = 0.01)
  expect_equal(peakX(1 / 3), 1.5, tolerance = 0.01)

  # A single genuine rise is unaffected.
  probOne <- 0.95 * stats::plogis((x - 1) / 0.05)
  derivOne <- c(diff(probOne) / diff(x), 0)
  one <- stimgate:::.getCpUnsLocDerivPeak(
    x, probOne, derivOne,
    alpha = 0.5, minRiseProb = 1 / 3
  )
  expect_equal(one$data$x[[one$index]], 1, tolerance = 0.01)

  # No rise reaching the level: no peak.
  low <- stimgate:::.getCpUnsLocDerivPeak(
    x, 0.2 * stats::plogis((x - 0.5) / 0.02),
    c(diff(0.2 * stats::plogis((x - 0.5) / 0.02)) / diff(x), 0),
    alpha = 0.5, minRiseProb = 1 / 3
  )
  expect_true(is.na(low$index))
  expect_identical(low$info$reason, "no_derivative_peak_rise_reached_min_probability")
  expect_identical(stimControl()$locMinRiseProb, 1 / 3)
  expect_error(stimControl(locMinRiseProb = 2), "locMinRiseProb")
})
