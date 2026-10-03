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
