test_that("refinement preserves marginal validity and interval diagnostics", {
  refine <- function(lower, gate) {
    .getCpPosTautString(
      ex = data.frame(A = 1:5), inc = rep(TRUE, 5), chnl = "A",
      cpOrig = gate,
      peakX = if (is.finite(lower)) 1 else NA_real_,
      windowWidth = if (is.finite(lower)) 3 else NA_real_, lower = lower
    )
  }
  invalid <- refine(NA_real_, 4)
  expect_identical(invalid$reason, "invalid_refinement_interval")
  expect_identical(invalid$lowerX, NA_real_)
  empty <- refine(2, 2)
  expect_identical(empty$reason, "empty_refinement_interval")
  expect_identical(empty$lowerX, 2)
  available <- refine(2, 4)
  expect_identical(available$reason, "too_few_other_cytokine_positive_cells")
  expect_identical(available$lowerX, 2)
})
