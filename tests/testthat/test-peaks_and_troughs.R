test_that(".getPeakMainLeftIdx handles NA when no local maxima exist", {
  # Monotonic or NA-bounded vector with no local maxima
  y_na <- c(1, 5, NA)
  expect_equal(stimgate:::.getPeakMainLeftIdx(y_na), 2L)

  # Tie preservation: should pick the last index of the maximum
  y_na_tie <- c(1, 5, 5, NA)
  expect_equal(stimgate:::.getPeakMainLeftIdx(y_na_tie), 3L)
})

test_that(".getPeakMainLeftIdx preserves tie behaviour on non-NA input without local maxima", {
  y_tie <- c(1, 5, 5)
  expect_equal(stimgate:::.getPeakMainLeftIdx(y_tie), 3L)
})
