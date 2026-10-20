# Tests for ordinary scalar bandwidth selection in R/bandwidth.R.

.makeTestVec <- function(seed = 42L, n = 500L) {
  set.seed(seed)
  c(
    stats::rnorm(round(n * 0.8), mean = 0, sd = 1),
    stats::rnorm(round(n * 0.2), mean = 4, sd = 0.5)
  )
}

test_that("bwCalcOne supports all six ordinary bandwidth selectors", {
  withr::local_preserve_seed()
  x <- .makeTestVec()
  for (method in c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")) {
    bw <- stimgate:::.bwCalcOne(x, bwMtd = method)
    expect_type(bw, "double")
    expect_length(bw, 1L)
    expect_true(is.finite(bw) && bw > 0, info = method)
    expected <- switch(method,
      nrd0 = stats::bw.nrd0(x),
      sj = stats::bw.SJ(x),
      ks::hpi(x, deriv.order = as.integer(sub("hpi", "", method)))
    )
    expect_equal(bw, as.numeric(expected), info = method)
    expect_no_error(stimControl(bwMtd = method))
  }
})

test_that("bwAdj scales ordinary bandwidth proportionally", {
  withr::local_preserve_seed()
  x <- .makeTestVec()
  bw1 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwAdj = 1)
  bw2 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwAdj = 2)
  expect_equal(bw2, bw1 * 2, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Edge cases — degenerate inputs
# ---------------------------------------------------------------------------

test_that("bwCalcOne returns NA for fewer than 2 observations", {
  bw <- stimgate:::.bwCalcOne(1, bwMtd = "nrd0")
  expect_true(is.na(bw))
})

test_that("bwCalcOne returns NA for all-identical values", {
  bw <- stimgate:::.bwCalcOne(rep(3.14, 50L), bwMtd = "nrd0")
  expect_true(is.na(bw))
})

test_that("bwCalcOne handles all-NA / non-finite input gracefully", {
  bw_na <- stimgate:::.bwCalcOne(rep(NA_real_, 30L), bwMtd = "nrd0")
  expect_true(is.na(bw_na))
  bw_inf <- stimgate:::.bwCalcOne(c(Inf, -Inf, NaN), bwMtd = "nrd0")
  expect_true(is.na(bw_inf))
})

test_that("bwCalcOne handles mixed finite and non-finite input", {
  withr::local_preserve_seed()
  x <- c(.makeTestVec(n = 200L), NA_real_, Inf, -Inf, NaN)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_true(is.finite(bw) && bw > 0)
})

# ---------------------------------------------------------------------------
# bwNcellMax / bwNcellMin sampling constraints
# ---------------------------------------------------------------------------

# bwNcellMax caps the number of cells used for bandwidth estimation.  With a
# large input (5000) and a tiny cap (20), the estimator runs on a 20-cell
# subsample.  The resulting bandwidth should differ from the uncapped result
# because the sample the estimator sees is different in size and composition.
test_that("bwCalcOne ordinary path bwNcellMax changes the result vs uncapped", {
  withr::local_preserve_seed()
  set.seed(1L)
  x <- stats::rnorm(5000L)
  set.seed(1L)
  bw_full <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  set.seed(1L)
  bw_capped <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwNcellMax = 20L)
  expect_true(is.finite(bw_full) && bw_full > 0)
  expect_true(is.finite(bw_capped) && bw_capped > 0)
  # A 20-cell subsample must produce a different bandwidth than the full 5000.
  expect_false(isTRUE(all.equal(as.numeric(bw_full), as.numeric(bw_capped))))
})

# bwNcellMin bootstraps up when the input is smaller than the requested
# minimum.  With 5 cells and bwNcellMin = 200, the estimator runs on a
# 200-cell bootstrap — which should give a different result than running on
# the original 5 cells (the latter would fail .bwCalcOneBase with n < 2
# unique values, but nrd0 can handle 5; the key is the size differs).
test_that("bwCalcOne ordinary path bwNcellMin bootstraps input up to minimum size", {
  withr::local_preserve_seed()
  set.seed(2L)
  x <- stats::rnorm(5L)
  # Without bwNcellMin the estimator sees 5 cells.
  set.seed(2L)
  bw_small <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  # With bwNcellMin = 200 the estimator sees a 200-cell bootstrap sample.
  set.seed(2L)
  bw_boosted <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwNcellMin = 200L)
  expect_true(is.finite(bw_small) && bw_small > 0)
  expect_true(is.finite(bw_boosted) && bw_boosted > 0)
  # The bootstrap changes the effective sample, so the bandwidth must differ.
  expect_false(isTRUE(all.equal(as.numeric(bw_small), as.numeric(bw_boosted))))
})

# ---------------------------------------------------------------------------
# Stability regression — protect a key invariant, not an exact number
# ---------------------------------------------------------------------------

test_that("bwCalcOne nrd0 bandwidth is in a plausible range for a known distribution", {
  withr::local_preserve_seed()
  # For rnorm(500, mean=0, sd=1), Silverman's rule-of-thumb gives ~ 0.2-0.5.
  # This is a loose sanity check, not an exact regression fixture.
  set.seed(42L)
  x <- stats::rnorm(500L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_gt(as.numeric(bw), 0.1)
  expect_lt(as.numeric(bw), 1.0)
})
