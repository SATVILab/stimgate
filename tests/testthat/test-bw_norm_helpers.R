# Tests for .bwCalcOne() in R/bw_norm_helpers.R
# These are package-level unit/regression tests that exercise the function
# directly, without depending on scripts/r analysis helpers.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# A simple bimodal-ish mixture: large "background" peak + small "signal" tail.
# n must be large enough to exercise the *Norm code paths (>= 20 unique values).
.makeTestVec <- function(seed = 42L, n = 500L) {
  set.seed(seed)
  c(
    stats::rnorm(round(n * 0.8), mean = 0, sd = 1),
    stats::rnorm(round(n * 0.2), mean = 4, sd = 0.5)
  )
}

# ---------------------------------------------------------------------------
# Ordinary (non-Norm) bandwidth paths
# ---------------------------------------------------------------------------

test_that("bwCalcOne returns finite positive scalar for nrd0", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_true(is.numeric(bw))
  expect_length(bw, 1L)
  expect_true(is.finite(bw) && bw > 0)
  expect_false(isTRUE(attr(bw, "adaptive")))
})

test_that("bwCalcOne returns finite positive scalar for sj", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "sj")
  expect_true(is.finite(bw) && bw > 0)
  expect_false(isTRUE(attr(bw, "adaptive")))
})

test_that("bwCalcOne returns finite positive scalar for hpi3", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "hpi3")
  expect_true(is.finite(bw) && bw > 0)
  expect_false(isTRUE(attr(bw, "adaptive")))
})

test_that("bwAdj scales ordinary bandwidth proportionally", {
  x <- .makeTestVec()
  bw1 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwAdj = 1)
  bw2 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwAdj = 2)
  expect_equal(as.numeric(bw2), as.numeric(bw1) * 2, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Norm bandwidth path (normMtd = "moments", the default)
# ---------------------------------------------------------------------------

test_that("bwCalcOne with nrd0Norm returns finite positive scalar", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", normMtd = "moments")
  expect_true(is.numeric(bw))
  expect_length(bw, 1L)
  expect_true(is.finite(bw) && bw > 0)
  expect_false(isTRUE(attr(bw, "adaptive")))
})

test_that("bwCalcOne with sjNorm returns finite positive scalar", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "sjNorm", normMtd = "moments")
  expect_true(is.finite(bw) && bw > 0)
})

test_that("bwCalcOne with hpi3Norm returns finite positive scalar", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "hpi3Norm", normMtd = "moments")
  expect_true(is.finite(bw) && bw > 0)
})

test_that("bwAdj scales Norm bandwidth proportionally", {
  x <- .makeTestVec()
  set.seed(11L)
  bw1 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", bwAdj = 1)
  set.seed(11L)
  bw2 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", bwAdj = 2)
  expect_equal(as.numeric(bw2), as.numeric(bw1) * 2, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Adaptive mode
# ---------------------------------------------------------------------------

test_that("bwCalcOne adaptive=TRUE with NormMtd returns list with adaptive attribute", {
  x <- .makeTestVec(n = 1000L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", adaptive = TRUE)
  expect_true(isTRUE(attr(bw, "adaptive")))
  expect_true(is.list(bw))
  expect_true(all(c("bin", "bw", "bwCore", "bwExtra") %in% names(bw)))
  expect_true(all(is.finite(bw$bwCore) & bw$bwCore > 0))
  expect_true(all(is.finite(bw$bwExtra) & bw$bwExtra > 0))
  expect_true(all(is.finite(bw$bw) & bw$bw > 0))
})

test_that("bwCalcOne adaptive=TRUE with non-Norm method returns scalar (not adaptive)", {
  x <- .makeTestVec()
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", adaptive = TRUE)
  # ordinary (non-Norm) path ignores adaptive flag and returns scalar
  expect_false(isTRUE(attr(bw, "adaptive")))
  expect_length(bw, 1L)
  expect_true(is.finite(bw) && bw > 0)
})

test_that("bwCalcOne adaptive=TRUE boxcox raises an error", {
  x <- .makeTestVec()
  expect_error(
    stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", adaptive = TRUE, normMtd = "boxcox"),
    regexp = "boxcox"
  )
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
  x <- c(.makeTestVec(n = 200L), NA_real_, Inf, -Inf, NaN)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_true(is.finite(bw) && bw > 0)
})

# ---------------------------------------------------------------------------
# Reproducibility with seeding
# ---------------------------------------------------------------------------

test_that("bwCalcOne Norm scalar result is reproducible with same seed", {
  x <- .makeTestVec(seed = 7L, n = 600L)
  set.seed(123L)
  bw1 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm")
  set.seed(123L)
  bw2 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm")
  expect_equal(as.numeric(bw1), as.numeric(bw2))
})

test_that("bwCalcOne adaptive Norm result is reproducible with same seed", {
  x <- .makeTestVec(seed = 99L, n = 1000L)
  set.seed(55L)
  bw1 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", adaptive = TRUE)
  set.seed(55L)
  bw2 <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", adaptive = TRUE)
  expect_equal(bw1$bwCore, bw2$bwCore)
  expect_equal(bw1$bwExtra, bw2$bwExtra)
  expect_equal(bw1$bw, bw2$bw)
})

# ---------------------------------------------------------------------------
# bwNcellMax / bwNcellMin sampling constraints
# ---------------------------------------------------------------------------

test_that("bwCalcOne ordinary path respects bwNcellMax (downsampled)", {
  set.seed(1L)
  x <- stats::rnorm(5000L)
  # cap at 100 cells: should still return a finite positive bandwidth
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwNcellMax = 100L)
  expect_true(is.finite(bw) && bw > 0)
})

test_that("bwCalcOne ordinary path respects bwNcellMin (bootstraps up)", {
  set.seed(2L)
  x <- stats::rnorm(5L)   # small – below the default bwNcellMin of 100
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0", bwNcellMin = 20L)
  expect_true(is.finite(bw) && bw > 0)
})

# ---------------------------------------------------------------------------
# Stability regression — protect a key invariant, not an exact number
# ---------------------------------------------------------------------------

test_that("bwCalcOne nrd0 bandwidth is in a plausible range for a known distribution", {
  # For rnorm(500, mean=0, sd=1), Silverman's rule-of-thumb gives ~ 0.2-0.5.
  # This is a loose sanity check, not an exact regression fixture.
  set.seed(42L)
  x <- stats::rnorm(500L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_gt(as.numeric(bw), 0.1)
  expect_lt(as.numeric(bw), 1.0)
})
