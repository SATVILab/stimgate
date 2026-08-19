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
  set.seed(11L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", normMtd = "moments")
  expect_true(is.numeric(bw))
  expect_length(bw, 1L)
  expect_true(is.finite(bw) && bw > 0)
  expect_false(isTRUE(attr(bw, "adaptive")))
})

test_that("bwCalcOne with sjNorm returns finite positive scalar", {
  x <- .makeTestVec()
  set.seed(11L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "sjNorm", normMtd = "moments")
  expect_true(is.finite(bw) && bw > 0)
})

test_that("bwCalcOne with hpi3Norm returns finite positive scalar", {
  x <- .makeTestVec()
  set.seed(11L)
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

# The moments-normalisation path builds a synthetic mixed distribution from a
# "core" background component and an "extra" signal-side component.  The size
# of the extra component is controlled by normExtraFrac.  When extra cells are
# suppressed (normExtraFrac = 0) the effective distribution fed to the bandwidth
# estimator is more tightly concentrated, so the returned bandwidth should
# differ from the default.  This test verifies that changing normExtraFrac
# actually changes the result, which would fail if the Norm path silently fell
# back to the ordinary selector for both calls.
test_that("bwCalcOne Norm moments path is sensitive to normExtraFrac", {
  x <- .makeTestVec(seed = 42L, n = 800L)
  set.seed(77L)
  bw_default <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", normMtd = "moments")
  set.seed(77L)
  bw_noextra <- stimgate:::.bwCalcOne(
    x,
    bwMtd = "nrd0Norm",
    normMtd = "moments",
    normExtraFrac = 0
  )
  expect_true(is.finite(bw_default) && bw_default > 0)
  expect_true(is.finite(bw_noextra) && bw_noextra > 0)
  # The two calls use different extra-component sizes, so the resulting
  # synthetic distribution — and therefore the bandwidth — must differ.
  expect_false(isTRUE(all.equal(as.numeric(bw_default), as.numeric(bw_noextra))))
})

# ---------------------------------------------------------------------------
# Adaptive mode
# ---------------------------------------------------------------------------

test_that("bwCalcOne adaptive=TRUE with NormMtd returns list with adaptive attribute", {
  x <- .makeTestVec(n = 1000L)
  set.seed(55L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0Norm", adaptive = TRUE)
  expect_true(isTRUE(attr(bw, "adaptive")))
  expect_true(is.list(bw))
  expect_true(all(c("bin", "bw", "bwCore", "bwExtra") %in% names(bw)))
  expect_true(all(is.finite(bw$bwCore) & bw$bwCore > 0))
  expect_true(all(is.finite(bw$bwExtra) & bw$bwExtra > 0))
  expect_true(all(is.finite(bw$bw) & bw$bw > 0))
})

# bwAdaptiveCore / bwAdaptiveExtra override the estimated core and extra
# component bandwidths.  The implementation substitutes them directly into the
# returned list, so the result must match exactly.
test_that("bwCalcOne adaptive mode honours manual bwAdaptiveCore and bwAdaptiveExtra", {
  x <- .makeTestVec(n = 1000L)
  manualCore <- 0.3
  manualExtra <- 0.8
  set.seed(55L)
  bw <- stimgate:::.bwCalcOne(
    x,
    bwMtd = "nrd0Norm",
    adaptive = TRUE,
    bwAdaptiveCore = manualCore,
    bwAdaptiveExtra = manualExtra
  )
  expect_true(isTRUE(attr(bw, "adaptive")))
  # The returned core and extra bandwidths must equal the manual values.
  expect_equal(bw$bwCore, manualCore, tolerance = 1e-8)
  expect_equal(bw$bwExtra, manualExtra, tolerance = 1e-8)
  # The per-bin blended bandwidth must be bounded by the two manual values.
  expect_true(all(bw$bw >= min(manualCore, manualExtra) - 1e-8))
  expect_true(all(bw$bw <= max(manualCore, manualExtra) + 1e-8))
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

# bwNcellMax caps the number of cells used for bandwidth estimation.  With a
# large input (5000) and a tiny cap (20), the estimator runs on a 20-cell
# subsample.  The resulting bandwidth should differ from the uncapped result
# because the sample the estimator sees is different in size and composition.
test_that("bwCalcOne ordinary path bwNcellMax changes the result vs uncapped", {
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
  # For rnorm(500, mean=0, sd=1), Silverman's rule-of-thumb gives ~ 0.2-0.5.
  # This is a loose sanity check, not an exact regression fixture.
  set.seed(42L)
  x <- stats::rnorm(500L)
  bw <- stimgate:::.bwCalcOne(x, bwMtd = "nrd0")
  expect_gt(as.numeric(bw), 0.1)
  expect_lt(as.numeric(bw), 1.0)
})
