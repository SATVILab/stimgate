pkg_ns <- asNamespace("stimgate")

.getCpUnsLocFilterAfterSmoothing <- get(
  ".getCpUnsLocFilterAfterSmoothing",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocDominanceBoundaryCurrent <- get(
  ".getCpUnsLocDominanceBoundaryCurrent",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocQualityBoundaryCurrent <- get(
  ".getCpUnsLocQualityBoundaryCurrent",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocAntimodeBoundaryCurrent <- get(
  ".getCpUnsLocAntimodeBoundaryCurrent",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocHighProbabilityReference <- get(
  ".getCpUnsLocHighProbabilityReference",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocTautStringExtremaCurrent <- get(
  ".getCpUnsLocTautStringExtremaCurrent",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocFilterMarginalBins <- get(
  ".getCpUnsLocFilterMarginalBins",
  envir = pkg_ns,
  mode = "function"
)

test_that("ordinary post-smoothing filter respects boundary hierarchy", {
  x_grid <- seq(0, 5, length.out = 100)
  stim_dens <- stats::dnorm(x_grid, mean = 3.5, sd = 0.8)
  unstim_dens <- stats::dnorm(x_grid, mean = 1.0, sd = 0.6)
  dens_comp <- data.frame(x = x_grid, stim = stim_dens, unstim = unstim_dens)

  x_vals <- seq(0.2, 4.8, length.out = 50)
  prob_vals <- 1 / (1 + exp(-(x_vals - 2.5) * 2))

  data_mod <- data.frame(
    IFNg = x_vals,
    probSmooth = prob_vals,
    pred = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "IFNg"
  attr(data_mod, "idxMod") <- seq_along(x_vals)
  attr(data_mod, "ind") <- 1L
  attr(data_mod, "minProbXPos") <- 0.5
  attr(data_mod, "locDensityBw") <- 0.4
  attr(data_mod, "locDensityComparison") <- dens_comp

  deriv_x <- seq(0.2, 4.8, length.out = 100)
  deriv_pred <- 1 / (1 + exp(-(deriv_x - 2.5) * 2))
  deriv_val <- 2 * deriv_pred * (1 - deriv_pred)
  attr(data_mod, "locProbDerivTbl") <- tibble::tibble(
    x = deriv_x,
    pred = deriv_pred,
    deriv = deriv_val
  )

  res <- .getCpUnsLocFilterAfterSmoothing(
    dataMod = data_mod,
    exTblStimNoMin = data_mod,
    exTblUnsBias = data_mod,
    cpMin = NULL,
    stage = "init",
    chnlSettings = list(),
    exTblStimOrig = data_mod,
    exTblUnsOrig = data_mod
  )

  expect_s3_class(res$dataMod, "data.frame")
  expect_true(is.null(res$cp))
  expect_true(is.list(res$info))

  info_final <- res$info$final
  expect_true(is.list(info_final))
  expect_true(is.finite(info_final$xClearInit))
  expect_true(is.finite(info_final$xClear))
  expect_true(is.finite(info_final$xQual))
  expect_true(is.finite(info_final$xBase))
  expect_true(is.finite(info_final$xSum))

  # Boundary hierarchy: xSum <= xBase <= xClear <= xClearInit
  expect_true(info_final$xClear <= info_final$xClearInit + 1e-6)
  expect_true(info_final$xBase <= info_final$xClear + 1e-6)
  expect_true(info_final$xSum <= info_final$xBase + 1e-6)

  # Distinguish from legacy global filter: no global filter applied
  expect_false(info_final$globalFilterApplied)
  expect_false(grepl("unavailable", res$info$marginal$trimReason))
  expect_equal(
    res$info$reason,
    "filtered_at_lowest_supported_post_smoothing_boundary"
  )

  # Retained expression values satisfy x >= xSum
  final_x <- attr(res$dataMod, "locFinalFilterX")
  expect_equal(final_x, info_final$xSum)
  expect_true(all(res$dataMod$IFNg >= final_x))
})

test_that("density dominance moves boundary left through contiguous region", {
  x_grid <- seq(0, 5, length.out = 100)
  # Stimulated density dominant for x >= 1.8, non-dominant for x < 1.8
  stim_dens <- stats::dnorm(x_grid, mean = 3.5, sd = 0.8)
  unstim_dens <- stats::dnorm(x_grid, mean = 1.0, sd = 0.6)
  dens_comp <- data.frame(x = x_grid, stim = stim_dens, unstim = unstim_dens)

  # 1. startX inside dominant region (x = 3.5)
  dom_valid <- .getCpUnsLocDominanceBoundaryCurrent(
    density = dens_comp,
    startX = 3.5,
    densityBw = 0.5,
    lowerBoundX = 0.5
  )

  expect_true(dom_valid$info$applied)
  expect_equal(
    dom_valid$info$reason,
    "identified_contiguous_density_dominance_boundary"
  )
  expect_true(is.finite(dom_valid$startX))
  expect_true(dom_valid$startX < 3.5)
  expect_true(dom_valid$startX >= 0.5)

  # 2. startX in non-dominant region (x = 0.8)
  dom_nondom <- .getCpUnsLocDominanceBoundaryCurrent(
    density = dens_comp,
    startX = 0.8,
    densityBw = 0.5,
    lowerBoundX = 0.5
  )

  expect_false(dom_nondom$info$applied)
  expect_equal(
    dom_nondom$info$reason,
    "density_not_dominant_at_clear_response_reference"
  )
  expect_true(is.na(dom_nondom$startX))

  # 3. Invalid inputs guard
  expect_true(is.na(.getCpUnsLocDominanceBoundaryCurrent(NULL, 3.0)$startX))
  expect_true(
    is.na(
      .getCpUnsLocDominanceBoundaryCurrent(
        density = dens_comp,
        startX = NA_real_
      )$startX
    )
  )
})

test_that("quality boundary respects preliminary lower bound", {
  x_vals <- seq(0, 5, length.out = 40)
  prob_vals <- 1 / (1 + exp(-(x_vals - 2.5) * 2))

  data_mod <- data.frame(
    TNF = x_vals,
    probSmooth = prob_vals,
    pred = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "TNF"
  attr(data_mod, "idxMod") <- seq_along(x_vals)

  qual_out <- .getCpUnsLocQualityBoundaryCurrent(
    dataMod = data_mod,
    chnlSettings = list(),
    probCol = "pred",
    xClear = 3.5,
    lowerBoundX = 1.5
  )

  expect_true(is.finite(qual_out$thresholdX))
  expect_true(qual_out$thresholdX >= 1.5)
  expect_equal(qual_out$info$preliminaryLowerBoundX, 1.5)
  expect_equal(qual_out$info$referenceBasis, "x_clear")
  expect_true(all(qual_out$dataMod$TNF >= qual_out$thresholdX))
})

test_that("antimode boundary moves xBase lower only when deep trough exists", {
  withr::local_preserve_seed()
  # Bimodal distribution with deep separation at x ~ 2.5
  set.seed(42)
  x_bimodal <- c(
    stats::rnorm(80, mean = 1.0, sd = 0.3),
    stats::rnorm(80, mean = 4.0, sd = 0.3)
  )
  x_bimodal <- sort(x_bimodal[x_bimodal >= 0 & x_bimodal <= 5])

  data_bimodal <- data.frame(
    CD8 = x_bimodal,
    probSmooth = seq(0, 1, length.out = length(x_bimodal)),
    stringsAsFactors = FALSE
  )
  attr(data_bimodal, "chnlCut") <- "CD8"
  attr(data_bimodal, "idxMod") <- seq_along(x_bimodal)
  attr(data_bimodal, "locDensityBw") <- 0.3

  # xBase = 3.5 (above trough, near mode of right component)
  antimode_valid <- .getCpUnsLocAntimodeBoundaryCurrent(
    dataMod = data_bimodal,
    chnlSettings = list(),
    xBase = 3.5,
    lowerBoundX = 0.2,
    heightFrac = 0.95
  )

  expect_true(antimode_valid$info$applied)
  expect_equal(
    antimode_valid$info$reason,
    "selected_rightmost_valid_antimode"
  )
  expect_true(is.finite(antimode_valid$thresholdX))
  expect_true(antimode_valid$thresholdX < 3.5)
  expect_true(antimode_valid$thresholdX > 1.5)

  # When xBase is below all troughs (xBase = 0.5), no antimode is below xBase
  antimode_low <- .getCpUnsLocAntimodeBoundaryCurrent(
    dataMod = data_bimodal,
    chnlSettings = list(),
    xBase = 0.5,
    lowerBoundX = 0.2,
    heightFrac = 0.95
  )

  expect_false(antimode_low$info$applied)
  expect_equal(
    antimode_low$info$reason,
    "no_antimode_below_supported_boundary"
  )
  expect_true(is.na(antimode_low$thresholdX))
})

test_that("high probability reference locates 85pct target point", {
  x_vals <- seq(0, 5, length.out = 50)
  prob_vals <- seq(0, 1, length.out = 50)

  data_mod <- data.frame(
    IL4 = x_vals,
    probSmooth = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "IL4"

  high_ref <- .getCpUnsLocHighProbabilityReference(
    dataMod = data_mod,
    probCol = "probSmooth",
    fraction = 0.85
  )

  expect_equal(high_ref$info$reason, "used_probability_85pct_reference")
  expect_true(is.finite(high_ref$thresholdX))
  # 85% of max (1.0) on linear prob over [0, 5] occurs at x = 4.25
  expect_equal(high_ref$thresholdX, 4.25)
  expect_equal(high_ref$info$peakProb, 1.0)
  expect_equal(high_ref$info$targetProb, 0.85)

  # Empty/short input returns NA
  data_empty <- data.frame(
    IL4 = numeric(0),
    probSmooth = numeric(0)
  )
  attr(data_empty, "chnlCut") <- "IL4"
  expect_true(
    is.na(
      .getCpUnsLocHighProbabilityReference(
        dataMod = data_empty,
        probCol = "probSmooth"
      )$thresholdX
    )
  )
})

test_that("taut string extrema helper identifies modes and antimodes", {
  # Synthetic piecewise-constant density with bimodal profile:
  # Mode 1 around x = 1.5, trough around x = 3.5, Mode 2 around x = 5.5
  x_pts <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)
  y_pts <- c(0.1, 0.8, 0.8, 0.2, 0.2, 0.9, 0.1)

  density_mock <- list(
    x = x_pts,
    y = y_pts,
    method = "taut_string"
  )

  extrema <- .getCpUnsLocTautStringExtremaCurrent(density_mock)
  expect_named(extrema, c("modes", "antimodes"))

  expect_identical(extrema, list(
    modes = data.frame(
      x = c(2, 5.5), height = c(0.8, 0.9), row.names = c("2", "4")
    ),
    antimodes = data.frame(x = 4, height = 0.2, row.names = "3")
  ))

  # Sorting and removal of non-finite pairs preserve plateau centres and labels.
  shuffled <- list(
    x = c(rev(x_pts), Inf, NA_real_),
    y = c(rev(y_pts), 0.5, 0.5),
    method = "taut_string"
  )
  expect_identical(.getCpUnsLocTautStringExtremaCurrent(shuffled), extrema)

  # Near-equal adjacent heights still form one plateau, using its first height.
  density_mock$y[[3L]] <- 0.8 + .Machine$double.eps
  expect_identical(.getCpUnsLocTautStringExtremaCurrent(density_mock), extrema)

  constant <- list(x = 1:5, y = rep(0.2, 5), method = "taut_string")
  empty <- data.frame(x = numeric(), height = numeric())
  expect_identical(
    .getCpUnsLocTautStringExtremaCurrent(constant),
    list(modes = empty, antimodes = empty)
  )

  # Invalid inputs return empty data frames
  empty_extrema <- .getCpUnsLocTautStringExtremaCurrent(NULL)
  expect_equal(nrow(empty_extrema$modes), 0L)
  expect_equal(nrow(empty_extrema$antimodes), 0L)
})

# Unit-spaced bins make scan boundaries and half-open span counts explicit.
.marginalExpression <- function(x) {
  out <- data.frame(marker = x)
  attr(out, "chnlCut") <- "marker"
  out
}

.marginalModel <- function(leftX, leftProb) {
  out <- .marginalExpression(c(leftX, seq(10.5, 48.5, by = 1)))
  out$pred <- c(leftProb, rep(1, 39))
  attr(out, "idxMod") <- seq_len(nrow(out))
  attr(out, "binVec") <- c(0, 49)
  out
}

test_that("an empty gap followed only by rejected cells leaves the cut at start", {
  # C57_M2-like: responding stim cells are above start, negatives below a gap,
  # and raw control tail cells lie in that gap. Empty bins provide no support.
  dm <- .marginalModel(c(3.5, 4.5, 5.5), rep(0.003, 3))
  out <- .getCpUnsLocFilterMarginalBins(
    dm, list(), "pred", startX = 10,
    exTblStimOrig = dm,
    exTblUnsOrig = .marginalExpression(c(0, 6.2, 7.1))
  )
  expect_equal(out$info$finalStartX, 10)
  expect_equal(out$dataMod$marker, seq(10.5, 48.5, by = 1))
  expect_equal(out$info$stopReason, "three_consecutive_rejections")
  expect_equal(sum(out$info$scanTbl$pending), 4)
  expect_false(any(out$info$scanTbl$accepted))
  expect_false(any(out$info$scanTbl$retained))
  expect_equal(nrow(out$info$trimTbl), 0L)
  expect_equal(out$info$trimReason, "no_acceptance_steps")
})

test_that("a non-empty acceptance retains a preceding empty gap only", {
  dm <- .marginalModel(c(0.5, 1.5, 3.5, 5.5, 6.5), c(1, 0, 0, 0, 1))
  out <- .getCpUnsLocFilterMarginalBins(
    dm, list(), "pred", startX = 10,
    exTblStimOrig = dm, exTblUnsOrig = .marginalExpression(c(0, 0))
  )
  expect_equal(out$info$finalStartX, 6)
  expect_equal(out$dataMod$marker, c(6.5, seq(10.5, 48.5, by = 1)))
  scan <- out$info$scanTbl
  expect_equal(scan$left[scan$accepted], 6)
  expect_true(all(scan$retained[scan$left >= 6]))
  # Empty bins after the last acceptance remain pending and are not retained.
  expect_false(any(scan$retained[scan$left < 6]))
  expect_false(any(scan$accepted[scan$pending]))
  # Pending bins between rejections do not reset their counter.
  expect_equal(sum(!scan$pending & !scan$accepted), 3)
  expect_equal(tail(scan$left, 1), 1)
  expect_equal(out$info$stopReason, "three_consecutive_rejections")
})

test_that("marginal trim undoes the last span and keeps the next positive span", {
  dm <- .marginalModel(c(1.5, 3.5, 5.5, 6.5, 8.5, 9.5), c(0, 0, 0, 1, 0, 1))
  # Extra raw cells are absent from dataMod, but enter counts and denominators.
  # Exact edges check [cutTo, cutFrom): 6 and 9 enter, 10 does not.
  stim <- .marginalExpression(c(dm$marker, rep(0, 100), 6, 9, 10))
  uns <- .marginalExpression(c(rep(0, 97), 6, 6.8, 8.9, 10))
  out <- .getCpUnsLocQualityBoundaryCurrent(
    dm, list(), "pred", xClear = 10,
    exTblStimOrig = stim, exTblUnsOrig = uns
  )
  expect_equal(out$thresholdX, 9)
  expect_equal(out$info$finalStartX, 9)
  expect_equal(out$dataMod$marker, c(9.5, seq(10.5, 48.5, by = 1)))
  trim <- out$info$trimTbl
  expect_equal(trim$step, 1:2)
  expect_equal(trim$cutFrom, c(10, 9))
  expect_equal(trim$cutTo, c(9, 6))
  expect_equal(trim$nStim, c(2, 3))
  expect_equal(trim$nUns, c(0, 3))
  expect_equal(trim$contribution, c(2, 3) / nrow(stim) - c(0, 3) / nrow(uns))
  expect_identical(trim$kept, c(TRUE, FALSE))
  expect_equal(out$info$nStepsTrimmed, 1L)
  expect_equal(out$info$trimReason, "positive_span_contribution")
  # Scan acceptance remains recorded even when trimming removes that bin.
  expect_equal(out$info$scanTbl$left[out$info$scanTbl$accepted], c(9, 6))
  expect_equal(out$info$scanTbl$left[out$info$scanTbl$retained], 9)
})

test_that("non-positive marginal spans can all be trimmed to the start", {
  dm <- .marginalModel(c(1.5, 3.5, 5.5, 6.5, 8.5, 9.5), c(0, 0, 0, 1, 0, 1))
  stim <- .marginalExpression(c(dm$marker, rep(0, 100)))
  for (uns in list(.marginalExpression(c(6.5, 8.5, 9.5)), stim)) {
    out <- .getCpUnsLocFilterMarginalBins(
      dm, list(), "pred", startX = 10,
      exTblStimOrig = stim, exTblUnsOrig = uns
    )
    expect_equal(out$info$finalStartX, 10)
    expect_equal(out$dataMod$marker, seq(10.5, 48.5, by = 1))
    expect_equal(out$info$nStepsTrimmed, 2L)
    expect_equal(out$info$trimReason, "trimmed_all_steps")
    expect_true(all(out$info$trimTbl$contribution <= 0))
    expect_false(any(out$info$trimTbl$kept))
    expect_false(any(out$info$scanTbl$retained))
  }
})

test_that("unavailable unstimulated expression skips trim but preserves scan steps", {
  dm <- .marginalModel(c(1.5, 3.5, 5.5, 6.5, 8.5, 9.5), c(0, 0, 0, 1, 0, 1))
  out <- .getCpUnsLocFilterMarginalBins(
    dm, list(), "pred", startX = 10,
    exTblStimOrig = dm, exTblUnsOrig = NULL
  )
  expect_equal(out$info$finalStartX, 6)
  expect_equal(out$info$nStepsTrimmed, 0L)
  expect_equal(out$info$trimReason, "unstimulated_expression_unavailable")
  expect_true(all(out$info$trimTbl$kept))
  expect_equal(out$info$trimTbl$nStim, c(1, 2))
  expect_true(all(is.na(out$info$trimTbl$nUns)))
  expect_true(all(is.na(out$info$trimTbl$contribution)))
  # The rejected bin at 8 and empty bin at 7 were pulled in by acceptance at 6.
  expect_true(all(out$info$scanTbl$retained[out$info$scanTbl$left %in% 6:9]))
})


test_that("acceptance retains two rejected non-empty bins between steps", {
  dm <- .marginalModel(c(1.5, 3.5, 5.5, 6.5, 7.5, 8.5, 9.5), c(0, 0, 0, 1, 0, 0, 1))
  out <- .getCpUnsLocFilterMarginalBins(dm, list(), "pred", startX = 10)
  expect_equal(out$info$finalStartX, 6)
  expect_equal(out$dataMod$marker, c(6.5, 7.5, 8.5, 9.5, seq(10.5, 48.5, by = 1)))
  scan <- out$info$scanTbl
  expect_false(any(scan$accepted[scan$left %in% 7:8]))
  expect_true(all(scan$retained[scan$left %in% 7:8]))
})

test_that("the shape-route marginal helper uses the same raw-span trim", {
  dm <- .marginalModel(c(1.5, 3.5, 5.5, 6.5, 8.5, 9.5), c(0, 0, 0, 1, 0, 1))
  out <- pkg_ns$.getCpUnsLocFilterMarginal(
    dataMod = dm, chnlSettings = list(), probCol = "pred",
    threshold = list(thresholdX = 10, info = list()),
    dominance = list(startX = NA_real_, info = list()),
    exTblStimOrig = dm, exTblUnsOrig = dm
  )
  expect_equal(out$info$finalStartX, 10)
  expect_equal(out$info$nStepsTrimmed, 2L)
  expect_equal(out$info$trimReason, "trimmed_all_steps")
})
