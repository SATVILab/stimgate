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
    chnlSettings = list()
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

  # Identifies 2 modes and 1 antimode
  expect_equal(nrow(extrema$modes), 2L)
  expect_equal(nrow(extrema$antimodes), 1L)

  # Antimode is located in the middle trough
  expect_true(extrema$antimodes$x[[1L]] > extrema$modes$x[[1L]])
  expect_true(extrema$antimodes$x[[1L]] < extrema$modes$x[[2L]])

  # Invalid inputs return empty data frames
  empty_extrema <- .getCpUnsLocTautStringExtremaCurrent(NULL)
  expect_equal(nrow(empty_extrema$modes), 0L)
  expect_equal(nrow(empty_extrema$antimodes), 0L)
})
