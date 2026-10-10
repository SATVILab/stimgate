#' Get the safe lower boundary from the full stimulated marginal distribution
#'
#' The left/main modal-complex peak is identified in the same way as in the
#' initial one-marker procedure. `windowWidth` is the span from the 5th
#' percentile of values below that peak to the peak itself. Cytokine-positive
#' refinement is not allowed at or below `peakX + windowWidth / 3`.
#'
#' @keywords internal
.getCytPosMarginalReference <- function(
  ex,
  chnl,
  bwMin = NA_real_
) {
  out <- list(
    peakX = NA_real_,
    windowWidth = NA_real_,
    lowerX = NA_real_,
    densityBw = NA_real_,
    reason = "marginal_reference_unavailable"
  )

  x <- suppressWarnings(
    as.numeric(ex[[chnl]])
  )
  x <- x[is.finite(x)]

  if (length(x) > 0L) {
    x <- x[x > min(x)]
  }

  if (
    length(x) < 5L ||
      length(unique(x)) < 3L
  ) {
    out$reason <- "too_few_marginal_values"
    return(out)
  }

  # Calculate the SJ bandwidth first so that the bwMin floor can be applied
  # before evaluating the density. This avoids evaluating the density twice
  # when the SJ bandwidth is below bwMin.
  bw <- try(
    suppressWarnings(
      stats::bw.SJ(x)
    ),
    silent = TRUE
  )

  if (
    inherits(bw, "try-error") ||
      length(bw) != 1L ||
      !is.finite(bw) ||
      bw <= 0
  ) {
    out$reason <- "marginal_density_failed"
    return(out)
  }

  bwMin <- bwMin[1L]

  if (
    is.finite(bwMin) &&
      bwMin > 0
  ) {
    bw <- max(bw, bwMin)
  }

  dens <- try(
    suppressWarnings(
      stats::density(
        x,
        bw = bw
      )
    ),
    silent = TRUE
  )

  if (inherits(dens, "try-error")) {
    out$reason <- "marginal_density_failed"
    return(out)
  }

  peakIdx <- .getPeakMainLeftIdx(
    dens$y
  )

  if (
    length(peakIdx) != 1L ||
      !is.finite(peakIdx)
  ) {
    out$reason <- "marginal_peak_unavailable"
    return(out)
  }

  peakX <- dens$x[peakIdx]

  if (!is.finite(peakX)) {
    out$reason <- "marginal_left_region_unavailable"
    return(out)
  }

  xLeft <- x[x < peakX]

  if (length(xLeft) < 2L) {
    out$reason <- "marginal_left_region_unavailable"
    return(out)
  }

  windowWidth <- abs(
    diff(
      stats::quantile(
        xLeft,
        probs = c(0.05, 1),
        na.rm = TRUE,
        names = FALSE
      )
    )
  )

  if (
    !is.finite(windowWidth) ||
      windowWidth <= 0
  ) {
    out$reason <- "marginal_window_width_unavailable"
    return(out)
  }

  out$peakX <- peakX
  out$windowWidth <- windowWidth
  out$lowerX <- peakX + windowWidth / 3
  out$densityBw <- dens$bw
  out$reason <- "marginal_reference_available"

  out
}

#' Locate internal troughs in a piecewise-constant taut-string density
#'
#' Flat troughs are represented by the midpoint of the complete flat interval.
#' Boundary runs are never treated as antimodes.
#'
#' @keywords internal
.getCytPosTautStringAntimodes <- function(density) {
  if (is.null(density)) {
    return(numeric(0L))
  }

  x <- suppressWarnings(as.numeric(density$x))
  y <- suppressWarnings(as.numeric(density$y))

  finite <- is.finite(x) & is.finite(y)
  x <- x[finite]
  y <- y[finite]

  if (length(x) != length(y) || length(y) < 3L) {
    return(numeric(0L))
  }

  change <- c(
    TRUE,
    !dplyr::near(
      y[-1L],
      y[-length(y)]
    )
  )

  if (sum(change) < 3L) {
    return(numeric(0L))
  }

  .getCpUnsLocPiecewiseConstantAntimodes(x, y)
}

#' Select a cytokine-positive refinement threshold from a taut-string density
#'
#' The leftmost internal antimode strictly between the marginal lower boundary
#' and the existing clustered gate is selected. If no such antimode exists, the
#' existing gate is retained by the caller.
#'
#' @keywords internal
.getCpPosTautString <- function(
  ex,
  inc,
  chnl,
  cpOrig,
  peakX,
  windowWidth,
  lower
) {
  out <- list(
    threshold = NA_real_,
    peakX = suppressWarnings(as.numeric(peakX))[1L],
    windowWidth = suppressWarnings(as.numeric(windowWidth))[1L],
    lowerX = NA_real_,
    gateOriginal = suppressWarnings(as.numeric(cpOrig))[1L],
    antimodes = numeric(0L),
    eligibleAntimodes = numeric(0L),
    nOtherCytPos = 0L,
    reason = "taut_string_threshold_unavailable"
  )

  if (
    !is.finite(out$peakX) ||
      !is.finite(out$windowWidth) ||
      out$windowWidth <= 0 ||
      !is.finite(out$gateOriginal)
  ) {
    out$reason <- "invalid_refinement_interval"
    return(out)
  }
  out$lowerX <- lower
  if (out$lowerX >= out$gateOriginal) {
    out$reason <- "empty_refinement_interval"
    return(out)
  }

  xAll <- suppressWarnings(as.numeric(ex[[chnl]]))
  finiteAll <- is.finite(xAll)
  if (!any(finiteAll)) {
    out$reason <- "no_finite_expression_values"
    return(out)
  }
  minX <- min(xAll[finiteAll])
  keep <- inc %in% TRUE & finiteAll & xAll > minX
  xPos <- xAll[keep]
  out$nOtherCytPos <- length(xPos)
  if (length(xPos) < 10L) {
    out$reason <- "too_few_other_cytokine_positive_cells"
    return(out)
  }

  density <- .getCpUnsLocAntimodeDensity(xPos)
  if (is.null(density)) {
    out$reason <- "taut_string_density_failed"
    return(out)
  }

  antimodes <- .getCytPosTautStringAntimodes(density)
  eligible <- antimodes[
    is.finite(antimodes) &
      antimodes > out$lowerX &
      antimodes < out$gateOriginal
  ]
  out$antimodes <- antimodes
  out$eligibleAntimodes <- eligible

  if (length(eligible) == 0L) {
    out$reason <- "no_internal_antimode_in_refinement_interval"
    return(out)
  }

  out$threshold <- min(eligible)
  out$reason <- "leftmost_eligible_taut_string_antimode"
  out
}


#' @keywords internal
.getCytPosGatesGateTblGet <- function(
  chnlVec,
  pop,
  pathProject,
  chnlLab
) {
  .debug("Getting gateTbl") # nolint
  purrr::map_df(chnlVec, function(chnlCurr) {
    .gatesGetPathAll(
      pathProject = pathProject,
      pop = pop,
      chnlCut = chnlCurr,
      init = TRUE
    ) |>
      readRDS() |>
      dplyr::mutate(chnl = chnlCurr, marker = chnlLab[chnlCurr]) |>
      dplyr::select(
        chnl, marker, gateName, batch, ind, gate,
        dplyr::any_of(c(
          "locGenerated", "locGeneratedDirect", "locSource", "locReason",
          "locResponder", "propBsEst", "locOwnFreq", "locShareLimit",
          "locShareProposed", "locShiftedPeakRef"
        ))
      ) # nolint
  })
}

#' Precompute base-threshold positivity for one sample
#' @keywords internal
.getCytPosBasePos <- function(ex, gateTblInd) {
  gateTblValid <- gateTblInd |>
    dplyr::filter(!is.na(.data$gate))

  chnl <- as.character(gateTblValid$chnl)
  gate <- suppressWarnings(as.numeric(gateTblValid$gate))

  posList <- stats::setNames(
    lapply(seq_along(chnl), function(i) {
      (ex[[chnl[[i]]]] > gate[[i]]) %in% TRUE
    }),
    chnl
  )

  nPos <- integer(nrow(ex))

  for (pos in posList) {
    nPos <- nPos + pos
  }

  list(
    pos = posList,
    nPos = nPos
  )
}
