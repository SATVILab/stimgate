.getLocalMaximaIdx <- function(y) {
  y <- suppressWarnings(as.numeric(y))
  if (length(y) < 3L) {
    return(integer(0L))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[
    is.finite(y[idx]) &
      is.finite(y[idx - 1L]) &
      is.finite(y[idx + 1L]) &
      y[idx] >= y[idx - 1L] &
      y[idx] > y[idx + 1L]
  ]
}

.getLocalMinimaIdx <- function(y) {
  y <- suppressWarnings(as.numeric(y))
  if (length(y) < 3L) {
    return(integer(0L))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[
    is.finite(y[idx]) &
      is.finite(y[idx - 1L]) &
      is.finite(y[idx + 1L]) &
      y[idx] <= y[idx - 1L] &
      y[idx] < y[idx + 1L]
  ]
}

#' Return the right-most peak belonging to the left/main modal complex.
#'
#' Peaks whose height is at least `peakMinRel * max(y)` are treated as
#' meaningful. If meaningful peaks are separated by a trough that is low relative
#' to both adjacent peaks and to the absolute peak, the first such trough ends
#' the left/main modal complex. Otherwise, shoulders and unresolved peaks are
#' allowed to belong to the same background complex.
#'
#' @keywords internal
.getPeakMainLeftIdx <- function(
  y,
  peakMinRel = 0.75,
  troughMaxRel = 0.75
) {
  y <- suppressWarnings(as.numeric(y))
  if (length(y) == 0L || all(!is.finite(y))) {
    return(integer(0L))
  }

  y <- pmax(y, 0)
  peakIdxAll <- .getLocalMaximaIdx(y)

  if (length(peakIdxAll) == 0L) {
    return(which.max(y))
  }
  if (length(peakIdxAll) == 1L) {
    return(peakIdxAll)
  }

  peakHeightMax <- max(y[peakIdxAll], na.rm = TRUE)
  if (!is.finite(peakHeightMax) || peakHeightMax <= 0) {
    return(which.max(y))
  }

  peakIdxMeaningful <- peakIdxAll[
    y[peakIdxAll] >= peakMinRel * peakHeightMax
  ]

  if (length(peakIdxMeaningful) == 0L) {
    return(which.max(y))
  }
  if (length(peakIdxMeaningful) == 1L) {
    return(peakIdxMeaningful)
  }

  nextTroughIdx <- .getPeakIdxNextTroughIdx(
    y = y,
    peakIdxMeaningful = peakIdxMeaningful,
    peakMinRel = troughMaxRel,
    peakHeightRef = peakHeightMax
  )

  if (length(nextTroughIdx) == 0L) {
    return(peakIdxMeaningful[length(peakIdxMeaningful)])
  }

  peakBefore <- peakIdxMeaningful[peakIdxMeaningful < nextTroughIdx]
  if (length(peakBefore) == 0L) {
    return(peakIdxMeaningful[1L])
  }

  max(peakBefore)
}

#' Return the first deep trough separating meaningful peaks.
#'
#' @keywords internal
.getPeakIdxNextTroughIdx <- function(
  y,
  peakIdxMeaningful,
  peakMinRel = 0.75,
  peakHeightRef = NULL
) {
  y <- suppressWarnings(as.numeric(y))
  y <- pmax(y, 0)
  peakIdxMeaningful <- sort(unique(as.integer(peakIdxMeaningful)))
  peakIdxMeaningful <- peakIdxMeaningful[
    is.finite(peakIdxMeaningful) &
      peakIdxMeaningful >= 1L &
      peakIdxMeaningful <= length(y)
  ]

  if (length(peakIdxMeaningful) < 2L) {
    return(integer(0L))
  }

  troughIdxAll <- .getLocalMinimaIdx(y)
  if (length(troughIdxAll) == 0L) {
    return(integer(0L))
  }

  if (is.null(peakHeightRef)) {
    peakHeightRef <- max(y[peakIdxMeaningful], na.rm = TRUE)
  }

  for (troughIdx in troughIdxAll) {
    leftPeak <- peakIdxMeaningful[peakIdxMeaningful < troughIdx]
    rightPeak <- peakIdxMeaningful[peakIdxMeaningful > troughIdx]

    if (length(leftPeak) == 0L || length(rightPeak) == 0L) {
      next
    }

    leftPeakIdx <- max(leftPeak)
    rightPeakIdx <- min(rightPeak)

    troughHeight <- y[troughIdx]
    leftPeakHeight <- y[leftPeakIdx]
    rightPeakHeight <- y[rightPeakIdx]

    lowEnoughAdjacent <-
      troughHeight <= peakMinRel * leftPeakHeight &&
      troughHeight <= peakMinRel * rightPeakHeight

    lowEnoughAbsolute <-
      is.finite(peakHeightRef) &&
      peakHeightRef > 0 &&
      troughHeight <= peakMinRel * peakHeightRef

    if (isTRUE(lowEnoughAdjacent) && isTRUE(lowEnoughAbsolute)) {
      return(troughIdx)
    }
  }

  integer(0L)
}
