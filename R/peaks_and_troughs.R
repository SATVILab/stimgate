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
    troughMaxRel = 0.75) {
  y <- suppressWarnings(as.numeric(y))
  if (length(y) == 0L || all(!is.finite(y))) {
    return(integer(0L))
  }

  y <- pmax(y, 0)
  peakIdxAll <- .getLocalMaximaIdx(y)

  if (length(peakIdxAll) == 0L) {
    out <- max(which(y == max(y, na.rm = TRUE)))
    return(out)
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

  # First deep trough separating meaningful peaks (peaks are sorted, unique and
  # in range, as they are a subset of `peakIdxAll`).
  nextTroughIdx <- integer(0L)
  for (troughIdx in .getLocalMaximaIdx(-y)) {
    leftPeak <- peakIdxMeaningful[peakIdxMeaningful < troughIdx]
    rightPeak <- peakIdxMeaningful[peakIdxMeaningful > troughIdx]

    if (length(leftPeak) == 0L || length(rightPeak) == 0L) {
      next
    }

    troughHeight <- y[troughIdx]
    lowEnoughAdjacent <-
      troughHeight <= troughMaxRel * y[max(leftPeak)] &&
        troughHeight <= troughMaxRel * y[min(rightPeak)]
    lowEnoughAbsolute <-
      is.finite(peakHeightMax) &&
        peakHeightMax > 0 &&
        troughHeight <= troughMaxRel * peakHeightMax

    if (isTRUE(lowEnoughAdjacent) && isTRUE(lowEnoughAbsolute)) {
      nextTroughIdx <- troughIdx
      break
    }
  }

  if (length(nextTroughIdx) == 0L) {
    return(peakIdxMeaningful[length(peakIdxMeaningful)])
  }

  peakBefore <- peakIdxMeaningful[peakIdxMeaningful < nextTroughIdx]
  if (length(peakBefore) == 0L) {
    return(peakIdxMeaningful[1L])
  }

  max(peakBefore)
}
