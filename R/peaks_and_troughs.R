.getLocalMaximaIdx <- function(y) {
  if (length(y) < 3L) {
    return(integer(0))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[y[idx] >= y[idx - 1L] & y[idx] > y[idx + 1L]]
}

.getLocalMinimaIdx <- function(y) {
  if (length(y) < 3L) {
    return(integer(0))
  }

  idx <- seq.int(2L, length(y) - 1L)
  idx[y[idx] <= y[idx - 1L] & y[idx] < y[idx + 1L]]
}

.getPeakMainLeftIdx <- function(y, peakMinRel = 0.75) {
  # get furthest left meaningful peak,
  # which is meaningful in the sense
  # that it can't truly be separated from
  # any peaks afterwards by being after a local minimum that is too low relative to the main peak
  # function assumes that y is on the ratio scale, i.e. that the minimum value is 0
  stopifnot(all(y >= 0, na.rm = TRUE))
  peakIdxAll <- .getLocalMaximaIdx(y)

  # return early if there aren't multiple peaks

  if (length(peakIdxAll) == 0L) {
    return(which.max(y))
  }
  if (length(peakIdxAll) == 1L) {
    return(peakIdxAll)
  }

  peakHeightMax <- max(y, na.rm = TRUE)
  peakIdxMeaningful <- peakIdxAll[
    y[peakIdxAll] >= peakMinRel * peakHeightMax
  ]

  # return early if there aren't multiple
  # meaningful peaks
  if (length(peakIdxMeaningful) == 0L) {
    return(which.max(y))
  }
  if (length(peakIdxMeaningful) == 1L) {
    return(peakIdxMeaningful)
  }

  # multiple peaks, choose right-most
  # peak that isn't after an antimode (local minimum) that is
  # too low relative to the main peak
  nextTroughIdx <- .getPeakIdxNextTrough(
    y = y,
    peakIdxMeaningful = peakIdxMeaningful,
    peakMinRel = peakMinRel
  )
  # if no such trough exists, return the right-most meaningful peak
  if (length(nextTroughIdx) == 0L) {
    return(peakIdxMeaningful[length(peakIdxMeaningful)])
  }

  # if there is a meaningful trough, return the right-most peak that is before it
  max(peakIdxMeaningful[peakIdxMeaningful < nextTroughIdx], na.rm = TRUE)
}

.getPeakIdxNextTroughIdx <- function(
  y,
  peakIdxMeaningful,
  peakMinRel = 0.75
) {
  stopifnot(all(y >= 0, na.rm = TRUE))
  # get lowest meaningful trough, which is meaningful
  # in the sense that it is the first trough
  # that is low relative to the adjacent peaks on either side. If there is no such trough, return integer(0)
  troughIdxAll <- .getLocalMinimaIdx(y)

  # return early if there aren't any troughs
  if (length(troughIdxAll) == 0L) {
    return(integer(0L))
  }

  # for each trough, look at the peaks on either side.
  # take a trough as "meaningful" if it is less than
  # peakMinRel * the height of the adjacent peaks on either side. If there are no meaningful troughs, return integer(0)
  troughIdxMeaningful <- NULL
  for (troughIdx in troughIdxAll) {
    leftPeakIdx <- min(
      peakIdxMeaningful[peakIdxMeaningful < troughIdx],
      na.rm = TRUE
    )
    rightPeakIdx <- max(
      peakIdxMeaningful[peakIdxMeaningful > troughIdx],
      na.rm = TRUE
    )
    if (is.finite(leftPeakIdx) && is.finite(rightPeakIdx)) {
      leftPeakHeight <- y[leftPeakIdx]
      rightPeakHeight <- y[rightPeakIdx]
      troughHeight <- y[troughIdx]
      if (
        troughHeight < 0.75 * leftPeakHeight ||
          troughHeight < 0.75 * rightPeakHeight
      ) {
        troughIdxMeaningful <- c(troughIdxMeaningful, troughIdx)
      }
    }
  }
  if (length(troughIdxMeaningful) == 0L) {
    return(integer(0L))
  }
  min(troughIdxAll, na.rm = TRUE)
}
