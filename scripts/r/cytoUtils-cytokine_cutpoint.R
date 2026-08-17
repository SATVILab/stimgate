#' @keywords internal
.cytokineCutpoint <- function(
  x,
  adjust = 1,
  numPeaks = 1,
  refPeak = 1,
  method = c("firstDeriv", "secondDeriv"),
  tol = 1e-2,
  side = "right",
  strict = TRUE,
  plot = FALSE,
  autoTol = FALSE,
  ...
) {
  method <- match.arg(method)

  peaks <- sort(.findPeaks(x, numPeaks = numPeaks, adjust = adjust)[, "x"])
  numPeaks <- length(peaks)

  if (refPeak > numPeaks) {
    msg <- paste(
      "The reference peak is larger than the number of peaks found.",
      "Setting the reference peak to 'numPeaks'..."
    )
    if (strict) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
    refPeak <- numPeaks
  }

  if (method == "firstDeriv") {
    # Calculate the first derivative
    derivOut <- .derivDensity(x = x, adjust = adjust, deriv = 1, ...)

    if (autoTol) {
      tol <- 0.01 * max(abs(derivOut$y))
    }

    if (side == "right") {
      derivValleys <- .findValleys(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust
      )
      derivValleys <- derivValleys[derivValleys > peaks[refPeak]]
      derivValleys <- sort(derivValleys)[1]

      # Safe NA Check
      if (is.na(derivValleys)) {
        cutpoint <- NA_real_
      } else {
        cutpointCandidates <- derivOut$x[
          derivOut$x > derivValleys & abs(derivOut$y) < tol
        ]
        cutpoint <- if (length(cutpointCandidates) > 0) {
          cutpointCandidates[1]
        } else {
          NA_real_
        }
      }
    } else if (side == "left") {
      derivOut$y <- -derivOut$y
      derivValleys <- .findValleys(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust
      )
      derivValleys <- derivValleys[derivValleys < peaks[refPeak]]
      derivValleys <- sort(derivValleys, decreasing = TRUE)[1]

      # Safe NA Check
      if (is.na(derivValleys)) {
        cutpoint <- NA_real_
      } else {
        cutpointCandidates <- derivOut$x[
          derivOut$x < derivValleys & abs(derivOut$y) < tol
        ]
        cutpoint <- if (length(cutpointCandidates) > 0) {
          cutpointCandidates[length(cutpointCandidates)]
        } else {
          NA_real_
        }
      }
    } else {
      stop("Unrecognized 'side' argument (was '", side, "').")
    }
  } else {
    # secondDeriv method
    derivOut <- .derivDensity(x = x, adjust = adjust, deriv = 2, ...)

    if (side == "right") {
      derivPeaks <- .findPeaks(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust
      )[, "x"]
      derivPeaks <- derivPeaks[derivPeaks > peaks[refPeak]]
      cutpoint <- sort(derivPeaks)[1]
    } else if (side == "left") {
      derivOut$y <- -derivOut$y
      derivPeaks <- .findPeaks(
        x = derivOut$x,
        y = derivOut$y,
        adjust = adjust
      )[, "x"]
      derivPeaks <- derivPeaks[derivPeaks < peaks[refPeak]]
      cutpoint <- sort(derivPeaks, decreasing = TRUE)[length(derivPeaks)]
    } else {
      stop("Unrecognized 'side' argument (was '", side, "').")
    }
  }

  cutpoint
}

#' @keywords internal
.derivDensity <- function(
  x,
  deriv = 1,
  bandwidth = NULL,
  adjust = 1,
  numPoints = 10000,
  ...
) {
  if (is.null(bandwidth)) {
    bandwidth <- ks::hpi(x, deriv.order = deriv)
  }
  kdeObj <- ks::kdde(
    x = x,
    deriv.order = deriv,
    h = bandwidth * adjust,
    gridsize = numPoints,
    ...
  )
  list(x = kdeObj$eval.points, y = kdeObj$estimate)
}
