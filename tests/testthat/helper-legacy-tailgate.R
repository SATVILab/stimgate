# Fallback helper definitions for legacy comparator functions used in tests
# when scripts/r/ is not available (e.g., in installed package R CMD check).

.findPeaksLegacyFallback <- function(
  x,
  y = NULL,
  numPeaks = NULL,
  adjust = 2,
  plot = FALSE,
  ...
) {
  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find peaks.")
    return(NA)
  }

  if (is.null(y)) {
    dens <- stats::density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  secondDeriv <- diff(sign(diff(dens$y)))
  whichMaxima <- which(secondDeriv == -2) + 1

  whichMaxima <- whichMaxima[
    findInterval(dens$x[whichMaxima], range(x)) == 1
  ]

  whichMaxima <- whichMaxima[order(dens$y[whichMaxima], decreasing = TRUE)]

  if (length(whichMaxima) > 0) {
    peaksX <- dens$x[whichMaxima]
    if (is.null(numPeaks) || numPeaks > length(peaksX)) {
      numPeaks <- length(peaksX)
    }
    peaks <- data.frame(
      x = peaksX[seq_len(numPeaks)],
      y = dens$y[whichMaxima][seq_len(numPeaks)]
    )
  } else {
    peaks <- data.frame(x = NA_real_, y = NA_real_)
  }

  if (plot) {
    graphics::plot(dens, main = paste("adjust =", adjust))
    graphics::points(peaks$x, peaks$y, col = "red")
  }

  peaks
}

.findValleysLegacyFallback <- function(
  x,
  y = NULL,
  numValleys = NULL,
  adjust = 2,
  ...
) {
  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find valleys.")
    return(NA)
  }

  if (is.null(y)) {
    dens <- stats::density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  secondDeriv <- diff(sign(diff(dens$y)))
  whichMinima <- which(secondDeriv == 2) + 1

  whichMinima <- whichMinima[
    findInterval(dens$x[whichMinima], range(x)) == 1
  ]

  whichMinima <- whichMinima[order(dens$y[whichMinima], decreasing = FALSE)]

  if (length(whichMinima) > 0) {
    valleys <- dens$x[whichMinima]
    if (is.null(numValleys) || numValleys > length(valleys)) {
      numValleys <- length(valleys)
    }
    valleys <- valleys[seq_len(numValleys)]
  } else {
    valleys <- NA
  }
  valleys
}

.derivDensityLegacyFallback <- function(
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

.cytokineCutpointLegacyFallback <- function(
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

.loadLegacyTailgateFallback <- function(env = parent.frame()) {
  assign(".findPeaks", .findPeaksLegacyFallback, envir = env)
  assign(".findValleys", .findValleysLegacyFallback, envir = env)
  assign(".derivDensity", .derivDensityLegacyFallback, envir = env)
  assign(".cytokineCutpoint", .cytokineCutpointLegacyFallback, envir = env)
}
