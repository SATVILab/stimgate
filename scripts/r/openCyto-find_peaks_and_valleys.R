#' @keywords internal
.findPeaks <- function(
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
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  secondDeriv <- diff(sign(diff(dens$y)))
  whichMaxima <- which(secondDeriv == -2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield peaks outside this range.  We remove
  # any such peaks.
  whichMaxima <- whichMaxima[
    findInterval(dens$x[whichMaxima], range(x)) == 1
  ]

  # Next, we sort the peaks in descending order based on the density heights.
  whichMaxima <- whichMaxima[order(dens$y[whichMaxima], decreasing = TRUE)]

  # Safely construct the return dataframe
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
    plot(dens, main = paste("adjust =", adjust))
    points(peaks$x, peaks$y, col = "red") # Fixed plotting syntax
  }

  peaks
}

#' @keywords internal
.findValleys <- function(x, y = NULL, numValleys = NULL, adjust = 2, ...) {
  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find valleys.")
    return(NA)
  }

  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  secondDeriv <- diff(sign(diff(dens$y)))
  whichMinima <- which(secondDeriv == 2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield valleys outside this range. We remove
  # any such valleys.
  whichMinima <- whichMinima[
    findInterval(dens$x[whichMinima], range(x)) == 1
  ]

  # Next, we sort the valleys in descending order based on the density heights.
  whichMinima <- whichMinima[order(dens$y[whichMinima], decreasing = FALSE)]

  # Returns the local minima. If there are none, we return 'NA' instead.
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
