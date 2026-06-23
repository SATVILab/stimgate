#' @internal
.find_peaks <- function(x, y = NULL, num_peaks = NULL, adjust = 2, plot = FALSE, ...) {
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
  second_deriv <- diff(sign(diff(dens$y)))
  which_maxima <- which(second_deriv == -2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield peaks outside this range.  We remove
  # any such peaks.
  which_maxima <- which_maxima[findInterval(dens$x[which_maxima], range(x)) == 1]

  # Next, we sort the peaks in descending order based on the density heights.
  which_maxima <- which_maxima[order(dens$y[which_maxima], decreasing = TRUE)]
  
  # Safely construct the return dataframe
  if (length(which_maxima) > 0) {
    peaks_x <- dens$x[which_maxima]
    if (is.null(num_peaks) || num_peaks > length(peaks_x)) {
      num_peaks <- length(peaks_x)
    }
    peaks <- data.frame(
      x = peaks_x[seq_len(num_peaks)], 
      y = dens$y[which_maxima][seq_len(num_peaks)]
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

#' @internal
.find_valleys <- function(x, y = NULL, num_valleys = NULL, adjust = 2, ...) {

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
  second_deriv <- diff(sign(diff(dens$y)))
  which_minima <- which(second_deriv == 2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield valleys outside this range. We remove
  # any such valleys.
  which_minima <- which_minima[findInterval(dens$x[which_minima], range(x)) == 1]

  # Next, we sort the valleys in descending order based on the density heights.
  which_minima <- which_minima[order(dens$y[which_minima], decreasing = FALSE)]

  # Returns the local minima. If there are none, we return 'NA' instead.
  if (length(which_minima) > 0) {
    valleys <- dens$x[which_minima]
    if (is.null(num_valleys) || num_valleys > length(valleys)) {
      num_valleys <- length(valleys)
    }
    valleys <- valleys[seq_len(num_valleys)]
  } else {
    valleys <- NA
  }
  valleys
}
