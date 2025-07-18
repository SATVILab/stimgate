#' Find local maxima (“peaks”) in a numeric vector via KDE
#'
#' @param x Numeric vector.
#' @param num_peaks Integer: how many of the highest peaks to return.
#' @param adjust Bandwidth multiplier, as in stats::density().
#' @return A data.frame with columns x, y (peak locations and heights).
#' @noRd
.find_peaks <- function(x, num_peaks = NULL, adjust = 2, ...) {
  x <- as.numeric(x)
  if (length(x) < 2) return(data.frame(x = NA_real_, y = NA_real_))

  dens <- stats::density(x, adjust = adjust, ...)
  # second discrete derivative:
  sd2 <- diff(sign(diff(dens$y)))
  idx <- which(sd2 == -2) + 1L
  idx <- idx[findInterval(dens$x[idx], range(x)) == 1L]
  heights <- dens$y[idx]
  if (length(idx)==0) return(data.frame(x=NA_real_, y=NA_real_))

  ord <- order(heights, decreasing = TRUE)
  if (!is.null(num_peaks)) ord <- head(ord, num_peaks)
  idx <- idx[ord]

  data.frame(
    x = dens$x[idx],
    y = dens$y[idx],
    row.names = NULL
  )
}

#' Find local minima (“valleys”) in a numeric vector via KDE
#'
#' @param x Numeric vector.
#' @param num_valleys Integer: how many of the lowest valleys to return.
#' @param adjust Bandwidth multiplier, as in stats::density().
#' @return Numeric vector of valley locations.
#' @noRd
.find_valleys <- function(x, num_valleys = NULL, adjust = 2, ...) {
  x <- as.numeric(x)
  if (length(x) < 2) return(NA_real_)

  dens <- stats::density(x, adjust = adjust, ...)
  sd2 <- diff(sign(diff(dens$y)))
  idx <- which(sd2 ==  2) + 1L
  idx <- idx[findInterval(dens$x[idx], range(x)) == 1L]
  vals  <- dens$x[idx]
  if (length(vals)==0) return(NA_real_)

  ord <- order(dens$y[idx], decreasing = FALSE)
  if (!is.null(num_valleys)) ord <- head(ord, num_valleys)
  vals[ord]
}
