
# “Adapted from cytoUtils by Mike Jiang, GPL-2-only”
#' Constructs a cutpoint for a flowFrame by using a derivative of the kernel
#' density estimate
#'
#' We determine a gating cutpoint using either the first or second derivative of
#' the kernel density estimate (KDE) of the \code{x}.
#'
#' By default, we compute the first derivative of the kernel density estimate. 
#' Next, we determine the lowest valley from the derivative, which corresponds to the
#' density's mode for cytokines. We then contruct a gating cutpoint as the value
#' less than the tolerance value \code{tol} in magnitude and is also greater
#' than the lowest valley.
#'
#' Alternatively, if the \code{method} is selected as \code{second_deriv}, we
#' select a cutpoint from the second derivative of the KDE. Specifically, we
#' choose the cutpoint as the largest peak of the second derivative of the KDE
#' density which is greater than the reference peak.
#'
#' @rdname gate_tail
#' @param x a \code{numeric} vector used as input data
#' @param num_peaks the number of peaks expected to see. This effectively removes
#' any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained. By
#' default, the peak farthest to the left is used.
#' @param strict \code{logical} when the actual number of peaks detected is less than \code{ref_peak}. 
#'                               an error is reported by default. But if \code{strict} is set to FALSE, then the reference peak will be reset to the peak of the far right.      
#' @param method the method used to select the cutpoint. See details.
#' @param tol the tolerance value
#' @param auto_tol when TRUE, it tries to set the tolerance automatically.
#' @param adjust the scaling adjustment applied to the bandwidth used in the
#' first derivative of the kernel density estimate
#' @param side character specifying the side of the gate, either
#'  \code{'right'} (default) or \code{'left'}
#' @param plot logical specifying whether to plot the peaks found
#' @param ... additional arguments passed to \code{.deriv_density}
#' @return the cutpoint along the x-axis
.cytokine_cutpoint <- function(x, num_peaks = 1, ref_peak = 1,
    method = c("first_deriv", "second_deriv"),
    tol = 1e-2, adjust = 1, side = "right", strict = TRUE, plot = FALSE, auto_tol = FALSE, ...) {
  
  method <- match.arg(method)
  peaks <- sort(.find_peaks(x, num_peaks = num_peaks, adjust = adjust, plot = plot)[, "x"])
  
  #update peak count since it can be less than num_peaks
  num_peaks <- length(peaks)
  
  if (ref_peak > num_peaks) {
    outFunc <- ifelse(strict, stop, warning)
    outFunc("The reference peak is larger than the number of peaks found.",
        "Setting the reference peak to 'num_peaks'...",
        call. = FALSE)
    ref_peak <- num_peaks
  }
  
  # TODO: Double-check that a cutpoint minimum found via 'first_deriv'
  # passes the second-derivative test.
  
  if (method == "first_deriv") {
    # Finds the deepest valleys from the kernel density and sorts them.
    # The number of valleys identified is determined by 'num_peaks'
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 1, ...)
    if(auto_tol){
      #Try to set the tolerance automatigically.
      tol <- 0.01*max(abs(deriv_out$y))
    }
    if (side == "right") {
      
      deriv_valleys <- with(deriv_out, .find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys)[1]
      cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[1]
      
    } else if (side == "left") {
      
      deriv_out$y <- -deriv_out$y
      deriv_valleys <- with(deriv_out, .find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys < peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys, decreasing=TRUE)[1]
      cutpoint <- with(deriv_out, x[x < deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[ length(cutpoint) ]
      
    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }
    
  } else {
    # The cutpoint is selected as the first peak from the second derivative
    # density which is to the right of the reference peak.
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 2, ...)
    
    if (side == "right") {
      deriv_peaks <- with(deriv_out, .find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks > peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks)[1]
    } else if (side == "left") {
      deriv_out$y <- -deriv_out$y
      deriv_peaks <- with(deriv_out, .find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks < peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks, decreasing=TRUE)[length(deriv_peaks)]
    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }
    
  }
  
  cutpoint
}

.deriv_density <- function(x, deriv = 1L, adjust = 1, n = 2048L, ...) {
  # 1) Compute a fine KDE on [min(x), max(x)]
  dens <- stats::density(x, adjust = adjust, n = n, ...)
  xx  <- dens$x
  yy  <- dens$y

  # 2) Approximate derivatives by repeated central differences
  for (k in seq_len(deriv)) {
    # central difference: f'(x_i) ≈ (y_{i+1} - y_{i-1}) / (x_{i+1} - x_{i-1})
    dx <- diff(xx, lag = 2)
    dy <- yy[-1] - yy[-length(yy)]
    # but dy is length n−1; for central we align:
    yy <- dy[-1] / dx
    xx <- xx[-c(1, length(xx))]
  }

  list(x = xx, y = yy)
}
