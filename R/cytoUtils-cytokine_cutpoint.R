# “Adapted from cytoUtils by Mike Jiang, GPL-2-only”
#' Constructs a cutpoint for a flowFrame by using a derivative of the kernel
#' density estimate
#'
#' We determine a gating cutpoint using either the first or second derivative of
#' the kernel density estimate (KDE) of the \code{x}.
#'
#' By default, we compute the first derivative of the kernel density estimate.
#' Next, we determine the lowest valley from the derivative, which corresponds
#' to the density's mode for cytokines. We then contruct a gating cutpoint as
#' the value less than the tolerance value \code{tol} in magnitude and is also
#' greater than the lowest valley.
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
#' @param num_peaks the number of peaks expected to see. This effectively
#' removes any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained. By
#' default, the peak farthest to the left is used.
#' @param strict \code{logical} when the actual number of peaks detected is
#' less than \code{ref_peak}, an error is reported by default. But if
#' \code{strict} is set to FALSE, then the reference peak will be reset to
#' the peak of the far right.
#' @param tol the tolerance value
#' @param auto_tol when TRUE, it tries to set the tolerance automatically.
#' @param adjust the scaling adjustment applied to the bandwidth used in the
#' first derivative of the kernel density estimate
#' @param side character specifying the side of the gate, either
#'  \code{'right'} (default) or \code{'left'}
#' @param plot logical specifying whether to plot the peaks found
#' @param ... additional arguments passed to \code{.deriv_density}
#' @return the cutpoint along the x-axis
#' @keywords internal
.cytokine_cutpoint <- function(x,
                               adjust = 1,
                               num_peaks = 1,
                               ref_peak = 1,
                               method = c("first_deriv", "second_deriv"),
                               tol = 1e-2,
                               side = "right",
                               strict = TRUE,
                               plot = FALSE,
                               auto_tol = FALSE,
                               ...) {
  method <- match.arg(method)
  peaks <- sort(.find_peaks(x, num_peaks = num_peaks, adjust = adjust)[, "x"])

  # update peak count since it can be less than num_peaks
  num_peaks <- length(peaks)

  if (ref_peak > num_peaks) {
    out_func <- ifelse(strict, stop, warning)
    out_func("The reference peak is larger than the number of peaks found.",
      "Setting the reference peak to 'num_peaks'...",
      call. = FALSE
    )
    ref_peak <- num_peaks
  }

  # Finds the deepest valleys from the kernel density and sorts them.
  # The number of valleys identified is determined by 'num_peaks'
  deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 1, ...)
  if (auto_tol) {
    # Try to set the tolerance automatigically.
    tol <- 0.01 * max(abs(deriv_out$y))
  }
  if (side == "right") {
    deriv_valleys <- with(deriv_out, .find_valleys(x = x, adjust = adjust))
    deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
    deriv_valleys <- sort(deriv_valleys)[1]
    cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
    cutpoint <- cutpoint[1]
  } else if (side == "left") {
    deriv_out$y <- -deriv_out$y
    deriv_valleys <- with(deriv_out, .find_valleys(x = x, adjust = adjust))
    deriv_valleys <- deriv_valleys[deriv_valleys < peaks[ref_peak]]
    deriv_valleys <- sort(deriv_valleys, decreasing = TRUE)[1]
    cutpoint <- with(deriv_out, x[x < deriv_valleys & abs(y) < tol])
    cutpoint <- cutpoint[length(cutpoint)]
  } else {
    stop("Unrecognized 'side' argument (was '", side, "'.")
  }

  cutpoint
}

#' @keywords internal
.deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
                           num_points = 10000, ...) {
  # 1. Bandwidth Selection
  # ks::hpi is the standard plugin selector (likely what you were using before)
  if (is.null(bandwidth)) {
    bandwidth <- ks::hpi(x, deriv.order = deriv)
  }

  # 2. Compute Derivative using ks::kdde
  # We MUST pass 'gridsize' to maintain the high resolution
  # your cutpoint logic needs
  kde_obj <- ks::kdde(
    x = x, deriv.order = deriv,
    h = bandwidth * adjust, gridsize = num_points,
    ...
  )

  # 3. Format Output
  # ks::kdde returns 'eval.points' (x) and 'estimate' (y)
  # We map them to 'x' and 'y' so cytokine_cutpoint doesn't break
  list(x = kde_obj$eval.points, y = kde_obj$estimate)
}
