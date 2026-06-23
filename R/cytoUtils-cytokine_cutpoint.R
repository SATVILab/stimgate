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
  num_peaks <- length(peaks)

  if (ref_peak > num_peaks) {
    msg <- paste("The reference peak is larger than the number of peaks found.",
                 "Setting the reference peak to 'num_peaks'...")
    if (strict) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
    ref_peak <- num_peaks
  }

  if (method == "first_deriv") {
    # Calculate the first derivative
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 1, ...)
    
    if (auto_tol) {
      tol <- 0.01 * max(abs(deriv_out$y))
    }
    
    if (side == "right") {
      deriv_valleys <- .find_valleys(x = deriv_out$x, y = deriv_out$y, adjust = adjust)
      deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys)[1]
      
      # Safe NA Check
      if (is.na(deriv_valleys)) {
        cutpoint <- NA_real_
      } else {
        cutpoint_candidates <- deriv_out$x[deriv_out$x > deriv_valleys & abs(deriv_out$y) < tol]
        cutpoint <- if (length(cutpoint_candidates) > 0) cutpoint_candidates[1] else NA_real_
      }
      
    } else if (side == "left") {
      deriv_out$y <- -deriv_out$y
      deriv_valleys <- .find_valleys(x = deriv_out$x, y = deriv_out$y, adjust = adjust)
      deriv_valleys <- deriv_valleys[deriv_valleys < peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys, decreasing = TRUE)[1]
      
      # Safe NA Check
      if (is.na(deriv_valleys)) {
        cutpoint <- NA_real_
      } else {
        cutpoint_candidates <- deriv_out$x[deriv_out$x < deriv_valleys & abs(deriv_out$y) < tol]
        cutpoint <- if (length(cutpoint_candidates) > 0) cutpoint_candidates[length(cutpoint_candidates)] else NA_real_
      }
    } else {
      stop("Unrecognized 'side' argument (was '", side, "').")
    }

  } else {
    # second_deriv method
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 2, ...)
    
    if (side == "right") {
      deriv_peaks <- .find_peaks(x = deriv_out$x, y = deriv_out$y, adjust = adjust)[, "x"]
      deriv_peaks <- deriv_peaks[deriv_peaks > peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks)[1]
      
    } else if (side == "left") {
      deriv_out$y <- -deriv_out$y
      deriv_peaks <- .find_peaks(x = deriv_out$x, y = deriv_out$y, adjust = adjust)[, "x"]
      deriv_peaks <- deriv_peaks[deriv_peaks < peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks, decreasing = TRUE)[length(deriv_peaks)]
      
    } else {
      stop("Unrecognized 'side' argument (was '", side, "').")
    }
  }

  cutpoint
}

#' @keywords internal
.deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
                           num_points = 10000, ...) {
  if (is.null(bandwidth)) {
    bandwidth <- ks::hpi(x, deriv.order = deriv)
  }
  kde_obj <- ks::kdde(
    x = x,  deriv.order = deriv,
    h = bandwidth * adjust, gridsize = num_points,
    ...
  )
  list(x = kde_obj$eval.points, y = kde_obj$estimate)
}

