#' @title Get log-likelihood at specified cutpoints
#'
#' @description
#' Calculates the log-likelihood of models where the cutpoint is the "break" point in a
#' piece-wise logistic regression model.
#'
#' @param high_ind_tbl tibble. Tibble containing a column \code{cut} with expression values for the marker
#' to be cut and a numeric columm \code{high} that has value \code{1} when the cell was high for at least one of the specified markers and 0 otherwise.
#' @param cps numeric vector. Points to treat as cutpoints when evaluating the log-likelihood.
#'
#' @return Numeric vector of log-likelihoods for each cutpoint in \code{cps}.
#'
#' @details
#' This function fits a model that is flat until the cutpoint and
#' can then jump and take on a linear term at the cutpoint.
.get_ll_scp <- function(high_ind_tbl, cps) {
  # get log-likelihoods
  vapply(cps, function(cp) {
    high_ind_tbl <- high_ind_tbl |>
      dplyr::mutate(
        diff_from_cut = pmax(cut - cp, 0),
        ind_x_high = as.numeric(cut > cp)
      )
    fit_alt <- glm(high ~ cut + ind_x_high + diff_from_cut,
      family = binomial, .data = high_ind_tbl
    )
    (AIC(fit_alt) - 4 * 2) / -2
  }, numeric(1))
}


#' @title Get cutpoint (helper function to main function)
#'
#' @description
#' Get the cutpoint for a given channel, based on being "high" for other channels. The
#' cutpoint is chosen as the cutpoint such that a piecewise function using it as the "break"
#' point in a logistic regression model has the highest possible log-likelihood.
#'
#' @param ll numeric vector. Numeric vector of log-likelihoods, where
#' the i-th element is the log-likelihood where the cutpoint is the
#' i-th element in \code{cps}.
#' @param cps numeric vector. Points that were treated as the cutpoints when evaluating
#' the log-likelihood in \code{ll}.
#'
#' @details
#' If two possible cutpoints have the same log-likelihood, then the
#' larger cutpoint is returned.
#'
#' @return A single number specifying the cutpoint. Note that if
.get_cp_scp <- function(ll, cps) {
  cp <- cps[which(ll == max(ll))]
  cp[length(cp)]
}

.test_for_changepoint <- function(high_ind_tbl, ll) {
  fit_h0 <- glm(high ~ 1, family = binomial, .data = high_ind_tbl)
  ll0 <- (AIC(fit_h0) - 1 * 2) / -2

  # get test stat

  R_vec <- 2 * (ll - ll0)
  G <- max(R_vec)

  # work out critical values
  n <- ceiling(floor(high_ind_tbl[["cut"]])) - ceiling(min(high_ind_tbl[["cut"]])) + 1
  n <- nrow(high_ind_tbl)
  critical_point_vec <- c(
    "bic" = 2 * log(n),
    "aic" = 2 * 2,
    "hannan-quinn" = 2 * 2 * log(log(n))
  )
  sig_vec <- (G > critical_point_vec)
  list(
    G = G, R = R_vec, critical_point = critical_point_vec,
    sig = sig_vec
  )
}

.get_ll_dcp <- function(ex, cut, cp_scp, max = NULL, n_break = 1) {
  ex_vec <- ex[[cut]][ex[[cut]] >= cp_scp]
  if (length(ex_vec) <= 10) {
    return(stats::setNames(NA, cp_scp))
  }
  hist_obj <- hist(ex_vec, breaks = 50, plot = FALSE) # changed it to 30 breaks now
  count_vec <- hist_obj$counts
  break_vec <- hist_obj$breaks
  midpt_vec <- purrr::map_dbl(seq_along(break_vec)[-1], function(i) (break_vec[i - 1] + break_vec[i]) / 2)
  init_mod_tbl <- tibble::tibble(ex = midpt_vec, count = count_vec)
  # cut_vec <- seq(ceiling(min(ex_vec)), floor(max(ex_vec)))
  cut_vec <- break_vec[-length(break_vec)]
  ll_vec <- purrr::map_dbl(cut_vec, function(cut) {
    if (n_break == 1) {
      final_mod_tbl <- init_mod_tbl |>
        dplyr::mutate(
          low_ind = as.numeric(ex < cut),
          diff_from_cut = abs(ex - cut),
          high_ind = as.numeric(ex >= cut),
          high_diff_from_cut = high_ind * diff_from_cut,
          low_diff_from_cut = low_ind * diff_from_cut
        )
      # fit <- lm(count ~ high_ind + high_diff_from_cut + low_diff_from_cut, .data = final_mod_tbl)
      fit <- glm(count ~ high_diff_from_cut + low_diff_from_cut, .data = final_mod_tbl, family = "poisson")
    } else if (n_break == 2) {
      final_mod_tbl <- init_mod_tbl |>
        dplyr::mutate(
          low_ind = as.numeric(ex < cut),
          diff_from_cut = abs(ex - cut),
          high_ind = as.numeric(ex >= cut),
          high_diff_from_cut = high_ind * diff_from_cut,
          low_diff_from_cut = low_ind * diff_from_cut
        )
      knots <- c((min(final_mod_tbl[["ex"]]) + cut) / 2, cut)
      fit <- glm(count ~ splines::bs(ex, knots = knots, degree = 1), .data = final_mod_tbl, family = "poisson")
    }

    # final_mod_tbl$count - predict(fit))^2)
    (AIC(fit) - 4 * 2) / -2
  })
  max <- ifelse(is.null(max), max(cut_vec), max)
  stats::setNames(ll_vec, as.character(cut_vec))[cut_vec <= max]
}

.get_cp_dcp <- function(ll, shift_prop = NULL, shift_abs = NULL) {
  if (is.na(ll[1]) & length(ll) == 1) {
    return(as.numeric(names(ll)))
  }
  if (FALSE) {
    # old - from when minimising sums of squares
    if (is.na(ll)) {
      return(as.numeric(names(ll)))
    }
    if (min(ll) == ll[length(ll)]) {
      return(as.numeric(names(ll)[length(ll)]))
    }
    min_min_ind <- min(which(ll == min(ll)))
    min_min_pt <- as.numeric(names(ll)[min_min_ind]) # minimum point at which minimum is reached
    min_min_ll <- ll[[min_min_ind]]
    max_min_ind <- max(which(ll == min(ll)))
    max_min_pt <- as.numeric(names(ll)[max_min_ind]) # minimum point at which minimum is reached
    max_min_ll <- ll[[max_min_ind]]
    if (is.null(shift_prop) & is.null(shift_abs)) {
      return(max_min_pt)
    }
    max_after_max_min_pt_ind <- min(which(ll == max(ll[max_min_ind:length(ll)])))
    max_after_max_min_pt_ll <- ll[max_after_max_min_pt_ind]
    within_ind <- ll > max_min_ll & (ll < (max_min_ll + (max_after_max_min_pt_ll - max_min_ll) * 0.4))
    max(as.numeric(names(ll))[within_ind])
  }
  max_max_ind <- max(which(ll == max(ll)))
  max_max_pt <- as.numeric(names(ll)[max_max_ind]) # minimum point at which minimum is reached
  max_max_ll <- ll[[max_max_ind]]
  max_max_pt
}

