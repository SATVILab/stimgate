#' @title Check that parameters for each marker for which a gate is required are complete
#'
#' @inheritParams gate # marker
#' @data_name 'gs_cytof' or `gs_proto`. Name of dataset to be gated in R environment.
#'
#' @return
.complete_marker_list <- function(marker,
                                  data_name,
                                  bias_uns,
                                  .data,
                                  pop_gate,
                                  cut,
                                  debug,
                                  bw_min,
                                  cp_min,
                                  ind_batch_list,
                                  ind_in_batch_gate,
                                  ind_in_batch_uns,
                                  ind_in_batch_lab_vec) {
  if (!str_detect_any(data_name, c( # nolint
    "gs_cytof", "gs_proto", "gs_cd8_base",
    "gs_cytof_acs"
  ))) {
    stop("data_name not recognised.")
  }
  purrr::map(marker, function(marker_curr) {
    # fill in optionally-specified parameters
    if ((!"fdr" %in% names(marker_curr)) ||
      is.null(marker_curr$fdr)) {
      marker_curr$fdr <- NULL
    }

    if (!"min_cell" %in% names(marker_curr)) {
      marker_curr$min_cell <- 1e2
    }

    if (!is.null(marker_curr[["high"]])) {
      marker_curr[["high"]] <- marker_curr[["high"]][
        names(marker_curr[["high"]]) != marker_curr[["cut"]]
      ]
    }

    if ((!"tol" %in% names(marker_curr)) ||
      is.null(marker_curr$tol)) {
      if (str_detect_any(data_name, c( # nolint
        "gs_cytof", "gs_proto", "gs_cd8_base",
        "gs_cytof_acs"
      ))) {
        marker_curr$tol <- 0.5e-8
      }
    }

    if ((!"gate_combn" %in% names(marker_curr)) |
      is.null(marker_curr$gate_combn)) {
      marker_curr$gate_combn <- NULL
    }
    marker_curr$gate_combn <- .get_gate_combn_list(
      gate_combn = marker_curr$gate_combn,
      fdr = marker_curr$fdr
    )
    if (!"pop_man_match_exact" %in% names(marker_curr)) {
      marker_curr$pop_man_match_exact <- TRUE
    }

    marker_curr$bias_uns <- .complete_marker_list_bias_uns(
      bias_uns = marker_curr$bias_uns,
      .data = .data,
      data_name = data_name,
      pop_gate = pop_gate,
      cut = marker_curr$cut,
      debug = debug,
      ind_batch_list = ind_batch_list,
      ind_in_batch_gate = ind_in_batch_gate,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec
    )

    marker_curr$bw_min <- .complete_marker_list_min_bw(
      bw_min = bw_min,
      bias_uns = max(marker_curr$bias_uns)
    )

    marker_curr$cp_min <- .complete_marker_list_cp_min(
      cp_min = cp_min,
      .data = .data,
      data_name = data_name,
      pop_gate = pop_gate,
      cut = marker_curr$cut,
      debug = debug,
      ind_batch_list = ind_batch_list,
      ind_in_batch_gate = ind_in_batch_gate,
      ind_in_batch_uns = ind_in_batch_uns,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec
    )

    marker_curr
  })
}

.complete_marker_list_bias_uns <- function(bias_uns,
                                           .data,
                                           data_name,
                                           pop_gate,
                                           cut,
                                           debug,
                                           ind_batch_list,
                                           ind_in_batch_gate,
                                           ind_in_batch_uns,
                                           ind_in_batch_lab_vec) {
  if (!is.null(bias_uns)) {
    return(bias_uns)
  }
  .debug(debug, "calculating bias_uns automatically") # nolint
  mean_range <- .complete_marker_list_bias_uns_get_mean_range(
    ind_batch_list = ind_batch_list,
    .data = .data,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    cut = cut,
    data_name = data_name
  )
  (mean_range / 40) |> signif(3)
}

.complete_marker_list_bias_uns_get_mean_range <- function(ind_batch_list,
                                                          .data,
                                                          ind_in_batch_gate,
                                                          ind_in_batch_uns,
                                                          ind_in_batch_lab_vec,
                                                          pop_gate,
                                                          cut,
                                                          data_name) {
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        data = .data, ind_batch = ind_batch_list[[i]],
        ind_in_batch_gate = ind_in_batch_gate,
        ind_in_batch_uns = ind_in_batch_uns,
        ind_in_batch_lab_vec = ind_in_batch_lab_vec,
        pop = pop_gate,
        cut = cut, high = NULL,
        data_name = data_name
      )
      purrr::map_dbl(ex_list, function(ex) {
        abs(diff(
          quantile(
            ex$cut[ex$cut > min(ex$cut)], c(0.99, 0.01)
          ),
          na.rm = TRUE
        ))[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


.complete_marker_list_min_bw <- function(bw_min, bias_uns) {
  if (!is.null(bw_min)) {
    return(bw_min)
  }
  bias_uns * 1.5
}

.complete_marker_list_cp_min <- function(cp_min,
                                         .data,
                                         data_name,
                                         pop_gate,
                                         cut,
                                         debug,
                                         ind_batch_list,
                                         ind_in_batch_gate,
                                         ind_in_batch_uns,
                                         ind_in_batch_lab_vec) {
  if (!is.null(cp_min)) {
    return(cp_min)
  }
  .debug(debug, "calculating cp_min automatically") # nolint
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        data = .data, ind_batch = ind_batch_list[[i]],
        ind_in_batch_gate = ind_in_batch_gate,
        ind_in_batch_uns = ind_in_batch_uns,
        ind_in_batch_lab_vec = ind_in_batch_lab_vec,
        pop = pop_gate,
        cut = cut, high = NULL,
        data_name = data_name
      )
      purrr::map_dbl(ex_list, function(ex) {
        median(ex$cut[ex$cut > min(ex$cut)], na.rm = TRUE)[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}


#' @title Get named vector specifying gate combination method for each cutpoint type
#'
#' @inheritParams get_cp # gate_batch_combn, fdr
# get gate_batch_combn vec
#'
#' @return Named character vector, where all elements together represent names of all cutpoints.
#' Each element name is name of a cutpoint, and corresponding value is name of gate combination
#' method for that cutpoint.
.get_gate_combn_list <- function(gate_combn, fdr) {
  # get all possible cp types
  cp_type_vec_full <- gate_combn |>
    unlist() |>
    unique()
  # cp_type_vec_full <- .get_full_cp_type_vec(fdr = fdr)

  # get gating method for each type of gate
  purrr::map(cp_type_vec_full, function(x) {
    gate_combn_vec <- c()
    for (i in seq_along(gate_combn)) {
      if (x %in% gate_combn[[i]]) {
        gate_combn_vec <- c(gate_combn_vec, names(gate_combn)[i])
      }
    }
    if (!length(gate_combn_vec)) {
      return(stats::setNames(list("no"), x))
    }
    list(gate_combn_vec) |> stats::setNames(x)
  }) |>
    flatten()
}

#' @title Get all cp type names
#'
#' @inheritParams get_cp #fdr
#'
#' @return Character vector, where each element is name of a
#' cutpoint, and all elements together represent names of all cutpoints.
.get_full_cp_type_vec <- function(fdr) {
  # Get cutpoint names for unstim-based cuts
  # cp_name_vec_uns <- .get_cp_uns_name_vec(fdr)
  # cp_name_vec_uns_root <- paste0(cp_name_vec_uns, 'root')

  # output all cutpoint names
  c(
    "man", "tg", "dcp", "midp", "scp",
    "uns", "unsr", "loc"
  )
}

#' @title Get name vec for unstim-based cuts
#'
#' @inheritParams get_cp # fdr
#'
#' @return Character vector.
.get_cp_uns_fdr_name_vec <- function(fdr) {
  # get cutpoint names
  cp_name_vec <- as.character(fdr * 100)
  len_max <- max(stringr::str_length(cp_name_vec))
  purrr::map_chr(cp_name_vec, function(cp_name) {
    len_curr <- stringr::str_length(cp_name)
    if (len_curr == len_max) {
      return(paste0("uns", cp_name))
    }
    cp_name <- paste0(rep("0", len_max - len_curr), cp_name)
    paste0(cp_name)
  })
}



#' @title Get stats for manual cutpoint
#'
#' @description
#' Get implied manual cutpoint for each sample with index in \code{ind_gate}
#' within original GatingSet
#'
#' @inheritParams get_cp data, cut, pop_man
#' @param ex_list list. List where each element is a dataframe corresponding to an expression matrix.
#' @param ind numeric vector. Each element represents the index of a sample within \code{data}
#' for which a gate is required.
#'
#' @return A named character vector with elements named \code{"count"},
#' \code{"freq"} and \code{"cp"} containg the count, frequency and implied cutpoint
# for the manually gated positive population.
.get_cp_man <- function(data, ind,
                        ex_list,
                        pop_man,
                        pop_man_match_exact,
                        gate_combn) {
  # get counts for each sample
  count_man_vec <- purrr::map_dbl(ind, function(ind_curr) {
    gh <- data[[ind_curr]]

    # get all data for specified node
    man_stats_tbl <- gh_pop_get_stats(gh, xml = FALSE) # was TRUE
    if (pop_man_match_exact) man_stats_tbl <- man_stats_tbl |> dplyr::filter(pop %in% pop_man)
    if (!pop_man_match_exact) {
      for (x in pop_man) {
        man_stats_tbl <- man_stats_tbl |> dplyr::filter(stringr::str_detect(pop, x))
      }
    }

    # return count
    man_stats_tbl[["count"]] |> sum()
  }) |>
    stats::setNames(ind)

  # cutpoints to be calculated as a group
  cp_list <- list()


  if ("prejoin" %in% gate_combn) {
    cp_list[["man_prejoin"]] <- .get_cp_man_prejoin(
      count = count_man_vec,
      ex_list = ex_list,
      ind = ind
    )
  }

  non_prejoin_combn_vec <- setdiff(gate_combn, "prejoin")

  if (length(non_prejoin_combn_vec) > 0) {
    cp_list <- cp_list |>
      append(.get_cp_man_non_prejoin(
        count = count_man_vec,
        ex_list = ex_list,
        ind = ind,
        gate_combn = non_prejoin_combn_vec
      ) |>
        purrr::flatten())
  }

  cp_list
}

#' @title Get the manual cutpoint based on pre-joining the samples
#'
#' @param count numeric vector. Vector of manually-gated counts across samples.
#' @param ind  numeric vector. Vector of indices of samples in original GatingSet.
#'
#' @return A named numeric vector, where the names are the sample indices in the
#' original GatingSet and the values are the corresponding manual cutpoints.
.get_cp_man_non_prejoin <- function(count, ex_list, ind, gate_combn) {
  purrr::map(gate_combn, function(gate_combn_curr) {
    cp_man_vec_init <- purrr::map_dbl(ind, function(ind_curr) {
      # get current count
      count_man <- count[as.character(ind_curr)]
      # get current expression tibble
      ex <- ex_list[[as.character(ind_curr)]] |> dplyr::arrange(desc(cut))
      # get the range to search over
      cps <- c(ex$cut[nrow(ex)], ex$cut[1])
      # return max observed value + a fraction of range if count is zero
      if (count_man == 0) {
        return(cps[2] + (cps[2] - cps[1]) / 200)
      }

      # get cutpoint for this sample
      (ex$cut[count_man] + ex$cut[count_man + 1]) / 2
    }) |>
      stats::setNames(ind)

    .combine_cp(cp = cp_man_vec_init, gate_combn = gate_combn_curr)
  }) |>
    stats::setNames(gate_combn)
}

#' @title Get the manual cutpoint based on pre-joining the samples
#'
#' @param count numeric vector. Vector of manually-gated counts across samples.
#' @param ind  numeric vector. Vector of indices of samples in original GatingSet.
#'
#' @return A named numeric vector, where the names are the sample indices in the
#' original GatingSet and the values are the corresponding manual cutpoints.
.get_cp_man_prejoin <- function(count, ex_list, ind) {
  # get overall count
  count_man <- sum(count)
  # get expression matrix of all the expression matrices joined together
  ex <- purrr::map(ind, function(ind_curr) ex_list[[as.character(ind_curr)]]) |>
    dplyr::bind_rows() |>
    dplyr::arrange(desc(cut))
  # get the range to search over
  cps <- c(ex$cut[nrow(ex)], ex$cut[1])
  # find the single cutpoint

  if (count_man == 0) {
    return(rep(cps[2] + (cps[2] - cps[1]) / 200, length(ind)) |>
      stats::setNames(ind))
  }
  cp_man <- (ex_stim$cut[count_pos] + ex_stim$cut[count_pos + 1]) / 2
  # repeat the cutpoints
  rep(cp_man, length(ind)) |>
    stats::setNames(ind)
}

#' @title COmbine gates across indices
#'
#' @param cp named numeric vector. Values are gates and names are indices within
#' original GatingSet to which gate applies.
#' @param gate_combn 'no', 'min', 'mean', 'trim20',  'median' or 'max'. Specifies
#' method of combining the gates in cp. 'no' and 'prejoin' do not apply common gates,
#' whereas each of the others do.
#'
#' @return named list. Each list corresponds to the gates for the set of indices
#' for a given gating method and gate grouping method. The name indicates
#' the gating method and the gate-grouping method, which are separated by an underscore.
#' The element is a named vector, where values are the gates for the individual samples and
#' the names are indices of the samples with the original GatingSet to which gates apply.
.combine_cp <- function(cp, gate_combn) {
  purrr::map(gate_combn, function(gate_combn_curr) {
    if (all(purrr::map_lgl(cp, is.na))) {
      return(stats::setNames(cp, names(cp)))
    }
    if (gate_combn_curr %in% c("no", "prejoin")) {
      return(cp)
    }
    if (gate_combn_curr == "min") {
      return(stats::setNames(
        rep(
          min(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gate_combn_curr == "mean") {
      return(stats::setNames(
        rep(
          mean(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gate_combn_curr == "trim20") {
      return(stats::setNames(
        rep(
          mean(cp, trim = 0.2, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gate_combn_curr == "median") {
      return(stats::setNames(
        rep(
          median(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gate_combn_curr == "max") {
      return(stats::setNames(
        rep(
          max(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
  }) |>
    stats::setNames(gate_combn)
}




#' @title Get the measurements for a given set of sample(s) for the marker to be cut on
#' as a single numeric vector
#'
#' @return A list where each element is a numeric vector.
.prepare_ex_list_with_bias_and_noise <- function(ex_list,
                                                 ind,
                                                 exc_min = TRUE,
                                                 bias = 0,
                                                 debug = FALSE,
                                                 noise_sd = NULL) {
  purrr::map(ind, function(ind_curr) {
    cut_tbl <- ex_list[[as.character(ind_curr)]] |>
      dplyr::mutate(
        sample = paste0(batch, "_", stim), # nolint
        expr = cut
      ) |>
      dplyr::select(sample, ind, ind_cell, expr) # nolint
    if (exc_min) {
      cut_tbl <- cut_tbl |> dplyr::filter(.data$expr > min(.data$expr)) # nolint
    }
    cut_tbl <- cut_tbl |> dplyr::mutate(expr = expr + bias) # nolint
    if (!is.null(noise_sd)) {
      cut_tbl <- purrr::map_df(1:5, function(i) {
        cut_tbl |>
          dplyr::mutate(expr = expr + rnorm(nrow(cut_tbl), sd = noise_sd)) # nolint
      })
    }
    cut_tbl
  }) |>
    stats::setNames(as.character(ind))
}



#' @title Replace zero p-values
#'
#' @description
#' Replaces zero p-values with equally spaced p-values from 0 (exclusive) to the
#' minimum positive p-value (exclusive).
#'
#' @param p numeric vector. Vector of p-values.
#'
#' @return Numeric vector.
.rep_zero_p_val <- function(p) {
  if (!sum(p)) {
    return("No p-values are non-zero. Will not work.")
  }
  p_zero_ind_vec <- which(p == 0)
  if (!length(p_zero_ind_vec)) {
    return(p)
  }
  min_non_zero_p <- min(p[p > 0])
  p_rep_vec <- seq(0, min_non_zero_p, length.out = length(p_zero_ind_vec) + 2)
  p[p_zero_ind_vec] <- p_rep_vec[-c(1, length(p_rep_vec))]
  p
}



#' @title Linearly interpolate
#'
#' @param val numeric. A single numeric value that represents the input
#' for which the value of the function needs to be interpolated at.
#' @param x,y numeric vector. Input and corresponding output values of function.
#'
#' @return Interpolated value of function at \code{val}.
.interp <- function(val, x, y) {
  x_low <- x[x <= val] |> max()
  if (x_low == val) {
    return(y[which(x == x_low)])
  }
  x_high <- x[x >= val] |> min()
  x_low_ind <- which(x == x_low)
  x_high_ind <- which(x == x_high)

  y[x_low_ind] + (val - x_low) * (y[x_high_ind] - y[x_low_ind]) / (x_high - x_low)
}

#' @title Get cytoUtils tailgate cutpoint
#'
#' @description
#' Returns the \code{cytoUtils::gate_tail} cutpoint.
#'
#'
#' @inheritParams .wins_ex
.get_cp_tg <- function(ex_list,
                       gate_combn,
                       cut,
                       ind_gate,
                       exc_min,
                       tol,
                       min_cell,
                       cp_min,
                       bw,
                       debug) {
  .debug(debug, "Getting tg cutpoint")
  cp_list <- list()

  if ("prejoin" %in% gate_combn) {
    .debug(debug, "prejoin")
    ex <- dplyr::bind_rows(ex_list)
    ex <- ex |> dplyr::filter(!is.na(cut))
    if (exc_min) ex <- ex |> dplyr::filter(cut > min(cut))
    if (nrow(ex) < max(min_cell, 5)) {
      cp_vec <- stats::setNames(rep(NA, length(ind_gate)), ind_gate)
    } else {
      dens <- density(ex[["cut"]])
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
      .ensure_cytoutils()
      cp <- cytoUtils:::.cytokine_cutpoint(
        x = ex[["cut"]], num_peaks = 1,
        ref_peak = 1, tol = tol, side = "right",
        strict = FALSE, adjust = adjust
      )
      # if(is.na(cp)) cp <- max(ex[['cut']]) + 0.005 * diff(range(ex[['cut']]))
      if (
        is.na(cp) || length(ex[["cut"]]) < min_cell
      ) {
        cp <- max(cp_min, max(ex[["cut"]]) +
          (max(ex[["cut"]]) - min(ex[["cut"]])) / 5)
      }
      cp_vec <- stats::setNames(rep(cp, length(ind_gate)), ind_gate)
    }

    cp_list <- cp_list |> append(list(prejoin = cp_vec))
  }

  # get cutpoint if group method is not prejoin
  non_prejoin_combn_vec <- setdiff(gate_combn, "prejoin")

  if (length(non_prejoin_combn_vec) > 0) {
    .debug(debug, "non-prejoin")
    cp_tg_vec <- purrr::map_dbl(ind_gate, function(ind) {
      .debug(debug, "ind", ind)
      # print(ind)
      ex <- ex_list[[as.character(ind)]]
      ex <- ex |> dplyr::filter(!is.na(cut))
      if (exc_min) ex <- ex |> dplyr::filter(cut > min(cut))
      if (nrow(ex) < max(min_cell, 5)) {
        return(NA)
      }
      if (length(ex[["cut"]]) < min_cell) {
        return(
          max(cp_min, max(ex[["cut"]]) +
            (max(ex[["cut"]]) - min(ex[["cut"]])) / 5)
        )
      }
      dens <- density(ex[["cut"]])
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
      .ensure_cytoutils() # nolint
      cytoUtils:::.cytokine_cutpoint(
        x = ex[["cut"]], num_peaks = 1,
        ref_peak = 1, tol = tol, side = "right",
        strict = FALSE, adjust = adjust
      )
    }) |>
      stats::setNames(ind_gate)

    .debug(debug, "combining thresholds") # nolint

    cp_tg_list <- .combine_cp(
      cp = cp_tg_vec,
      gate_combn = gate_combn
    ) |>
      stats::setNames(gate_combn)

    cp_list <- cp_list |> append(cp_tg_list)
  }
  .debug(debug, "Done tg cutpoint")

  cp_list
}



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
      family = binomial, data = high_ind_tbl
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
  fit_h0 <- glm(high ~ 1, family = binomial, data = high_ind_tbl)
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
      # fit <- lm(count ~ high_ind + high_diff_from_cut + low_diff_from_cut, data = final_mod_tbl)
      fit <- glm(count ~ high_diff_from_cut + low_diff_from_cut, data = final_mod_tbl, family = "poisson")
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
      fit <- glm(count ~ splines::bs(ex, knots = knots, degree = 1), data = final_mod_tbl, family = "poisson")
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





#' @title Get counts and frequencies for across populations and samples, with cutpoints set separately for groups of samples
#'
#' @description
#' Gets cell counts and frequencies for samples using pre-specified cutpoints.
#' The cutpoints are to be calculated before this function is applied using
#' \code{get_cp}, and are extracted
#' from the list returned by \code{get_cp}. This means that, for example, that a
#' single cutpoint can be calculated using a group of samples, and then we
#' can apply that cutpoint on the constituent samples to find counts and frequencies
#' on a per-sample (rather than grouped) basis.
#'
#' @inheritParams .get_cp
#' @param data GatingSet. Each element corresponds to a single sample to which the cutpoint(s) in
#' a specified \code{cp_obj} output must be applied.
#' @param pop character vector. Each element fully specifies the parent population of cells which must be
#' divided into positive/negative cells based on their individual values for the variable specified by \code{cut}
#' and the cutpoint(s) specified by \code{cp_obj}.
#' @param ind list. A list, where each element is a numeric vector specifying the indices of the
#' GatingHierarchy objects in \code{data} to get the counts and frequency for the
#' cutpoints specified by \code{cp} for.
#' @param cp_obj list. A list, where each element is a list returned by \code{get_cp}.
#'
#' @details
#' Note that the f-score assumes that the smaller count of cells is a strict subset of
#' the larger count of cells. This holds when the difference between manual and automatic
#' gating is due to differences in the placement of a single univariate cutpoint,
#' but does not hold in general.
#'
#' Note that this differs from .get_cp_stats_tbl_pop_samples_elem in that
#' it may apply different cutpoints to different sets of samples within \code{data}. \code{.get_cp_stats_tbl_pop_samples_elem}
#' applies the single set of cutpoints it receives to all specified samples.
#'
#' @return
#' @return A tibble with columns fcs, pop, type, cp, count and freq, with corresponding
#' values being the fcs file gated on, parent population, cutpoint type,
#' cutpoint value, count of positive cells and frequency of positive cells.
.get_cp_stats_tbl_pop_samples <- function(data, ind, pop, wins, cut, high, cp_obj) {
  purrr::map_df(seq_along(cp_obj), function(i) {
    # get cp_obj that provides the cutpoints
    cp_obj_curr <- cp_obj[[i]]

    # get initial cp_vec
    cp_vec <- stats::setNames(cp_obj_curr[["cp_stats"]][["cp"]], cp_obj_curr[["cp_stats"]][["type"]])

    # get statistics (count and frequency) from different pop(s) from different samples
    cp_stats_tbl_pop_samples <- .get_cp_stats_tbl_pop_samples_elem(
      data = data,
      ind = ind[[i]],
      pop = pop,
      wins = wins,
      cut = cut,
      high = high,
      cp = cp_vec
    )
  })
}

#' @title Get mid-probability cut
#'
#'
.get_cp_pwmid <- function(high_ind_tbl, cp_scp) {
  # get table to model - all values above changepoint
  mod_tbl <- high_ind_tbl |>
    dplyr::filter(cut >= cp_scp)

  # check if too little data to model with here
  if (nrow(mod_tbl) < 5 || length(unique(high_ind_tbl$high)) == 1) {
    print("cp_pwmid set to cp_scp + 5 due to too few obs above cp_sc or too few postive obs")
    return(cp_scp + 5)
  }

  # model the probability of positivity above this point
  if (nrow(mod_tbl) < 40) {
    fit_pw <- glm(high ~ cut,
      family = binomial, data = mod_tbl
    )
  } else {
    fit_pw <- glm(high ~ splines::ns(cut, df = 3),
      family = binomial, data = mod_tbl
    )
  }

  # get predictions over a range of cut values
  pred_tbl <- tibble::tibble(cut = seq(min(mod_tbl$cut), max(mod_tbl$cut)))
  pred_vec <- predict(fit_pw, pred_tbl, type = "response")
  pred_tbl <- pred_tbl |> dplyr::mutate(pred = pred_vec)

  # find prediction in between middle and max
  min <- mean(
    high_ind_tbl |> dplyr::filter(cut < cp_scp) |> dplyr::pull("high")
  )
  max <- max(pred_tbl$pred)
  mid_prob <- mean(c(max, min))

  # find the cutpoint that minimises the difference
  # objective function
  opt_func <- function(x) {
    pred_tbl <- tibble::tibble(cut = x)
    pred_vec <- predict(fit_pw, pred_tbl, type = "response")
    (pred_vec - mid_prob)^2
  }
  # initial search parameter
  init_par <- mean(c(cp_scp, max(high_ind_tbl$cut)))
  # optimisation
  optim(init_par,
    fn = opt_func,
    method = "Brent", lower = cp_scp,
    upper = max(high_ind_tbl$cut)
  )$par
}







#' @title Get axis labels from annotated data frame
#'
#' @description Get axis labels for cut and high channels
#' using the annotated data frame.
#'
#' @param data GatingHierarchy.
#' Code written to take a GatingHierarchy at present, but this can easily be extended.
#' @param inheritParams get_cp
#'
#' @return A named vector, where the names are the channel names and the
#' values are the corresponding marker (common) names.
.get_labs <- function(data, cut, high = NULL) {
  force(data)
  adf_data <- flowWorkspace::gh_pop_get_data(data) |>
    flowCore::parameters() |>
    flowCore::pData()

  if (!is.null(high)) {
    cut_lab <- adf_data[["desc"]][[which(adf_data$name == cut)]] |>
      stats::setNames(cut)
    high_lab_vec <- purrr::map(seq_along(high), function(i) {
      high_lab <- adf_data[["desc"]][[
        which(adf_data$name == names(high)[i])
      ]] |>
        stats::setNames(names(high)[i])
    }) |>
      unlist()
    return(c(cut_lab, high_lab_vec))
  }

  purrr::map_chr(cut, function(cut_curr) {
    adf_data[["desc"]][[which(adf_data$name == cut_curr)]]
  }) |>
    stats::setNames(cut)
}


#' @title Saves objects within a list
#'
#' @description
#' Save the parameters as a CSV file and as an RDS object.
#'
#' @param dir character. Directory in which to save parameters. Note that the parameters will be saved to a sub-folder named "params" inside this directory.
#' @param obj_name character. Name to give to saved CSV and RDS files.
#' @param list_to_save named list. List to be saved.
#' @param tbl_elems_to_save character vector. Each element is the name of an element within \code{list_to_save} to be saved as part of a table. If \code{NULL}, then all elements are saved. If no elements in \code{tbl_elems_to_save} are found in \code{list_to_save}, then no dataframe is saved.
#' @param list_elems_to_save character vector. Each element is the name of an element within \code{list_to_save} to be saved as a R-readable list. If \code{NULL}, then all elements are saved. If no elements in \code{list_elems_to_save} are found in \code{list_to_save}, then no list is saved.
#'
#' @return Returns \code{TRUE} invisibly.
.save_list <- function(dir, obj_name, list_to_save,
                       tbl_elems_to_save, list_elems_to_save) {
  # ===============================
  # Preparation
  # ===============================

  # check that list elems to save and tbl elems to save are provided
  if (missing(tbl_elems_to_save)) {
    stop("tbl_elems_to_save must be specified")
  }
  if (missing(list_elems_to_save)) {
    stop("list_elems_to_save must be specified")
  }

  # get names of elems to save, if need be
  if (is.null(tbl_elems_to_save)) tbl_elems_to_save <- names(list_to_save)
  if (is.null(list_elems_to_save)) list_elems_to_save <- names(list_to_save)

  # ===============================
  # Create dataframe to save
  # ===============================

  # get length that each dataframe column will need to be
  max_n <- purrr::map_dbl(tbl_elems_to_save, function(x) {
    length(list_to_save[[x]])
  }) |>
    max()

  # create dataframe to save
  save_tbl <- tibble::tibble(x = rep(NA, max_n))
  # loop over parameters to save
  for (x in tbl_elems_to_save) {
    # if the parameter name to save is actually in params
    if (x %in% names(list_to_save)) {
      # get param element
      list_elem <- list_to_save[[x]]

      # if the param element has names
      if (!is.null(names(list_elem))) {
        # create column name to save under
        col_name <- paste0(x, "_name")
        # get column values, and extend to required length if need be
        col_val <- names(list_elem)
        if (length(col_val) < max_n) col_val <- c(col_val, rep("", max_n - length(col_val)))
        # add column to table to save
        save_tbl <- save_tbl |> dplyr::mutate(!!ensym(col_name) := col_val)
      }

      # if the param element has names
      # create column name to save under
      col_name <- paste0(x, "_val")
      # get column values, and extend to required length if need be
      col_val <- stats::setNames(list_elem, NULL)
      if (length(col_val) < max_n) col_val <- c(col_val, rep("", max_n - length(col_val)))
      # add column to table to save
      save_tbl <- save_tbl |> dplyr::mutate(!!ensym(col_name) := col_val)
    }
  }

  # remove placeholder column
  save_tbl <- save_tbl |> dplyr::select(-x)

  # ===============================
  # Create list to save as an R object
  # ===============================

  list_elems_to_save_ind_vec <- list_elems_to_save %in% names(list_to_save)
  list_elems_to_save_found_vec <- list_elems_to_save[list_elems_to_save_ind_vec]
  save_list <- purrr::map(list_elems_to_save_found_vec, function(x) list_to_save[[x]])

  # ===============================
  # Output
  # ===============================

  # create directory to save to, if need be
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  # save the csv of parameter values created above
  if (ncol(save_tbl) > 0) {
    write.csv(save_tbl, file.path(dir, paste0(obj_name, ".csv")),
      row.names = FALSE
    )
  }
  # save the parameter object in an R-readable format
  if (length(save_list) > 0) saveRDS(save_list, file.path(dir, obj_name))

  # return invisibly
  invisible(TRUE)
}
