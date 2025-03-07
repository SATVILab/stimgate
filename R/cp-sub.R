#' @title Check that parameters for each marker for which a gate is required are complete
#'
#' @inheritParams gate # marker
#' @data_name 'gs_cytof' or `gs_proto`. Name of dataset to be gated in R environment.
#'
#' @return

.complete_marker_list <- function(marker,
                                  bias_uns,
                                  .data,
                                  pop_gate,
                                  ind_batch_list,
                                  bw_min,
                                  cp_min,
                                  min_cell,
                                  tol,
                                  max_pos_prob_x,
                                  gate_combn,
                                  marker_settings,
                                  debug) {
  marker_settings_common <- list(
    bias_uns = bias_uns, cp_min = cp_min, bw_min = bw_min,
    min_cell = min_cell, tol = tol, gate_combn = gate_combn,
    max_pos_prob_x = max_pos_prob_x
  )
  purrr::map(marker, function(marker_curr) {
   .complete_marker_list_ind(
      marker = marker_curr,
      marker_settings_common = marker_settings_common,
      marker_settings_spec = list(cut = marker_curr) |>
        append(marker_settings[[marker_curr]]),
      .data = .data,
      pop_gate = pop_gate,
      debug = debug,
      ind_batch_list = ind_batch_list
    )
  })
}

.complete_marker_list_ind <- function(marker_settings_common,
                                      marker_settings_spec,
                                      marker,
                                      .data,
                                      pop_gate,
                                      debug,
                                      ind_batch_list) {

  marker_settings <- .complete_marker_list_add_common(
    marker_settings_common = marker_settings_common,
    marker_settings = marker_settings_spec
  )

  marker_settings$bias_uns <- .complete_marker_list_bias_uns(
    bias_uns = marker_settings$bias_uns,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    debug = debug,
    ind_batch_list = ind_batch_list
  )

  marker_settings$bw_min <- .complete_marker_list_min_bw(
    bw_min = marker_settings$bw_min,
    bias_uns = max(marker_settings$bias_uns)
  )

  marker_settings$cp_min <- .complete_marker_list_cp_min(
    cp_min = marker_settings$cp_min,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = marker,
    debug = debug,
    ind_batch_list = ind_batch_list
  )

  marker_settings
}

.complete_marker_list_add_common <- function(marker_settings_common,
                                             marker_settings) {
  marker_settings |>
    append(marker_settings_common[
      setdiff(names(marker_settings_common), names(marker_settings))
    ])
}

.complete_marker_list_bias_uns <- function(bias_uns,
                                           .data,
                                           pop_gate,
                                           chnl_cut,
                                           debug,
                                           ind_batch_list) {
  if (!is.null(bias_uns)) {
    return(bias_uns)
  }
  .debug(debug, "calculating bias_uns automatically") # nolint
  mean_range <- .complete_marker_list_bias_uns_get_mean_range(
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut
  )
  (mean_range / 12) |> signif(3)
}

.complete_marker_list_bias_uns_get_mean_range <- function(ind_batch_list,
                                                          .data,
                                                          pop_gate,
                                                          chnl_cut) {
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i]
      )
      purrr::map_dbl(ex_list, function(ex) {
        abs(
          diff(quantile(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], c(0.99, 0.01)),
            na.rm = TRUE
          )
        )[[1]]
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
  bias_uns * 2.25
}

.complete_marker_list_cp_min <- function(cp_min,
                                         .data,
                                         pop_gate,
                                         chnl_cut,
                                         debug,
                                         ind_batch_list) {
  if (!is.null(cp_min)) {
    return(cp_min)
  }
  .debug(debug, "calculating cp_min automatically") # nolint
  purrr::map(
    seq_len(min(2, length(ind_batch_list))),
    function(i) {
      ex_list <- .get_ex_list( # nolint
        .data = .data,
        ind_batch = ind_batch_list[[i]],
        pop = pop_gate,
        chnl_cut,
        batch = names(ind_batch_list)[i]
      )
      purrr::map_dbl(ex_list, function(ex) {
        median(.get_cut(ex)[.get_cut(ex) > min(.get_cut(ex))], na.rm = TRUE)[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
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
        expr = chnl_cut
      ) |>
      dplyr::select(ind, ind_cell, expr, everything()) # nolint
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
                       chnl_cut,
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
    ex <- ex |> dplyr::filter(!is.na(chnl_cut))
    if (exc_min) ex <- ex |> dplyr::filter(chnl_cut > min(chnl_cut))
    if (nrow(ex) < max(min_cell, 5)) {
      cp_vec <- stats::setNames(rep(NA, length(ind_gate)), ind_gate)
    } else {
      dens <- density(.get_cut(ex))
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
      .ensure_cytoutils()
      cp <- cytoUtils:::.cytokine_cutpoint(
        x = .get_cut(ex), num_peaks = 1,
        ref_peak = 1, tol = tol, side = "right",
        strict = FALSE, adjust = adjust
      )
      if (
        is.na(cp) || length(.get_cut(ex)) < min_cell
      ) {
        cp <- max(cp_min, max(.get_cut(ex)) +
          (max(.get_cut(ex)) - min(.get_cut(ex))) / 5)
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
      ex <- ex |> dplyr::filter(!is.na(chnl_cut))
      if (exc_min) ex <- ex |> dplyr::filter(chnl_cut > min(chnl_cut))
      if (nrow(ex) < max(min_cell, 5)) {
        return(NA)
      }
      if (length(.get_cut(ex)) < min_cell) {
        return(
          max(cp_min, max(.get_cut(ex)) +
            (max(.get_cut(ex)) - min(.get_cut(ex))) / 5)
        )
      }
      dens <- density(.get_cut(ex))
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
      .ensure_cytoutils() # nolint
      cytoUtils:::.cytokine_cutpoint(
        x = .get_cut(ex), num_peaks = 1,
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
#' @param .data GatingSet. Each element corresponds to a single sample to which the cutpoint(s) in
#' a specified \code{cp_obj} output must be applied.
#' @param pop character vector. Each element fully specifies the parent population of cells which must be
#' divided into positive/negative cells based on their individual values for the variable specified by \code{cut}
#' and the cutpoint(s) specified by \code{cp_obj}.
#' @param ind list. A list, where each element is a numeric vector specifying the indices of the
#' GatingHierarchy objects in \code{.data} to get the counts and frequency for the
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
#' it may apply different cutpoints to different sets of samples within \code{.data}. \code{.get_cp_stats_tbl_pop_samples_elem}
#' applies the single set of cutpoints it receives to all specified samples.
#'
#' @return
#' @return A tibble with columns fcs, pop, type, cp, count and freq, with corresponding
#' values being the fcs file gated on, parent population, cutpoint type,
#' cutpoint value, count of positive cells and frequency of positive cells.
.get_cp_stats_tbl_pop_samples <- function(.data, ind, pop, wins, chnl_cut, high, cp_obj) {
  purrr::map_df(seq_along(cp_obj), function(i) {
    # get cp_obj that provides the cutpoints
    cp_obj_curr <- cp_obj[[i]]

    # get initial cp_vec
    cp_vec <- stats::setNames(cp_obj_curr[["cp_stats"]][["cp"]], cp_obj_curr[["cp_stats"]][["type"]])

    # get statistics (count and frequency) from different pop(s) from different samples
    cp_stats_tbl_pop_samples <- .get_cp_stats_tbl_pop_samples_elem(
      .data = .data,
      ind = ind[[i]],
      pop = pop,
      wins = wins,
      chnl_cut,
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
    dplyr::filter(chnl_cut >= cp_scp)

  # check if too little .data to model with here
  if (nrow(mod_tbl) < 5 || length(unique(high_ind_tbl$high)) == 1) {
    print("cp_pwmid set to cp_scp + 5 due to too few obs above cp_sc or too few postive obs")
    return(cp_scp + 5)
  }

  # model the probability of positivity above this point
  if (nrow(mod_tbl) < 40) {
    fit_pw <- glm(high ~ chnl_cut,
      family = binomial, .data = mod_tbl
    )
  } else {
    fit_pw <- glm(high ~ splines::ns(chnl_cut, df = 3),
      family = binomial, .data = mod_tbl
    )
  }

  # get predictions over a range of cut values
  pred_tbl <- tibble::tibble(chnl_cut = seq(min(mod_tbl$chnl_cut), max(mod_tbl$chnl_cut)))
  pred_vec <- predict(fit_pw, pred_tbl, type = "response")
  pred_tbl <- pred_tbl |> dplyr::mutate(pred = pred_vec)

  # find prediction in between middle and max
  min <- mean(
    high_ind_tbl |> dplyr::filter(chnl_cut < cp_scp) |> dplyr::pull("high")
  )
  max <- max(pred_tbl$pred)
  mid_prob <- mean(c(max, min))

  # find the cutpoint that minimises the difference
  # objective function
  opt_func <- function(x) {
    pred_tbl <- tibble::tibble(chnl_cut = x)
    pred_vec <- predict(fit_pw, pred_tbl, type = "response")
    (pred_vec - mid_prob)^2
  }
  # initial search parameter
  init_par <- mean(c(cp_scp, max(high_ind_tbl$chnl_cut)))
  # optimisation
  optim(init_par,
    fn = opt_func,
    method = "Brent", lower = cp_scp,
    upper = max(high_ind_tbl$chnl_cut)
  )$par
}







#' @title Get axis labels from annotated .data frame
#'
#' @description Get axis labels for cut and high channels
#' using the annotated .data frame.
#'
#' @param .data GatingHierarchy.
#' Code written to take a GatingHierarchy at present, but this can easily be extended.
#' @param inheritParams get_cp
#'
#' @return A named vector, where the names are the channel names and the
#' values are the corresponding marker (common) names.
.get_labs <- function(.data, chnl_cut, high = NULL) {
  force(.data)
  adf_data <- flowWorkspace::gh_pop_get_data(.data) |>
    flowCore::parameters() |>
    flowCore::pData()

  if (!is.null(high)) {
    cut_lab <- adf_data[["desc"]][[which(adf_data$name == chnl_cut)]] |>
      stats::setNames(chnl_cut)
    return(cut_lab)
  }

  purrr::map_chr(chnl_cut, function(cut_curr) {
    adf_data[["desc"]][[which(adf_data$name == cut_curr)]]
  }) |>
    stats::setNames(chnl_cut)
}

