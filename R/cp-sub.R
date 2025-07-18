
#' @title Get the measurements for a given set of sample(s) for the marker to be cut on
#' as a single numeric vector
#'
#' @return A list where each element is a numeric vector.
.prepare_ex_list_with_bias_and_noise <- function(ex_list,
                                                 ind,
                                                 exc_min,
                                                 bias = 0,
                                                 .debug = FALSE,
                                                 noise_sd = NULL) {
  purrr::map(ind, function(ind_curr) {
    cut_tbl <- ex_list[[as.character(ind_curr)]]
    if (exc_min) {
      n_row_init <- nrow(cut_tbl)
      cut_tbl <- cut_tbl[
        .get_cut(cut_tbl) > min(.get_cut(cut_tbl)),
      ] # nolint
      n_row_fin <- nrow(cut_tbl)
      attr(cut_tbl, "prob_g_min") <- n_row_fin / n_row_init
    }
    cut_tbl[[attr(cut_tbl, "chnl_cut")]] <- .get_cut(cut_tbl) + bias # nolint
    if (!is.null(noise_sd)) {
      cut_tbl <- cut_tbl[[attr(cut_tbl, "chnl_cut")]] +
        rnorm(nrow(cut_tbl), sd = noise_sd) # nolint
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

.get_cp_tg <- function(ex_list,
                       gate_combn,
                       chnl_cut,
                       exc_min,
                       tol,
                       min_cell,
                       cp_min,
                       bw,
                       .debug) {
  # get cytoUtils tailgate cutpoint
  .debug_msg(.debug, "Getting tg cutpoint")
  cp_list <- list()

  if ("prejoin" %in% gate_combn) {
    .debug_msg(.debug, "prejoin")
    ind_gate <- names(ex_list)[-length(ex_list)]
    ex <- dplyr::bind_rows(ex_list[ind_gate])
    ex <- ex[!is.na(.get_cut(ex)), ]
    if (exc_min) ex <- ex[.get_cut(ex) > min(.get_cut(ex)), ]
    if (nrow(ex) < max(min_cell, 5)) {
      cp_vec <- stats::setNames(rep(NA, length(ind_gate)), ind_gate)
    } else {
      dens <- density(.get_cut(ex))
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
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
    .debug_msg(.debug, "non-prejoin")
    ind_gate <- names(ex_list)[-length(ex_list)]
    cp_tg_vec <- purrr::map_dbl(ind_gate, function(ind) {
      .debug_msg(.debug, "ind", ind)
      # print(ind)
      ex <- ex_list[[as.character(ind)]]
      ex <- ex[!is.na(.get_cut(ex)), ]
      if (exc_min) ex <- ex[.get_cut(ex) > min(.get_cut(ex)), ]
      if (nrow(ex) < max(min_cell, 5)) {
        return(
          max(cp_min, max(.get_cut(ex)) +
            (max(.get_cut(ex)) - min(.get_cut(ex))) / 5)
        )
      }
      dens <- density(.get_cut(ex))
      adjust <- ifelse(dens$bw < bw, bw / dens$bw, 1)
      cytoUtils:::.cytokine_cutpoint(
        x = .get_cut(ex), num_peaks = 1,
        ref_peak = 1, tol = tol * 1e3, side = "right",
        strict = FALSE, adjust = adjust
      )
    }) |>
      stats::setNames(ind_gate)

    .debug_msg(.debug, "combining thresholds") # nolint

    cp_tg_list <- .combine_cp(
      cp = cp_tg_vec,
      gate_combn = gate_combn
    ) |>
      stats::setNames(gate_combn)

    cp_list <- cp_list |> append(cp_tg_list)
  }
  .debug_msg(.debug, "Done tg cutpoint")

  cp_list
}



#' @title Get mid-probability cut
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


.combine_cp <- function(cp, gate_combn) {
  purrr::map(gate_combn, function(gate_combn_curr) {
    if (all(purrr::map_lgl(cp, is.na))) {
      return(stats::setNames(cp, names(cp)))
    }
    if (is.null(gate_combn_curr) || gate_combn_curr %in% c("no", "prejoin")) {
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
