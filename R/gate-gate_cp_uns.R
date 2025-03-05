# get middle pred value
#' @title Get unstim-based cutpoints
#'
#' @inheritParams get_cp # ind, ind_in_batch_uns, cut, fdr
#' @inheritParams .get_cp_man # ex_list
#' @params ind_gate, ind_uns. Numeric vectors. Numeric vectors where
#' each element specifies the index of samples in \code{.data}
#' for which gates are to be calculated (ind_gate) or is the unstim (ind_uns) in this batch.
#' @return
.get_cp_uns <- function(ex_list,
                        ind_gate,
                        ind_uns,
                        .data,
                        fdr,
                        gate_combn,
                        bias_uns = 5,
                        noise_sd = 2,
                        cp_min,
                        min_cell) {
  # get cutpoints for each level of bias
  .get_cp_uns_bias(
    ex_list = ex_list, ind_gate = ind_gate,
    ind_uns = ind_uns, .data = .data, bias_uns = bias_uns,
    noise_sd = noise_sd, cp_min = cp_min, gate_combn = gate_combn,
    fdr = fdr, min_cell = min_cell
  )
}

#' @title Get the unstim-based FDR-method cutpoint for
.get_cp_uns_bias <- function(ex_list, ind_gate, ind_uns,
                             .data, bias_uns, noise_sd, cp_min,
                             gate_combn, fdr, min_cell) {
  # get ecdf of uns
  purrr::map(bias_uns, function(bias) {
    # get ecdf of uns sample
    cut_vec_uns <- .prepare_ex_list_with_bias_and_noise( # nolint
      ex_list = ex_list, ind = ind_uns,
      exc_min = TRUE, bias = bias,
      noise_sd = noise_sd
    )[[1]]
    ecdf_uns <- ecdf(cut_vec_uns)

    # get gates for given level of bias across gate combination methods
    cp_uns_gate_combn_list <- .get_cp_uns_gate_combn(
      ex_list = ex_list,
      ind_uns = ind_uns, ind_gate = ind_gate,
      exc_min = TRUE, gate_combn = gate_combn,
      ecdf_uns = ecdf_uns, cp_min = cp_min, fdr = fdr,
      min_cell = min_cell
    )
    names(cp_uns_gate_combn_list) <-
      paste0(names(cp_uns_gate_combn_list), "b", bias)
    cp_uns_gate_combn_list
  }) |>
    purrr::flatten()
}

#' @title Get the unstim-based FDR-method cutpoint for
.get_cp_uns_gate_combn <- function(ex_list = ex_list,
                                   ind_uns = ind_uns, ind_gate = ind_gate,
                                   exc_min = TRUE, gate_combn,
                                   ecdf_uns, cp_min, fdr, min_cell) {
  # get cutpoints for prejoin gate combination method
  if ("prejoin" %in% gate_combn) {
    # get marker expression for stim samples, join and then sort into descending order
    cut_list_stim <- .prepare_ex_list_with_bias_and_noise( # nolint
      ex_list = ex_list, ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE, bias = 0,
      noise_sd = NULL
    ) |>
      unlist() |>
      sort() |>
      rev() |>
      list()

    # get cutpoints for gate combn method for a range of fdr's
    cp_uns_list_prejoin <- .get_cp_uns_fdr(
      cut_stim = cut_list_stim,
      ecdf_uns = ecdf_uns,
      fdr = fdr,
      cp_min = cp_min,
      ind_gate = ind_gate,
      ind_uns = ind_uns, min_cell = min_cell
    ) |>
      purrr::map(function(x) list("prejoin" = x))
  } else {
    cp_uns_list_prejoin <- list()
  }

  # get cutpoint if group method is not only prejoin
  non_prejoin_combn_vec <- setdiff(gate_combn, "prejoin")

  if (length(non_prejoin_combn_vec) > 0) {
    # get marker expression for stim samples, and sort into descending order
    cut_list_stim <- .prepare_ex_list_with_bias_and_noise( # nolint
      ex_list = ex_list, ind = setdiff(ind_gate, ind_uns),
      exc_min = TRUE, bias = 0,
      noise_sd = NULL
    ) |>
      purrr::map(function(x) {
        x |>
          sort() |>
          rev()
      })

    cp_uns_list_nonjoin <- .get_cp_uns_fdr(
      cut_stim = cut_list_stim,
      ecdf_uns = ecdf_uns,
      fdr = fdr,
      cp_min = cp_min,
      ind_gate = ind_gate,
      ind_uns = ind_uns,
      min_cell = min_cell
    )

    cp_uns_list_nonjoin <- purrr::map(cp_uns_list_nonjoin, function(x) {
      .combine_cp(cp = x, gate_combn = non_prejoin_combn_vec) # nolint
    }) |>
      stats::setNames(names(cp_uns_list_nonjoin))
  } else {
    cp_uns_list_nonjoin <- list()
  }

  combined_list <- cp_uns_list_prejoin |> append(cp_uns_list_nonjoin)
  cp_uns_list <- purrr::map(unique(names(combined_list)), function(x) {
    cp_uns_list_prejoin[[x]] |>
      append(cp_uns_list_nonjoin[[x]])
  }) |>
    stats::setNames(unique(names(combined_list)))

  cp_uns_list
}

#' @title Get cutpoints for a range of FDRs across sample(s)
#'
#' @inheritParams gate cp_min, fdr
#' @inheritParams .get_cp_uns ind_gate
#' @param ecdf_uns function. Empirical CDF for the unstim stimulation sample
#' for this blood sample.
#' @param cut_stim list. Each element must be a numeric vector
#' that corresponds to the readings for the marker for a sample for which a cutpoint
#' is to be  calculated. Each such numeric vector must be ordered
#' in descending order.
#'
#' @return A named list, where each element is a named vector
#' referring to the cutpoints for all the samples for which
#' gates are required for this batch for a given FDR.
.get_cp_uns_fdr <- function(cut_stim, ecdf_uns, fdr, cp_min,
                            ind_gate, ind_uns, min_cell) {
  # get list of q-values
  q_list <- purrr::map(cut_stim, function(cut_vec_stim) {
    p_vec <- .rep_zero_p_val(1 - ecdf_uns(cut_vec_stim)) # nolint
    p.adjust(p_vec, method = "BH")
  })

  # get cutpoint for each fdr
  purrr::map(fdr, function(fdr_curr) {
    # get cutpoint for each sample
    cp_uns_fdr <- .get_cp_uns_sample(
      q_list = q_list, cut_stim = cut_stim,
      fdr = fdr_curr, cp_min = cp_min, min_cell = min_cell
    )

    # Repeat cutpoint if samples were prejoined
    ind_stim <- setdiff(ind_gate, ind_uns)
    if (length(cp_uns_fdr) == 1 && ind_stim != 1) {
      cp_uns_vec_fdr <- stats::setNames(rep(cp_uns_fdr, length(ind_stim)), ind_stim)
    } else {
      cp_uns_vec_fdr <- stats::setNames(cp_uns_fdr, ind_stim)
    } # name gate indices if not prejoined

    # if unstim is in the samples to be gated, then
    # set the cutpoint for unstim as the mean of all the cutpoints
    if (ind_uns %in% ind_gate) cp_uns_vec_fdr <- c(cp_uns_vec_fdr, stats::setNames(mean(cp_uns_vec_fdr), ind_uns))

    cp_uns_vec_fdr
  }) |>
    stats::setNames(.get_cp_uns_fdr_name_vec(fdr))
}

#' @title Get cutpoint for a range of samples given the q-value and fdr
#'
#' @description Calculate the cutpoint for each sample in a batch at a given FDR.
#'
#' @param q_list list. List where each element is a numeric vector specifying the
#' q-values for a given sample.
#' @param cut_stim list. List where each element are the marker expression readings
#' of the marker to be cut on for the cells in a sample. Note that the i-th element
#' in \code{cut_stim} must correspond to the i-th element in \code{q_list}, i.e.
#' must be related to the same marker in same cell population from the same blood sample and stimulation.
#' @param fdr numeric. A value between 0 and 1 specifying the false discovery rate
#' the sample should be cut at.
#'
#' @return Numeric vector. A cutpoint for each sample.
.get_cp_uns_sample <- function(q_list, cut_stim, fdr, cp_min, min_cell) {
  purrr::map_dbl(seq_along(q_list), function(i) {
    # get vector of q-values
    q_vec <- q_list[[i]]

    # get expression matrix
    stim_vec <- cut_stim[[i]]

    # get count and frequency .data
    count_pos <- sum(q_vec <= fdr)

    range_vec <- c(min(stim_vec), max(stim_vec))

    # return a value slightly higher than last value if count is zero or if all values are selected
    if (count_pos == 0 || count_pos == length(q_vec) || length(cut_stim) < min_cell) {
      return(max(cp_min, range_vec[2] + (range_vec[2] - range_vec[1]) / 5))
    }

    max((stim_vec[count_pos] + stim_vec[count_pos + 1]) / 2, cp_min)
  })
}
