#' @title Calculate the local fdr-based cut
.get_cp_uns_loc <- function(ex_list, ind_gate, ind_uns, gate_combn,
                            pop_root = NULL, data, bias_uns = 0,
                            noise_sd = NULL, min_bw = 80,
                            cp_min, min_cell, params, plot,
                            path_project, debug = FALSE) {
  # get cutpoints for each level of bias
  .get_cp_uns_loc_bias( # nolint
    ex_list = ex_list, ind_gate = ind_gate, ind_uns = ind_uns,
    data = data, bias_uns = bias_uns, noise_sd = noise_sd,
    cp_min = cp_min, gate_combn = gate_combn, min_bw = min_bw,
    min_cell = min_cell, gate_tbl = params$gate_tbl,
    gate_name_curr = params$gate_name_curr, cut = params$cut,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates, plot = plot,
    path_project = path_project, debug = debug
  )
}



#' @title Get the unstim-based local fdr-method cutpoint for each level of bias
.get_cp_uns_loc_bias <- function(ex_list, ind_gate, ind_uns, data,
                                 bias_uns, noise_sd, cp_min,
                                 gate_combn, min_bw, min_cell,
                                 gate_tbl, gate_name_curr, cut,
                                 calc_cyt_pos_gates, plot, path_project,
                                 debug) {
  # get ecdf of uns
  purrr::map(bias_uns, function(bias) {
    .debug(debug, "bias_uns", bias) # nolint

    ex_list_prep <- .get_cp_uns_loc_bias_data_prep( # nolint
      ex_list = ex_list, ind_gate = ind_gate, ind_uns = ind_uns,
      bias = bias, noise_sd = noise_sd, debug = debug
    )

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cp_uns_gate_combn_obj <- .get_cp_uns_loc_gate_combn( # nolint
      ex_list_orig = ex_list_prep[["ex_list_orig"]],
      ex_list_no_min = ex_list_prep[["ex_list_no_min"]],
      ex_tbl_uns_bias = ex_list_prep[["ex_tbl_uns_bias"]],
      ind_uns = ind_uns,
      ind_gate = ind_gate,
      exc_min = TRUE,
      gate_combn = gate_combn,
      cp_min = cp_min,
      min_bw = min_bw,
      min_cell = min_cell,
      gate_tbl = gate_tbl,
      gate_name_curr = gate_name_curr,
      cut = cut,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )

    # extract and add bias label to gates
    # ---------------------------------------
    .get_cp_uns_loc_gate_label(cp_uns_gate_combn_obj, bias)
  }) |>
    purrr::flatten()
}

.get_cp_uns_loc_gate_label <- function(cp_uns_gate_combn_obj, bias) {
  cp_uns_gate_combn_list <- cp_uns_gate_combn_obj[["cp_uns"]]
  names(cp_uns_gate_combn_list) <-
    paste0(names(cp_uns_gate_combn_list), "b", bias)
  cp_uns_gate_combn_list
}

.get_cp_uns_loc_bias_data_prep <- function(ex_list,
                                           ind_gate,
                                           ind_uns,
                                           bias,
                                           noise_sd,
                                           debug) {
  # rename `cut` column` to `expr`
  # -------------------------------------
  ex_list_orig <- .get_cut_list( # nolint
    ex_list = ex_list, ind = union(ind_gate, ind_uns), exc_min = FALSE,
    bias = 0, noise_sd = NULL, debug = debug
  )

  # separate stim and uns samples, rename
  # `cut` column to `expr` and exclude min
  # values
  # -------------------------------------
  ex_list_no_min <- .get_cut_list( # nolint
    ex_list = ex_list, ind = union(ind_gate, ind_uns), exc_min = TRUE,
    bias = 0, noise_sd = NULL, debug = debug
  )

  # adjust expression in unstim,
  # applying bias, excluding the min val
  # and /or adding noise
  # -------------------------------------
  ex_tbl_uns_bias <- .get_cut_list( # nolint
    ex_list = ex_list, ind = ind_uns, exc_min = TRUE,
    bias = bias, debug = debug, noise_sd = NULL
  )[[1]]

  list(
    "ex_list_orig" = ex_list_orig,
    "ex_list_no_min" = ex_list_no_min,
    "ex_tbl_uns_bias" = ex_tbl_uns_bias
  )
}


#' @title Get the unstim-based local fdr-method
#' cutpoint for a given bias across gate combination methods
.get_cp_uns_loc_gate_combn <- function(ex_list_orig,
                                       ex_list_no_min,
                                       ex_tbl_uns_bias,
                                       ind_uns,
                                       ind_gate,
                                       exc_min = TRUE,
                                       gate_combn,
                                       cp_min,
                                       min_bw,
                                       min_cell,
                                       gate_tbl,
                                       gate_name_curr,
                                       cut,
                                       calc_cyt_pos_gates,
                                       plot,
                                       bias,
                                       path_project,
                                       debug = FALSE) {
  .debug(debug, "getting gate_combn") # nolint

  # get cutpoints for prejoin gate combination method
  cp_uns_list_prejoin <- .get_cp_uns_loc_gate_combn_prejoin( # nolint
    gate_combn = gate_combn, ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig, ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min, ind_uns = ind_uns, ind_gate = ind_gate,
    min_bw = min_bw, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates, plot = plot,
    bias = bias, path_project = path_project, debug = debug
  )

  # get cutpoint if group method is not only prejoin
  cp_uns_list_prejoin_non <- .get_cp_uns_loc_gate_combn_prejoin_non(
    non_prejoin_combn = setdiff(gate_combn, "prejoin"),
    ex_list_no_min_stim = ex_list_no_min,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min, ind_uns = ind_uns, ind_gate = ind_gate,
    min_bw = min_bw, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates, plot = plot,
    bias = bias, path_project = path_project, debug = debug
  )

  # merge above two lists
  cp_uns_list <- .get_cp_uns_loc_gate_combn_merge(
    cp_uns_list_prejoin, cp_uns_list_prejoin_non, debug
  )

  list(
    "cp_uns" = list("loc" = cp_uns_list),
    "p_list" = list()
  )
}

.get_cp_uns_loc_gate_combn_merge <- function(cp_uns_list_prejoin,
                                             cp_uns_list_prejoin_non,
                                             debug) {
  .debug(debug, "done getting gate_combn") # nolint

  combined_list <- cp_uns_list_prejoin |>
    append(cp_uns_list_prejoin_non)
  purrr::map(
    unique(names(combined_list)),
    function(x) {
      .debug(debug, "cutpoint name", paste0(x, collapse = "-")) # nolint
      cp_uns_list_prejoin[[x]] |>
        append(cp_uns_list_prejoin_non[[x]])
    }
  ) |>
    stats::setNames(unique(names(combined_list)))
}

# --------------------------------
# gate using prejoined data
# --------------------------------
.get_cp_uns_loc_gate_combn_prejoin <- function(gate_combn,
                                               ex_list_no_min,
                                               ex_list_orig,
                                               ex_tbl_uns_bias,
                                               cp_min,
                                               ind_uns,
                                               ind_gate,
                                               min_bw,
                                               min_cell,
                                               gate_tbl,
                                               gate_name_curr,
                                               cut,
                                               calc_cyt_pos_gates,
                                               plot,
                                               bias,
                                               path_project,
                                               debug) {
  if (!"prejoin" %in% gate_combn) {
    return(.get_cp_uns_loc_gate_combn_prejoin_not())
  }
  .get_cp_uns_loc_gate_combn_prejoin_actual(
    ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    ind_uns = ind_uns,
    ind_gate = ind_gate,
    min_bw = min_bw,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    plot = plot,
    bias = bias,
    path_project = path_project,
    debug = debug
  )
}

.get_cp_uns_loc_gate_combn_prejoin_not <- function() {
  list("cp" = list(), "p_list" = list())
}

.get_cp_uns_loc_gate_combn_prejoin_actual <- function(ex_list_no_min,
                                                      ex_list_orig,
                                                      ex_tbl_uns_bias,
                                                      cp_min,
                                                      ind_uns,
                                                      ind_gate,
                                                      min_bw,
                                                      min_cell,
                                                      gate_tbl,
                                                      gate_name_curr,
                                                      cut,
                                                      calc_cyt_pos_gates,
                                                      plot,
                                                      bias,
                                                      path_project,
                                                      debug) {
  .debug(debug, "prejoin") # nolint

  # get marker expression for stim samples,
  # join and then sort into descending order
  ex_list_no_min_stim <- ex_list_no_min[names(ex_list_no_min) != ind_uns] |>
    purrr::map(function(x) x |> dplyr::arrange(dplyr::desc(expr))) # nolint

  # get cutpoints for gate combn method for a range of fdr's
  .get_cp_uns_loc_sample(
    ex_list_orig = ex_list_orig,
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    ind_uns = ind_uns,
    ind_gate = ind_gate,
    min_bw = min_bw,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    plot = plot,
    bias = bias,
    path_project = path_project
  ) |>
    purrr::map(function(x) list("prejoin" = x))
}

# --------------------------------
# gate each sample individually
# --------------------------------
.get_cp_uns_loc_gate_combn_prejoin_non <- function(non_prejoin_combn,
                                                   ex_list_no_min_stim,
                                                   ex_list_orig,
                                                   ex_tbl_uns_bias,
                                                   cp_min,
                                                   ind_uns,
                                                   ind_gate,
                                                   min_bw,
                                                   min_cell,
                                                   gate_tbl,
                                                   gate_name_curr,
                                                   cut,
                                                   calc_cyt_pos_gates,
                                                   plot,
                                                   bias,
                                                   path_project,
                                                   debug) {
  if (length(non_prejoin_combn) == 0L) {
    return(
      .get_cp_uns_loc_gate_combn_prejoin_non_not()
    )
  }
  .get_cp_uns_loc_gate_combn_prejoin_non_actual(
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    ind_uns = ind_uns,
    ind_gate = ind_gate,
    min_bw = min_bw,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    plot = plot,
    bias = bias,
    non_prejoin_combn = non_prejoin_combn,
    path_project = path_project,
    debug = debug
  )
}
.get_cp_uns_loc_gate_combn_prejoin_non_not <- function() {
  list("cp" = list(), "p_list" = list())
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual <- function(ex_list_no_min_stim,
                                                          ex_list_orig,
                                                          ex_tbl_uns_bias,
                                                          cp_min,
                                                          ind_uns,
                                                          ind_gate,
                                                          min_bw,
                                                          min_cell,
                                                          gate_tbl,
                                                          gate_name_curr,
                                                          cut,
                                                          calc_cyt_pos_gates,
                                                          plot,
                                                          bias,
                                                          path_project,
                                                          debug,
                                                          non_prejoin_combn) {
  .debug(debug, "non-prejoin") # nolint
  cp_uns_list_nonjoin <- .get_cp_uns_loc_sample(
    ex_list_orig =
      .get_cp_uns_loc_gate_combn_prejoin_non_actual_prep(
        ex_list_no_min_stim
      ),
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min, ind_uns = ind_uns, ind_gate = ind_gate,
    min_bw = min_bw, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, cut = cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates, plot = plot,
    bias = bias, path_project = path_project
  )
  p_list <-
    .get_cp_uns_loc_gate_combn_prejoin_non_actual_plot_org(
      cp_uns_list_nonjoin, debug, ex_list_orig, path_project
    )
  cp_uns_list_nonjoin <- .get_cp_uns_loc_gate_combn_prejoin_non_actual_combn(
    debug, cp_uns_list_nonjoin, non_prejoin_combn
  )
  list("cp" = cp_uns_list_nonjoin, "p_list" = p_list)
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual_prep <- function(ex_list_no_min_stim) { # nolint
  # get marker expression for stim samples,
  # removing minimum values
  # and sorting into descending order
  ex_list_no_min_stim |>
    purrr::map(function(x) x |> dplyr::arrange(desc(expr))) # nolint
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual_plot_org <- function(cp_uns_list_nonjoin, # nolint
                                                                   debug,
                                                                   ex_list_orig,
                                                                   path_project) { # nolint
  # get list of plots organised
  # ---------------------------
  .debug(debug, "Organising plots") # nolint

  indices_name_vec <- paste0(
    names(cp_uns_list_nonjoin[["loc"]]),
    collapse = "_"
  )
  sample_name_vec <- purrr::map_chr(
    names(cp_uns_list_nonjoin[["p_list"]]),
    function(ind) {
      paste0(
        ex_list_orig[[ind]]$batch_sh[1],
        "_", ex_list_orig[[ind]]$stim[1]
      )
    }
  )

  p_list_sample_level <- stats::setNames(
    cp_uns_list_nonjoin[["p_list"]],
    sample_name_vec
  )

  stats::setNames(
    list(p_list_sample_level), indices_name_vec
  )
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual_combn <- function(debug,
                                                                cp_uns_list_nonjoin, # nolint
                                                                non_prejoin_combn_vec) { # nolint
  # get list of cutpoints combined in the appropriate way
  # ---------------------------
  .debug(debug, "Combining cutpoints") # nolint
  .combine_cp( # nolint
    cp = cp_uns_list_nonjoin[["loc"]],
    gate_combn = non_prejoin_combn_vec
  )
}

# ------------------------------------------
# get cutpoints for a range of samples, and then individual samples
# ------------------------------------------

#' @title Get cutpoint for a range of samples given the q-value and fdr
#'
#' @description Calculate the cutpoint for each
#' sample in a batch at a given FDR.
#'
#' @param cut_stim list.
#' List where each element are the marker expression readings
#' of the marker to be cut on for the cells in a sample.
#' Note that the i-th element
#' in \code{cut_stim} must correspond to the i-th element in
#' \code{q_list}, i.e.
#' must be related to the same marker in same cell population
#' from the same blood sample and stimulation.
#' @param fdr numeric. A value between 0 and 1 specifying the
#' false discovery rate
#' the sample should be cut at.
#'
#' @return Numeric vector. A cutpoint for each sample.
.get_cp_uns_loc_sample <- function(ex_list_orig,
                                   ex_list_no_min_stim,
                                   ex_tbl_uns_bias,
                                   cp_min,
                                   min_bw,
                                   ind_uns,
                                   ind_gate,
                                   min_cell,
                                   gate_tbl,
                                   gate_name_curr,
                                   cut,
                                   calc_cyt_pos_gates,
                                   plot = TRUE,
                                   bias,
                                   path_project,
                                   debug = FALSE) {
  .debug(debug, "getting loc gate at sample level") # nolint

  # get cutpoints for each sample
  cp_uns_loc_obj_list <- purrr::map(seq_along(ex_list_no_min_stim), function(i) { # nolint
    .debug(debug, "sample", i) # nolint

    # return early if there are too few cells
    too_few_cells_lgl <- .get_cp_uns_loc_sample_check_cell_number(
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      min_cell = min_cell, ex_tbl_uns_bias = ex_tbl_uns_bias
    )
    if (too_few_cells_lgl) {
      return(.get_cp_uns_loc_sample_out_cell_number(debug))
    }

    # remove any cytokine-positive cells from unstim using gates from
    # sample for which single-positive gates are required
    ex_tbl_uns_bias <- .get_cp_uns_loc_sample_uns_rm_cyt_pos(
      debug = debug,
      ex_tbl_uns_orig = ex_list_orig[[as.character(ind_uns)]],
      gate_tbl = gate_tbl,
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      gate_name = gate_name_curr,
      cut = cut,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      bias = bias,
      ex_tbl_uns_bias = ex_tbl_uns_bias
    )

    cp_uns_loc_obj <- .get_cp_uns_loc_ind( # nolint
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      ex_tbl_stim_orig = ex_list_orig[[as.character(ind_gate)[i]]],
      ex_tbl_uns_orig = ex_list_orig[[as.character(ind_uns)]],
      min_bw = min_bw,
      cp_min = cp_min,
      min_cell = min_cell,
      plot = plot,
      bias = bias,
      path_project = path_project,
      debug = debug
    )
  })

  # name sample
  # ------------------
  .get_cp_uns_loc_output(
    debug = debug, cp_uns_loc_obj_list = cp_uns_loc_obj_list,
    ind_gate = ind_gate, ind_uns = ind_uns,
    ind_stim = seq_along(ex_list_no_min_stim)
  )
}



.get_cp_uns_loc_sample_check_cell_number <- function(ex_tbl_stim_no_min,
                                                     min_cell,
                                                     ex_tbl_uns_bias) {
  nrow(ex_tbl_stim_no_min) < min_cell ||
    nrow(ex_tbl_uns_bias) < min_cell
}

.get_cp_uns_loc_sample_out_cell_number <- function(debug) {
  .debug(debug, "Too few cells") # nolint
  p_list <- .get_cp_uns_loc_p_list_empty()
  list(p = NA, p_list = p_list)
}

.get_cp_uns_loc_sample_uns_rm_cyt_pos <- function(debug,
                                                  ex_tbl_uns_orig,
                                                  gate_tbl,
                                                  ex_tbl_stim_no_min,
                                                  gate_name,
                                                  cut,
                                                  calc_cyt_pos_gates,
                                                  bias,
                                                  ex_tbl_uns_bias) {
  if (is.null(gate_tbl)) {
    return(ex_tbl_uns_bias)
  }
  .debug(debug, "Removing cytokine-positive cells from unstim") # nolint

  # first filter
  gate_tbl_gn_ind <- gate_tbl |>
    # filter to use gates from cut_stim
    dplyr::filter(
      ind == ex_tbl_stim_no_min$ind[1], # nolint
      .data$gate_name == .env$gate_name # nolint
    )

  pos_ind_vec_but_single_pos_curr <-
    .get_pos_ind_but_single_pos_for_one_cyt( # nolint
      ex = ex_tbl_uns_orig,
      gate_tbl = gate_tbl_gn_ind,
      chnl_single_exc = cut,
      chnl = NULL,
      gate_type_cyt_pos = ifelse(calc_cyt_pos_gates, "cyt", "base"),
      gate_type_single_pos = "base"
    )

  ex_tbl_uns_orig <- ex_tbl_uns_orig[
    !pos_ind_vec_but_single_pos_curr, ,
    drop = FALSE
  ]

  # re-apply bias, noise and exclude minimum after
  # removing cytokine-positive cells
  .get_cut_list( # nolint
    ex_list = stats::setNames(list(ex_tbl_uns_orig), ex_tbl_uns_orig$ind[1]),
    ind = ex_tbl_uns_orig$ind[1],
    exc_min = TRUE,
    bias = bias,
    noise_sd = NULL
  )[[1]]
}

.get_cp_uns_loc_p_list_empty <- function() {
  lapply(1:3, function(x) ggplot2::ggplot()) |>
    stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))
}

.get_cp_uns_loc_ind_orig <- function(ex_tbl_stim_no_min, ex_tbl_uns_bias) {
  max_dens_x <- .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min)
  list(stim = ex_tbl_stim_no_min, uns = ex_tbl_uns_bias, max_x = max_dens_x)
}
.get_cp_uns_loc_ind <- function(ex_tbl_uns_bias, # was cut_unstim
                                ex_tbl_stim_no_min, # was cut_stim
                                min_bw,
                                cp_min,
                                min_cell,
                                # replacement of params
                                ex_tbl_stim_orig,
                                ex_tbl_uns_orig,
                                plot = TRUE,
                                prob_min = 0.1,
                                bias,
                                path_project,
                                debug = FALSE) {
  .debug(debug, "getting loc gate for single sample") # nolint

  # estimate densities for stim and unstim over stim range
  if (.get_cp_uns_loc_check_early(ex_tbl_stim_no_min, min_cell, cp_min)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_no_min, debug, "Too few cells"
    ))
  }

  # stop expr being higher than max_x to prevent
  # really far away values creating modes
  ex_tbl_stim_threshold <- .get_cp_uns_loc_set_max_expr(
    ex_tbl_stim_no_min, .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min)
  )
  ex_tbl_uns_threshold <- .get_cp_uns_loc_set_max_expr(
    ex_tbl_uns_bias, .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min)
  )

  # get smoothed probabilities
  data_mod <- .get_cp_uns_loc_get_prob(
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_stim_threshold = ex_tbl_stim_threshold,
    ex_tbl_uns_threshold = ex_tbl_uns_threshold,
    debug = debug, min_bw = min_bw, cp_min = cp_min,
    ex_tbl_uns_orig = ex_tbl_uns_orig, ex_tbl_stim_orig = ex_tbl_stim_orig
  )

  # get threshold
  cp <- .get_cp_uns_loc_get_cp(data_mod, ex_tbl_stim_no_min)

  list(cp = cp, p_list = .get_cp_uns_loc_p_list_empty())
}

# initial checks
# ---------------------
.get_cp_uns_loc_ind_check_n_cell <- function(ex_tbl_stim_no_min, min_cell) {
  nrow(ex_tbl_stim_no_min) < min_cell
}

.get_cp_uns_loc_ind_max_dens_x <- function(ex_tbl_stim_no_min) {
  max(ex_tbl_stim_no_min$expr) - 0.05 * (diff(range(ex_tbl_stim_no_min$expr)))
}

.get_cp_uns_loc_ind_check_max_x <- function(ex_tbl_stim_no_min, cp_min) {
  .get_cp_uns_loc_ind_max_dens_x(ex_tbl_stim_no_min) <= cp_min
}

.get_cp_uns_loc_check_early <- function(ex_tbl_stim_no_min, min_cell, cp_min) {
  .get_cp_uns_loc_ind_check_n_cell(ex_tbl_stim_no_min, min_cell) ||
    .get_cp_uns_loc_ind_check_max_x(ex_tbl_stim_no_min, cp_min)
}

.get_cp_uns_loc_ind_check_out <- function(cp_min,
                                          ex_tbl_stim_no_min,
                                          debug,
                                          msg) {
  .debug(debug, msg) # nolint
  list(
    cp = get_cp_uns_loc_ind_cp_non_loc(cp_min, ex_tbl_stim_no_min), # nolint
    p_list = .get_cp_uns_loc_p_list_empty()
  )
}

.get_cp_uns_loc_ind_cp_non_loc <- function(cp_min, ex_tbl_stim_no_min) {
  # get the threshold that is returned
  # automatically when we decide
  # we cannot apply this algorithm
  # (perhaps due to too few cells)
  max(
    cp_min,
    ex_tbl_stim_no_min$expr +
      (max(ex_tbl_stim_no_min$expr) -
        min(ex_tbl_stim_no_min$expr)) / 5
  )
}

.get_cp_uns_loc_set_max_expr <- function(.data, max_x) {
  .data |> dplyr::mutate(expr = pmin(.data$expr, max_x))
}

# get dens_tbl_raw
# -----------------------
.get_cp_uns_loc_get_dens_raw <- function(ex_tbl_stim_no_min,
                                         ex_tbl_uns_threshold,
                                         debug,
                                         min_bw) {
  .debug(debug, "Calculating densities") # nolint

  dens_list <- .get_cp_uns_loc_get_dens_raw_densities(
    ex_tbl_stim_no_min, ex_tbl_uns_threshold, debug, min_bw
  )

  # put raw densities into table
  get_cp_uns_loc_get_dens_raw_tabulate(dens_list = dens_list) # nolint
}

.get_cp_uns_loc_get_dens_raw_densities <- function(ex_tbl_stim_no_min,
                                                   ex_tbl_uns_threshold,
                                                   debug,
                                                   min_bw) {
  bw <- .get_cp_uns_loc_get_dens_raw_densities_bw(
    ex_tbl_stim_no_min, ex_tbl_uns_threshold, min_bw
  )
  dens_stim <- .get_cp_uns_loc_get_dens_raw_densities_stim(
    ex_tbl_stim_no_min, bw
  )
  dens_uns <- .get_cp_uns_loc_get_dens_raw_densities_uns(
    ex_tbl_uns_threshold, dens_stim, bw
  )
  list(stim = dens_stim, uns = dens_uns, bw = bw)
}

.get_cp_uns_loc_get_dens_raw_densities_bw <- function(ex_tbl_stim_no_min,
                                                      ex_tbl_uns_threshold,
                                                      min_bw) {
  bw_stim <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    ex_tbl_stim_no_min$expr, min_bw
  )
  bw_uns <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    ex_tbl_uns_threshold$expr, min_bw
  )
  max(bw_stim, min_bw, bw_uns)
}

.get_cp_uns_loc_get_dens_raw_densities_bw_init <- function(.data, min_bw) {
  bw_calc <- try(density(.data, bw = "SJ")$bw)
  if (inherits(bw_calc, "try-error")) min_bw else bw_calc
}

.get_cp_uns_loc_get_dens_raw_densities_stim <- function(ex_tbl_stim_no_min,
                                                        bw) {
  density(ex_tbl_stim_no_min$expr, bw = bw)
}
.get_cp_uns_loc_get_dens_raw_densities_uns <- function(ex_tbl_uns_threshold,
                                                       dens_stim,
                                                       bw) {
  density(
    ex_tbl_uns_threshold$expr,
    from = min(dens_stim$x), to = max(dens_stim$x), bw = bw
  )
}

.get_cp_uns_loc_get_dens_raw_tabulate <- function(stim_x,
                                                  stim_y,
                                                  uns_x,
                                                  uns_y) {
  dens_tbl_raw_stim <- tibble::tibble(x_stim = stim_x, y_stim = stim_y)
  dens_tbl_raw_wide <- .get_cp_uns_loc_get_dens_raw_tabulate_uns_interp(
    dens_tbl_raw_stim, uns_x, uns_y
  )
  .get_cp_uns_loc_get_dens_raw_tabulate_format(dens_tbl_raw_wide)
}

.get_cp_uns_loc_get_dens_raw_tabulate_uns_interp <- function(.data,
                                                             uns_x,
                                                             uns_y) {
  .data |>
    dplyr::mutate(y_uns = purrr::map_dbl(
      x_stim, # nolint
      function(marker) .interp(val = marker, x = uns_x, y = uns_y) # nolint
    ))
}

.get_cp_uns_loc_get_dens_raw_tabulate_format <- function(.data) {
  .data |>
    tidyr::pivot_longer(y_stim:y_uns, names_to = "stim", values_to = "dens") |> # nolint
    dplyr::mutate(stim = ifelse(stim == "y_stim", "yes", "no")) # nolint
}

#
.get_cp_uns_loc_get_prob_tbl <- function(dens_tbl_raw,
                                         debug,
                                         cp_min,
                                         ex_stim) {
  .debug(debug, "Normalising probabilities") # nolint

  # calculate raw and
  # normed probability based on densities for density measurements
  prob_tbl <- .get_cp_uns_loc_get_prob_tbl_init(dens_tbl_raw, cp_min) # nolint

  # choose probabilities to right of largest peak and
  # with sufficient evidence of a response ito probabilities
  # to be worth taking time smoothing over
  prob_tbl_pos <- .get_cp_uns_loc_prob_tbl_filter(ex_stim, prob_tbl, debug)
  list(all = prob_tbl, pos = prob_tbl_pos)
}

.get_cp_uns_loc_get_prob_tbl_init <- function(dens_tbl_raw, cp_min) {
  dens_tbl_raw |>
    tidyr::pivot_wider(
      id_cols = x_stim, names_from = stim, values_from = dens # nolint
    ) |>
    dplyr::mutate(
      prob_stim = 1 - no / yes, # nolint
      prob_stim = ifelse(yes == 0 & no == 0, 0, prob_stim), # nolint
      prob_stim_norm = pmin(1, prob_stim), # nolint
      prob_stim_norm = pmax(0, prob_stim_norm) # nolint
    ) |>
    dplyr::filter(x_stim > cp_min)
}

.get_cp_uns_loc_prob_tbl_filter <- function(ex_stim,
                                            prob_tbl,
                                            debug) {
  .debug(debug, "Filtering before smoothing") # nolint
  # I don't know why this is ex_stim orig, which
  # excludes the minimum. Shouldn't it be
  # ex_stim? (2024 June 14)
  # changing it to ex_stim_orig.

  # find highest peak (assumed to be the left-most one)
  density_exc_min <- density(ex_stim)
  dens_tbl <- tibble::tibble(x = density_exc_min$x, y = density_exc_min$y)
  peak <- dens_tbl |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint

  prob_tbl <- prob_tbl |>
    dplyr::filter(x_stim > peak + 0.02 * diff(range(x_stim))) # nolint

  # get range of values for which we'd want to calculate
  # probability:
  # - those cells for which the probability
  # of responding from their position onwards
  # is 0.025 or more
  prob_tbl |>
    dplyr::mutate(
      ge10 = prob_stim_norm >= 0.025, ge10 = cumsum(ge10) > 0 # nolint
    ) |>
    dplyr::filter(ge10) |>
    dplyr::select(-ge10)
}

.get_cp_uns_loc_get_min_prob_x <- function(prob_tbl_pos) {
  min(prob_tbl_pos$x_stim)
}

.get_cp_uns_loc_check_response <- function(prob_tbl_pos, ex_stim_orig) {
  nrow(prob_tbl_pos) == 0 ||
    max(ex_stim_orig) < .get_cp_uns_loc_get_min_prob_x(prob_tbl_pos)
}

.get_cp_uns_loc_get_data_mod <- function(ex_tbl_stim_no_min,
                                         ex_tbl_stim_orig,
                                         ex_tbl_uns_orig,
                                         prob_tbl_list) {
  if (.get_cp_uns_loc_check_response(prob_tbl_list$pos, ex_tbl_stim_orig)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_orig, debug, "No responding cells" # nolint
    ))
  }
  margin <- get_cp_uns_loc_get_data_mod_margin(
    ex_tbl_stim_orig, ex_tbl_uns_orig
  )

  ex_tbl_stim_no_min |>
    dplyr::filter(expr > # nolint
      (min(.get_cp_uns_loc_get_min_prob_x(prob_tbl_list$pos) - margin))) |>
    dplyr:::mutate(prob_smooth = expr)
}

get_cp_uns_loc_get_data_mod_margin <- function(ex_tbl_stim_orig,
                                               ex_tbl_uns_orig) {
  abs(max(diff(ex_tbl_stim_orig$expr), diff(ex_tbl_uns_orig$expr))) * 0.05
}

# smooth
# ---------------------
.get_cp_uns_loc_get_prob_smooth <- function(data_mod) {
  # enough cells to bother smoothing
  if (!.get_cp_uns_loc_get_prob_smooth_check_n_cell(data_mod)) {
    return(.get_cp_uns_loc_get_prob_smooth_check_n_cell_out(data_mod))
  }

  # get predictions after smoothing
  pred_vec <- .get_cp_uns_loc_get_prob_smooth_actual(data_mod, debug)
  data_mod |> dplyr::mutate(pred = pred_vec)
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell <- function(data_mod) {
  nrow(data_mod) < 10
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell_out <- function(data_mod) {
  data_mod |> dplyr::mutate(pred = prob_smooth - 1e-4) # nolint
}

.get_cp_uns_loc_get_prob_smooth_actual <- function(data_mod, debug) {
  fit_1 <- .get_cp_uns_loc_get_prob_smooth_actual_first(data_mod, debug)
  .get_cp_uns_loc_get_prob_smooth_actual_first_response(
    fit_1, data_mod, debug
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first <- function(data_mod, debug) {
  .debug(debug, "Smoothing I") # nolint
  try(
    scam::scam(
      prob_smooth ~ s(expr, bs = "mpi"), # nolint
      family = "binomial",
      data = data_mod |>
        dplyr::mutate(
          prob_smooth = pmin(prob_smooth, 0.999),
          prob_smooth = pmax(prob_smooth, 0.001)
        ),
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
      )
    ),
    silent = TRUE
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first_response <- function(fit,
                                                                  data_mod,
                                                                  debug) {
  # return predictions if success
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit, data_mod)) {
    .debug(debug, "Smoothed") # nolint
    return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit, data_mod
    )$pred)
  }
  # fit again if not a success
  .get_cp_uns_loc_get_prob_smooth_actual_first_response_failure(
    debug, data_mod
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_check <- function(fit, data_mod) {
  if (inherits(fit, "try-error")) {
    return(FALSE)
  }
  out_list <- .get_cp_uns_loc_get_prob_smooth_actual_response_success(
    fit, data_mod
  )
  !(all(out_list$pred > 0.99) || out_list$mean_abs_error > 0.3) # nolint
}

.get_cp_uns_loc_get_prob_smooth_actual_response_success <- function(fit, # nolint
                                                                    data_mod) { # nolint
  pred_vec <- predict(fit, type = "response")
  mean_abs_error <- mean(abs(pred_vec - data_mod$prob_smooth))
  list("pred" = pred_vec, "mean_abs_error" = mean_abs_error)
}

.get_cp_uns_loc_get_prob_smooth_actual_first_response_failure <- function(debug, # nolint
                                                                          data_mod) { # nolint
  fit_2 <- .get_cp_uns_loc_get_prob_smooth_actual_second(data_mod, debug)
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit_2, data_mod)) {
    .debug(debug, "Smoothed") # nolint
    return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit_2, data_mod
    )$pred)
  }
  .get_cp_uns_loc_get_prob_smooth_actual_third(data_mod, debug)
}

.get_cp_uns_loc_get_prob_smooth_actual_second <- function(data_mod, debug) {
  .debug(debug, "Smoothing II") # nolint
  try(
    scam::scam(
      prob_smooth ~ s(expr, bs = "micv"),
      family = "binomial",
      data = data_mod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.01
      )
    ),
    silent = TRUE
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_third <- function(data_mod) {
  .debug(debug, "Failed to smooth") # nolint
  data_mod$prob_smooth - 0.0001
}

# get probabilities
.get_cp_uns_loc_get_prob <- function(ex_tbl_stim_no_min,
                                     ex_tbl_stim_threshod,
                                     ex_tbl_uns_threshold,
                                     debug,
                                     min_bw,
                                     cp_min,
                                     ex_tbl_uns_orig,
                                     ex_tbl_stim_orig) {
  # get raw densities
  dens_tbl_raw <- .get_cp_uns_loc_get_dens_raw(
    ex_tbl_stim_no_min, ex_tbl_uns_orig, debug, min_bw
  )

  # get probabilities
  prob_tbl_list <- .get_cp_uns_loc_get_prob_tbl(
    dens_tbl_raw, debug, cp_min, ex_tbl_stim_no_min$expr
  )

  # get data to smooth over
  data_mod <- .get_cp_uns_loc_get_data_mod(
    ex_tbl_stim_no_min, ex_tbl_stim_orig, ex_tbl_uns_orig, prob_tbl_list
  )

  # smooth
  .get_cp_uns_loc_get_prob_smooth(data_mod)
}

# get cp
.get_cp_uns_loc_get_cp <- function(data_mod,
                                   ex_tbl_stim_orig,
                                   debug) {
  if (!inherits(data_mod, "data.frame")) {
    return(data_mod)
  }

  data_threshold <- .get_cp_uns_loc_get_cp_data_threshold(
    data_mod, ex_tbl_stim_orig
  )
  cp <- .get_cp_uns_loc_get_cp_actual(data_threshold)
  .debug(debug, "Completed loc gate for single sample") # nolint
  cp
}

.get_cp_uns_loc_get_cp_data_threshold <- function(data_mod,
                                                  ex_tbl_stim_orig,
                                                  ex_tbl_stim_no_min,
                                                  ex_tbl_uns_bias,
                                                  ex_tbl_uns_no_min,
                                                  bias) {
  data_count <- .get_cp_uns_loc_get_cp_data_threshold_count(data_mod)
  prob_bs_est <- .get_cp_uns_loc_get_cp_data_threshold_prop_bs_est(
    data_count, ex_tbl_stim_orig
  )
  .get_cp_uns_loc_get_cp_data_threshold_actual(
    data_count, prob_bs_est, ex_tbl_stim_orig,
    ex_tbl_stim_no_min, ex_tbl_uns_bias, ex_tbl_uns_no_min, bias # nolint
  )
}

.get_cp_uns_loc_get_cp_data_threshold_count <- function(data_mod) {
  data_mod |>
    dplyr::filter(expr >= min(.data$expr)) |> # nolint
    dplyr::arrange(expr) |>
    dplyr::mutate(n_row = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > prob_smooth) != n_row) |> # nolint
    dplyr::select(-n_row)
}

.get_cp_uns_loc_get_cp_data_threshold_prop_bs_est <- function(data_count,
                                                              ex_tbl_stim_no_min) { # nolint
  sum(data_count$pred) / nrow(ex_tbl_stim_no_min)
}

.get_cp_uns_loc_get_cp_data_threshold_actual <- function(data_count,
                                                         prop_bs_est,
                                                         ex_tbl_stim_orig,
                                                         ex_tbl_uns_bias,
                                                         ex_tbl_uns_orig,
                                                         bias) {
  data_count |>
    dplyr::arrange(desc(expr)) |> # nolint
    dplyr::mutate(count_stim = seq_len(dplyr::n())) |>
    dplyr::mutate(
      prop_stim = count_stim / nrow(ex_tbl_stim_orig), # nolint
      prop_uns = purrr::map_dbl(expr, function(x) {
        # work out proportion greater than the threshold,
        # adding back the bias that was removed
        # and dividing by the number of unstim cells
        # (before any were removed due to being cytokine-positive).
        # TODO: think about how to handle:
        # - bias
        # - cytokine-positive cells having been removed
        sum(ex_tbl_uns_orig$expr >= x) / nrow(ex_tbl_uns_orig)
      }),
      prop_bs = prop_stim - prop_uns, # nolint
      prop_bs_diff = prop_bs - prop_bs_est # nolint
    )
}

.get_cp_uns_loc_get_cp_actual <- function(data_threshold) {
  data_threshold |>
    dplyr::filter(abs(prop_bs_diff) == min(abs(prop_bs_diff))) |> # nolint
    dplyr::slice(1) |>
    dplyr::pull(expr) # nolint
}

.get_cp_uns_loc_output <- function(debug,
                                   cp_uns_loc_obj_list,
                                   ind_gate,
                                   ind_uns,
                                   ind_stim) {
  cp_vec <- .get_cp_uns_loc_sample_cp_rep(
    debug, cp_uns_loc_obj_list, ind_gate, ind_uns
  )
  p_list <- purrr::map(cp_uns_loc_obj_list, function(x) x$p_list) |>
    stats::setNames(ind_stim)
  .debug(debug, "done getting loc gate at sample level") # nolint
  # collate plots
  list(
    "loc" = cp_vec,
    "p_list" = p_list
  )
}

.get_cp_uns_loc_sample_cp_rep <- function(debug,
                                          cp_uns_loc_obj_list,
                                          ind_gate,
                                          ind_uns) {
  .debug(debug, "Possibly re-using calculated cutpoints") # nolint

  # extract vector of cutpoints
  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # stimulation indices
  ind_stim <- setdiff(ind_gate, ind_uns)

  cp_vec <- purrr::map_dbl(
    cp_uns_loc_obj_list,
    function(cp_uns_loc_obj) {
      cp_uns_loc_obj$cp
    }
  )

  # Repeat cutpoint if it was prejoined
  if (length(cp_vec) == 1 && ind_stim != 1) {
    cp_vec <- stats::setNames(
      rep(cp_vec, length(ind_stim)), ind_stim
    )
  } else {
    # name gate indices if not prejoined
    cp_vec <- stats::setNames(cp_vec, ind_stim)
  }

  # add unstim if unstim gate required
  if (ind_uns %in% ind_gate) {
    if (!all(purrr::map_lgl(cp_vec, is.na))) {
      cp_vec <- c(cp_vec, stats::setNames(mean(cp_vec, na.rm = TRUE), ind_uns))
    } else {
      cp_vec <- c(cp_vec, NA)
    }
  }

  cp_vec
}
