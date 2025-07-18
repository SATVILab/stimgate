# Calculate the local fdr-based cut
.get_cp_uns_loc <- function(ex_list, gate_combn,
                            .data, bias_uns = 0, exc_min,
                            noise_sd = NULL, bw_min = 80,
                            cp_min, min_cell, params,
                            path_project, .debug = FALSE) {
  # get cutpoints for each level of bias
  .get_cp_uns_loc_bias( # nolint
    ex_list = ex_list,
    .data = .data, bias_uns = bias_uns, exc_min = exc_min, noise_sd = noise_sd,
    cp_min = cp_min, gate_combn = gate_combn, bw_min = bw_min,
    min_cell = min_cell, gate_tbl = params$gate_tbl,
    gate_name_curr = params$gate_name_curr, chnl_cut = params$chnl_cut,
    calc_cyt_pos_gates = params$calc_cyt_pos_gates,
    path_project = path_project, .debug = .debug
  )
}



# Get the unstim-based local fdr-method cutpoint for each level of bias
.get_cp_uns_loc_bias <- function(ex_list, .data,
                                 bias_uns, exc_min, noise_sd, cp_min,
                                 gate_combn, bw_min, min_cell,
                                 gate_tbl, gate_name_curr, chnl_cut,
                                 calc_cyt_pos_gates, path_project,
                                 .debug) {
  # get ecdf of uns
  purrr::map(bias_uns, function(bias) {
    .debug_msg(.debug, "bias_uns", bias) # nolint

    ex_list_prep <- .prepare_data_with_bias_and_noise( # nolint
      ex_list = ex_list,
      bias = bias, noise_sd = noise_sd, exc_min = exc_min, .debug = .debug
    )

    # get gates for given level of bias across gate combination methods
    # --------------------------------------
    cp_uns_gate_combn_obj <- .get_cp_uns_loc_gate_combn( # nolint
      ex_list_orig = ex_list_prep[["ex_list_orig"]],
      ex_list_no_min = ex_list_prep[["ex_list_no_min"]],
      ex_tbl_uns_bias = ex_list_prep[["ex_tbl_uns_bias"]],
      gate_combn = gate_combn,
      cp_min = cp_min,
      bw_min = bw_min,
      min_cell = min_cell,
      gate_tbl = gate_tbl,
      gate_name_curr = gate_name_curr,
      chnl_cut,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      bias = bias,
      exc_min = exc_min,
      path_project = path_project,
      .debug = .debug
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
    paste0(names(cp_uns_gate_combn_list))
  cp_uns_gate_combn_list
}

.prepare_data_with_bias_and_noise <- function(ex_list,
                                              bias,
                                              exc_min,
                                              noise_sd,
                                              .debug) {
  # rename `cut` column` to `expr`
  # -------------------------------------
  ex_list_orig <- .prepare_ex_list_with_bias_and_noise( # nolint
    # exc_min should always be FALSE here, as we're trying
    # to keep all the original data
    ex_list = ex_list, ind = names(ex_list), exc_min = FALSE,
    bias = 0, noise_sd = NULL, .debug = .debug
  ) |>
    .arrange_samples_by_desc_expr()

  # separate stim and uns samples, rename
  # `cut` column to `expr` and exclude min
  # values
  # -------------------------------------

  ex_list_no_min <- .prepare_ex_list_with_bias_and_noise( # nolint
    # exc_min will be whatever the user wanted it to be.
    ex_list = ex_list, ind = names(ex_list), exc_min = exc_min,
    bias = 0, noise_sd = NULL, .debug = .debug
  ) |>
    .arrange_samples_by_desc_expr()

  # adjust expression in unstim,
  # applying bias, excluding the min val
  # and /or adding noise
  # -------------------------------------
  ex_list_bias <- .prepare_ex_list_with_bias_and_noise( # nolint
    ex_list = ex_list, ind = names(ex_list)[length(ex_list)], exc_min = exc_min,
    bias = bias, .debug = .debug, noise_sd = NULL
  ) |>
    .arrange_samples_by_desc_expr()
  ex_tbl_uns_bias <- ex_list_bias[[1]]

  list(
    "ex_list_orig" = ex_list_orig,
    "ex_list_no_min" = ex_list_no_min,
    "ex_tbl_uns_bias" = ex_tbl_uns_bias
  )
}


# Get the unstim-based local fdr-method
# cutpoint for a given bias across gate combination methods
.get_cp_uns_loc_gate_combn <- function(ex_list_orig,
                                       ex_list_no_min,
                                       ex_tbl_uns_bias,
                                       exc_min,
                                       gate_combn,
                                       cp_min,
                                       bw_min,
                                       min_cell,
                                       gate_tbl,
                                       gate_name_curr,
                                       chnl_cut,
                                       calc_cyt_pos_gates,
                                       bias,
                                       path_project,
                                       .debug = FALSE) {
  .debug_msg(.debug, "getting gate_combn") # nolint

  # get cutpoints for prejoin gate combination method
  cp_uns_list_prejoin <- .get_cp_uns_loc_gate_combn_prejoin( # nolint
    gate_combn = gate_combn, ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig, ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    bw_min = bw_min, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias, path_project = path_project, .debug = .debug
  )

  # get cutpoint for non-prejoin grouping methods
  cp_uns_list_prejoin_non <- .get_cp_uns_loc_gate_combn_prejoin_non(
    non_prejoin_combn = setdiff(gate_combn, "prejoin"),
    ex_list_no_min_stim = ex_list_no_min[-length(ex_list_no_min)],
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    bw_min = bw_min, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias, exc_min = exc_min, path_project = path_project, .debug = .debug
  )

  # merge above two lists
  cp_uns_list <- .get_cp_uns_loc_gate_combn_merge(
    cp_uns_list_prejoin, cp_uns_list_prejoin_non, .debug
  )

  list(
    "cp_uns" = list("loc" = cp_uns_list),
    "p_list" = list()
  )
}

.get_cp_uns_loc_gate_combn_merge <- function(cp_uns_list_prejoin,
                                             cp_uns_list_prejoin_non,
                                             .debug) {
  .debug_msg(.debug, "done getting gate_combn") # nolint

  combined_list <- cp_uns_list_prejoin |>
    append(cp_uns_list_prejoin_non)
  purrr::map(
    unique(names(combined_list)),
    function(x) {
      .debug_msg(.debug, "cutpoint name", paste0(x, collapse = "-")) # nolint
      cp_uns_list_prejoin[[x]] |>
        append(cp_uns_list_prejoin_non[[x]])
    }
  ) |>
    stats::setNames(unique(names(combined_list)))
}

# --------------------------------
# gate using prejoined .data
# --------------------------------
.get_cp_uns_loc_gate_combn_prejoin <- function(gate_combn,
                                               ex_list_no_min,
                                               ex_list_orig,
                                               ex_tbl_uns_bias,
                                               cp_min,
                                               bw_min,
                                               min_cell,
                                               gate_tbl,
                                               gate_name_curr,
                                               chnl_cut,
                                               calc_cyt_pos_gates,
                                               bias,
                                               path_project,
                                               .debug) {
  if (!"prejoin" %in% gate_combn) {
    return(.get_cp_uns_loc_gate_combn_prejoin_not())
  }
  .get_cp_uns_loc_gate_combn_prejoin_actual(
    ex_list_no_min = ex_list_no_min,
    ex_list_orig = ex_list_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    bw_min = bw_min,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias,
    path_project = path_project,
    .debug = .debug
  )
}

.get_cp_uns_loc_gate_combn_prejoin_not <- function() {
  list("cp" = list(), "p_list" = list())
}

.prepare_data_for_prejoin <- function(ex_list_orig,
                                      ex_list_no_min) {
  ex_list_no_min <- .prepare_data_for_prejoin_ind(
    ex_list = ex_list_no_min
  )
  ex_list_orig <- .prepare_data_for_prejoin_ind(
    ex_list = ex_list_orig
  )

  list(
    "ex_list_orig" = ex_list_orig,
    "ex_list_no_min" = ex_list_no_min
  )
}

.prepare_data_for_prejoin_ind <- function(ex_list) {
  ind_uns <- names(ex_list)[length(ex_list)]
  ex_tbl_stim <- ex_list[seq_len(length(ex_list) - 1)] |>
    dplyr::bind_rows()
  ex_tbl_stim <- ex_tbl_stim[order(.get_cut(ex_tbl_stim)), ]
  list(
    ex_tbl_stim,
    ex_list[[length(ex_list)]]
  ) |>
    stats::setNames(c(names(ex_list)[1], ind_uns))
}

.get_cp_uns_loc_gate_combn_prejoin_actual <- function(ex_list_no_min,
                                                      ex_list_orig,
                                                      ex_tbl_uns_bias,
                                                      cp_min,
                                                      bw_min,
                                                      min_cell,
                                                      gate_tbl,
                                                      gate_name_curr,
                                                      chnl_cut,
                                                      calc_cyt_pos_gates,
                                                      bias,
                                                      path_project,
                                                      .debug) {
  .debug_msg(.debug, "prejoin") # nolint

  # get marker expression for stim samples,
  # join and then sort into descending order
  ex_list_prejoin <- .prepare_data_for_prejoin(
    ex_list_orig = ex_list_orig, ex_list_no_min = ex_list_no_min
  )

  # get cutpoints for gate combn method for a range of fdr's
  .get_cp_uns_loc_sample(
    ex_list_orig = ex_list_prejoin[["ex_list_orig"]],
    ex_list_no_min_stim = ex_list_prejoin[["ex_list_no_min"]],
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    bw_min = bw_min,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias,
    exc_min = exc_min,
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
                                                   bw_min,
                                                   min_cell,
                                                   gate_tbl,
                                                   gate_name_curr,
                                                   chnl_cut,
                                                   calc_cyt_pos_gates,
                                                   bias,
                                                   exc_min,
                                                   path_project,
                                                   .debug) {
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
    bw_min = bw_min,
    min_cell = min_cell,
    gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr,
    chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias,
    exc_min = exc_min,
    non_prejoin_combn = non_prejoin_combn,
    path_project = path_project,
    .debug = .debug
  )
}
.get_cp_uns_loc_gate_combn_prejoin_non_not <- function() {
  list("cp" = list(), "p_list" = list())
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual <- function(ex_list_no_min_stim,
                                                          ex_list_orig,
                                                          ex_tbl_uns_bias,
                                                          cp_min,
                                                          bw_min,
                                                          min_cell,
                                                          gate_tbl,
                                                          gate_name_curr,
                                                          chnl_cut,
                                                          calc_cyt_pos_gates,
                                                          bias,
                                                          exc_min,
                                                          path_project,
                                                          .debug,
                                                          non_prejoin_combn) {
  .debug_msg(.debug, "non-prejoin") # nolint
  cp_uns_list_nonjoin <- .get_cp_uns_loc_sample(
    ex_list_orig = ex_list_orig,
    ex_list_no_min_stim = ex_list_no_min_stim,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    cp_min = cp_min,
    bw_min = bw_min, min_cell = min_cell, gate_tbl = gate_tbl,
    gate_name_curr = gate_name_curr, chnl_cut,
    calc_cyt_pos_gates = calc_cyt_pos_gates,
    bias = bias, exc_min = exc_min, path_project = path_project
  )

  cp_uns_list_nonjoin <- .get_cp_uns_loc_gate_combn_prejoin_non_actual_combn(
    .debug, cp_uns_list_nonjoin, non_prejoin_combn
  )
  list("cp" = cp_uns_list_nonjoin, "p_list" = list())
}

.arrange_samples_by_desc_expr <- function(ex_list) {
  # arrange in descending order of expression
  ex_list |>
    purrr::map(function(x)  x[order(.get_cut(x)), ]) # nolint
}

.get_cp_uns_loc_gate_combn_prejoin_non_actual_combn <- function(.debug,
                                                                cp_uns_list_nonjoin, # nolint
                                                                non_prejoin_combn_vec) { # nolint
  # get list of cutpoints combined in the appropriate way
  # ---------------------------
  .debug_msg(.debug, "Combining cutpoints") # nolint
  .combine_cp( # nolint
    cp = cp_uns_list_nonjoin[["loc"]],
    gate_combn = non_prejoin_combn_vec
  )
}

# ------------------------------------------
# get cutpoints for a range of samples, and then individual samples
# ------------------------------------------

# Get cutpoint for a range of samples given the q-value and fdr
# 
# Calculate the cutpoint for each
# sample in a batch at a given FDR.
# 
# @param cut_stim list.
# List where each element are the marker expression readings
# of the marker to be cut on for the cells in a sample.
# Note that the i-th element
# in cut_stim must correspond to the i-th element in
# q_list, i.e.
# must be related to the same marker in same cell population
# from the same blood sample and stimulation.
# @param fdr numeric. A value between 0 and 1 specifying the
# false discovery rate
# the sample should be cut at.
# @return Numeric vector. A cutpoint for each sample.
.get_cp_uns_loc_sample <- function(ex_list_orig,
                                   ex_list_no_min_stim,
                                   ex_tbl_uns_bias,
                                   cp_min,
                                   bw_min,
                                   min_cell,
                                   gate_tbl,
                                   gate_name_curr,
                                   chnl_cut,
                                   calc_cyt_pos_gates,
                                   bias,
                                   exc_min,
                                   path_project,
                                   .debug = FALSE) {
  .debug_msg(.debug, "getting loc gate at sample level") # nolint

  # get cutpoints for each sample
  cp_uns_loc_obj_list <- purrr::map(seq_along(ex_list_no_min_stim), function(i) { # nolint
    .debug_msg(.debug, "sample", i) # nolint

    # return early if there are too few cells
    too_few_cells_lgl <- .get_cp_uns_loc_sample_check_cell_number(
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      min_cell = min_cell, ex_tbl_uns_bias = ex_tbl_uns_bias
    )
    if (too_few_cells_lgl) {
      return(.get_cp_uns_loc_ind_check_out(
          cp_min, ex_list_no_min_stim[[i]],
          ex_tbl_uns_bias, .debug, "Too few cells"
        ))
    }

    # remove any cytokine-positive cells from unstim using gates from
    # sample for which single-positive gates are required
    ex_tbl_uns_bias <- .get_cp_uns_loc_sample_uns_rm_cyt_pos(
      .debug = .debug,
      ex_tbl_uns_orig = ex_list_orig[[length(ex_list_orig)]],
      gate_tbl = gate_tbl,
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      gate_name = gate_name_curr,
      chnl_cut,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      bias = bias,
      exc_min = exc_min,
      ex_tbl_uns_bias = ex_tbl_uns_bias
    )

    .get_cp_uns_loc_ind( # nolint
      ex_tbl_stim_no_min = ex_list_no_min_stim[[i]],
      ex_tbl_uns_bias = ex_tbl_uns_bias,
      ex_tbl_stim_orig = ex_list_orig[[i]],
      ex_tbl_uns_orig = ex_list_orig[[length(ex_list_orig)]],
      bw_min = bw_min,
      cp_min = cp_min,
      min_cell = min_cell,
      bias = bias,
      path_project = path_project,
      .debug = .debug
    )
  }) |>
    stats::setNames(names(ex_list_no_min_stim))

  # name sample
  # ------------------
  .get_cp_uns_loc_output(
    .debug = .debug, cp_uns_loc_obj_list = cp_uns_loc_obj_list,
    ind_uns = names(ex_list_orig)[length(ex_list_orig)],
    ind_stim = names(ex_list_no_min_stim)
  )
}

.get_cp_uns_loc_sample_check_cell_number <- function(ex_tbl_stim_no_min,
                                                     min_cell,
                                                     ex_tbl_uns_bias) {
  nrow(ex_tbl_stim_no_min) < min_cell ||
    nrow(ex_tbl_uns_bias) < min_cell
}

.get_cp_uns_loc_sample_uns_rm_cyt_pos <- function(.debug,
                                                  ex_tbl_uns_orig,
                                                  gate_tbl,
                                                  ex_tbl_stim_no_min,
                                                  gate_name,
                                                  chnl_cut,
                                                  calc_cyt_pos_gates,
                                                  bias,
                                                  exc_min,
                                                  ex_tbl_uns_bias) {
  if (is.null(gate_tbl)) {
    return(ex_tbl_uns_bias)
  }
  .debug_msg(.debug, "Removing cytokine-positive cells from unstim") # nolint

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
      chnl_single_exc = chnl_cut,
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
  .prepare_ex_list_with_bias_and_noise( # nolint
    ex_list = stats::setNames(
      list(ex_tbl_uns_orig), attr(ex_tbl_uns_orig, "ind_uns")
      ),
    ind = attr(ex_tbl_uns_orig, "ind_uns"),
    exc_min = exc_min,
    bias = bias,
    noise_sd = NULL
  )[[1]]
}

.get_cp_uns_loc_p_list_empty <- function() {
  lapply(1:3, function(x) list()) |>
    stats::setNames(c("p_loc_dens", "p_loc_prob", "p_loc_ctb"))
}

.get_cp_uns_loc_ind <- function(ex_tbl_uns_bias, # was cut_unstim
                                ex_tbl_stim_no_min, # was cut_stim
                                bw_min,
                                cp_min,
                                min_cell,
                                # replacement of params
                                ex_tbl_stim_orig,
                                ex_tbl_uns_orig,
                                plot = TRUE,
                                prob_min = 0.1,
                                bias,
                                path_project,
                                .debug = FALSE) {
  .debug_msg(.debug, "getting loc gate for single sample") # nolint

  # estimate densities for stim and unstim over stim range
  if (.get_cp_uns_loc_check_early(ex_tbl_stim_no_min, min_cell, cp_min)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_no_min, ex_tbl_uns_bias, .debug, "Too few cells"
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
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    .debug = .debug, bw_min = bw_min, cp_min = cp_min + bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig
  )

  # get threshold
  .get_cp_uns_loc_get_cp(
    data_mod = data_mod,
    .debug = .debug,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    bias = bias,
    cp_min = cp_min
  )
}

# initial checks
# ---------------------
.get_cp_uns_loc_ind_check_n_cell <- function(ex_tbl_stim_no_min, min_cell) {
  nrow(ex_tbl_stim_no_min) < min_cell
}

.get_cp_uns_loc_ind_max_dens_x <- function(ex_tbl_stim_no_min) {
  max(.get_cut(ex_tbl_stim_no_min)) - 0.05 * (diff(range(.get_cut(ex_tbl_stim_no_min))))
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
                                          ex_tbl_uns_bias,
                                          .debug,
                                          msg) {
  .debug_msg(.debug, msg) # nolint
  list(
    cp = .get_cp_uns_loc_ind_cp_non_loc(
      cp_min, ex_tbl_stim_no_min, ex_tbl_uns_bias
      ), # nolint
    p_list = .get_cp_uns_loc_p_list_empty()
  )
}

.get_cp_uns_loc_ind_cp_non_loc <- function(cp_min,
                                           ex_tbl_stim_no_min,
                                           ex_tbl_uns_bias) {
  # get the threshold that is returned
  # automatically when we decide
  # we cannot apply this algorithm
  # (perhaps due to too few cells)
  range_vec_stim <- range(.get_cut(ex_tbl_stim_no_min))
  range_vec_uns <- range(.get_cut(ex_tbl_uns_bias))
  range_len <- max(diff(range_vec_stim), diff(range_vec_uns))
  max(
    cp_min,
    range_vec_stim[[2]] + range_len / 5,
    range_vec_uns[[2]] + range_len / 3
  )
}

.get_cp_uns_loc_set_max_expr <- function(.data, max_x) {
  .data[, attr(.data, "chnl_cut")] <- pmin(.get_cut(.data), max_x)
  .data
}

# get probabilities
.get_cp_uns_loc_get_prob <- function(ex_tbl_stim_no_min,
                                     ex_tbl_stim_threshold,
                                     ex_tbl_uns_threshold,
                                     ex_tbl_uns_bias,
                                     .debug,
                                     bw_min,
                                     cp_min,
                                     ex_tbl_uns_orig) {
  # get raw densities
  dens_tbl_raw <- .get_cp_uns_loc_get_dens_raw(
    ex_tbl_stim_threshold, ex_tbl_uns_threshold, .debug, bw_min
  )

  # get probabilities
  prob_tbl_list <- .get_cp_uns_loc_get_prob_tbl(
    dens_tbl_raw, .debug, cp_min, .get_cut(ex_tbl_stim_threshold),
    .get_cut(ex_tbl_uns_threshold)
  )

  # get .data to smooth over
  data_mod <- .get_cp_uns_loc_get_data_mod(
    ex_tbl_stim_threshold, ex_tbl_stim_no_min, ex_tbl_uns_threshold,
    ex_tbl_uns_bias, prob_tbl_list, cp_min, .debug
  )

  # smooth
  .get_cp_uns_loc_get_prob_smooth(data_mod, .debug)
}

# get dens_tbl_raw
# -----------------------
.get_cp_uns_loc_get_dens_raw <- function(ex_tbl_stim_threshold,
                                         ex_tbl_uns_threshold,
                                         .debug,
                                         bw_min) {
  .debug_msg(.debug, "Calculating densities") # nolint

  dens_list <- .get_cp_uns_loc_get_dens_raw_densities(
    ex_tbl_stim_threshold, ex_tbl_uns_threshold, .debug, bw_min
  )

  # put raw densities into table
  .get_cp_uns_loc_get_dens_raw_tabulate(
    stim_x = dens_list$stim$x,
    stim_y = dens_list$stim$y,
    uns_x = dens_list$uns$x,
    uns_y = dens_list$uns$y
  ) # nolint
}

.get_cp_uns_loc_get_dens_raw_densities <- function(ex_tbl_stim_threshold,
                                                   ex_tbl_uns_threshold,
                                                   .debug,
                                                   bw_min) {
  bw <- .get_cp_uns_loc_get_dens_raw_densities_bw(
    ex_tbl_stim_threshold, ex_tbl_uns_threshold, bw_min
  )
  dens_stim <- .get_cp_uns_loc_get_dens_raw_densities_stim(
    ex_tbl_stim_threshold, bw
  )
  dens_uns <- .get_cp_uns_loc_get_dens_raw_densities_uns(
    ex_tbl_uns_threshold, dens_stim, bw
  )
  list(stim = dens_stim, uns = dens_uns, bw = bw)
}

.get_cp_uns_loc_get_dens_raw_densities_bw <- function(ex_tbl_stim_threshold,
                                                      ex_tbl_uns_threshold,
                                                      bw_min) {
  bw_stim <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    .get_cut(ex_tbl_stim_threshold), bw_min
  )
  bw_uns <- .get_cp_uns_loc_get_dens_raw_densities_bw_init(
    .get_cut(ex_tbl_uns_threshold), bw_min
  )
  bw_init <- max(bw_stim, bw_min, bw_uns)
  # don't allow bandwidth to be too large, as otherwise
  # we can push up the negative pop's right side
  # by smoothing the positive pop into it.
  min(bw_min * 1.5, bw_init)
}

.get_cp_uns_loc_get_dens_raw_densities_bw_init <- function(.data, bw_min) {
  bw_calc <- try(density(.data, bw = "SJ")$bw)
  if (inherits(bw_calc, "try-error")) bw_min else bw_calc
}

.get_cp_uns_loc_get_dens_raw_densities_stim <- function(ex_tbl_stim_threshold,
                                                        bw) {
  dens_obj <- density(.get_cut(ex_tbl_stim_threshold), bw = bw)
  if (is.null(attr(ex_tbl_stim_threshold, "prob_g_min"))) {
    return(dens_obj)
  }
  dens_obj$y <- dens_obj$y * attr(ex_tbl_stim_threshold, "prob_g_min")
  dens_obj
}
.get_cp_uns_loc_get_dens_raw_densities_uns <- function(ex_tbl_uns_threshold,
                                                       dens_stim,
                                                       bw) {
  dens_obj <- density(
    .get_cut(ex_tbl_uns_threshold),
    from = min(dens_stim$x), to = max(dens_stim$x), bw = bw
  )
  if (is.null(attr(ex_tbl_uns_threshold, "prob_g_min"))) {
    return(dens_obj)
  }
  dens_obj$y <- dens_obj$y * attr(ex_tbl_uns_threshold, "prob_g_min")
  dens_obj
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

# get probabilities
# -------------------


.get_cp_uns_loc_get_prob_tbl <- function(dens_tbl_raw,
                                         .debug,
                                         cp_min,
                                         ex_vec_stim_threshold,
                                         ex_vec_uns_threshold) {
  .debug_msg(.debug, "Normalising probabilities") # nolint

  # calculate raw and
  # normed probability based on densities for density measurements
  prob_tbl <- .get_cp_uns_loc_get_prob_tbl_init(dens_tbl_raw, cp_min) # nolint

  # choose probabilities to right of largest peak and
  # with sufficient evidence of a response ito probabilities
  # to be worth taking time smoothing over
  prob_tbl_pos <- .get_cp_uns_loc_prob_tbl_filter(
    ex_vec_stim_threshold, ex_vec_uns_threshold, prob_tbl, .debug
  )
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

.get_cp_uns_loc_prob_tbl_filter <- function(ex_vec_stim_threshold,
                                            ex_vec_uns_threshold,
                                            prob_tbl,
                                            .debug) {
  .debug_msg(.debug, "Filtering before smoothing") # nolint
  # I don't know why this is ex_stim orig, which
  # excludes the minimum. Shouldn't it be
  # ex_stim? (2024 June 14)
  # changing it to ex_tbl_stim_orig.

  # find highest peak (assumed to be the left-most one)
  density_exc_min_stim <- density(ex_vec_stim_threshold)
  dens_tbl_stim <- tibble::tibble(
    x = density_exc_min_stim$x, y = density_exc_min_stim$y
    )
  peak_stim <- dens_tbl_stim |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  density_exc_min_uns <- density(ex_vec_uns_threshold)
  dens_tbl_uns <- tibble::tibble(
    x = density_exc_min_uns$x, y = density_exc_min_uns$y
    )
  peak_uns <- dens_tbl_uns |>
    dplyr::filter(y == max(y)) |> # nolint
    dplyr::pull("x") # nolint
  peak <- max(peak_stim, peak_uns)

  window_width <- 0.15 * diff(quantile(prob_tbl$x_stim, c(0.05, 0.5)))

  # move to the right of the peak, to handle peak misalignment

  prob_tbl <- prob_tbl |>
    dplyr::filter(x_stim > peak + window_width) # nolint

  # filter out observations where no observation
  # to the left of it had a prob greater than 0.025
  prob_tbl |>
    dplyr::filter(
      cumsum(prob_stim_norm >= 0.025) > 0
    )

  # get range of values for which we'd want to calculate
  # probability:
  # - those cells for which the probability
  # of responding from their position onwards
  # is 0.025 or more, within the window to their right
  prob_tbl <- prob_tbl |>
    dplyr::mutate(
      minor_response_ind = prob_stim_norm >= 0.025,
      moderate_response_ind = prob_stim_norm >= 0.075,
      # ge10 = cumsum(minor_response_ind) > 0,
      n_remaining = dplyr::n() - seq_len(dplyr::n()) + 1
      # ge10 = (cumsum(ge10) / n_remaining) > 0.25 # nolint
    )
  prob_tbl |>
    dplyr::filter(
      # clear out all the near-negatives
      cumsum(minor_response_ind) > 0 # nolint
    ) |>
    # detect places of consistent non-response
    dplyr::mutate(
      prob_larger_count = purrr::map_int(x_stim, function(x) {
        sum(prob_tbl$moderate_response_ind[prob_tbl$x_stim >= x]) 
      }),
      prob_larger_prop = prob_larger_count / n_remaining # nolint
    ) |>
    dplyr::filter(
      prob_larger_prop > 0.25 # nolint
    ) |>
    dplyr::select(
      -c(
        prob_larger_prop, minor_response_ind, moderate_response_ind,
        n_remaining, prob_larger_count
        )
    )
  
}

.get_cp_uns_loc_get_min_prob_x <- function(prob_tbl_pos) {
  min(prob_tbl_pos$x_stim)
}

.get_cp_uns_loc_check_response <- function(prob_tbl_pos, ex_tbl_stim_orig) {
  nrow(prob_tbl_pos) == 0 ||
    max(.get_cut(ex_tbl_stim_orig)) < .get_cp_uns_loc_get_min_prob_x(prob_tbl_pos)
}

.get_cp_uns_loc_get_data_mod <- function(ex_tbl_stim_threshold,
                                         ex_tbl_stim_no_min,
                                         ex_tbl_uns_threshold,
                                         ex_tbl_uns_bias,
                                         prob_tbl_list,
                                         cp_min,
                                         .debug) {
  if (.get_cp_uns_loc_check_response(prob_tbl_list$pos, ex_tbl_stim_no_min)) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_no_min,
      ex_tbl_uns_bias, .debug, "No responding cells" # nolint
    ))
  }
  margin <- get_cp_uns_loc_get_data_mod_margin(
    ex_tbl_stim_no_min, ex_tbl_uns_threshold
  )

  data_mod <- ex_tbl_stim_threshold
  data_mod <- data_mod[
    .get_cut(data_mod) >=
      (min(.get_cp_uns_loc_get_min_prob_x(prob_tbl_list$pos) - margin)),
  ]
  if (nrow(data_mod) == 0L) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_no_min, ex_tbl_uns_bias,
      .debug, "No responding cells" # nolint
    ))
  }
  prob_vec <- approx(
      x = prob_tbl_list$pos$x_stim,
      y = prob_tbl_list$pos$prob_stim_norm,
      xout = data_mod[[1]],
      method = "linear",
      f = 0.5,    # midpoint rule → nearest
      rule = 2    # outside range: use end‐points
    )$y

  data_mod |>
    dplyr::mutate(prob_smooth = prob_vec)
}

get_cp_uns_loc_get_data_mod_margin <- function(ex_tbl_stim_no_min,
                                               ex_tbl_uns_no_min) {
  abs(max(diff(.get_cut(ex_tbl_stim_no_min)), diff(.get_cut(ex_tbl_uns_no_min)))) * 0.05
}

# smooth
# ---------------------
.get_cp_uns_loc_get_prob_smooth <- function(data_mod, .debug) {
  # enough cells to bother smoothing
  if (!.get_cp_uns_loc_get_prob_smooth_check_n_cell(data_mod)) {
    return(.get_cp_uns_loc_get_prob_smooth_check_n_cell_out(data_mod))
  }

  # get predictions after smoothing
  pred_vec <- .get_cp_uns_loc_get_prob_smooth_actual(data_mod, .debug)
  data_mod |> dplyr::mutate(pred = pred_vec)
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell <- function(data_mod) {
  # it's already the cutpoint, possibly, so check
  # if it's a dataframe first and then that
  # there are enough cells
  is.data.frame(data_mod) && nrow(data_mod) > 10
}

.get_cp_uns_loc_get_prob_smooth_check_n_cell_out <- function(data_mod) {
  if (is.data.frame(data_mod)) {
    data_mod |> dplyr::mutate(pred = prob_smooth - 1e-4) # nolint
  } else {
    # just return the list with the cutpoint
    data_mod
  }
}

.get_cp_uns_loc_get_prob_smooth_actual <- function(data_mod, .debug) {
  fit_1 <- .get_cp_uns_loc_get_prob_smooth_actual_first(data_mod, .debug)
  .get_cp_uns_loc_get_prob_smooth_actual_first_response(
    fit_1, data_mod, .debug
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first <- function(data_mod, .debug) {
  .debug_msg(.debug, "Smoothing I") # nolint
  try({
    fml <- as.formula(paste0(
      "prob_smooth ~ s(x_stim), bs = 'mpi')"
    ))
    scam::scam(
      fml, # nolint
      family = "binomial",
      .data = data_mod |>
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
    )
  },
  silent = TRUE
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_first_response <- function(fit,
                                                                  data_mod,
                                                                  .debug) {
  # return predictions if success
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit, data_mod)) {
    .debug_msg(.debug, "Smoothed") # nolint
    return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit, data_mod
    )$pred)
  }
  # fit again if not a success
  .get_cp_uns_loc_get_prob_smooth_actual_first_response_failure(
    .debug, data_mod
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

.get_cp_uns_loc_get_prob_smooth_actual_first_response_failure <- function(.debug, # nolint
                                                                          data_mod) { # nolint
  fit_2 <- .get_cp_uns_loc_get_prob_smooth_actual_second(data_mod, .debug)
  if (.get_cp_uns_loc_get_prob_smooth_actual_check(fit_2, data_mod)) {
    .debug_msg(.debug, "Smoothed") # nolint
    return(.get_cp_uns_loc_get_prob_smooth_actual_response_success(
      fit_2, data_mod
    )$pred)
  }
  .get_cp_uns_loc_get_prob_smooth_actual_third(data_mod, .debug)
}

.get_cp_uns_loc_get_prob_smooth_actual_second <- function(data_mod, .debug) {
  .debug_msg(.debug, "Smoothing II") # nolint
  try({
    fml <- as.formula(paste0(
      "prob_smooth ~ s(x_stim), bs = 'micv')"
    ))
    scam::scam(
      fml,
      family = "binomial",
      .data = data_mod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.01
      )
    )},
    silent = TRUE
  )
}

.get_cp_uns_loc_get_prob_smooth_actual_third <- function(data_mod, .debug) {
  .debug_msg(.debug, "Failed to smooth") # nolint
  data_mod$prob_smooth - 0.0001
}






# get cp
.get_cp_uns_loc_get_cp <- function(data_mod,
                                   ex_tbl_stim_orig,
                                   .debug,
                                   ex_tbl_stim_no_min,
                                   ex_tbl_uns_orig,
                                   ex_tbl_uns_bias,
                                   bias,
                                   cp_min) {
  if (!is.data.frame(data_mod)) {
    return(data_mod)
  }

  data_threshold <- .get_cp_uns_loc_get_cp_data_threshold(
    data_mod = data_mod, ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_stim_no_min = ex_tbl_stim_no_min,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    bias = bias
  )
  cp <- .get_cp_uns_loc_get_cp_actual(
    data_threshold, ex_tbl_stim_no_min, ex_tbl_uns_bias, cp_min, .debug
  )
  .debug_msg(.debug, "Completed loc gate for single sample") # nolint
  list("cp" = cp, "p_list" = .get_cp_uns_loc_p_list_empty())
}

.get_cp_uns_loc_get_cp_data_threshold <- function(data_mod,
                                                  ex_tbl_stim_orig,
                                                  ex_tbl_stim_no_min,
                                                  ex_tbl_uns_bias,
                                                  ex_tbl_uns_orig,
                                                  bias) {
  data_count <- .get_cp_uns_loc_get_cp_data_threshold_count(data_mod)
  prob_bs_est <- .get_cp_uns_loc_get_cp_data_threshold_prop_bs_est(
    data_count = data_count, ex_tbl_stim_orig = ex_tbl_stim_orig
  )
  .get_cp_uns_loc_get_cp_data_threshold_actual(
    data_count = data_count,
    prop_bs_est = prob_bs_est,
    ex_tbl_stim_orig = ex_tbl_stim_orig,
    ex_tbl_uns_bias = ex_tbl_uns_bias,
    ex_tbl_uns_orig = ex_tbl_uns_orig,
    bias = bias
  )
}

.get_cp_uns_loc_get_cp_data_threshold_count <- function(data_mod) {
  if (nrow(data_mod) == 1L) {
    min_val <- min(.get_cut(data_mod)) - 1
  } else {
    min_val <- min(.get_cut(data_mod))
  }
  data_mod <- data_mod[.get_cut(data_mod) > min_val, ]
  data_mod <- data_mod[order(.get_cut(data_mod)), ]
  data_mod |>
    dplyr::mutate(n_row = seq_len(dplyr::n())) |>
    dplyr::filter(cumsum(pred > prob_smooth) != n_row) |> # nolint
    dplyr::select(-n_row)
}

.get_cp_uns_loc_get_cp_data_threshold_prop_bs_est <- function(data_count,
                                                              ex_tbl_stim_orig) { # nolint
  sum(data_count$pred) / nrow(ex_tbl_stim_orig)
}

.get_cp_uns_loc_get_cp_data_threshold_actual <- function(data_count,
                                                         prop_bs_est,
                                                         ex_tbl_stim_orig,
                                                         ex_tbl_uns_bias,
                                                         ex_tbl_uns_orig,
                                                         bias) {
  data_count <- data_count[order(.get_cut(data_count)), ]
  # data_count <- data_count |>
  #   dplyr::mutate(count_stim = seq_len(dplyr::n())) |>
  #   dplyr::mutate(
  #     prop_stim = count_stim / nrow(ex_tbl_stim_orig)
  #   )
  prop_stim_vec <- purrr::map_dbl(.get_cut(data_count), function(x) {
    sum(.get_cut(ex_tbl_stim_orig) >= x) / nrow(ex_tbl_stim_orig)
  })
  prop_uns_vec <- purrr::map_dbl(.get_cut(data_count), function(x) {
    sum(.get_cut(ex_tbl_uns_orig) >= x) / nrow(ex_tbl_uns_orig)
  })
  data_count |>
    dplyr::mutate(
      prop_stim = prop_stim_vec,
      prop_uns = prop_uns_vec,
    ) |>
    dplyr::mutate(
      prop_bs = prop_stim - prop_uns, # nolint
      prop_bs_diff = prop_bs - prop_bs_est # nolint
    )
}

.get_cp_uns_loc_get_cp_actual <- function(data_threshold,
                                          ex_tbl_stim_no_min,
                                          ex_tbl_uns_bias,
                                          cp_min,
                                          .debug) {
  if (nrow(data_threshold) == 0L) {
    return(.get_cp_uns_loc_ind_check_out(
      cp_min, ex_tbl_stim_no_min, ex_tbl_uns_bias, .debug,
      "Too few responding cells"
    )[["cp"]])
  }
  data_threshold <- data_threshold |>
    dplyr::filter(abs(prop_bs_diff) == min(abs(prop_bs_diff))) |> # nolint
    dplyr::slice(1) |>
    .get_cut()
}

.get_cp_uns_loc_output <- function(.debug,
                                   cp_uns_loc_obj_list,
                                   ind_uns,
                                   ind_stim) {
  cp_vec <- .get_cp_uns_loc_sample_cp_rep(
    .debug = .debug,
    cp_uns_loc_obj_list = cp_uns_loc_obj_list,
    ind_uns = ind_uns,
    ind_stim = ind_stim
  )
  .debug_msg(.debug, "done getting loc gate at sample level") # nolint
  # collate plots
  list(
    "loc" = cp_vec,
    "p_list" = list()
  )
}

.get_cp_uns_loc_sample_cp_rep <- function(.debug,
                                          cp_uns_loc_obj_list,
                                          ind_uns,
                                          ind_stim) {
  .debug_msg(.debug, "Possibly re-using calculated cutpoints") # nolint

  # extract vector of cutpoints
  cp_vec <- purrr::map_dbl(cp_uns_loc_obj_list, ~ .x[["cp"]])

  # Repeat cutpoint if it was prejoined
  if (length(cp_vec) != length(ind_stim)) {
    cp_vec <- stats::setNames(
      rep(cp_vec, length(ind_stim)), ind_stim
    )
  } else {
    # name gate indices if not prejoined
    cp_vec <- stats::setNames(cp_vec, ind_stim)
  }

  # add unstim if unstim gate required
  if (!all(purrr::map_lgl(cp_vec, is.na))) {
    cp_vec <- c(cp_vec, stats::setNames(mean(cp_vec, na.rm = TRUE), ind_uns))
  } else {
    cp_vec <- c(cp_vec, NA)
  }

  cp_vec
}
