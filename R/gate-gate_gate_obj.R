#' @title Get gates for each sample
#' within each batch
#' for a given gating population and marker
#' @inheritParams gate # data, batch_size, ind_in_batch_uns, ind_in_batch_gate,
#' ind_in_batch_lab_vec,
#' @param pop_gate character
#' Population for which separate gates are to be
#' calculated.
#' @param cut character.
#' Name of channel to get gate for.
#' @param high named numeric vector.
#' Names are channel names and values are
#' thresholds above which a value for a given cell in this channel is deemed
#' 'high'.
#' @param cps_scp numeric vector. Values at which the
#' changepoint for probability of positivity in the cut
#' channel is to be checked. Default is 90.
#' @param fdr numeric vector. False discovery rates at
#' which the unstim-based gate is to be calculated.
#' Default is \code{c(0.3, 0.2, 0.1)}.
#' @param tol numeric. Tolerance value for the
#' \code{cytoUtils:::.cytokine_cutpoint} method.
#' Default is 0.5e-8 for CyTOF and flow.
#' @param gate_combn named list. Named list where names
#' are one of 'pre_join', 'mean', 'median', 'trim20' or
#' 'min', and elements are character vectors of 'scp',
#' 'dcp', 'tg', 'midp', 'uns#' and 'uns#r' (where # are
#' FDRs expressed as percentages). Each element therefore
#' specifies the method of combining gates from individual
#' samples within a group for a subset of the automated
#' gating methods. \code{'grp'} means to gate jointly,
#' whereas all of the others do what they sound like. If
#' not specified (i.e. \code{NULL}), then all gates are
#' performed individually on each sample. Default is
#' \code{NULL}.
.get_gate_obj <- function(data,
                          ind_batch_list,
                          ind_in_batch_gate,
                          ind_in_batch_uns,
                          ind_in_batch_lab_vec,
                          pop_gate = pop_gate_curr,
                          data_name,
                          cut,
                          high,
                          cps_scp,
                          gate_combn,
                          pop_man_sub,
                          pop_man_match_exact,
                          tol,
                          fdr,
                          noise_sd,
                          bias_uns,
                          min_bw,
                          cp_min,
                          boot_n = NULL,
                          boot_sd = NULL,
                          min_cell,
                          plot,
                          max_pos_prob_x,
                          gate_quant,
                          tol_ctrl,
                          tol_gate,
                          pop_sub,
                          gate_tbl = NULL,
                          tol_gate_single,
                          calc_cyt_pos_gates,
                          path_project,
                          debug) {
  # print progress
  .debug(debug, "pop_gate: ", pop_gate) # nolint

  # Preparation
  # ----------------
  pop_man_vec <- .get_pop_man_vec( # nolint
    pop_man_sub, pop_man_match_exact, pop_gate
  )

  # Parameters list
  # ----------------

  # named
  chnl_lab_vec <- .get_labs( # nolint
    data = data[[ind_batch_list[[1]][1]]],
    cut = cut,
    high = high
  )

  # get parameters
  params <- list(
    data = data,
    data_name = data_name,
    ind_batch_list = ind_batch_list,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    cut = cut, high = high,
    cps_scp = cps_scp,
    gate_combn = gate_combn,
    pop_man_sub = pop_man_sub,
    pop_man_match_exact = pop_man_match_exact,
    tol = tol, fdr = fdr,
    chnl_lab = chnl_lab_vec,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    min_bw = min_bw,
    cp_min = cp_min,
    boot_sd = boot_sd,
    boot_n = boot_n, min_cell = min_cell,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    gate_tbl = gate_tbl,
    calc_cyt_pos_gates = calc_cyt_pos_gates
  )

  # delete locb gates
  .get_gate_obj_delete_old_gates( # nolint
    params$bias_uns, params$data_name
  )

  # get indices of
  ind_uns_vec <- .get_gate_obj_ind_uns_vec_get( # nolint
    ind_batch_list, ind_in_batch_uns
  )

  # Initial gates
  # ----------------
  gate_tbl <- .get_gate_obj_pre_adj_gates_gate( # nolint
    ind_batch_list = ind_batch_list,
    data = data,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    cut = cut,
    high = high,
    cps_scp = cps_scp,
    gate_combn = gate_combn,
    pop_man = pop_man_vec,
    pop_man_match_exact = pop_man_match_exact,
    tol = tol,
    fdr = fdr,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    min_bw = min_bw,
    cp_min = cp_min,
    boot_n = boot_n,
    boot_sd = boot_sd,
    min_cell = min_cell,
    data_name = data_name,
    params = params,
    plot = plot,
    debug = debug
  )

  gate_tbl <- gate_tbl |>
    dplyr::filter(!ind %in% ind_uns_vec) # nolint

  # Get adjusted and/or clustered gates
  # ----------------

  print("getting clustered and/or controlled gates")

  # =================================================
  # Cluster-based or tail-gate controlled gating
  # =================================================

  # =============================
  # For all cells
  # =============================
  .get_gate_obj_get_adj_gates( # nolint
    gate_tbl_params = params$gate_tbl,
    gate_tbl = gate_tbl,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    gate_quant = gate_quant,
    pop_sub = pop_sub,
    params = params,
    cut = cut,
    min_bw = min_bw,
    data = data,
    path_project = path_project,
    debug = debug,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    ind_batch_list = ind_batch_list,
    data_name = data_name,
    ind_in_batch_uns = ind_in_batch_uns
  )
}
