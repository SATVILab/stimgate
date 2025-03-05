
.get_gate_obj <- function(.data,
                          ind_batch_list,
                          pop_gate,
                          cut,
                          gate_combn,
                          tol,
                          noise_sd,
                          bias_uns,
                          bw_min,
                          cp_min,
                          min_cell,
                          plot,
                          max_pos_prob_x,
                          gate_quant,
                          tol_ctrl,
                          tol_gate,
                          gate_tbl = NULL,
                          tol_gate_single,
                          calc_cyt_pos_gates,
                          path_project,
                          debug) {
  # print progress
  .debug(debug, "pop_gate: ", pop_gate) # nolint


  # Parameters list
  # ----------------

  # named
  chnl_lab_vec <- .get_labs( # nolint
    .data = .data[[ind_batch_list[[1]][1]]],
    cut = cut
  )

  # get parameters
  params <- list(
    .data = .data,
    data_name = data_name,
    ind_batch_list = ind_batch_list,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    cut = cut,
    gate_combn = gate_combn,
    pop_man_sub = pop_man_sub,
    pop_man_match_exact = pop_man_match_exact,
    tol = tol,
    chnl_lab = chnl_lab_vec,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    bw_min = bw_min,
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
  ind_uns_vec <- ind_batch_list[length(ind_batch_list)]

  # Initial gates
  # ----------------
  gate_tbl <- .get_gate_obj_pre_adj_gates_gate( # nolint
    ind_batch_list = ind_batch_list,
    .data = .data,
    ind_in_batch_gate = ind_in_batch_gate,
    ind_in_batch_uns = ind_in_batch_uns,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    cut = cut,
    gate_combn = gate_combn,
    tol = tol,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    bw_min = bw_min,
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
    params = params,
    cut = cut,
    bw_min = bw_min,
    .data = .data,
    path_project = path_project,
    debug = debug,
    ind_in_batch_lab_vec = ind_in_batch_lab_vec,
    pop_gate = pop_gate,
    ind_batch_list = ind_batch_list,
    data_name = data_name,
    ind_in_batch_uns = ind_in_batch_uns
  )
}
