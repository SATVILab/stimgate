
.gate_marker <- function(.data,
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
    .data = .data[[ind_batch_list[[1]]]],
    cut = cut
  )

  # get parameters
  params <- list(
    .data = .data,
    ind_batch_list = ind_batch_list,
    pop_gate = pop_gate,
    cut = cut,
    gate_combn = gate_combn,
    tol = tol,
    chnl_lab = chnl_lab_vec,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    bw_min = bw_min,
    cp_min = cp_min,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_ctrl = tol_ctrl,
    tol_gate = tol_gate,
    tol_gate_single = tol_gate_single,
    gate_tbl = gate_tbl,
    calc_cyt_pos_gates = calc_cyt_pos_gates
  )

  # delete locb gates
  .gate_marker_delete_old_gates()

  # Initial gates
  # ----------------
  gate_tbl <- .gate_marker_pre_adj_gates_gate( # nolint
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    cut = cut,
    gate_combn = gate_combn,
    tol = tol,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    bw_min = bw_min,
    cp_min = cp_min,
    min_cell = min_cell,
    params = params,
    plot = plot,
    debug = debug
  )

  gate_tbl <- gate_tbl |>
    dplyr::filter(!ind %in% ind_batch_list[length(ind_batch_list)]) # nolint

  # Get adjusted and/or clustered gates
  # ----------------

  print("getting clustered and/or controlled gates")

  # =================================================
  # Cluster-based or tail-gate controlled gating
  # =================================================

  # =============================
  # For all cells
  # =============================

  .gate_marker_get_adj_gates( # nolint
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
    pop_gate = pop_gate,
    ind_batch_list = ind_batch_list
  )
}
