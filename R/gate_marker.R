.gate_marker <- function(.data,
                         ind_batch_list,
                         pop_gate,
                         chnl_cut,
                         gate_combn,
                         tol,
                         noise_sd,
                         bias_uns,
                         exc_min,
                         bw_min,
                         cp_min,
                         min_cell,
                         max_pos_prob_x,
                         gate_quant,
                         tol_clust,
                         gate_tbl = NULL,
                         tol_gate_single,
                         calc_cyt_pos_gates,
                         path_project,
                         .debug) {
  # message progress
  .debug_msg(.debug, "pop_gate: ", pop_gate) # nolint

  # Parameters list
  # ----------------

  # named
  chnl_lab_vec <- .get_labs( # nolint
    .data = .data[[ind_batch_list[[1]]]],
    chnl_cut = chnl_cut
  )

  # get parameters
  params <- list(
    ind_batch_list = ind_batch_list,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut,
    gate_combn = gate_combn,
    chnl_lab = chnl_lab_vec,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    bw_min = bw_min,
    cp_min = cp_min,
    max_pos_prob_x = max_pos_prob_x,
    gate_quant = gate_quant,
    tol_clust = tol_clust,
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
    chnl_cut,
    gate_combn = gate_combn,
    tol_clust = tol_clust,
    noise_sd = noise_sd,
    bias_uns = bias_uns,
    exc_min = exc_min,
    bw_min = bw_min,
    cp_min = cp_min,
    min_cell = min_cell,
    params = params,
    .debug = .debug,
    path_project = path_project
  )

  gate_tbl <- gate_tbl |>
    dplyr::filter(
      !ind %in% vapply(
        ind_batch_list, function(x) as.character(x[length(x)]), character(1)
      )
    ) # nolint

  # Get adjusted and/or clustered gates
  # ----------------

  message("getting clustered and/or controlled gates")

  # =================================================
  # Cluster-based or tail-gate controlled gating
  # =================================================

  # =============================
  # For all cells
  # =============================

  .gate_marker_get_adj_gates( # nolint
    gate_tbl_params = params$gate_tbl,
    gate_tbl = gate_tbl,
    tol_clust = tol_clust,
    gate_quant = gate_quant,
    params = params,
    chnl_cut,
    bw_min = bw_min,
    .data = .data,
    path_project = path_project,
    .debug = .debug,
    pop_gate = pop_gate,
    ind_batch_list = ind_batch_list
  )
}
