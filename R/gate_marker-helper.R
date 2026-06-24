#' @keywords internal
.gate_chnl <- function(
  .data,
  ind_batch_list,
  pop_gate,
  chnl_cut,
  chnl_settings,
  gate_tbl = NULL,
  tol_gate_single = NULL,
  calc_cyt_pos_gates = NULL,
  path_project,
  stage
) {
  # message progress
  .debug("pop_gate: ", pop_gate) # nolint

  # Parameters list
  # ----------------

  # named
  chnl_lab_vec <- .get_labs(
    # Get first list element, then first item to cleanly evaluate
    .data = .data[[ind_batch_list[[1]][1]]], 
    chnl_cut = chnl_cut
  )

  # delete locb gates
  .gate_chnl_delete_old_gates()

  # Initial gates
  # ----------------
  gate_tbl <- .gate_chnl_pre_adj_gates_gate(
    ind_batch_list = ind_batch_list,
    .data = .data,
    pop_gate = pop_gate,
    chnl_cut = chnl_cut,
    chnl_settings = chnl_settings,
    stage = stage,
    path_project = path_project
  )

  gate_tbl <- gate_tbl |>
    dplyr::filter(
      !as.character(ind) %in%
        vapply(
          ind_batch_list,
          function(x) as.character(x[length(x)]),
          character(1)
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

  .gate_chnl_get_adj_gates(
    gate_tbl_params = chnl_settings$gate_tbl,
    gate_tbl = gate_tbl,
    tol_clust = chnl_settings$tol_clust,
    gate_quant = chnl_settings$gate_quant,
    .data = .data,
    chnl_cut = chnl_cut,
    bw_min = chnl_settings$bw_min,
    path_project = path_project,
    ind_batch_list = ind_batch_list,
    pop_gate = pop_gate,
    stage = stage,
    chnl_lab = chnl_lab_vec,
    calc_cyt_pos_gates = calc_cyt_pos_gates %||% chnl_settings$calc_cyt_pos_gates
  )
}
