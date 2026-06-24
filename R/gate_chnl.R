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
    # nolint
    .data = .data[[ind_batch_list[[1]]]],
    chnl_cut = chnl_cut
  )

  # delete locb gates
  .gate_chnl_delete_old_gates()

  # Initial gates
  # ----------------
  gate_tbl <- .gate_chnl_pre_adj_gates_gate(
    # nolint
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
    # nolint
    gate_tbl = gate_tbl,
    gate_tbl_params = chnl_settings$gate_tbl,
    chnl_settings = chnl_settings,
    chnl_cut = chnl_cut,
    .data = .data,
    path_project = path_project,
    pop_gate = pop_gate,
    ind_batch_list = ind_batch_list,
    stage = stage
  )
}
