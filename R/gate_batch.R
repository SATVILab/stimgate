# Get cutpoints for a single batch
#' @keywords internal
.gate_batch <- function(
  .data,
  ind_batch,
  chnl_settings,
  batch,
  stage,
  path_project
) {
  # get list of dataframes
  ex_list <- .get_ex_list(
    # nolint
    .data = .data,
    ind_batch = ind_batch,
    pop = chnl_settings$pop_gate,
    chnl_cut = chnl_settings$chnl_cut,
    batch = batch,
    path_project = path_project
  )

  if (is.null(chnl_settings$gate_tbl)) {
    .gate_batch_all(
      ind_batch = ind_batch,
      batch = batch,
      ex_list = ex_list,
      .data = .data,
      chnl_settings = chnl_settings,
      stage = stage,
      path_project = path_project
    )
  } else {
    .gate_batch_single(
      ind_batch = ind_batch,
      batch = batch,
      ex_list = ex_list,
      .data = .data,
      chnl_settings = chnl_settings,
      calc_cyt_pos_gates = calc_cyt_pos_gates,
      stage = stage,
      path_project = path_project
    )
  }
}
