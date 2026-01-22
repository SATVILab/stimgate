# Get cutpoints for a single batch
#' @keywords internal
.gate_batch <- function(.data,
                        ind_batch,
                        pop_gate,
                        chnl_cut,
                        gate_combn,
                        tol_clust,
                        noise_sd,
                        bias_uns,
                        exc_min,
                        bw_min,
                        cp_min,
                        min_cell,
                        params,
                        batch,
                        stage,
                        path_project) {
  # get list of dataframes
  ex_list <- .get_ex_list( # nolint
    .data = .data,
    ind_batch = ind_batch,
    pop = pop_gate,
    chnl_cut,
    batch = batch,
    path_project = path_project
  )
  if (is.null(params$gate_tbl)) {
    .gate_batch_all(
      ind_batch, batch, ex_list, gate_combn, .data, noise_sd,
      bias_uns, exc_min, cp_min, min_cell, tol_clust,
      bw_min, params, stage, path_project
    )
  } else {
    .gate_batch_single(
      .debug, ind_batch, batch,
      ex_list, .data, noise_sd,
      bias_uns, exc_min, cp_min, min_cell, chnl_cut, tol_clust, bw_min, params,
      stage, path_project
    )
  }
}
