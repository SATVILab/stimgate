#' @param data GatingSet.
#' Original GatingSet algorithm was applied to, with the same name.
#' @param pop_gate character.
#' Name of population to gate.
#' @param cut character vector.
#' Channels to gate on.
#' @param ind_in_batch_lab_vec character.
#' I-th element specifies name of i-th sample within a batch.
#' @param gate_name character. Name of gate.
#' @param ind_in_batch_gate numeric vector.
#' Specifies indices for each batch to gate.
#' @param ind_in_batch_uns numeric.
#' Specifies which index is unstim within each batch.
#'
#' @return
#' Gate table with gates for each sample for each marker.
#' @export
get_gate_tbl <- function(data,
                         pop_gate,
                         cut,
                         ind_in_batch_lab_vec,
                         gate_name,
                         ind_in_batch_gate,
                         ind_in_batch_uns,
                         path_project) {

  purrr::map_df(cut, function(cut_curr) {
    # get stats tbl
    gate_stats <- readRDS(file.path(path_project, cut_curr, "gate_stats_tbl"))

    gate_stats <- gate_stats |>
      dplyr::mutate(cut = cut_curr, marker = chnl_lab_vec[cut_curr]) |>
      dplyr::select(cut, marker, everything())

    gate_stats
  })
}
