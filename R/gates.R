#' @title Get gates
#' 
#' @description
#' Get all the gates for each of the markers gated.
#' 
#' @param path_project character. Path to the project directory.
#' 
#' @return
#' Gate table with gates for each sample for each marker.
#' @export
get_gate_tbl <- function(path_project) {
  dir_vec <- list.dirs(path_project, recursive = FALSE) |>
    basename()
  purrr::map_df(dir_vec, function(dir_curr) {
    # get stats tbl
    gate_tbl <- readRDS(file.path(path_project, dir_curr, "gate_tbl.rds"))

    gate_stats <- gate_stats |>
      dplyr::mutate(chnl_cut = dir_curr) |>
      dplyr::select(chnl_cut, everything())

    gate_stats
  })
}
