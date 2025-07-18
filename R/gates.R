#' @title Get gates
#'
#' @description Get all the gates for each of the markers gated.
#'
#' @param path_project character. Path to the project directory.
#'
#' @return Gate table with gates for each sample for each marker.
#' @examples
#' \dontrun{
#' # Get gate table from project directory
#' gate_tbl <- get_gate_tbl("/path/to/project")
#' }
#' @export
get_gate_tbl <- function(path_project) {
  dir_vec <- list.dirs(path_project, recursive = FALSE) |>
    basename()
  purrr::map_df(dir_vec, function(dir_curr) {
    # get stats tbl
    readRDS(file.path(path_project, dir_curr, "gate_tbl.rds")) |>
      dplyr::filter(chnl == dir_curr)
  })
}
