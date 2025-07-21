#' @title Get gates
#'
#' @description Get all the gates for each of the markers gated.
#'
#' @param path_project character. Path to the project directory.
#'
#' @return Gate table with gates for each sample for each marker.
#' @examples{
#' # Get example dataset
#' example_data <- get_example_data()
#' gs <- flowWorkspace::load_gs(example_data$path_gs)
#'
#' # Run the stimgate pipeline
#' path_project <- stimgate_gate(
#'   path_project = file.path(tempdir(), "stimgate_example"),
#'   .data = gs,
#'   batch_list = example_data$batch_list,
#'   marker = example_data$marker,
#'   pop_gate = "root"
#' )
#'
#' # Get statistics for the identified gates
#' gates <- get_gate_tbl(path_project)
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
