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
#'   path_project = file.path(tempdir(), "get_gate_example"),
#'   .data = gs,
#'   batch_list = example_data$batch_list,
#'   marker = example_data$marker,
#'   pop_gate = "root"
#' )
#'
#' # Get statistics for the identified gates
#' gates <- stimgate_gate_get(path_project)
#' }
#' @export
stimgate_gate_get <- function(path_project,
                              pop = NULL,
                              marker = NULL,
                              chnl = NULL) {
  pop <- pop %||% .gate_get_pop(path_project)
  purrr::map_df(pop, function(pop_curr) {
    if (!is.null(marker)) {
      marker_lab <- stimgate_meta_read_marker_lab(path_project)
      chnl <- marker_lab[marker] |>
        stats::setNames(NULL)
    } else {
      chnl %||% .gate_get_chnl(path_project, pop_curr)
    }
    purrr::map_df(chnl, function(chnl_curr) {
      marker_curr <- stimgate_meta_read_chnl_lab(path_project)[chnl_curr] |>
        stats::setNames(NULL)
      .gates_get_path_all(path_project, pop_curr, chnl_curr, FALSE) |>
        readRDS() |>
        dplyr::filter(chnl == chnl_curr) |>
        dplyr::mutate(marker = marker_curr) |>
        dplyr::mutate(pop = pop_curr) |>
        dplyr::select(pop, dplyr::everything())
    })
  })
}

#' @keywords internal
.gate_get_pop <- function(path_project) {
  path_dir <- file.path(path_project, "gates")
  if (!dir.exists(path_dir)) {
    return(character(0))
  }
  dir_vec <- list.dirs(path_dir, full.names = FALSE, recursive = FALSE)
  pop_vec <- unique(sub("^pop_(.*)$", "\\1", dir_vec))
  pop_vec
}

#' @keywords internal
.gate_get_chnl <- function(path_project, pop) {
  path_dir <- file.path(path_project, "gates", paste0("pop_", pop))
  if (!dir.exists(path_dir)) {
    return(character(0))
  }
  dir_vec <- list.dirs(path_dir, full.names = FALSE, recursive = FALSE)
  chnl_vec <- unique(sub("^chnl_(.*)$", "\\1", dir_vec))
  chnl_vec
}

#' @keywords internal
.gates_get_path_all <- function(path_project,
                                pop,
                                chnl,
                                init) {
  file.path(
    path_project,
    "gates",
    paste0("pop_", pop),
    paste0("chnl_", chnl),
    "all",
    if (init) "gate_tbl_init.rds" else "gate_tbl.rds"
  )
}
