#' @title Get gates
#'
#' @description Get all the gates for each of the markers gated.
#'
#' @param path_project character. Path to the project directory.
#'
#' @return Gate table with gates for each sample for each marker.
#' @examples
#' \dontrun{
#' # Load example data
#' library(flowWorkspace)
#' library(HDCytoData)
#' 
#' # Get flowSet (using HDCytoData example)
#' fs <- HDCytoData::Bodenmiller_BCR_XL_flowSet()
#' 
#' # Set up batch structure (unstim/stim pairs)
#' batch_list <- list(
#'   batch1 = c(1, 9),   # sample 1 = unstim, sample 9 = stim
#'   batch2 = c(2, 10),  # sample 2 = unstim, sample 10 = stim  
#'   batch3 = c(3, 11)   # sample 3 = unstim, sample 11 = stim
#' )
#' 
#' # Create GatingSet
#' gs <- GatingSet(fs)
#' 
#' # Define markers to gate
#' markers <- c("pS6", "pNFkB", "pp38")
#' 
#' # Create temporary project directory
#' path_project <- tempfile("stimgate_project")
#' 
#' # Run stimgate_gate first to generate results
#' stimgate_gate(
#'   path_project = path_project,
#'   .data = gs,
#'   batch_list = batch_list,
#'   marker = markers
#' )
#' 
#' # Now get the gate table
#' gate_tbl <- stimgate_gate_get(path_project)
#' print(gate_tbl)
#' }
#' @export
stimgate_gate_get <- function(path_project) {
  dir_vec <- list.dirs(path_project, recursive = FALSE) |>
    basename()
  purrr::map_df(dir_vec, function(dir_curr) {
    # get stats tbl
    readRDS(file.path(path_project, dir_curr, "gate_tbl.rds")) |>
      dplyr::filter(chnl == dir_curr)
  })
}
