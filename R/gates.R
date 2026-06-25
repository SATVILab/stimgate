#' @title Get gates
#'
#' @description Get all the gates for each of the markers gated.
#'
#' @param pathProject character. Path to the project directory.
#' @param pop character. Optional population name(s) to filter gates by. Default is NULL (all populations).
#' @param marker character. Optional marker name(s) to filter gates by. Default is NULL (all markers).
#' @param chnl character. Optional channel name(s) to filter gates by. Default is NULL (all channels).
#'
#' @return Gate table with gates for each sample for each marker.
#' @examples{
#' # Get example dataset
#' exampleData <- get_example_data()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#'
#' # Run the stimgate pipeline
#' pathProject <- stimgate_gate(
#'   pathProject = file.path(tempdir(), "getGateExample"),
#'   .data = gs,
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker,
#'   popGate = "root"
#' )
#'
#' # Get statistics for the identified gates
#' gates <- getStimGates(pathProject)
#' }
#' @export
getStimGates <- function(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL
) {
  pop <- pop %|c|% .gateGetPop(pathProject)
  
  purrr::map_df(pop, function(popCurr) {
    chnlVec <- if (!is.null(marker)) {
      markerLab <- stimgate_meta_read_marker_lab(pathProject)
      markerLab[marker] |> stats::setNames(NULL)
    } else {
      chnl %|c|% .gateGetChnl(pathProject, popCurr)
    }
    
    purrr::map_df(chnlVec, function(chnlCurr) {
      markerCurr <- stimgate_meta_read_chnl_lab(pathProject)[chnlCurr] |>
        stats::setNames(NULL)
        
      .gatesGetPathAll(
        pathProject = pathProject, 
        pop = popCurr, 
        chnlCut = chnlCurr, 
        init = FALSE
      ) |>
        readRDS() |>
        dplyr::filter(chnl == chnlCurr) |>
        dplyr::mutate(marker = markerCurr) |>
        dplyr::mutate(pop = popCurr) |>
        dplyr::select(pop, dplyr::everything())
    })
  })
}

#' @keywords internal
.gateGetPop <- function(pathProject) {
  pathDir <- file.path(pathProject, "gates")
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  dirVec <- list.dirs(pathDir, full.names = FALSE, recursive = FALSE)
  popVec <- unique(sub("^pop_(.*)$", "\\1", dirVec))
  popVec
}

#' @keywords internal
.gateGetChnl <- function(pathProject, pop) {
  pathDir <- file.path(pathProject, "gates", paste0("pop_", pop))
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  dirVec <- list.dirs(pathDir, full.names = FALSE, recursive = FALSE)
  chnlVec <- unique(sub("^chnl_(.*)$", "\\1", dirVec))
  chnlVec
}

#' @keywords internal
.gatesGetPathAll <- function(pathProject, pop, chnlCut, init) {
  file.path(
    pathProject,
    "gates",
    paste0("pop_", pop),
    paste0("chnl_", chnlCut),
    "all",
    if (init) "gateTblInit.rds" else "gateTbl.rds"
  )
}
