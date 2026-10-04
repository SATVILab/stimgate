#' @keywords internal
.getStats <- function(
  gateTbl = NULL,
  chnl = NULL,
  filterOtherCytPos = FALSE,
  combn = TRUE,
  gateTypeCytPosFilter = "base",
  gateTypeCytPosCalc,
  popGate,
  chnlLab = NULL,
  .data,
  save = FALSE,
  indBatchList,
  saveGateTbl = FALSE,
  gateName = NULL,
  tolClust = NULL,
  pathProject
) {
  # prep
  # ---------------
  chnlLab <- chnlLab %||% .getLabs( # nolint: object_usage_linter.
    .data = .data[[1]], chnlCut = chnl
  )

  gateTbl <- .getStatsGateTblGet(
    gateTbl = gateTbl,
    chnlLab = chnlLab,
    pathProject = pathProject,
    popGate = popGate,
    gateName = gateName,
    tolClust = tolClust
  )

  chnl <- chnl %||% unique(gateTbl$chnl)
  gateName <- gateName %||% unique(gateTbl$gateName)

  if ((!filterOtherCytPos) && combn) {
    nChnl <- length(chnl)
    combnMatList <- .getStatsCombnMatListGet(
      nChnl = nChnl
    )
    cytCombnVecList <- .getStatsCytCombnVecListGet(
      combnMatList = combnMatList,
      chnl = chnl
    )
  } else {
    combnMatList <- NULL
    cytCombnVecList <- NULL
  }

  .getStatsGateTblSave(
    gateTbl = gateTbl,
    pathProject = pathProject,
    popGate = popGate,
    chnlLab = chnlLab,
    chnl = chnl,
    save = saveGateTbl
  )

  statTbl <- .getStatsOverall(
    indBatchList = indBatchList,
    gateTbl = gateTbl,
    chnl = chnl,
    combn = combn,
    cytCombnVecList = cytCombnVecList,
    gateTypeCytPosCalc = gateTypeCytPosCalc,
    gateTypeCytPosFilter = gateTypeCytPosFilter,
    popGate = popGate,
    .data = .data,
    chnlLab = chnlLab,
    filterOtherCytPos = filterOtherCytPos,
    combnMatList = combnMatList,
    gateName = gateName,
    pathProject = pathProject
  )

  # save it
  .statsSave(
    save = save,
    statTbl = statTbl,
    pathProject = pathProject
  )
}

#' @title Get gating statistics
#'
#' @description Read and return gating statistics computed during gating.
#'
#' @param pathProject character. Path to the project directory.
#' @return A data frame with gating statistics.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   pathProject = file.path(tempdir(), "getStimStatsExample"),
#'   .data = gs,
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker,
#'   popGate = "root"
#' )
#'
#' # Get gating statistics
#' statTbl <- getStimStats(pathProject)
#' @export
getStimStats <- function(pathProject) {
  pathStatsPartial <- file.path(pathProject, "gateStats")
  if (file.exists(paste0(pathStatsPartial, ".rds"))) {
    readRDS(paste0(pathStatsPartial, ".rds"))
  } else if (file.exists(paste0(pathStatsPartial, ".csv"))) {
    utils::read.csv(paste0(pathStatsPartial, ".csv"))
  } else {
    stop(
      "No stats file found"
    )
  }
}
