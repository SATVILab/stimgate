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

#' @title Read gating statistics
#' @description Read cell counts and background-subtracted frequencies saved
#'   by [gateStim()].
#' @param pathProject character Project directory from [gateStim()].
#' @return A tibble (or data.frame when read from CSV) with sample and gate
#'   identifiers, `countStim`, `countUns`, `nCellStim`, `nCellUns`, proportions
#'   `propStim`, `propUns`, `propBs`, and percentages `freqStim`, `freqUns`,
#'   `freqBs`. Background subtraction is stimulated minus unstimulated.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' getStimStats(pathProject)
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
