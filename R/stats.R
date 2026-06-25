#' @keywords internal
.getStats <- function(
  gateTbl = NULL,
  chnl = NULL,
  filterOtherCytPos = FALSE,
  combn = TRUE,
  gateTypeCytPosFilter = "base",
  gateTypeSinglePosFilter = "base",
  gateTypeCytPosCalc,
  gateTypeSinglePosCalc,
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
  chnlLab <- .getStatsChnlLabGet(
    chnlLab = chnlLab,
    .data = .data,
    chnl = chnl
  )

  gateTbl <- .getStatsGateTblGet(
    gateTbl = gateTbl,
    chnlLab = chnlLab,
    pathProject = pathProject,
    popGate = popGate,
    gateName = gateName,
    tolClust = tolClust
  )

  chnl <- .getStatsChnlGet(
    chnl = chnl,
    gateTbl = gateTbl
  )

  gateName <- .getStatsGateNameGet(
    gateName = gateName,
    gateTbl = gateTbl
  )

  if ((!filterOtherCytPos) && combn) {
    nChnl <- length(chnl)
    combnMatList <- .getStatsCombnMatListGet(
      nChnl = nChnl,
      nPos = 2
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
    gateTypeSinglePosCalc = gateTypeSinglePosCalc,
    gateTypeCytPosFilter = gateTypeCytPosFilter,
    gateTypeSinglePosFilter = gateTypeSinglePosFilter,
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

#' @keywords internal
.readGateStats <- function(statsSaveOutput) {
  if (inherits(statsSaveOutput, "data.frame")) {
    return(statsSaveOutput)
  }
  if (!inherits(statsSaveOutput, "character")) {
    stop(
      "statsSaveOutput must be a character string if not a data.frame."
    )
  }
  pathStats <- file.path(statsSaveOutput, "gate_stats.rds")
  gateStatsTbl <- readRDS(pathStats)
  if ("ind" %in% colnames(gateStatsTbl)) {
    gateStatsTbl[, "ind"] <- as.character(gateStatsTbl[["ind"]])
  }
  if ("batch" %in% colnames(gateStatsTbl)) {
    gateStatsTbl[, "batch"] <- as.character(gateStatsTbl[["batch"]])
  }
  gateStatsTbl
}

#' @title Get gating statistics
#' @param pathProject character. Path to the project directory.
#' @return A data frame with gating statistics.
#' @export
getStimGates <- function(pathProject) {
  pathStatsPartial <- file.path(pathProject, "gate_stats")
  if (file.exists(paste0(pathStatsPartial, ".rds"))) {
    readRDS(paste0(pathStatsPartial, ".rds"))
  } else if (file.exists(paste0(pathStatsPartial, ".csv"))) {
    read.csv(paste0(pathStatsPartial, ".csv"))
  } else {
    stop(
      "No stats file found"
    )
  }
}