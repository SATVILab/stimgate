# Get cutpoints using clustering approach
# Clusters thresholds from similar distributions to identify optimal cutpoints
#' @keywords internal
.getCpCluster <- function(
  .data,
  gateTbl,
  gateStatsTbl,
  gateTblCtrl,
  chnlSettings,
  stage,
  pathProject,
  control = list(),
  filterOtherCytPos,
  calcCytPosGates,
  indBatchList
) {
  # ==================================
  # PREPARATION
  # ==================================
  .browse(seq_along(.data)) # nolint

  stageChnl <- file.path(stage,  chnlSettings$chnlCut)

  # control
  control <- .get_cp_cluster_control_update(control) # nolint

  # statistics
  gateStatsTbl <- .get_cp_cluster_gate_stats_tbl_update(
    gateStatsTbl
  )

  # cp
  cpMin <- .getCpClusterCpGetMin(gateTbl, gateTblCtrl) # nolint
  maxCp <- .getCpClusterCpGetMax(gateTbl, gateTblCtrl) # nolint

  propBsByCpTblObj <- .getCpClusterPropBsByCpTblObj(
    .data = .data,
    gateTbl = gateTbl,
    indBatchList = indBatchList,
    calcCytPosGates = calcCytPosGates,
    cpMin = cpMin,
    maxCp = maxCp,
    gateStatsTbl = gateStatsTbl,
    filterOtherCytPos = filterOtherCytPos,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )

  propBsByCpTbl <- propBsByCpTblObj[["propBsByCpTbl"]]
  exprMax <- propBsByCpTblObj[["exprMax"]]
  exprMin <- propBsByCpTblObj[["exprMin"]]

  .intSave("all", stageChnl, pathProject, propBsByCpTbl)

  densTbl <- .getCpClusterDensTblGet(
    indBatchList = indBatchList,
    .data = .data,
    filterOtherCytPos = filterOtherCytPos,
    calcCytPosGates = calcCytPosGates,
    chnlSettings = chnlSettings,
    exprMin = exprMin,
    exprMax = exprMax,
    gateTbl = gateTbl,
    control = control,
    pathProject = pathProject
  )

  .intSave("all", stageChnl, pathProject, densTbl)

  nClus <- .getCpClusterNClus(
    densTbl
  )
  clusVec <- .getCpClusterClus(
    densTbl,
    nClus
  )

  .intSave("all", stageChnl, pathProject, nClus, clusVec)

  densTbl <- densTbl |>
    dplyr::mutate(grp = as.character(clusVec))

  densTblGrpLab <- densTbl |>
    dplyr::group_by(grp, ind) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup()

  grpIndLabVec <- densTblGrpLab$grp |>
    stats::setNames(densTblGrpLab$ind)

  if (FALSE) {
    .getCpClusterPlotCheck1(densTbl)
  }

  # bind onto the summary statistics the group labels
  densTblJoin <- densTbl |>
    dplyr::group_by(batch, ind) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup()
  densTblJoin <- densTblJoin[, c("ind", "grp")]
  propBsByCpTbl <- propBsByCpTbl |>
    dplyr::left_join(densTblJoin, by = "ind")

  # =====================================
  # CALCULATE THRESHOLDS
  # =====================================

  cpGrpLabVec <- .getCpClusterCpGrpLabVecGet(
    propBsByCpTbl = propBsByCpTbl,
    exprMax = exprMax,
    exprMin = exprMin
  )

  gateSummStatTbl <- .getCpClusterGateSummStatTblGet(
    gateTbl = gateTbl,
    chnlCut = chnlSettings$chnlCut,
    grpIndLabVec = grpIndLabVec
  )

  # calculate thresholds
  # ---------------------------

  # calculate cpJoin for each ind
  cpTbl <- .getCpClusterCpJoinGet(
    propBsByCpTbl = propBsByCpTbl,
    cpGrpLabVec = cpGrpLabVec
  )

  gateTblChnl <- .getCpClusterGateTblChnlGet(
    gateTbl,
    chnl
  )

  # add other tables with useful information for creating new gates
  cpTbl <- .getCpClusterCpTblAddInfo(
    cpTbl = cpTbl,
    gateStatsTbl = gateStatsTbl,
    gateSummStatTbl = gateSummStatTbl,
    gateTblCtrl = gateTblCtrl,
    gateTblChnl = gateTblChnl
  )

  cpTbl <- .getCpClusterCpTblAddOrigQuantMin(
    cpTbl = cpTbl
  )

  # calculate for each individual cpLse (less than 0.01 standard errors)
  cpTbl <- .getCpClusterCpJoinLseGet(
    cpTbl = cpTbl
  )

  # add tail-gate-based thresholds
  cpTbl <- .getCpClusterCpJoinTgGet(
    cpTbl = cpTbl
  )

  # calculate cp where you don't go all the way to the new cp,
  # but only halfway from original cp (if original cp higher)
  cpTbl <- cpTbl |>
    .getCpClusterCpLseOrigMean() # nolint

  cpTbl <- cpTbl |>
    .getCpClusterCpJoinTgOrigMean() # nolint

  # calculate cp that is the minimum of lseOrigMean and tgOrig
  cpTbl <- cpTbl |>
    .getCpClusterCpJoinLseOrigMeanTg() # nolint

  # filter at cp just above cpJoinLseOrigMeanTg, in order
  # to get the propBsCpDiff closest to it
  cpTbl <- cpTbl |>
    .getCpClusterCpFilterAboveCpJoinLseOrigMeanTg() # nolint

  # if no gate found above, then set it to base threshold OR
  # impute based on group.
  cpTbl <- .getCpClusterCpImputeMissingBatch(
    cpTbl = cpTbl,
    chnlCut = chnlSettings$chnlCut,
    gateTbl = gateTbl,
    densTbl = densTbl
  )

  cpTbl <- .getCpClusterImputeMissingInd(
    cpTbl = cpTbl,
    chnlCut = chnlSettings$chnlCut,
    gateTbl = gateTbl,
    densTbl = densTbl
  )

  cpTbl <- .getCpClusterImputeMissingFinal(
    cpTbl = cpTbl,
    chnlCut = chnlSettings$chnlCut,
    gateTbl = gateTbl,
    densTbl = densTbl
  )

  cpTbl <- .getCpClusterImputeMissingFinalBatch(
    cpTbl = cpTbl,
    chnlCut = chnlSettings$chnlCut,
    gateTbl = gateTbl,
    densTbl = densTbl
  )

  .intSave("all", stageChnl, pathProject, cpTbl)

  # =========================
  # OUTPUT
  # =========================
  cpTbl
}
