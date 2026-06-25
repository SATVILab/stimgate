#' @keywords internal
.getCpClusterControlUpdate <- function(control) {
  if (!"minThresholdFrac" %in% names(control)) {
    control[["minThresholdFrac"]] <- 0.8
  }
  if (!"minThresholdQuant" %in% names(control)) {
    control[["minThresholdQuant"]] <- 0.1
  }
  control
}

#' @keywords internal
.getCpClusterGateStatsTblUpdate <- function(gateStatsTbl) {
  .debug("Updating gate statistics table") # nolint
  gateStatsTbl |>
    dplyr::mutate(
      propStimPos = pmax(countStim, 1) / nCellStim, # nolint
      propUnsPos = pmax(countUns, 1) / nCellUns, # nolint
      propStimSd = sqrt(
        propStimPos * (1 - propStimPos) / nCellStim # nolint
      ),
      propUnsSd = sqrt(
        propUnsPos * (1 - propUnsPos) / nCellUns # nolint
      ),
      propBsSd = sqrt(
        propStimSd^2 + propUnsSd^2 # nolint
      )
    )
}

#' @keywords internal
.getCpClusterCpGetMin <- function(gateTbl, gateTblCtrl) {
  suppressWarnings(
    min(gateTbl$gateCyt, gateTbl$gate, gateTblCtrl$gate, na.rm = TRUE)
  )
}

#' @keywords internal
.getCpClusterCpGetMax <- function(gateTbl, gateTblCtrl) {
  suppressWarnings(
    max(gateTbl$gate, gateTbl$gateSingle, gateTblCtrl$gate, na.rm = TRUE)
  )
}

#' @keywords internal
.getCpClusterPropBsByCpTblObj <- function(
  .data,
  gateTbl,
  indBatchList,
  calcCytPosGates,
  cpMin,
  maxCp,
  gateStatsTbl,
  filterOtherCytPos,
  pathProject,
  chnlSettings
) {
  .debug("Getting propBsByCpTbl object") # nolint

  dataListObj <- .getPropBsByCpTblDataList(
    .data = .data,
    gateTbl = gateTbl,
    indBatchList = indBatchList,
    chnlSettings = chnlSettings,
    calcCytPosGates = calcCytPosGates,
    maxCp = maxCp,
    filterOtherCytPos = filterOtherCytPos,
    cpMin = cpMin,
    pathProject = pathProject
  )

  propBsByCpTbl <- .getPropBsByCpTblActual(
    dataList = dataListObj[["dataList"]],
    cpMin = cpMin,
    maxCp = maxCp,
    gateStatsTbl = gateStatsTbl,
    indBatchList = indBatchList
  )

  list(
    "propBsByCpTbl" = propBsByCpTbl,
    "exprMax" = dataListObj[["exprMax"]],
    "exprMin" = dataListObj[["exprMin"]]
  )
}


#' @keywords internal
.getPropBsByCpTblDataList <- function(
  .data,
  gateTbl,
  indBatchList,
  chnlSettings,
  calcCytPosGates,
  maxCp,
  filterOtherCytPos,
  cpMin,
  pathProject
) {
  .debug("Getting .data list") # nolint
  dataList <- .getPropBsByCpTblDataListInit(
    indBatchList = indBatchList,
    .data = .data,
    chnlSettings = chnlSettings,
    filterOtherCytPos = filterOtherCytPos,
    calcCytPosGates = calcCytPosGates,
    cpMin = cpMin,
    pathProject = pathProject
  )
  .getPropBsByCpTblDataListFinal(dataList, maxCp)
}

#' @keywords internal
.getPropBsByCpTblDataListInit <- function(
  indBatchList,
  .data,
  chnlSettings,
  filterOtherCytPos,
  calcCytPosGates,
  cpMin,
  pathProject
) {
  purrr::map(seq_along(indBatchList), function(i) {
    indBatch <- indBatchList[[i]]
    exList <- .getExList(
      .data = .data,
      indBatch = indBatch,
      pop = chnlSettings$popGate,
      chnlCut = chnlSettings$chnlCut,
      batch = names(indBatchList)[i],
      pathProject = pathProject
    )

    minMaxVec <- .getPropBsByCpTblDataListMinmax(exList)

    exListFilter <- .getPropBsByCpTblDataListFilterCytPos(
      filterOtherCytPos = filterOtherCytPos,
      gateTbl = gateTbl,
      exList = exList,
      calcCytPosGates = calcCytPosGates
    )

    outTbl <- .getPropBsByCpTblDataListFilterAboveMin(
      exListFilter = exListFilter,
      cpMin = cpMin
    )

    list(
      "outTbl" = outTbl,
      "exprMin" = minMaxVec[[1]],
      "exprMax" = minMaxVec[[2]]
    )
  })
}

#' @keywords internal
.getPropBsByCpTblDataListMinmax <- function(exList) {
  exprRangeTbl <- purrr::map_df(
    seq_along(exList),
    function(i) {
      ex <- exList[[i]]
      if (nrow(ex) <= 5) {
        return(NULL)
      }
      quantVec <- quantile(.getCut(ex), c(0.0025, 0.999))
      tibble::tibble(
        lb = quantVec[[1]],
        ub = 3 * quantVec[[2]]
      )
    }
  )

  exprMin <- quantile(exprRangeTbl[["lb"]], 0.0025)
  exprMax <- max(exprRangeTbl[["ub"]])

  c("min" = exprMin, "max" = exprMax)
}

#' @keywords internal
.getPropBsByCpTblDataListFilterCytPos <- function(
  filterOtherCytPos,
  exList,
  gateTbl,
  calcCytPosGates
) {
  if (!filterOtherCytPos) {
    return(exList)
  }
  purrr::map(seq_along(exList), function(i) {
    if (i == 1) {
      return(exList[[i]])
    }
    gateTblInd <- gateTbl |>
      dplyr::filter(ind == attr(exList[[i]], "ind")) # nolint

    posIndVecButSinglePosCurr <-
      .get_pos_ind_but_single_pos_for_one_cyt(
        ex = exList[[i]],
        gateTbl = gateTblInd,
        chnlSingleExc = attr(exList[[i]], "chnlCut"),
        chnl = NULL,
        gateTypeCytPos = ifelse(calcCytPosGates, "cyt", "base"),
        gateTypeSinglePos = "base"
      )
    exList[[i]][!posIndVecButSinglePosCurr, , drop = FALSE]
  }) |>
    stats::setNames(names(exList))
}

#' @keywords internal
.getPropBsByCpTblDataListFilterAboveMin <- function(
  exListFilter,
  cpMin
) {
  exListFilter |>
    purrr::map(function(x) {
      attr(x, "nCell") <- nrow(x)
      xOut <- x |>
        dplyr::filter(.getCut(x) >= min(.env$cpMin, max(.getCut(x)))) # nolint
      if (nrow(xOut) == 0) {
        allCols <- colnames(x)
        batchIdx <- which(allCols == "batch")
        stimIdx <- which(allCols == "stim")
        selIdx <- seq(batchIdx, stimIdx)
        xOut <- x[1, selIdx, drop = FALSE]
        xAdd <- x[1, setdiff(seq_along(x), selIdx)]
        xAdd[1, ] <- NA
        xOut <- xOut |>
          dplyr::bind_cols(xAdd)
      }
      xOut
    })
}

#' @keywords internal
.getPropBsByCpTblDataListFinal <- function(dataList, maxCp) {
  exprMinVec <- vapply(dataList, function(x) x$exprMin, numeric(1))
  exprMaxVec <- vapply(dataList, function(x) x$exprMax, numeric(1))
  exprMin <- min(exprMinVec, na.rm = TRUE)
  exprMax <- max(
    max(exprMaxVec, na.rm = TRUE),
    maxCp + 0.2 * (max(exprMaxVec, na.rm = TRUE) - exprMin)
  )
  dataList <- lapply(dataList, function(x) x$outTbl) |>
    purrr::flatten()
  list(
    "dataList" = dataList,
    "exprMin" = exprMin,
    "exprMax" = exprMax
  )
}


#' @keywords internal
.getPropBsByCpTblActual <- function(
  dataList,
  cpMin,
  maxCp,
  gateStatsTbl,
  indBatchList
) {
  .debug("Getting propBsByCpTbl") # nolint
  cpParList <- .getPropBsByCpTblActualPrep(cpMin, maxCp)
  purrr::map(seq_along(dataList), function(i) {
    .getPropBsByCpTblActualInd(
      i = i,
      dataList = dataList,
      cpParList = cpParList,
      gateStatsTbl = gateStatsTbl,
      indBatchList = indBatchList
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()
}

#' @keywords internal
.getPropBsByCpTblActualInd <- function(
  i,
  cpParList,
  gateStatsTbl,
  indBatchList,
  dataList
) {
  .getPropBsByCpTblActualProgress(i, dataList)
  exList <- .getPropBsByCpTblActualExGet(
    dataList = dataList,
    i = i,
    indBatchList = indBatchList
  )
  if (is.null(exList)) {
    return(NULL)
  }
  .getPropBsByCpTblInd(
    exStim = exList$stim,
    exUns = exList$uns,
    cpSeq = cpParList[["seq"]],
    gateStatsTbl = gateStatsTbl
  )
}

#' @keywords internal
.getPropBsByCpTblActualPrep <- function(cpMin, cpMax) {
  cpRange <- c(cpMin, cpMax)
  cpSeqVec <- seq(cpRange[1], cpRange[2], length.out = 1e2)
  list("range" = cpRange, "seq" = cpSeqVec)
}

#' @keywords internal
.getPropBsByCpTblInd <- function(
  exStim,
  exUns,
  cpSeq,
  gateStatsTbl
) {
  parList <- .getPropBsByCpTblIndPrep(
    gateStatsTbl = gateStatsTbl,
    exStim = exStim,
    exUns = exUns,
    cpSeq = cpSeq
  )

  .getPropBsByCpTblIndInit(exStim, exUns, parList, cpSeq) |>
    .getPropBsByCpTblIndCalc(
      nCellStim = attr(exStim, "nCell"),
      nCellUns = attr(exUns, "nCell")
    )
}


#' @keywords internal
.getPropBsByCpTblIndPrep <- function(
  gateStatsTbl,
  exStim,
  exUns,
  cpSeq
) {
  gateStatsTblCurr <- gateStatsTbl |>
    dplyr::filter(.data$ind == attr(exStim, "ind")) # nolint
  countStimVec <- rep(NA, length(cpSeq))
  countUnsVec <- rep(NA, length(cpSeq))
  for (i in seq_along(cpSeq)) {
    cp <- cpSeq[i]
    countStimVec[i] <- sum(.getCut(exStim) > cp) # nolint
    countUnsVec[i] <- sum(.getCut(exUns) > cp) # nolint
  }
  propBsSd <- gateStatsTblCurr$propBsSd
  propBsOrig <- gateStatsTblCurr$propBs
  list(
    "gateStats" = gateStatsTblCurr,
    "countStim" = countStimVec,
    "countUns" = countUnsVec,
    "bsSd" = propBsSd,
    "bsOrig" = propBsOrig
  )
}

#' @keywords internal
.getPropBsByCpTblIndInit <- function(exStim, exUns, parList, cpSeq) {
  tibble::tibble(
    ind = attr(exStim, "ind"),
    propBsOrig = parList[["bsOrig"]],
    propBsSd = parList[["bsSd"]],
    cp = cpSeq,
    maxExpr = max(.getCut(exStim), .getCut(exUns)), # nolint
    countStimCp = parList[["countStim"]],
    countUnsCp = parList[["countUns"]]
  )
}

#' @keywords internal
.getPropBsByCpTblIndCalc <- function(.data, nCellStim, nCellUns) {
  .data |>
    dplyr::mutate(
      propStimCp = countStimCp / nCellStim, # nolint
      propUnsCp = countUnsCp / nCellUns # nolint
    ) |>
    dplyr::mutate(
      propBsCp = propStimCp - propUnsCp, # nolint
      propBsCpDiff = propBsCp - propBsOrig, # nolint
      propBsCpDiffSd = propBsCpDiff / propBsSd, # nolint
      propStimPosCp = pmax(countStimCp, 1) / nCellStim,
      propUnsPosCp = pmax(countUnsCp, 1) / nCellUns,
      propStimSdCp = sqrt(
        propStimPosCp * (1 - propStimPosCp) / nCellStim # nolint
      ),
      propUnsSdCp = sqrt(
        propUnsPosCp * (1 - propUnsPosCp) / nCellUns # nolint
      ),
      propBsSdCp = sqrt(
        propStimSdCp^2 + propUnsSdCp^2 # nolint
      ),
      propBsCpDiffSdMax = propBsCpDiff /
        pmax(propBsSd, propBsSdCp) # nolint
    )
}

#' @keywords internal
.getPropBsByCpTblActualProgress <- function(i, dataList) {
  if (i %% 20 == 0) {
    .debug(paste0("Processing ", i, " of ", length(dataList))) # nolint
  }
}

#' @keywords internal
.getPropBsByCpTblActualExGet <- function(dataList, i, indBatchList) {
  ex <- dataList[[i]]
  if (.getPropBsByCpReturnEarlyStim(ex)) {
    return(NULL)
  }
  exUns <- dataList[[attr(ex, "indUns")]]

  if (.getPropBsByCpReturnEarlyUns(exUns)) {
    return(NULL)
  }

  list("stim" = ex, "uns" = exUns)
}

#' @keywords internal
.getPropBsByCpReturnEarlyStim <- function(ex) {
  is.na(.getCut(ex)[1]) || attr(ex, "isUns")
}

#' @keywords internal
.getPropBsByCpReturnEarlyUns <- function(ex) {
  is.na(.getCut(ex)[1])
}

#' @keywords internal
.getCpClusterDensTblGetBatchPrep <- function(indBatch) {
  .debug(paste0("Processing batch ", indBatch)) # nolint
}


#' @keywords internal
.getCpClusterDensTblGet <- function(
  indBatchList,
  .data,
  filterOtherCytPos,
  calcCytPosGates,
  exprMin,
  exprMax,
  gateTbl,
  control,
  pathProject,
  chnlSettings
) {
  .debug("Getting density table") # nolint
  minThreshold <- .getCpClusterDensTblGetMinThreshold(
    gateTbl = gateTbl,
    control = control
  )
  densTbl <- purrr::map_df(seq_along(indBatchList), function(i) {
    indBatch <- indBatchList[[i]]
    exList <- .getCpClusterDensTblGetBatchPrepExList(
      .data = .data,
      indBatch = indBatch,
      popGate = chnlSettings$popGate,
      chnlCut = chnlSettings$chnlCut,
      filterOtherCytPos = filterOtherCytPos,
      gateTbl = gateTbl,
      calcCytPosGates = calcCytPosGates,
      control = control,
      batch = names(indBatchList)[i],
      pathProject = pathProject
    )

    purrr::map_df(exList, function(x) {
      .getCpClusterDensTblGetActualInd(
        exprVec = .getCut(x),
        batch = attr(x, "batch"),
        ind = attr(x, "ind"),
        minThreshold = minThreshold,
        chnlCut = chnlSettings$chnlCut,
        exprMin = exprMin,
        exprMax = exprMax,
        bw = chnlSettings$bwCluster
      )
    })
  }) |>
    dplyr::filter(!is.na(x1)) # nolint
  notAllNaVecInd <- purrr::map_lgl(
    seq_len(ncol(densTbl)),
    function(i) !all(is.na(densTbl[[i]]))
  )
  densTbl[, notAllNaVecInd]
}


#' @keywords internal
.getCpClusterDensTblGetActualInd <- function(
  exprVec,
  batch,
  ind,
  minThreshold,
  chnlCut,
  exprMin,
  exprMax,
  bw
) {
  if (.getCpClusterDensTblGetActualIndEarlyReturnCheck(exprVec)) {
    return(.getCpClusterDensTblGetActualIndEarlyReturn(
      batch = batch,
      ind = ind
    ))
  }
  dens <- density(
    exprVec[exprVec > min(exprVec)],
    from = exprMin,
    to = exprMax,
    bw = bw
  )
  tibble::tibble(
    batch = batch[[1]],
    ind = ind[[1]],
    y = dens[["y"]],
    x = dens[["x"]],
    xInd = paste0("x", seq_len(length(dens[["y"]])))
  ) |>
    dplyr::filter(x <= minThreshold) |> # nolint
    dplyr::select(-x) |>
    dplyr::mutate(y = y / sum(y)) |>
    tidyr::pivot_wider(names_from = xInd, values_from = y) # nolint
}

.getCpClusterDensTblGetActualIndEarlyReturnCheck <-
  function(exprVec) {
    length(exprVec[exprVec > min(exprVec)]) < 3
  }

#' @keywords internal
.getCpClusterDensTblGetActualIndEarlyReturn <- function(batch, ind) {
  tibble::tibble(
    batch = batch[1],
    ind = ind[1],
    y = rep(NA, 512),
    x = paste0("x", seq.int(from = 1, to = 512))
  ) |>
    tidyr::pivot_wider(names_from = x, values_from = y) # nolint
}


#' @keywords internal
.getCpClusterDensTblGetMinThreshold <- function(gateTbl, control) {
  minThresholdGateQuant <- quantile(
    gateTbl$gate,
    control$minThresholdQuant,
    na.rm = TRUE
  )
  control$minThresholdFrac * minThresholdGateQuant
}

#' @keywords internal
.getCpClusterDensTblGetBatchPrepExList <- function(
  .data,
  indBatch,
  batch,
  popGate,
  chnlCut,
  filterOtherCytPos,
  gateTbl,
  calcCytPosGates,
  control,
  pathProject
) {
  exList <- .getExList(
    .data = .data,
    indBatch = indBatch,
    pop = popGate,
    chnlCut = chnlCut,
    batch = batch,
    pathProject = pathProject
  )

  if (!filterOtherCytPos) {
    return(exList[-1])
  }

  .getCpClusterDensTblGetBatchPrepExListFilter(
    exList = exList,
    chnlCut = chnlCut,
    gateTbl = gateTbl,
    calcCytPosGates = calcCytPosGates
  )
}

#' @keywords internal
.getCpClusterDensTblGetBatchPrepExListFilter <- function(
  exList,
  chnlCut,
  gateTbl,
  calcCytPosGates
) {
  .debug("Filtering other cytokine positive cells") # nolint
  exListFilter <- purrr::map(seq_along(exList), function(i) {
    if (i == 1) {
      return(exList[[i]])
    }
    .getCpClusterDensTblGetBatchPrepExListFilterInd(
      exTbl = exList[[i]],
      gateTbl = gateTbl,
      chnlCut = chnlCut,
      calcCytPosGates = calcCytPosGates
    )
  }) |>
    stats::setNames(names(exList))

  exListFilter[-1]
}

#' @keywords internal
.getCpClusterDensTblGetBatchPrepExListFilterInd <- function(
  exTbl,
  gateTbl,
  chnlCut,
  calcCytPosGates
) {
  posIndVecButSinglePosCurr <-
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = exTbl,
      gateTbl = gateTbl[gateTbl[["ind"]] == attr(exTbl, "ind"), ],
      chnlSingleExc = chnlCut,
      chnl = NULL,
      gateTypeCytPos = if (calcCytPosGates) "cyt" else "base",
      gateTypeSinglePos = "base"
    )
  exTbl[!posIndVecButSinglePosCurr, , drop = FALSE]
}

#' @keywords internal
.getCpClusterNClus <- function(densTbl) {
  maxCluster <- min(6, nrow(densTbl) / 3) |>
    floor() |>
    max(1)
  if (maxCluster == 1L) {
    return(1L)
  }
  for (i in seq_len(ncol(densTbl))) {
    if (any(is.na(densTbl[[i]]))) {
      message(i)
    }
  }
  densMat <- densTbl[, grepl("^x\\d+", colnames(densTbl))] |>
    as.matrix()
  clusGapObj <- cluster::clusGap(
    densMat,
    FUNcluster = kmeans,
    K.max = maxCluster,
    B = 50
  )
  cluster::maxSE(
    clusGapObj$Tab[, "gap"],
    clusGapObj$Tab[, "SE.sim"]
  )
}

#' @keywords internal
.getCpClusterPlotCheck1 <- function(densTbl) {
  densPlot <- densTbl |>
    tidyr::pivot_longer(
      names_to = "x",
      values_to = "y",
      x1:x512 # nolint
    ) |>
    dplyr::group_by(grp, ind) |> # nolint
    dplyr::mutate(x = xVec) |> # nolint
    dplyr::ungroup()

  ggplot(
    densPlot,
    aes(x, y, col = grp, group = ind) # nolint
  ) +
    geom_line() # nolint
}

#' @keywords internal
.getCpClusterPlotCheck1lse <- function(propBsByCpTbl) {
  dataPlot <- propBsByCpTbl |>
    dplyr::group_by(grp, cp) |> # nolint
    dplyr::summarise(
      propL1se = sum(propBsCpDiffSd < 1) / dplyr::n() # nolint
    ) |>
    dplyr::ungroup()

  ggplot(dataPlot, aes(x = cp, y = propL1se, col = grp)) + # nolint
    geom_line() + # nolint
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = FALSE
    )
}

#' @keywords internal
.getCpClusterClus <- function(densTbl, nClus) {
  densMat <- densTbl[, grepl("^x\\d+", colnames(densTbl))] |>
    as.matrix()
  stats::kmeans(densMat, centers = nClus)$cluster
}

#' @keywords internal
.getCpClusterDataModPre <- function(propBsByCpTbl) {
  dataModPre <- propBsByCpTbl |>
    tidyr::pivot_longer(
      -c(ind:propBsCpDiffSdMax), # nolint
      names_to = "grp",
      values_to = "grpLevel"
    )

  dataModFilterVec <- dataModPre |>
    dplyr::group_by(ind, grp, grpLevel) |> # nolint
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::group_by(grp, grpLevel) |>
    dplyr::arrange(grp, grpLevel, ind) |>
    dplyr::summarise(indVec = paste0(ind, collapse = "_")) |>
    dplyr::ungroup() |>
    dplyr::group_by(indVec) |> # nolint
    dplyr::slice(1) |>
    dplyr::mutate(grpGrpLevel = paste0(grp, grpLevel)) |>
    dplyr::pull("grpGrpLevel")

  dataMod <- dataModPre |>
    dplyr::filter(
      paste0(grp, grpLevel) %in% dataModFilterVec # nolint
    ) |>
    dplyr::group_by(grp, grpLevel, cp) |> # nolint
    dplyr::summarise(
      propL1se = sum(propBsCpDiffSd < 1, na.rm = TRUE) / # nolint
        sum(!is.na(propBsCpDiffSd))
    ) |>
    dplyr::ungroup()

  dataMod |>
    dplyr::group_by(grp, grpLevel) |> # nolint
    dplyr::filter(sum(is.na(propL1se)) != dplyr::n()) |> # nolint
    dplyr::ungroup()
}

#' @keywords internal
.getCpClusterCpGrpLabVecGet <- function(
  propBsByCpTbl,
  exprMin,
  exprMax,
  .debug = FALSE
) {
  dataMod <- .getCpClusterDataModPre(
    propBsByCpTbl = propBsByCpTbl
  )
  purrr::map(unique(dataMod$grp), function(nGrpCurr) {
    .debug(paste0("Processing cluster ", nGrpCurr)) # nolint
    dataModCurr <- dataMod |>
      dplyr::filter(.data$grp == nGrpCurr) # nolint
    purrr::map(
      unique(dataModCurr$grpLevel),
      function(grpLevel) {
        dataModCurrGrp <- dataModCurr |>
          dplyr::filter(grpLevel == .env$grpLevel) # nolint
        dataModCurrGrpNotNa <- dataModCurrGrp |>
          dplyr::filter(!is.na(propL1se)) # nolint
        if (nrow(dataModCurrGrpNotNa) <= 5) {
          return(max(dataModCurrGrp$cp, na.rm = TRUE))
        }
        fit <- try(
          suppressWarnings(
            mgcv::gam(
              propL1se ~ s(cp, bs = "cs"),
              .data = dataModCurrGrpNotNa,
              family = mgcv::betar(link = "logit"),
              maxit = 20
            )
          ),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          return(max(dataModCurrGrp$cp, na.rm = TRUE))
        }
        dataModCurrCrp <- dataModCurrGrp |> # nolint
          dplyr::mutate(
            pred = predict(
              fit,
              dataModCurrGrp,
              type = "response"
            )
          )
        cpRange <- range(dataModCurrGrp$cp)
        minCpPermitted <- dataModCurrGrp$cp[
          dataModCurrGrp$propL1se > 0.5
        ] |>
          min(na.rm = TRUE)
        if (abs(minCpPermitted) == Inf) {
          return(max(dataModCurrGrp$cp, na.rm = TRUE))
        }
        dataPred <- tibble::tibble(
          cp = seq(minCpPermitted, cpRange[2], length.out = 1e5)
        )
        dataPred <- dataPred |>
          dplyr::mutate(
            pred = predict(fit, newdata = dataPred, type = "response")
          )
        dataPredDer <- dataPred[-1, ] |>
          dplyr::mutate(
            der = (pred - dataPred$pred[-nrow(dataPred)]) / # nolint
              (cp - dataPred$cp[-nrow(dataPred)]) # nolint
          )
        maxDer <- 0.1 / diff(c(exprMax, exprMin))
        cpDer <- dataPredDer |>
          dplyr::filter(der > 0) |> # nolint
          dplyr::filter(der < maxDer) |>
          dplyr::pull("cp") |>
          min()
        cpCdf <- dataPredDer |>
          dplyr::filter(pred > 0.85) |> # nolint
          dplyr::pull("cp") |>
          min()
        min(cpDer, cpCdf)
      }
    ) |>
      stats::setNames(paste0(nGrpCurr, unique(dataModCurr$grpLevel)))
  }) |>
    purrr::flatten() |>
    unlist()
}

#' @keywords internal
.getCpClusterGateSummStatTblGet <- function(
  gateTbl,
  chnlCut,
  grpIndLabVec
) {
  .debug(
    "Getting quantiles of original gates per clustered observations" # nolint
  )
  if ("chnl" %in% names(gateTbl)) {
    gateTbl <- gateTbl |>
      dplyr::filter(.data$chnl == .env$chnlCut) # nolint
  }
  gateTbl |>
    dplyr::mutate(grp = grpIndLabVec[as.character(ind)]) |> # nolint
    dplyr::group_by(grp) |> # nolint
    dplyr::summarise(
      gate05 = quantile(gate, 0.05, na.rm = TRUE), # nolint
      gate10 = quantile(gate, 0.1, na.rm = TRUE), # nolint
      gate25 = quantile(gate, 0.25, na.rm = TRUE), # nolint
      gate50 = median(gate, na.rm = TRUE), # nolint
      gate75 = quantile(gate, 0.75, na.rm = TRUE), # nolint
      gate90 = quantile(gate, 0.9, na.rm = TRUE), # nolint
      gate95 = quantile(gate, 0.95, na.rm = TRUE), # nolint
      .groups = "drop"
    )
}

#' @keywords internal
.getCpClusterCpJoinGet <- function(propBsByCpTbl, cpGrpLabVec) {
  .debug("Getting cpJoin") # nolint
  propBsByCpTbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      cpJoin = cpGrpLabVec[paste0("grp", grp)] # nolint
    ) |>
    dplyr::ungroup()
}

#' @keywords internal
.getCpClusterGateTblChnlGet <- function(gateTbl, chnl) {
  if ("chnl" %in% colnames(gateTbl)) {
    gateTblChnl <- gateTbl |>
      dplyr::filter(.data$chnl == .env$chnl) # nolint
  } else {
    gateTblChnl <- gateTbl
  }
  gateTblChnl
}

#' @keywords internal
.getCpClusterCpTblAddInfo <- function(
  cpTbl,
  gateStatsTbl,
  gateSummStatTbl,
  gateTblCtrl,
  gateTblChnl
) {
  .debug("Adding information to cp table") # nolint
  cpTbl |>
    dplyr::left_join(
      gateTblChnl |>
        dplyr::select(gate, ind) |> # nolint
        dplyr::rename(cpOrig = gate),
      by = "ind"
    ) |>
    dplyr::left_join(
      gateStatsTbl |>
        dplyr::select(ind, freqBs, freqStim) |> # nolint
        dplyr::rename(
          freqOrig = freqStim,
          freqBsOrig = freqBs
        ),
      by = "ind"
    ) |>
    dplyr::left_join(gateSummStatTbl, by = "grp") |>
    dplyr::left_join(
      gateTblCtrl |>
        dplyr::select(ind, gate) |>
        dplyr::filter(!is.na(gate)) |>
        dplyr::rename(cpTgCtrl = gate),
      by = "ind"
    )
}

#' @keywords internal
.getCpClusterCpTblAddOrigQuantMin <- function(cpTbl) {
  .debug("Adding original and minimum quantile threshold") # nolint
  cpTbl |>
    dplyr::mutate(
      cpOrigQuantMin = pmax(pmin(cpOrig, maxExpr), gate05) # nolint
    )
}

#' @keywords internal
.getCpClusterCpJoinLseGet <- function(cpTbl) {
  .debug("Getting cpJoinLse") # nolint
  cpTbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      lseOrig = propBsCpDiffSd < 0.01, # nolint
      cpJoinLse = min(cp[lseOrig & cp >= cpJoin]), # nolint
      cpJoinLseOrig = pmin(cpJoinLse, cpOrigQuantMin) # nolint
    ) |>
    dplyr::ungroup()
}

#' @keywords internal
.getCpClusterCpJoinTgGet <- function(cpTbl) {
  .debug("Getting cpJoinTg") # nolint
  cpTbl |>
    dplyr::group_by(ind) |> # nolint
    dplyr::mutate(
      cpJoinTg = min(cp[
        cp >= cpJoin &
          cp >= cpTgCtrl & # nolint
          propBsCpDiffSd <= 2 # nolint
      ]),
      cpJoinTg = ifelse(
        is.na(cpJoinTg),
        cpJoinLse,
        cpJoinTg # nolint
      ),
      cpJoinTgOrig = pmin(cpJoinTg, cpOrigQuantMin) # nolint
    ) |>
    dplyr::ungroup()
}

#' @keywords internal
.getCpClusterCpLseOrigMean <- function(cpTbl) {
  .debug("Getting cpJoinLseOrigMean") # nolint
  cpTbl |>
    dplyr::mutate(
      cpJoinLseOrigMean = pmin(
        cpOrigQuantMin, # nolint
        cpJoinLseOrig + # nolint
          (cpOrigQuantMin - cpJoinLseOrig) / 2
      )
    )
}

#' @keywords internal
.getCpClusterCpJoinTgOrigMean <- function(cpTbl) {
  .debug("Getting cpJoinTgOrigMean") # nolint
  cpTbl |>
    dplyr::mutate(
      cpJoinTgOrigMean = pmin(
        cpOrigQuantMin, # nolint
        cpJoinTgOrig + (cpOrigQuantMin - cpJoinTgOrig) / 2 # nolint
      )
    )
}

#' @keywords internal
.getCpClusterCpJoinLseOrigMeanTg <- function(cpTbl) {
  .debug("Getting cpJoinLseOrigMean_tg") # nolint
  cpTbl |>
    dplyr::mutate(
      cpJoinLseOrigMeanTg = pmin(
        cpJoinLseOrigMean,
        cpJoinTgOrig # nolint
      )
    )
}

.getCpClusterCpFilterAboveCpJoinLseOrigMeanTg <-
  function(cpTbl) {
    .debug(
      "Filtering above cpJoinLseOrigMean_tg"
    ) # nolint
    cpTbl |>
      dplyr::filter(cp >= cpJoinLseOrigMeanTg) |> # nolint
      dplyr::group_by(ind) |> # nolint
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(
        grp,
        ind,
        cpOrigQuantMin, # nolint
        cpJoin,
        cpJoinLse,
        cpJoinLseOrig, # nolint
        cpJoinLseOrigMean, # nolint
        cpJoinTgOrig,
        cpJoinTgOrigMean, # nolint
        cpJoinLseOrigMeanTg,
        propBsOrig,
        propBsCpDiff, # nolint
        propBsCpDiffSd,
        propBsCp # nolint
      ) |>
      dplyr::ungroup() |>
      dplyr::group_by(grp) |>
      dplyr::arrange(dplyr::desc(cpOrigQuantMin)) |>
      dplyr::ungroup()
  }

#' @keywords internal
.getCpClusterCpImputeMissingBatch <- function(
  cpTbl,
  chnlCut,
  gateTbl,
  densTbl
) {
  .debug(
    "considering imputing missing thresholds by batch"
  ) # nolint
  indWithMissingGates <- setdiff(
    gateTbl$ind,
    cpTbl$ind[!is.na(cpTbl$cpJoinTgOrig)]
  )
  indWithMissingGates <- indWithMissingGates[
    !is.na(indWithMissingGates)
  ]
  if (length(indWithMissingGates) == 0L) {
    .debug("no missing thresholds by batch") # nolint
    return(cpTbl)
  }
  .debug("imputing missing thresholds by batch") # nolint
  for (indCurr in indWithMissingGates) {
    batch <- gateTbl |>
      dplyr::filter(ind == indCurr) |> # nolint
      dplyr::pull("batch")

    gateTblInd <- gateTbl |>
      dplyr::filter(ind == indCurr) # nolint
    batch <- gateTblInd |>
      dplyr::pull("batch")

    gateTblBatchIndVec <- gateTbl |>
      dplyr::filter(batch == .env$batch) |> # nolint
      dplyr::pull("ind")
    gateTblBatchIndVec <- gateTblBatchIndVec[
      !is.na(gateTblBatchIndVec)
    ]

    cpTblBatch <- cpTbl |>
      dplyr::filter(ind %in% gateTblBatchIndVec) # nolint

    densTblInd <- densTbl |>
      dplyr::filter(ind == indCurr) # nolint

    cpVecImp <- cpTblBatch$cpJoinTgOrig[
      !is.na(cpTblBatch$cpJoinTgOrig)
    ]

    if (length(cpVecImp) == 0) {
      next
    }

    gateImpute <- .combineCp(
      stats::setNames(cpVecImp, paste0("a", seq_along(cpVecImp))),
      gateTbl$gateCombn[1]
    )[[1]][[1]]

    cpTblAdd <- cpTbl[1, ] |>
      tidyr::pivot_longer(cols = -grp) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = ifelse(
          nrow(densTblInd) > 0,
          densTblInd$grp[1],
          NA
        ),
        ind = indCurr,
        cpJoinTgOrig = gateImpute
      )

    cpTbl <- cpTbl |> dplyr::filter(ind != indCurr) # nolint
    cpTbl <- cpTbl |> dplyr::bind_rows(cpTblAdd)
  }
  cpTbl <- cpTbl |> dplyr::arrange(ind) # nolint
  cpTbl
}

#' @keywords internal
.getCpClusterImputeMissingInd <- function(
  cpTbl,
  chnlCut,
  gateTbl,
  densTbl
) {
  .debug("considering imputing missing thresholds individually") # nolint
  indWithMissingGates <- setdiff(
    gateTbl$ind,
    cpTbl$ind[!is.na(cpTbl$cpJoinTgOrig)]
  )
  indWithMissingGates <- indWithMissingGates[
    !is.na(indWithMissingGates)
  ]
  if (length(indWithMissingGates) == 0L) {
    .debug("no missing thresholds individually") # nolint
    return(cpTbl)
  }
  .debug("imputing missing thresholds individually") # nolint
  for (indCurr in indWithMissingGates) {
    densTblInd <- densTbl |>
      dplyr::select(ind, grp) |> # nolint
      dplyr::filter(ind == indCurr)
    if (nrow(densTblInd) == 0) {
      next
    }
    gateImpute <- cpTbl |>
      dplyr::filter(grp == densTblInd$grp[1]) |> # nolint
      dplyr::pull("cpJoinTgOrig") |>
      quantile(0.75)

    cpTblAdd <- cpTbl[1, ] |>
      tidyr::pivot_longer(cols = -c(ind, grp)) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = densTblInd$grp[1],
        ind = densTblInd$ind[1],
        cpJoinTgOrig = gateImpute
      )

    cpTbl <- cpTbl |> dplyr::filter(ind != indCurr) # nolint
    cpTbl <- cpTbl |> dplyr::bind_rows(cpTblAdd)
  }
  cpTbl |> dplyr::arrange(ind) # nolint
}

#' @keywords internal
.getCpClusterImputeMissingFinal <- function(
  cpTbl,
  chnlCut,
  gateTbl,
  densTbl
) {
  .debug("considering imputing missing thresholds finally") # nolint
  indWithMissingGates <- setdiff(
    gateTbl$ind,
    cpTbl$ind[!is.na(cpTbl$cpJoinTgOrig)]
  )
  indWithMissingGates <- indWithMissingGates[
    !is.na(indWithMissingGates)
  ]
  if (length(indWithMissingGates) == 0L) {
    .debug("no missing thresholds finally") # nolint
    return(cpTbl)
  }
  .debug("imputing missing thresholds finally") # nolint
  for (indCurr in indWithMissingGates) {
    densTblInd <- densTbl |>
      dplyr::select(ind, grp) |> # nolint
      dplyr::filter(ind == indCurr)
    if (nrow(densTblInd) == 0) {
      next
    }
    gateImpute <- cpTbl |>
      dplyr::filter(grp == densTblInd$grp[1]) |> # nolint
      dplyr::pull("cpJoinTgOrig") |>
      quantile(0.75)

    cpTblAdd <- cpTbl[1, ] |>
      tidyr::pivot_longer(cols = -c(ind, grp)) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = densTblInd$grp[1],
        ind = densTblInd$ind[1],
        cpJoinTgOrig = gateImpute
      )

    cpTbl <- cpTbl |> dplyr::filter(ind != indCurr) # nolint
    cpTbl <- cpTbl |> dplyr::bind_rows(cpTblAdd)
  }
  cpTbl |> dplyr::arrange(ind) # nolint
}

#' @keywords internal
.getCpClusterImputeMissingFinalBatch <- function(
  cpTbl,
  chnlCut,
  gateTbl,
  densTbl
) {
  .debug(
    "considering imputing missing thresholds finally by batch"
  )
  indWithMissingGates <- setdiff(
    gateTbl$ind,
    cpTbl$ind[!is.na(cpTbl$cpJoinTgOrig)]
  )
  indWithMissingGates <- indWithMissingGates[
    !is.na(indWithMissingGates)
  ]
  if (length(indWithMissingGates) == 0L) {
    .debug("no missing thresholds finally by batch") # nolint
    return(cpTbl)
  }
  .debug("imputing missing thresholds finally by batch") # nolint
  for (indCurr in indWithMissingGates) {
    gateTblInd <- gateTbl |> dplyr::filter(ind == indCurr) # nolint
    batch <- gateTblInd |>
      dplyr::pull("batch")

    gateTblBatchIndVec <- gateTbl |>
      dplyr::filter(batch == .env$batch) |> # nolint
      dplyr::pull("ind")
    gateTblBatchIndVec <- gateTblBatchIndVec[
      !is.na(gateTblBatchIndVec)
    ]

    cpTblBatch <- cpTbl |>
      dplyr::filter(ind %in% gateTblBatchIndVec) # nolint

    densTblInd <- densTbl |>
      dplyr::filter(ind == indCurr) # nolint

    cpVecImp <- cpTblBatch$cpJoinTgOrig[
      !is.na(cpTblBatch$cpJoinTgOrig)
    ]

    if (length(cpVecImp) == 0) {
      next
    }

    gateImpute <- .combineCp(
      stats::setNames(cpVecImp, paste0("a", seq_along(cpVecImp))),
      gateTbl$gateCombn[1]
    )[[1]][[1]]

    cpTblAdd <- cpTbl[1, ] |>
      tidyr::pivot_longer(cols = -c(ind, grp)) |> # nolint
      dplyr::mutate(value = NA_real_) |>
      tidyr::pivot_wider(id_cols = grp) |>
      dplyr::mutate(
        grp = ifelse(nrow(densTblInd) > 0, densTblInd$grp[1], NA),
        ind = indCurr,
        cpJoinTgOrig = gateImpute
      )

    cpTbl <- cpTbl |> dplyr::filter(ind != indCurr) # nolint
    cpTbl <- cpTbl |> dplyr::bind_rows(cpTblAdd)
  }
  cpTbl <- cpTbl |> dplyr::arrange(ind) # nolint
}
