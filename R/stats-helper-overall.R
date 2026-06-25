#' @keywords internal
.getStatsOverall <- function(
  indBatchList,
  .data,
  popGate,
  gateTbl,
  gateName,
  chnl,
  chnlLab,
  filterOtherCytPos,
  combn,
  gateTypeCytPosFilter,
  gateTypeSinglePosFilter,
  gateTypeSinglePosCalc,
  gateTypeCytPosCalc,
  combnMatList,
  cytCombnVecList,
  pathProject
) {
  statTbl <- purrr::map_df(
    seq_along(indBatchList),
    function(i) {
      .getStatsOverallProgress(
        indBatchList = indBatchList,
        i = i,
        combn = combn,
        filterOtherCytPos = filterOtherCytPos
      )
      .getStatsBatch(
        indBatch = indBatchList[[i]],
        batch = names(indBatchList)[i],
        .data = .data,
        popGate = popGate,
        gateTbl = gateTbl,
        chnl = chnl,
        filterOtherCytPos = filterOtherCytPos,
        combn = combn,
        gateTypeCytPosFilter = gateTypeCytPosFilter,
        gateTypeSinglePosFilter = gateTypeSinglePosFilter,
        gateTypeSinglePosCalc = gateTypeSinglePosCalc,
        gateTypeCytPosCalc = gateTypeCytPosCalc,
        combnMatList = combnMatList,
        cytCombnVecList = cytCombnVecList,
        gateName = gateName,
        pathProject = pathProject
      )
    }
  )

  statTbl <- statTbl |>
    dplyr::mutate(
      propStim = countStim / nCellStim, # nolint
      propUns = countUns / nCellUns, # nolint
      propBs = propStim - propUns, # nolint
      freqStim = propStim * 1e2, # nolint
      freqUns = propUns * 1e2, # nolint
      freqBs = freqStim - freqUns # nolint
    )

  statTbl <- .getStatsUpdateCombnN(
    combn = combn,
    statTbl = statTbl,
    chnlCut = chnl[[1]],
    chnlLab = chnlLab
  )

  if ("ind" %in% colnames(statTbl)) {
    statTbl[, "ind"] <- as.character(statTbl[["ind"]])
  }

  .getStatsLabel(
    statTbl = statTbl
  )
}

#' @keywords internal
.getStatsOverallProgress <- function(
  indBatchList,
  i,
  combn,
  filterOtherCytPos
) {
  indBatch <- indBatchList[[i]]
  .debug(
    "indBatch: ",
    paste0(indBatch, collapse = "-")
  )
  if (i %% 10 == 0 || i == length(indBatchList)) {
    if (combn && !filterOtherCytPos) {
      txt <- paste0("batch ", i, " of ", length(indBatchList))
      message(txt)
    }
  }
}

#' @keywords internal
.getStatsBatch <- function(
  indBatch,
  batch,
  .data,
  popGate,
  gateTbl,
  chnl,
  filterOtherCytPos,
  combn,
  gateTypeCytPosFilter,
  gateTypeSinglePosFilter,
  gateTypeSinglePosCalc,
  gateTypeCytPosCalc,
  combnMatList,
  cytCombnVecList,
  gateName,
  pathProject
) {
  .debug("Getting gate stats for a batch") # nolint
  .debug("indBatch: ", paste0(indBatch, collapse = "-")) # nolint

  exList <- .getExList(
    .data = .data,
    indBatch = indBatch,
    batch = batch,
    pop = popGate,
    chnlCut = unique(gateTbl$chnl),
    pathProject = pathProject
  )

  purrr::map_df(gateName, function(gn) {
    .getStatsBatchGn(
      gn = gn,
      exList = exList,
      gateTbl = gateTbl,
      chnl = chnl,
      filterOtherCytPos = filterOtherCytPos,
      gateTypeCytPosFilter = gateTypeCytPosFilter,
      gateTypeSinglePosFilter = gateTypeSinglePosFilter,
      gateTypeSinglePosCalc = gateTypeSinglePosCalc,
      gateTypeCytPosCalc = gateTypeCytPosCalc,
      combn = combn,
      combnMatList = combnMatList,
      cytCombnVecList = cytCombnVecList,
      indBatch = indBatch
    )
  })
}

#' @keywords internal
.getStatsBatchGn <- function(
  gn,
  exList,
  gateTbl,
  chnl,
  filterOtherCytPos,
  gateTypeCytPosFilter,
  gateTypeSinglePosFilter,
  gateTypeSinglePosCalc,
  gateTypeCytPosCalc,
  combn,
  combnMatList,
  cytCombnVecList,
  indBatch
) {
  .debug("gate name: ", gn) # nolint
  gateTblGn <- gateTbl |> dplyr::filter(gateName == gn) # nolint
  if (filterOtherCytPos || !combn) {
    statTblGn <- .getStatsBatchGnFilterOrNonCombn(
      exList = exList,
      indBatch = indBatch,
      gateTblGn = gateTblGn,
      gn = gn,
      chnl = chnl,
      filterOtherCytPos = filterOtherCytPos,
      gateTypeSinglePosCalc = gateTypeSinglePosCalc,
      gateTypeCytPosFilter = gateTypeCytPosFilter,
      gateTypeSinglePosFilter = gateTypeSinglePosFilter
    )
    return(statTblGn)
  }

  .getStatsBatchGnCombnLoopInd(
    exList = exList,
    gateTblGn = gateTblGn,
    gn = gn,
    chnl = chnl,
    combnMatList = combnMatList,
    cytCombnVecList = cytCombnVecList,
    gateTypeCytPosCalc = gateTypeCytPosCalc,
    gateTypeSinglePosCalc = gateTypeSinglePosCalc
  )
}

#' @keywords internal
.getStatsBatchGnCombnLoopInd <- function(
  exList,
  gateTblGn,
  gn,
  chnl,
  combnMatList,
  cytCombnVecList,
  gateTypeCytPosCalc,
  gateTypeSinglePosCalc
) {
  exListStim <- exList[-length(exList)]
  exUns <- exList[[length(exList)]]
  nCellUns <- nrow(exUns) # nolint
  
  purrr::map_df(seq_along(exListStim), function(i) {
    .debug("i: ", i) # nolint
    ex <- exListStim[[i]]
    gateTblGnInd <- gateTblGn |>
      dplyr::filter(ind == attr(ex, "ind")) # nolint
    combnTbl <- purrr::map_df(names(combnMatList), function(j) {
      .getStatsBatchGnCombn(
        j = j,
        ex = ex,
        exUns = exUns,
        gateTblGnInd = gateTblGnInd,
        gn = gn,
        chnl = chnl,
        combnMatList = combnMatList,
        cytCombnVecList = cytCombnVecList,
        gateTypeCytPosCalc = gateTypeCytPosCalc,
        gateTypeSinglePosCalc = gateTypeSinglePosCalc
      )
    }) |>
      dplyr::mutate(
        nCellStim = nrow(ex),
        nCellUns = .env$nCellUns # nolint
      )
    combnTbl |> .getStatsBatchGnCombnNeg(chnl)
  })
}

#' @keywords internal
.getStatsBatchGnCombn <- function(
  j,
  ex,
  exUns,
  gateTblGnInd,
  gn,
  chnl,
  combnMatList,
  cytCombnVecList,
  gateTypeCytPosCalc,
  gateTypeSinglePosCalc
) {
  .debug("number of cytokines positive: ", j) # nolint
  combnMat <- combnMatList[[j]]
  cytCombn <- cytCombnVecList[[j]]
  statTblGnInd <- tibble::tibble(
    ind = attr(ex, "ind"),
    gateName = gn,
    cytCombn = cytCombn,
    countStim = NA_integer_,
    nCellStim = NA_integer_,
    countUns = NA_integer_,
    nCellUns = NA_integer_
  )

  for (i in seq_len(nrow(statTblGnInd))) {
    .debug("i: ", i) # nolint
    chnlPos <- chnl[combnMat[i, , drop = TRUE]]
    chnlNeg <- chnl[
      setdiff(seq_along(chnl), combnMat[i, , drop = TRUE])
    ]
    statTblGnInd[i, "countStim"] <- sum(
      .getPosIndCytCombn(
        ex = ex,
        gateTbl = gateTblGnInd,
        chnlPos = chnlPos,
        chnlNeg = chnlNeg,
        chnlAlt = NULL,
        gateTypeCytPos = gateTypeCytPosCalc,
        gateTypeSinglePos = gateTypeSinglePosCalc
      )
    )
    statTblGnInd[i, "countUns"] <- sum(
      .getPosIndCytCombn(
        ex = exUns,
        gateTbl = gateTblGnInd |>
          dplyr::mutate(ind = attr(exUns, "ind")),
        chnlPos = chnlPos,
        chnlNeg = chnlNeg,
        chnlAlt = NULL,
        gateTypeCytPos = gateTypeCytPosCalc,
        gateTypeSinglePos = gateTypeSinglePosCalc
      )
    )
  }
  statTblGnInd
}

#' @keywords internal
.getStatsBatchGnCombnNeg <- function(.data, chnl) {
  allNegRow <- .data |>
    dplyr::mutate(cytCombn = paste0(paste0(chnl, collapse = "~-~"), "~-~")) |>
    dplyr::group_by(ind, cytCombn, gateName) |>
    dplyr::summarise(
      countStim = nCellStim[[1]] - sum(countStim),
      nCellStim = nCellStim[[1]],
      countUns = nCellUns[[1]] - sum(countUns),
      nCellUns = nCellUns[[1]],
      .groups = "drop"
    )
  .data |> dplyr::bind_rows(allNegRow)
}

#' @keywords internal
.getStatsBatchGnFilterOrNonCombn <- function(
  exList,
  indBatch,
  gateTblGn,
  gn,
  chnl,
  filterOtherCytPos,
  gateTypeSinglePosCalc,
  gateTypeCytPosFilter,
  gateTypeSinglePosFilter
) {
  .debug("filtering or not working out combinations") # nolint
  purrr::map_df(chnl, function(chnlCurr) {
    .debug("chnlCurr: ", chnlCurr) # nolint

    if (filterOtherCytPos) {
      exList <- .getStatsBatchGnFilterOrNonCombnFilter(
        exList = exList,
        gateTblGn = gateTblGn,
        chnlCurr = chnlCurr,
        gateTypeCytPosFilter = gateTypeCytPosFilter,
        gateTypeSinglePosFilter = gateTypeSinglePosFilter
      )
    }

    statTblGnInd <- tibble::tibble(
      ind = indBatch[-length(indBatch)],
      gateName = gn,
      chnl = chnlCurr,
      countStim = NA,
      nCellStim = NA,
      countUns = NA,
      nCellUns = NA
    )
    for (j in seq_len(nrow(statTblGnInd))) {
      .debug("j: ", j) # nolint
      ex <- exList[[j]]
      gateTblGnInd <- gateTblGn |>
        dplyr::filter(ind == attr(ex, "ind")) # nolint
      nothingToGate <- nrow(ex) == 0 ||
        nrow(gateTblGnInd) == 0 ||
        all(is.na(ex[[chnlCurr]]))
      if (nothingToGate) {
        .debug("filling in NAs") # nolint
        statTblGnInd[j, "countStim"] <- NA_integer_
        statTblGnInd[j, "nCellStim"] <- min(
          sum((!is.na(ex[[chnlCurr]])) & (!is.nan(ex[[chnlCurr]])))
        )
        statTblGnInd[j, "countUns"] <- NA_integer_
        statTblGnInd[j, "nCellUns"] <- nrow(exList[[length(exList)]])
        next
      }
      cnVec <- colnames(gateTblGnInd)
      gateColInd <- switch(
        gateTypeSinglePosCalc,
        "base" = which(cnVec == "gate"),
        "single" = which(cnVec == "gateSingle"),
        stop(paste0(
          "gateTypeSinglePosCalc value of ",
          gateTypeSinglePosCalc,
          " not either 'base' or 'single'."
        ))
      )
      gateGnIndChnl <- gateTblGnInd[[gateColInd]][gateTblGnInd$chnl == chnlCurr]
      statTblGnInd[j, "countStim"] <- sum(ex[[chnlCurr]] > gateGnIndChnl)
      statTblGnInd[j, "nCellStim"] <- nrow(ex)
      
      exUns <- exList[[length(exList)]]
      if (filterOtherCytPos) {
        posIndVecButSinglePosCurr <- .getPosIndButSinglePosForOneCyt(
          ex = exUns |> dplyr::mutate(isUns = FALSE),
          gateTbl = gateTblGnInd |> dplyr::mutate(ind = attr(ex, "ind")),
          chnlSingleExc = chnlCurr,
          chnl = NULL,
          gateTypeCytPos = gateTypeCytPosFilter,
          gateTypeSinglePos = gateTypeSinglePosFilter
        )
        exUns <- exUns[!posIndVecButSinglePosCurr, , drop = FALSE]
      }
      statTblGnInd[j, "countUns"] <- sum(exUns[[chnlCurr]] > gateGnIndChnl)
      statTblGnInd[j, "nCellUns"] <- nrow(exUns)
    }
    statTblGnInd
  })
}

#' @keywords internal
.getStatsBatchGnFilterOrNonCombnFilter <- function(
  exList,
  gateTblGn,
  chnlCurr,
  gateTypeCytPosFilter,
  gateTypeSinglePosFilter
) {
  .debug("filtering other cyt pos") # nolint

  purrr::map(seq_along(exList), function(i) {
    .debug("i: ", i) # nolint

    returnEarly <- .getStatsBatchGnFilterOrNonCombnFilterCheckEarly(
      i = i,
      exList = exList,
      chnlCurr = chnlCurr,
      gateTblGn = gateTblGn
    )
    if (returnEarly) {
      return(exList[[i]])
    }

    posIndVecButSinglePosCurr <- .getPosIndButSinglePosForOneCyt(
      ex = exList[[i]],
      gateTbl = gateTblGn |>
        dplyr::filter(ind == attr(exList[[i]], "ind")), # nolint
      chnlSingleExc = chnlCurr,
      chnl = NULL,
      gateTypeCytPos = gateTypeCytPosFilter,
      gateTypeSinglePos = gateTypeSinglePosFilter
    )
    exList[[i]][!posIndVecButSinglePosCurr, , drop = FALSE]
  }) |>
    stats::setNames(names(exList))
}

#' @keywords internal
.getStatsBatchGnFilterOrNonCombnFilterCheckEarly <- function(
  i,
  exList,
  chnlCurr,
  gateTblGn
) {
  if (i == length(exList)) {
    return(TRUE)
  }
  gateTblGnInd <- gateTblGn |>
    dplyr::filter(ind == attr(exList[[i]], "ind")) # nolint

  all(is.na(exList[[i]][[chnlCurr]])) ||
    nrow(gateTblGnInd) == 0 ||
    nrow(exList[[i]]) == 0
}

#' @keywords internal
.getStatsUpdateCombnN <- function(combn, statTbl, chnlCut, chnlLab) {
  if (combn) {
    return(statTbl)
  }
  if ((!"chnl" %in% colnames(statTbl))) {
    statTbl <- statTbl |>
      dplyr::mutate(
        chnl = chnlCut,
        marker = chnlLab[chnlCut]
      )
  }
  if ((!"marker" %in% colnames(statTbl))) {
    statTbl <- statTbl |>
      dplyr::mutate(marker = chnlLab[.data$chnl]) # nolint
  }
  statTbl
}

#' @keywords internal
.getStatsLabel <- function(statTbl) {
  cnVecOrder <- c(
    "gateName",
    "chnl",
    "marker",
    "ind"
  )
  cnVecOrderCurr <- cnVecOrder[cnVecOrder %in% colnames(statTbl)]
  statTbl |>
    dplyr::select(dplyr::all_of(cnVecOrderCurr), dplyr::everything()) # nolint
}
