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
    gateTypeCytPosCalc,
    combnMatList,
    cytCombnVecList,
    pathProject) {
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
    filterOtherCytPos) {
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
    gateTypeCytPosCalc,
    combnMatList,
    cytCombnVecList,
    gateName,
    pathProject) {
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
    gateTypeCytPosCalc,
    combn,
    combnMatList,
    cytCombnVecList,
    indBatch) {
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
      gateTypeCytPosFilter = gateTypeCytPosFilter
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
    gateTypeCytPosCalc = gateTypeCytPosCalc
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
    gateTypeCytPosCalc) {
  exListStim <- exList[-1]
  exUns <- exList[[1]]
  nCellUns <- nrow(exUns) # nolint

  purrr::map_df(seq_along(exListStim), function(i) {
    .debug("i: ", i) # nolint

    ex <- exListStim[[i]]

    gateTblGnInd <- gateTblGn |>
      dplyr::filter(
        ind == attr(ex, "ind") # nolint
      )

    # The unstimulated cells are classified using the gates belonging
    # to the corresponding stimulated condition.
    gateTblGnIndUns <- gateTblGnInd |>
      dplyr::mutate(
        ind = attr(exUns, "ind")
      )

    # Calculate context-dependent positivity once for the stimulated sample.
    posCacheStim <- .getPosIndCache(
      ex = ex,
      gateTbl = gateTblGnInd,
      chnl = chnl
    )

    posByChnlStim <- .getPosIndByChnl(
      ex = ex,
      gateTbl = gateTblGnInd,
      chnl = chnl,
      gateTypeCytPos = gateTypeCytPosCalc,
      posCache = posCacheStim
    )

    # And once for the corresponding unstimulated sample.
    posCacheUns <- .getPosIndCache(
      ex = exUns,
      gateTbl = gateTblGnIndUns,
      chnl = chnl
    )

    posByChnlUns <- .getPosIndByChnl(
      ex = exUns,
      gateTbl = gateTblGnIndUns,
      chnl = chnl,
      gateTypeCytPos = gateTypeCytPosCalc,
      posCache = posCacheUns
    )

    combnTbl <- purrr::map_df(
      names(combnMatList),
      function(j) {
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
          posByChnlStim = posByChnlStim,
          posByChnlUns = posByChnlUns
        )
      }
    ) |>
      dplyr::mutate(
        nCellStim = nrow(ex),
        nCellUns = .env$nCellUns # nolint
      )

    combnTbl |>
      .getStatsBatchGnCombnNeg(chnl)
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
    posByChnlStim,
    posByChnlUns) {
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

    chnlPos <- chnl[
      combnMat[i, , drop = TRUE]
    ]

    chnlNeg <- chnl[
      setdiff(
        seq_along(chnl),
        combnMat[i, , drop = TRUE]
      )
    ]

    statTblGnInd[i, "countStim"] <- sum(
      .getPosIndCytCombn(
        ex = ex,
        gateTbl = gateTblGnInd,
        chnlPos = chnlPos,
        chnlNeg = chnlNeg,
        chnlAlt = NULL,
        gateTypeCytPos = gateTypeCytPosCalc,
        posByChnl = posByChnlStim
      )
    )

    statTblGnInd[i, "countUns"] <- sum(
      .getPosIndCytCombn(
        ex = exUns,
        gateTbl = gateTblGnInd,
        chnlPos = chnlPos,
        chnlNeg = chnlNeg,
        chnlAlt = NULL,
        gateTypeCytPos = gateTypeCytPosCalc,
        posByChnl = posByChnlUns
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
.getStatsBatchGnFilterMasks <- function(
    exList,
    gateTblGn,
    chnl,
    gateTypeCytPosFilter) {
  exUns <- exList[[1]]

  purrr::map(
    exList[-1],
    function(ex) {
      gateTblGnInd <- gateTblGn |>
        dplyr::filter(
          ind == attr(ex, "ind") # nolint
        )

      if (
        nrow(ex) == 0L ||
          nrow(gateTblGnInd) == 0L
      ) {
        return(
          list(
            stim = NULL,
            uns = NULL
          )
        )
      }

      list(
        stim = .getPosIndButSinglePosByChnl(
          ex = ex,
          gateTbl = gateTblGnInd,
          chnl = chnl,
          gateTypeCytPos = gateTypeCytPosFilter
        ),
        uns = .getPosIndButSinglePosByChnl(
          ex = exUns,
          gateTbl = gateTblGnInd,
          chnl = chnl,
          gateTypeCytPos = gateTypeCytPosFilter
        )
      )
    }
  )
}

#' @keywords internal
.getStatsBatchGnFilterOrNonCombn <- function(
    exList,
    indBatch,
    gateTblGn,
    gn,
    chnl,
    filterOtherCytPos,
    gateTypeCytPosFilter) {
  .debug("filtering or not working out combinations") # nolint

  exUns <- exList[[1]]

  filterMasks <- if (filterOtherCytPos) {
    .getStatsBatchGnFilterMasks(
      exList = exList,
      gateTblGn = gateTblGn,
      chnl = chnl,
      gateTypeCytPosFilter = gateTypeCytPosFilter
    )
  } else {
    NULL
  }

  purrr::map_df(
    chnl,
    function(chnlCurr) {
      .debug("chnlCurr: ", chnlCurr) # nolint

      statTblGnInd <- tibble::tibble(
        ind = indBatch[-1],
        gateName = gn,
        chnl = chnlCurr,
        countStim = NA,
        nCellStim = NA,
        countUns = NA,
        nCellUns = NA
      )

      for (j in seq_len(nrow(statTblGnInd))) {
        .debug("j: ", j) # nolint

        ex <- exList[[j + 1]]

        gateTblGnInd <- gateTblGn |>
          dplyr::filter(
            ind == attr(ex, "ind") # nolint
          )

        xStim <- ex[[chnlCurr]]

        if (
          filterOtherCytPos &&
            !is.null(filterMasks[[j]]$stim)
        ) {
          xStim <- xStim[
            !filterMasks[[j]]$stim[[chnlCurr]]
          ]
        }

        nothingToGate <-
          length(xStim) == 0L ||
            nrow(gateTblGnInd) == 0L ||
            all(is.na(xStim))

        if (nothingToGate) {
          .debug("filling in NAs") # nolint

          statTblGnInd[j, "countStim"] <- NA_integer_

          statTblGnInd[j, "nCellStim"] <- sum(
            (!is.na(xStim)) &
              (!is.nan(xStim))
          )

          statTblGnInd[j, "countUns"] <- NA_integer_
          statTblGnInd[j, "nCellUns"] <- nrow(exUns)

          next
        }

        gateGnIndChnl <- gateTblGnInd$gate[
          gateTblGnInd$chnl == chnlCurr
        ]

        statTblGnInd[j, "countStim"] <- sum(
          xStim > gateGnIndChnl
        )

        statTblGnInd[j, "nCellStim"] <- length(
          xStim
        )

        xUns <- exUns[[chnlCurr]]

        if (filterOtherCytPos) {
          xUns <- xUns[
            !filterMasks[[j]]$uns[[chnlCurr]]
          ]
        }

        statTblGnInd[j, "countUns"] <- sum(
          xUns > gateGnIndChnl
        )

        statTblGnInd[j, "nCellUns"] <- length(
          xUns
        )
      }

      statTblGnInd
    }
  )
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
