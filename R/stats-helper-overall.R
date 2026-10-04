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

  if (!combn) {
    statTbl <- statTbl |>
      dplyr::mutate(marker = chnlLab[.data$chnl]) # nolint
  }

  if ("ind" %in% colnames(statTbl)) {
    statTbl[, "ind"] <- as.character(statTbl[["ind"]])
  }

  statTbl |>
    dplyr::select(
      dplyr::any_of(c("gateName", "chnl", "marker", "ind")),
      dplyr::everything()
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
  gateTypeCytPosCalc,
  combnMatList,
  cytCombnVecList,
  gateName,
  pathProject
) {
  .debug("Getting gate stats for a batch") # nolint
  .debug("indBatch: ", paste0(indBatch, collapse = "-")) # nolint

  if (combn && !filterOtherCytPos) {
    return(.getStatsBatchCombn(
      indBatch = indBatch, batch = batch, .data = .data, popGate = popGate,
      gateTbl = gateTbl, gateName = gateName, chnl = chnl,
      combnMatList = combnMatList, cytCombnVecList = cytCombnVecList,
      gateTypeCytPosCalc = gateTypeCytPosCalc, pathProject = pathProject
    ))
  }

  exList <- .getExList(
    .data = .data,
    indBatch = indBatch,
    batch = batch,
    pop = popGate,
    chnlCut = unique(gateTbl$chnl),
    pathProject = pathProject
  )

  purrr::map_df(gateName, function(gn) {
    .debug("gate name: ", gn) # nolint
    .getStatsBatchGnFilterOrNonCombn(
      exList = exList,
      indBatch = indBatch,
      gateTblGn = gateTbl |> dplyr::filter(gateName == gn), # nolint
      gn = gn,
      chnl = chnl,
      filterOtherCytPos = filterOtherCytPos,
      gateTypeCytPosFilter = gateTypeCytPosFilter
    )
  })
}

#' @keywords internal
.getStatsBatchCombn <- function(
  indBatch,
  batch,
  .data,
  popGate,
  gateTbl,
  gateName,
  chnl,
  combnMatList,
  cytCombnVecList,
  gateTypeCytPosCalc,
  pathProject
) {
  if (length(chnl) > 30L) {
    stop("Combination statistics support at most 30 channels.")
  }
  # Preserve combination-size order, then the original combn() row order.
  bitIndex <- unlist(lapply(combnMatList, function(mat) {
    as.integer(rowSums(2^(mat - 1L)))
  }), use.names = FALSE)
  cytCombn <- c(
    unlist(cytCombnVecList, use.names = FALSE),
    paste0(chnl, "~-~", collapse = "")
  )
  bitIndex <- c(bitIndex, 0L)
  # The old batch loader obtained tube sizes from the gated expression channels.
  nCellChnl <- if (nrow(gateTbl) > 0L) unique(gateTbl$chnl)[[1]] else chnl[[1]]

  purrr::map_df(gateName, function(gn) {
    .debug("gate name: ", gn)
    gateTblGn <- gateTbl |> dplyr::filter(gateName == gn) # nolint
    purrr::map_df(seq_along(indBatch[-1]), function(i) {
      .debug("i: ", i)
      ind <- indBatch[[i + 1L]]
      gates <- gateTblGn |> dplyr::filter(.data$ind == .env$ind)
      stim <- .getStatsCombnTube(
        .data, ind, indBatch[[1]], batch, popGate, pathProject,
        chnl, nCellChnl, gates, gateTypeCytPosCalc, combnMatList, bitIndex
      )
      # Re-read raw unstim expression for each stim sample's gates. Extra disk
      # reads keep memory bounded to one tube's logical comparisons and codes,
      # rather than retaining a full unstim double-expression table per batch.
      uns <- .getStatsCombnTube(
        .data, indBatch[[1]], indBatch[[1]], batch, popGate, pathProject,
        chnl, nCellChnl, gates, gateTypeCytPosCalc, combnMatList, bitIndex
      )
      tibble::tibble(
        ind = as.character(ind), gateName = gn, cytCombn = cytCombn,
        countStim = stim$count, nCellStim = stim$n,
        countUns = uns$count, nCellUns = uns$n
      )
    })
  })
}

#' @keywords internal
.getStatsCombnTube <- function(
  .data,
  ind,
  indUns,
  batch,
  popGate,
  pathProject,
  chnl,
  nCellChnl,
  gateTbl,
  gateTypeCytPos,
  combnMatList,
  bitIndex
) {
  gateTypeCytPos <- match.arg(gateTypeCytPos, c("base", "cyt"))
  .readChnl <- function(chnlCurr) {
    .getEx(
      .data = if (is.null(.data)) NULL else .data[[ind]],
      pop = popGate, chnlCut = chnlCurr, ind = ind, indUns = indUns,
      batch = batch, pathProject = pathProject
    )
  }
  ex <- .readChnl(nCellChnl)
  n <- nrow(ex)
  rm(ex)
  if (nrow(gateTbl) == 0L) {
    return(list(count = rep(NA_integer_, length(bitIndex)), n = n))
  }

  posCache <- NULL
  code <- integer(n)
  # Retain only logical comparisons for cyt+ context and the NA fallback;
  # discard each channel's double expression immediately after comparison.
  for (k in seq_along(chnl)) {
    ex <- if (chnl[[k]] %in% gateTbl$chnl) {
      .readChnl(chnl[[k]])
    } else {
      tibble::tibble(.rows = n)
    }
    posCache <- .getPosIndCache(ex, gateTbl, chnl[[k]], posCache)
    if (gateTypeCytPos == "base") {
      code <- code + as.integer(posCache$base[[chnl[[k]]]]) *
        bitwShiftL(1L, k - 1L)
    }
    rm(ex)
  }
  # Only nrow(ex) is used when the logical cache is already complete.
  ex <- tibble::tibble(.rows = n)
  posByChnl <- .getPosIndByChnl(
    ex, gateTbl, chnl, gateTypeCytPos, posCache
  )
  if (any(vapply(posByChnl, anyNA, logical(1)))) {
    .debug("Combination statistics: using NA-preserving Reduce fallback")
    # Logical AND/OR can resolve some NA cells to FALSE. Tabulation cannot
    # reproduce that per-combination behaviour, so use the original helper.
    count <- unlist(lapply(combnMatList, function(mat) {
      vapply(seq_len(nrow(mat)), function(i) {
        chnlPos <- chnl[mat[i, , drop = TRUE]]
        as.integer(sum(.getPosIndCytCombn(
          ex = ex, gateTbl = gateTbl, chnlPos = chnlPos,
          chnlNeg = setdiff(chnl, chnlPos), gateTypeCytPos = gateTypeCytPos,
          posByChnl = posByChnl
        )))
      }, integer(1))
    }), use.names = FALSE)
    return(list(count = c(count, as.integer(n - sum(count))), n = n))
  }
  if (gateTypeCytPos == "cyt") {
    for (k in seq_along(chnl)) {
      code <- code + as.integer(posByChnl[[chnl[[k]]]]) *
        bitwShiftL(1L, k - 1L)
    }
  }
  rm(posCache, posByChnl, ex)
  count <- tabulate(code + 1L, nbins = 2^length(chnl))
  list(count = count[bitIndex + 1L], n = n)
}

#' @keywords internal
.getStatsBatchGnFilterMasks <- function(
  exList,
  gateTblGn,
  chnl,
  gateTypeCytPosFilter
) {
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
  gateTypeCytPosFilter
) {
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

          statTblGnInd[j, "nCellStim"] <- sum(!is.na(xStim))

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
