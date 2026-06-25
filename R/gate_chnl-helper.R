#' @keywords internal
.gateChnlDeleteOldGates <- function() {
  dirSave <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dirSave)) {
    return(invisible(FALSE))
  }
  unlink(dirSave, recursive = TRUE)
  invisible(TRUE)
}

# Get gates for each sample within each batch
#' @keywords internal
.gateChnlPreAdjGatesGate <- function(
  indBatchList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  message("getting pre-adjustment gates")
  purrr::map_df(seq_along(indBatchList), function(i) {
    .debug("indBatchList", i) # nolint

    # message progress
    if (i %% 50 == 0 || i == length(indBatchList)) {
      txt <- paste0("batch ", i, " of ", length(indBatchList))
      message(txt)
    }
    .gateBatch(
      # nolint
      .data = .data,
      indBatch = indBatchList[[i]],
      batch = names(indBatchList)[i],
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject
    ) |>
      dplyr::select(
        gateName,
        gateType,
        gateCombn, # nolint
        batch,
        ind,
        gate,
        everything() # nolint
      )
  })
}

#' @keywords internal
.gateChnlGetAdjGates <- function(
  gateTbl,
  gateTblParams,
  chnlSettings,
  .data,
  stage,
  pathProject,
  indBatchList,
  calcCytPosGates
) {
  if (stage == "init") {
    .gateChnlGetAdjGatesAll(
      # nolint
      gateTbl = gateTbl,
      .data = .data,
      pathProject = pathProject,
      stage = stage,
      indBatchList = indBatchList,
      chnlSettings = chnlSettings,
      calcCytPosGates = calcCytPosGates
    )
  } else if (stage == "single") {
    .gateChnlGateAdjGatesSingle(
      gateTbl = gateTbl,
      gateTblParams = gateTblParams,
      .data = .data,
      calcCytPosGates = TRUE,
      pathProject = pathProject,
      indBatchList = indBatchList,
      stage = stage,
      chnlSettings = chnlSettings
    )
  } else {
    stop("stage not recognized")
  }
}

#' @keywords internal
.gateChnlGetAdjGatesAll <- function(
  gateTbl,
  .data,
  pathProject,
  stage,
  indBatchList,
  chnlSettings,
  calcCytPosGates
) {
  if (!is.null(chnlSettings$tolClust)) {
    gateTblTgGate <- gateTbl |>
      dplyr::filter(gateUse == "tgClust") # nolint
  } else {
    gateTblTgGate <- NULL
  }

  gateTbl <- gateTbl |> dplyr::filter(gateUse == "gate") # nolint
  gateTbl <- gateTbl |> dplyr::select(-gateUse) # nolint

  # =========================
  # Cluster-based gating
  # =========================

  if (!is.null(chnlSettings$tolClust)) {
    pathDirStats <- .getStats(
      # nolint
      gateTbl = gateTbl |> dplyr::mutate(chnl = chnlSettings$chnlCut),
      gateName = NULL,
      chnl = chnlSettings$chnlCut,
      .data = .data,
      filterOtherCytPos = FALSE,
      combn = FALSE,
      gateTypeSinglePosCalc = "base",
      pathProject = pathProject,
      indBatchList = indBatchList,
      popGate = chnlSettings$popGate,
      tolClust = chnlSettings$tolClust
    )
    gateStatsTbl <- pathDirStats |>
      .readGateStats() # nolint

    gateTblCluster <- purrr::map_df(
      unique(gateTbl$gateName),
      function(gn) {
        gateTblCluster <- .getCpCluster(
          # nolint
          .data = .data,
          gateTbl = gateTbl |>
            dplyr::filter(gateName == gn), # nolint
          gateStatsTbl = gateStatsTbl |>
            dplyr::filter(gateName == gn),
          gateTblCtrl = gateTblTgGate,
          chnlSettings = chnlSettings,
          filterOtherCytPos = FALSE,
          stage = stage,
          pathProject = pathProject,
          calcCytPosGates = calcCytPosGates,
          indBatchList = indBatchList
        )

        gateTblCluster |>
          dplyr::select(ind, cpJoinTgOrig) |> # nolint
          dplyr::rename(gate = cpJoinTgOrig) |>
          dplyr::left_join(
            gateTbl |>
              dplyr::filter(gateName == gn) |> # nolint
              dplyr::select(
                gateName,
                gateType,
                gateCombn, # nolint
                batch,
                ind # nolint
              ),
            by = c("ind")
          ) |>
          dplyr::select(
            gateName,
            gateType,
            gateCombn,
            batch,
            ind,
            gate # nolint
          ) |>
          dplyr::mutate(
            gateCombn = paste0(gateCombn, "Clust"),
            gateName = paste0(gateType, gateCombn)
          )
      }
    )
    gateTbl <- gateTbl |>
      dplyr::bind_rows(gateTblCluster)
  }
  # Output
  # ------------------

  list(
    gateTbl = gateTbl
  )
}

#' @keywords internal
.gateChnlGateAdjGatesSingle <- function(
  gateTbl,
  gateTblParams,
  .data,
  calcCytPosGates,
  pathProject,
  indBatchList,
  stage,
  chnlSettings
) {
  # get gates
  gateTblSingle <- gateTbl

  # merge
  gateTbl <- .gateChnlGateAdjGatesSingleMerge(
    # nolint
    gateTblSingle = gateTblSingle,
    gateTblParams = gateTblParams,
    chnlCut = chnlSettings$chnlCut
  )

  # get stats table (if needed)
  gateStatsTbl <- .gateChnlGateAdjGatesSingleStatsTblGet(
    gateTbl = gateTbl,
    .data = .data,
    chnlSettings = chnlSettings,
    calcCytPosGates = calcCytPosGates,
    pathProject = pathProject,
    indBatchList = indBatchList
  )

  gateTblOut <- .gateChnlGateAdjGatesSingleOutGet(
    gateTbl = gateTbl,
    gateStatsTbl = gateStatsTbl,
    gateTblSingle = gateTblSingle,
    gateQuant = gateQuant,
    .data = .data,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )

  list(
    gateTbl = gateTblOut
  )
}

#' @keywords internal
.gateChnlGateAdjGatesSingleMerge <- function(
  gateTblSingle,
  gateTblParams,
  chnlCut
) {
  gateTblParams |>
    dplyr::left_join(
      gateTblSingle |>
        dplyr::mutate(chnl = chnlCut) |>
        dplyr::rename(gateSingle = gate) |> # nolint
        dplyr::filter(gateUse == "gate") |> # nolint
        dplyr::select(
          chnl,
          gateName,
          ind,
          gateSingle # nolint
        ),
      by = c("chnl", "ind", "gateName")
    ) |>
    dplyr::mutate(
      gateType = purrr::map_chr(
        gateName,
        function(gn) {
          stringr::str_remove(gn, "Adj") |>
          stringr::str_remove("Clust")
        }
      ),
      gateCombn = gateName |>
        stringr::str_remove("Adj") |>
        stringr::str_remove("Clust") |>
        stringr::str_remove(gateType) |> # nolint
        stringr::str_remove("_")
    ) |>
    dplyr::select(
      chnl,
      marker,
      gateName, # nolint
      gateType,
      gateCombn,
      everything() # nolint
    ) |>
    dplyr::mutate(
      gateSingle = ifelse(is.na(gateSingle), gate, gateSingle) # nolint
    )
        }

#' @keywords internal
.gateChnlGateAdjGatesSingleStatsTblGet <- function(
  gateTbl,
  chnlSettings,
  .data,
  calcCytPosGates,
  pathProject,
  indBatchList
) {
  gateNameVec <- unique(gateTbl$gateName)
  if (!.gateChnlGateAdjGatesSingleStatsTblGetCheck(gateNameVec)) {
    return(NULL)
  }

  gateNameVecClust <- gateNameVec[
    stringr::str_detect(gateNameVec, "Clust")
  ]
  gateNameVecAdj <- gateNameVec[
    stringr::str_detect(gateNameVec, "Adj")
  ]

  .getStats(
    # nolint
    gateTbl = gateTbl |>
      dplyr::filter(
        gateName %in% c(gateNameVecClust, gateNameVecAdj) # nolint
      ),
    gateName = NULL,
    chnl = chnlSettings$chnlCut,
    filterOtherCytPos = TRUE,
    gateTypeCytPosFilter = if (calcCytPosGates) "cyt" else "base",
    .data = .data,
    gateTypeSinglePosFilter = "base",
    gateTypeSinglePosCalc = "base",
    combn = FALSE,
    pathProject = pathProject,
    indBatchList = indBatchList,
    popGate = chnlSettings$chnlCutpopGate
  )
}

#' @keywords internal
.gateChnlGateAdjGatesSingleStatsTblGetCheck <- function(
  gateNameVec
) {
  # nolint
  anyClustInd <- any(
    purrr::map_lgl(
      gateNameVec,
      function(gn) stringr::str_detect(gn, "Clust")
    )
  )
  anyAdjInd <- any(
    purrr::map_lgl(
      gateNameVec,
      function(gn) stringr::str_detect(gn, "Adj")
    )
  )
  anyClustInd || anyAdjInd
}

#' @keywords internal
.gateChnlGateAdjGatesSingleOutGet <- function(
  gateTbl,
  gateStatsTbl,
  gateTblSingle,
  gateQuant,
  .data,
  stage,
  pathProject,
  chnlSettings
) {
  # get tail-gate gates
  gateTblCtrlClust <- gateTblSingle |>
    dplyr::filter(gateUse == "tgClust") # nolint

  # get control gates
  gateTblCtrlCtrl <- gateTblSingle |> # nolint
    dplyr::filter(gateUse == "ctrl") # nolint

  # get gate names
  gateNameVec <- unique(gateTbl$gateName)
  purrr::map_df(gateNameVec, function(gn) {
    .gateChnlGateAdjGatesSingleOutGetGn(
      gn = gn,
      gateTbl = gateTbl,
      gateStatsTbl = gateStatsTbl,
      gateTblCtrlClust = gateTblCtrlClust,
      gateTblCtrlCtrl = gateTblCtrlCtrl,
      gateTblSingleGn = gateTblSingle |>
        dplyr::filter(gateName == gn), # nolint
      gateQuant = gateQuant,
      .data = .data,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    )
  }) |>
    dplyr::mutate(gateSingle = pmax(gate, gateSingle)) # nolint
}

#' @keywords internal
.gateChnlGateAdjGatesSingleOutGetGn <- function(
  gn,
  gateTbl,
  gateStatsTbl,
  gateTblCtrlClust,
  gateTblCtrlCtrl,
  gateTblSingleGn,
  gateQuant,
  .data,
  stage,
  pathProject,
  chnlSettings
) {
  .debug("gateNameVec", gn) # nolint
  gateTblGn <- gateTbl |>
    dplyr::filter(gateName == gn) # nolint

  adjInd <- stringr::str_detect(gn, "Adj")

  clustInd <- stringr::str_detect(gn, "Clust")
  if (!clustInd && !adjInd) {
    return(
      gateTblGn |>
        dplyr::filter(chnl == chnlSettings$chnlCut) |> # nolint
        dplyr::select(
          chnl,
          marker,
          gateName,
          gateType, # nolint
          gateCombn,
          batch,
          ind,
          gate,
          gateCyt,
          gateSingle # nolint
        )
    )
  }

  # =========================
  # Tail-gate controlled gating
  # =========================

  if (adjInd) {
    gateStatsTblGn <- gateStatsTbl |>
      dplyr::filter(gateName == gn) # nolint
    gateTblCtrlCtrl <- gateTblCtrlCtrl |>
      dplyr::filter(gateName == gn) # nolint

    gateTblGn2 <- .getCpAdjTbl(
      # nolint
      gateStatsTbl = gateStatsTblGn,
      gateQuant = gateQuant,
      gateTblCtrl = gateTblCtrlCtrl
    )
    gateTblGn2 <- gateTblGn2 |>
      dplyr::left_join(
        gateTblSingleGn,
        by = c("ind")
      ) |>
      dplyr::rename(gateSingle = gate) # nolint

    gateTblSingleGn <- gateTblGn2
    return(gateTblSingleGn)
  }

  # =========================
  # Cluster-based gating
  # =========================

  if (clustInd) {
    gateStatsTblGn <- gateStatsTbl |>
      dplyr::filter(gateName == gn) # nolint
    gateTblCtrlClustGn <- gateTblCtrlClust |>
      dplyr::filter(gateName == gn) # nolint

    gateTblClusterGn <- .getCpCluster(
      # nolint
      .data = .data,
      gateTbl = gateTblGn,
      gateStatsTbl = gateStatsTblGn,
      gateTblCtrl = gateTblCtrlClustGn,
      filterOtherCytPos = TRUE,
      stage = stage,
      pathProject = pathProject
    )

    gateTblClusterGn <- gateTblClusterGn |>
      dplyr::select(ind, cpJoinTgOrig) |> # nolint
      dplyr::rename(gateSingle = cpJoinTgOrig) |>
      dplyr::left_join(
        gateTbl |>
          dplyr::filter(
            gateName == gn, # nolint
            chnl == chnlSettings$chnlCut # nolint
          ) |>
          dplyr::select(
            gateName,
            gateType,
            gateCombn, # nolint
            batch,
            ind,
            gate,
            gateCyt # nolint
          ),
        by = c("ind")
      ) |>
      dplyr::mutate(
        chnl = chnlSettings$chnlCut,
        marker = chnlSettings$marker
      ) |>
      dplyr::select(
        chnl,
        marker,
        gateName,
        gateType,
        gateCombn, # nolint
        batch,
        ind,
        gate,
        gateCyt,
        gateSingle # nolint
      ) |>
      dplyr::mutate(
        gateCombn = paste0(gateCombn, "Clust"),
        gateName = paste0(gateType, gateCombn)
      )
  }

  gateTblClusterGn
}

# Placeholder function for adjustment table - this may need proper implementation
#' @keywords internal
.getCpAdjTbl <- function(gateStatsTbl, gateQuant, gateTblCtrl) {
  uniqueInd <- unique(gateStatsTbl$ind)
  tibble::tibble(
    ind = uniqueInd,
    gate = rep(0.5, length(uniqueInd)) # placeholder gate values
  )
}
