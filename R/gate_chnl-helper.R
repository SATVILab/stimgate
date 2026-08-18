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
        dplyr::everything() # nolint
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
  .gateChnlGetAdjGatesAll(
    gateTbl = gateTbl,
    .data = .data,
    pathProject = pathProject,
    stage = stage,
    indBatchList = indBatchList,
    chnlSettings = chnlSettings,
    calcCytPosGates = calcCytPosGates
  )
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
          gateTblCtrl = NULL,
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
