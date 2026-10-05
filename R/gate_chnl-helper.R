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
    )
  })
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
  gateTbl <- gateTbl |>
    dplyr::filter(gateUse == "gate") |> # nolint
    dplyr::select(-gateUse) # nolint

  # =========================
  # Cluster-based gating
  # =========================

  if (isTRUE(chnlSettings$clusterGates)) {
    # share one expression lookup across gate names; with a single gate name
    # .getCpCluster() builds it itself
    exLookup <- if (length(unique(gateTbl$gateName)) > 1L) {
      .getCpClusterLocExLookup(
        .data = .data,
        indBatchList = indBatchList,
        chnlSettings = chnlSettings,
        filterOtherCytPos = FALSE,
        calcCytPosGates = calcCytPosGates,
        gateTbl = gateTbl,
        pathProject = pathProject
      )
    }
    gateTblCluster <- purrr::map_df(
      unique(gateTbl$gateName),
      function(gn) {
        gateTblCluster <- .getCpCluster(
          # nolint
          .data = .data,
          gateTbl = gateTbl |>
            dplyr::filter(gateName == gn), # nolint
          chnlSettings = chnlSettings,
          filterOtherCytPos = FALSE,
          stage = stage,
          pathProject = pathProject,
          calcCytPosGates = calcCytPosGates,
          indBatchList = indBatchList,
          exLookup = exLookup
        )

        gateTblCluster |>
          dplyr::select(
            ind, gate = cpJoinTgOrig,
            dplyr::any_of(c(
              "locGenerated", "locGeneratedDirect", "locSource", "locReason"
            ))
          ) |> # nolint
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
          dplyr::relocate(gateName, gateType, gateCombn, batch) |> # nolint
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
