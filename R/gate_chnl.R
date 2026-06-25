#' @keywords internal
.gateChnl <- function(
  .data,
  indBatchList,
  chnlSettings,
  gateTbl = NULL,
  tolGateSingle = NULL,
  calcCytPosGates,
  pathProject,
  stage
) {
  # message progress
  .debug("popGate: ", popGate) # nolint

  # Parameters list
  # ----------------

  # named
  chnlLabVec <- .getLabs(
    # nolint
    .data = .data[[indBatchList[[1]]]],
    chnlCut = chnlSettings$chnlCut
  )

  # delete locb gates
  .gateChnlDeleteOldGates()

  # Initial gates
  # ----------------
  gateTbl <- .gateChnlPreAdjGatesGate(
    # nolint
    indBatchList = indBatchList,
    .data = .data,
    chnlSettings = chnlSettings,
    stage = stage,
    pathProject = pathProject
  )

  gateTbl <- gateTbl |>
    dplyr::filter(
      !as.character(ind) %in%
        vapply(
          indBatchList,
          function(x) as.character(x[1]),
          character(1)
        )
    ) # nolint

  # Get adjusted and/or clustered gates
  # ----------------

  message("getting clustered and/or controlled gates")

  # =================================================
  # Cluster-based or tail-gate controlled gating
  # =================================================

  # =============================
  # For all cells
  # =============================

  .gateChnlGetAdjGates(
    # nolint
    gateTbl = gateTbl,
    gateTblParams = chnlSettings$gateTbl,
    chnlSettings = chnlSettings,
    .data = .data,
    pathProject = pathProject,
    indBatchList = indBatchList,
    stage = stage,
    calcCytPosGates = calcCytPosGates
  )
}
