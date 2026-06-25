# Get cutpoints for a single batch
#' @keywords internal
.gateBatch <- function(
  .data,
  indBatch,
  chnlSettings,
  batch,
  stage,
  pathProject
) {
  # get list of dataframes
  exList <- .getExList(
    # nolint
    .data = .data,
    indBatch = indBatch,
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    batch = batch,
    pathProject = pathProject
  )

  if (is.null(chnlSettings$gateTbl)) {
    .gateBatchAll(
      indBatch = indBatch,
      batch = batch,
      exList = exList,
      .data = .data,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject
    )
  } else {
    .gateBatchSingle(
      indBatch = indBatch,
      batch = batch,
      exList = exList,
      .data = .data,
      chnlSettings = chnlSettings,
      calcCytPosGates = calcCytPosGates,
      stage = stage,
      pathProject = pathProject
    )
  }
}
