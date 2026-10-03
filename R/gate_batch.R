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

  .debug("chnlSettings$gateTbl is NULL") # nolint
  .debug(
    "gating ",
    paste0(indBatch, collapse = "-") # nolint
  )

  # create bare list
  gateList <- .getCpUnsLoc(
    exList = exList,
    .data = .data,
    chnlSettings = chnlSettings,
    stage = stage,
    pathProject = pathProject
  )

  .gateBatchTbl(gateList, attr(exList[[1]], "batch")) # nolint
}
