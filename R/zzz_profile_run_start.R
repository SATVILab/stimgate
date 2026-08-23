# Start debug profiling at the first project write ---------------------------
#
# gateStim() calls .saveMetaData() once after input verification and before
# channel setup or gating. Wrapping that call gives each debug run a fresh
# profile directory without changing the public gateStim() signature.

.profileOriginalSaveMetaData <- .saveMetaData

.saveMetaData <- function(.data, batchList, pathProject) {
  if (.profileEnabled()) {
    .profileInit(pathProject, reset = TRUE)
  }
  .profileOriginalSaveMetaData(
    .data = .data,
    batchList = batchList,
    pathProject = pathProject
  )
}
