# Initial channel gating uses only expression vectors cached on disk. Keep the
# worker entry point in the namespace so its closure cannot capture a GatingSet.

#' @keywords internal
.gateMapChnl <- function(
  chnlSettings,
  .data,
  indBatchList,
  pathProject,
  parallel = FALSE
) {
  if (parallel) {
    .gateCacheChnl(.data, indBatchList, chnlSettings, pathProject)
    pathProject <- normalizePath(pathProject, winslash = "/", mustWork = TRUE)
  }
  if (parallel && length(chnlSettings) > 1L) {
    future.apply::future_lapply(
      chnlSettings,
      .gateInitChnlWorker,
      indBatchList = indBatchList,
      pathProject = pathProject,
      env = Sys.getenv(
        c("STIMGATE_DEBUG", "STIMGATE_INTERMEDIATE"),
        unset = NA
      ),
      future.seed = TRUE
    )
  } else {
    lapply(
      chnlSettings,
      .gateInitChnl,
      .data = .data,
      indBatchList = indBatchList,
      pathProject = pathProject
    )
  }
}

#' @keywords internal
.gateCacheChnl <- function(.data, indBatchList, chnlSettings, pathProject) {
  pops <- unique(vapply(chnlSettings, function(x) x$popGate, character(1)))
  for (ind in unique(unlist(indBatchList, use.names = FALSE))) {
    for (pop in pops) {
      chnl <- unique(vapply(
        Filter(function(x) x$popGate == pop, chnlSettings),
        function(x) x$chnlCut,
        character(1)
      ))
      missing <- chnl[!vapply(chnl, function(x) {
        file.exists(.getExChnlPath(x, ind, pop, pathProject))
      }, logical(1))]
      if (length(missing) == 0L) {
        next
      }
      if (is.null(.data)) {
        stop(
          "Incomplete expression cache for sample ", ind,
          ", population ", pop, "."
        )
      }
      # Read each sample/population once, then save every missing channel.
      fr <- flowWorkspace::gh_pop_get_data(.data[[ind]], y = pop)
      ex <- flowCore::exprs(fr)
      for (chnlCurr in missing) {
        pathChnl <- .getExChnlPath(chnlCurr, ind, pop, pathProject)
        dir.create(dirname(pathChnl), recursive = TRUE, showWarnings = FALSE)
        saveRDS(ex[, chnlCurr], pathChnl)
      }
    }
  }
  invisible(NULL)
}

#' @keywords internal
.gateInitChnl <- function(chnlSettings, .data, indBatchList, pathProject) {
  message(paste0("chnl: ", chnlSettings$chnlCut))
  gateObj <- .gateChnl(
    .data = .data,
    indBatchList = indBatchList,
    chnlSettings = chnlSettings,
    pathProject = pathProject,
    stage = "init",
    calcCytPosGates = FALSE
  )
  pathSave <- .gatesGetPathAll(
    pathProject = pathProject,
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    init = TRUE
  )
  dir.create(dirname(pathSave), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gateObj$gateTbl, pathSave)
}

#' @keywords internal
.gateInitChnlWorker <- function(chnlSettings, indBatchList, pathProject, env) {
  oldEnv <- Sys.getenv(names(env), unset = NA)
  oldDebug <- as.list(.debugState)
  oldProfile <- as.list(.profileState)
  on.exit({
    Sys.unsetenv(names(oldEnv)[is.na(oldEnv)])
    do.call(Sys.setenv, as.list(oldEnv[!is.na(oldEnv)]))
    list2env(oldDebug, envir = .debugState)
    list2env(oldProfile, envir = .profileState)
  }, add = TRUE)
  Sys.unsetenv(names(env)[is.na(env)])
  do.call(Sys.setenv, as.list(env[!is.na(env)]))
  if (.profileEnabled()) {
    # Attach without resetting the parent's debug/profile directories.
    .profileAttach(pathProject)
    .profileState$context <- .profileContextDefault()
    .debugState$file <- file.path(pathProject, "debug", "debug.txt")
    .debugState$initialized <- TRUE
  }
  .gateInitChnl(chnlSettings, NULL, indBatchList, pathProject)
}
