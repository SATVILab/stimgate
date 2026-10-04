#' @title Identify cytokine-positive cells through automated gating
#'
#' @description
#' Identify cells responding to stimulation by comparing cytokine expression in
#' stimulated samples with unstimulated controls from the same donor/batch.
#' Saves gates and background-subtracted statistics to a project directory.
#'
#' @param pathProject character. Path to project directory where all results will be saved.
#'   This directory will contain subdirectories for each marker with gate tables,
#'   statistics, and plots. The directory will be created if it doesn't exist.
#' @param .data GatingSet, flowSet, cytoset, flowFrame, cytoframe, character,
#'   list or data.frame Cytometry samples, FCS paths/directory, numeric matrices
#'   or data frames per sample, or a long data frame with a `sample` column.
#' @param batchList list. List where each element contains integer indices or character names of samples belonging to the same batch/donor. The first index per element is the unstimulated control sample, e.g. if `batchList = list(c(3, 1, 2), c(6, 4, 5))`, then indices 3 and 6 correspond to the unstimulated samples for batches 1 and 2, respectively. If `batchList` is named, e.g. `list(pid1 = c(3, 1, 2), pid2 = c(6, 4, 5))`, then these names will be used for batch identification. Each element needs the unstimulated sample and at least one stimulated sample. An unstimulated sample may be shared across batches (first in each), but a stimulated sample may belong to only one batch.
#' @param chnl character vector. Channel names to gate on. Specify either
#'   `chnl` or `marker`. Default is NULL.
#' @param marker character vector. Alternative way to specify markers to gate on.
#'   When provided, this is used instead of chnl to determine which markers to analyze.
#'   Default is NULL.
#' @param popGate character vector. Population(s) within which to perform gating.
#'   Default is "root" to gate on all cells. Can specify other populations like
#'   "CD3+" or "CD4+" if these gates already exist in the GatingSet.
#' @param biasUns numeric. Bias adjustment for unstimulated samples to account for
#'   background cytokine production. When NULL (default), 1/4 of `bwFallback` is used
#'   (scaled by `biasUnsFactor`). Positive values shift the unstimulated distribution higher,
#'   making gates more conservative. Default is NULL.
#' @param bw numeric. Specify the bandwith for density estimation. When NULL (default), bandwidth is estimated automatically. A bandwidth may also be set per marker through `markerControl`. Default is `NULL`.
#' @param control stimControl Tuning settings from [stimControl()].
#'   Most users do not need to change these. Default: `stimControl()`.
#' @param markerControl list or NULL. Named per-marker overrides, keyed by
#'   marker label or channel name, for example
#'   `list(IL2 = list(bw = 0.12, biasUns = 0))`. Settings from [stimControl()]
#'   (except the global-only `locEnforceShapeThreshold` and `calcCytPosGates`),
#'   plus `bw`, `biasUns` and `popGate`, can be overridden. Default: NULL.
#' @return character. Returns the path to the project directory where all results
#'   have been saved. The directory structure created includes:
#'   \itemize{
#'     \item \code{pathProject/[markerName]/}: Directory for each marker containing:
#'     \item \code{gateTblInit.rds}: Initial gate table with preliminary gates
#'     \item \code{gateTbl.rds}: Final refined gate table
#'     \item \code{stats/}: Directory containing statistics files
#'     \item \code{plots/}: Directory containing visualization plots (if generated)
#'   }
#'
#' @param parallel logical If TRUE, gate channels in parallel during the initial gating stage using the active future::plan(). Default: FALSE.
#' @details
#' Initial local-FDR gates compare each stimulated sample with its batch's
#' unstimulated control. Thresholds may be shared across similar distributions,
#' then refined using cells positive for another cytokine. Results and statistics
#' are saved to `pathProject` for [getStimGates()], [getStimStats()], and [plotStim()].
#' Use [stimControl()] for tuning and `markerControl` for per-marker overrides.
#'
#' Inputs are normalised to a GatingSet; only GatingSet inputs can contain
#' populations other than "root", including in `markerControl`. StimGate applies
#' no transformation (such as arcsinh or logicle): provide FCS/matrix data on
#' the scale you want gated, as with a GatingSet. FCS files are read using
#' `flowWorkspace::load_cytoset_from_fcs()` with its default reader behaviour.
#' A directory is searched non-recursively for case-insensitive `.fcs` filenames,
#' sorted with `sort()`; a file vector preserves its supplied order. FCS sample
#' names are basenames. List names are sample names, or default to sample1,
#' sample2, etc. Channel columns must be numeric with matching names; their
#' order is aligned to the first sample and marker labels equal channel names.
#' Long data frames use observed factor-level order or first appearance of
#' `sample`. The first sample in each batch is always the unstimulated control.
#'
#' To gate channels in parallel, set `parallel = TRUE` and select a future
#' plan, for example `future::plan(future::multisession, workers = 4)`.
#' The default `parallel = FALSE` runs sequentially regardless of the active
#' plan. Only the initial per-channel gating stage is parallel; subsequent
#' cytokine-positive gating and statistics remain sequential. Workers read
#' expression data from the project disk cache rather than a GatingSet, so
#' the project directory must be accessible to all workers.
#' With `parallel = TRUE`, RNG-dependent subsampling uses parallel-safe
#' L'Ecuyer streams (`future.seed = TRUE`). Results are reproducible for a
#' given `set.seed()` and independent of the chosen non-sequential plan,
#' but may differ slightly from a sequential run.
#'
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- file.path(tempdir(), "demonstration")
#'
#' # Run gating
#' gateStim(
#'   .data = gs,
#'   pathProject = pathProject,
#'   popGate = "root",
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker
#' )
#'
#' # Customise tuning and override the bandwidth for the first marker
#' gateStim(
#'   pathProject = file.path(tempdir(), "custom-gating"),
#'   .data = gs,
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker,
#'   bw = 0.1,
#'   control = stimControl(bwAdj = 1.5, clusterGates = FALSE),
#'   markerControl = stats::setNames(
#'     list(list(bw = 0.12, biasUns = 0)), exampleData$marker[1]
#'   )
#' )
#'
#' # Use in-memory matrices, with channel names also serving as marker labels
#' matrices <- lapply(seq_along(gs), function(i) {
#'   flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]], y = "root"))
#' })
#' gateStim(
#'   pathProject = file.path(tempdir(), "matrix-gating"), .data = matrices,
#'   batchList = exampleData$batchList, chnl = exampleData$chnl,
#'   control = stimControl(calcCytPosGates = FALSE)
#' )
#'
#' # Create plots
#' if (requireNamespace("hexbin", quietly = TRUE)) {
#'   plots <- plotStim(
#'     ind = exampleData$batchList[[1]], # indices in `gs` to plot
#'     .data = gs, # GatingSet
#'     pathProject = pathProject,
#'     marker = exampleData$marker,
#'     grid = TRUE
#'   )
#' }
#' @export
gateStim <- function(
  pathProject,
  .data,
  batchList,
  marker = NULL,
  chnl = NULL,
  popGate = "root",
  biasUns = NULL,
  bw = NULL,
  control = stimControl(),
  markerControl = NULL,
  parallel = FALSE
) {
  isGatingSet <- inherits(.data, "GatingSet")
  .data <- .asStimGatingSet(.data, popGate)
  if (Sys.getenv("STIMGATE_DEBUG") == "") {
    Sys.setenv("STIMGATE_DEBUG" = "FALSE")
    on.exit(Sys.unsetenv("STIMGATE_DEBUG"), add = TRUE)
  }

  mustDebug <- tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")

  if (
    mustDebug &&
      is.character(pathProject) &&
      length(pathProject) == 1L &&
      !is.na(pathProject) &&
      nzchar(pathProject)
  ) {
    .debugInit(pathProject)
    pathDebug <- file.path(pathProject, "debug", "debug.txt")
    message(paste0("Saving debug output to ", pathDebug))
    .profileInit(pathProject, reset = TRUE)
  }

  runSuccess <- FALSE
  on.exit(
    {
      if (.profileEnabled() && isTRUE(.profileState$initialized)) {
        tryCatch(
          {
            .profileFinishRun(
              pathProject = pathProject,
              status = if (runSuccess) "completed" else "failed"
            )
          },
          error = function(e) {
            .profileStateReset()
            invisible(NULL)
          }
        )
      }
      .debugStateReset()
    },
    add = TRUE
  )

  # Verify global function inputs before reading anything from `control`, so an
  # invalid control object is reported as such.
  .verifyGateInputs(
    pathProject = pathProject,
    .data = .data,
    batchList = batchList,
    popGate = popGate,
    chnl = chnl,
    marker = marker,
    biasUns = biasUns,
    bw = bw,
    control = control,
    markerControl = markerControl
  )

  if (!isGatingSet) {
    .checkStimInputPop(unlist(lapply(markerControl, function(x) x$popGate)))
  }
  batchList <- .resolveBatchList(batchList, .data)

  calcCytPosGates <- control$calcCytPosGates

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE")
  }
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Install the 'future.apply' package to use `parallel = TRUE`.")
  }

  .verifyBatchList(batchList)
  if (is.null(names(batchList))) {
    batchList <- batchList |>
      stats::setNames(paste0("batch", seq_along(batchList)))
  }

  # get unspecified levels in marker elements
  .saveMetaData(.data, batchList, pathProject)
  chnl <- .extractChnl(chnl, marker, pathProject)

  chnlSettingsRaw <- .resolveMarkerControl(
    markerControl = markerControl,
    chnl = chnl,
    chnlLab = stimgateMetaReadChnlLab(pathProject)
  )
  chnlSettingsCache <- lapply(chnl, function(chnlCurr) {
    popCurr <- chnlSettingsRaw[[chnlCurr]]$popGate
    list(
      chnlCut = chnlCurr,
      popGate = if (!is.null(popCurr)) popCurr else popGate
    )
  })
  runPops <- unique(
    vapply(chnlSettingsCache, function(x) x$popGate, character(1))
  )
  .gateInvalidateRunPopulations(pathProject = pathProject, pops = runPops)
  .gateCacheChnl(
    .data = .data,
    indBatchList = batchList,
    chnlSettings = chnlSettingsCache,
    pathProject = pathProject
  )

  chnlSettings <- .completeChnlSettings(
    chnl = chnl,
    markerControl = markerControl,
    control = control,
    biasUns = biasUns,
    bw = bw,
    .data = .data,
    popGate = popGate,
    indBatchList = batchList,
    pathProject = pathProject
  )

  # inital gates
  .gateInit(
    chnlSettings = chnlSettings,
    .data = .data,
    indBatchList = batchList,
    pathProject = pathProject,
    parallel = parallel
  )

  # cytokine-positive gates
  gateTbl <- .gateCytPos(
    chnlSettings = chnlSettings,
    indBatchList = batchList,
    .data = .data,
    calcCytPos = calcCytPosGates,
    stage = "cytPos",
    pathProject = pathProject
  )

  message("getting cyt combn frequencies")

  .gateStats(
    .data = .data,
    gateTbl = gateTbl,
    calcCytPosGates = calcCytPosGates,
    chnlSettings = chnlSettings,
    indBatchList = batchList,
    pathProject = pathProject
  )

  runSuccess <- TRUE
  pathProject
}

#' @keywords internal
.gateInit <- function(
  chnlSettings,
  .data,
  indBatchList,
  pathProject,
  parallel = FALSE
) {
  message("getting base gates")

  .gateMapChnl(chnlSettings, .data, indBatchList, pathProject, parallel)
  invisible(chnlSettings)
}

#' @keywords internal
.gateStats <- function(
  .data,
  gateTbl = NULL,
  calcCytPosGates,
  chnlSettings,
  indBatchList,
  pathProject
) {
  force(.data)
  .getStats(
    gateTbl = gateTbl,
    filterOtherCytPos = FALSE,
    combn = TRUE,
    gateTypeCytPosFilter = if (calcCytPosGates) "cyt" else "base",
    gateTypeCytPosCalc = if (calcCytPosGates) "cyt" else "base",
    save = TRUE,
    popGate = chnlSettings[[1]]$popGate,
    chnl = purrr::map_chr(chnlSettings, function(x) x$chnlCut),
    indBatchList = indBatchList,
    .data = .data,
    saveGateTbl = TRUE,
    pathProject = pathProject
  )
}

#' @keywords internal
.gateInvalidateRunPopulations <- function(pathProject, pops) {
  for (pop in pops) {
    pathExPop <- dirname(.getExChnlPathDir("dummy", pop, pathProject))
    if (dir.exists(pathExPop)) {
      unlink(pathExPop, recursive = TRUE)
    }
    pathGatePop <- dirname(dirname(dirname(
      .gatesGetPathAll(pathProject, pop, "dummy", FALSE)
    )))
    if (dir.exists(pathGatePop)) {
      unlink(pathGatePop, recursive = TRUE)
    }
  }
  invisible(NULL)
}
