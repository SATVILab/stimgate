#' @title Gate cells responding to stimulation
#' @description Compare cytokine expression in stimulated samples with an
#'   unstimulated control from the same donor or batch. Save gates and
#'   background-subtracted statistics to a project directory.
#' @param pathProject character Directory for results; created if needed.
#' @param .data GatingSet, flowSet, cytoset, flowFrame, cytoframe, character,
#'   list or data.frame Cytometry samples: a GatingSet or other flowCore/
#'   flowWorkspace object, FCS file paths or one FCS directory, a list of
#'   numeric matrices or data frames (cells by channels), or one data frame
#'   with a `sample` column. See Details.
#' @param batchList list Samples grouped by donor or batch, as indices or
#'   sample names, with the unstimulated control first in each vector, e.g.
#'   `list(donor1 = c(3, 1, 2))`. List names identify batches.
#' @param marker character vector or NULL Marker labels to gate; supply
#'   either `marker` or `chnl`. Default: NULL.
#' @param chnl character vector or NULL Channel names to gate. Default: NULL.
#' @param popGate character Population(s) already present in `.data`; only
#'   GatingSets have populations other than "root". Default: "root" (all cells).
#' @param biasUns numeric or NULL Upward shift of unstimulated expression.
#'   NULL uses one quarter of `bwFallback`, scaled by `biasUnsFactor`.
#'   Positive shifts make gating more conservative. Default: NULL.
#' @param bw numeric or NULL Fixed density bandwidth; NULL estimates it
#'   automatically. Per-marker values go in `markerControl`. Default: NULL.
#' @param control stimControl Tuning settings from [stimControl()].
#'   Default: stimControl().
#' @param markerControl list or NULL Overrides keyed by marker label or channel,
#'   e.g. `list(IL2 = list(bw = 0.12, biasUns = 0))`. Accepts [stimControl()]
#'   settings except `locEnforceShapeThreshold` and `calcCytPosGates`, plus
#'   `bw`, `biasUns` and `popGate`. Default: NULL.
#' @param parallel logical Use the active [future::plan()] for initial channel
#'   gating. Default: FALSE.
#' @details
#' Thresholds can be shared across similar distributions, then refined using
#' cells positive for another cytokine. Read results with [getStimGates()],
#' [getStimStats()] and [getStimExpr()]; inspect them with [plotStim()].
#'
#' **Input data.** Non-GatingSet inputs are converted to a GatingSet with only
#' the "root" population. StimGate does not transform data (e.g. arcsinh), so
#' supply values on the scale to gate. Sample order, which `batchList` indices
#' refer to, is:
#' * FCS directory: `.fcs` files (any case, not recursive) in sorted order;
#'   a vector of paths keeps its order. Sample names are file basenames.
#' * List of matrices/data frames: list order; names are sample names
#'   (default `sample1`, `sample2`, ...). Columns must have the same unique
#'   names in every sample and serve as both channels and markers.
#' * Data frame with `sample`: factor-level order, otherwise first appearance.
#'
#' Pass the same data, in the same order, to [plotStim()], [getStimExpr()]
#' and [writeStimFCS()].
#'
#' For parallel gating, set `parallel = TRUE` and choose a future plan, e.g.
#' `future::plan(future::multisession, workers = 4)`. All workers must be able
#' to access `pathProject`. Later stages run sequentially. Set a seed for
#' reproducible parallel subsampling; results may differ from sequential runs.
#' @return A character string: `pathProject`. Gates are saved under `gates/`,
#'   expression under `sampleData/`, settings under `metaData/`, and statistics
#'   in `gateStats.rds` and `gateStats.csv`.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' getStimGates(pathProject)
#'
#' # Disable gate sharing and fix the first marker's bandwidth
#' gateStim(
#'   tempfile("custom_gating_"), gs, exampleData$batchList,
#'   marker = exampleData$marker, control = stimControl(clusterGates = FALSE),
#'   markerControl = stats::setNames(
#'     list(list(bw = 0.12, biasUns = 0)), exampleData$marker[1]
#'   )
#' )
#'
#' # Gate in-memory matrices; column names act as channels and markers
#' matrices <- lapply(seq_along(gs), function(i) {
#'   flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
#' })
#' gateStim(
#'   tempfile("matrix_gating_"), matrices, exampleData$batchList,
#'   chnl = exampleData$chnl, control = stimControl(calcCytPosGates = FALSE)
#' )
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
  sampleNames <- flowWorkspace::sampleNames(.data)
  batchList <- lapply(batchList, function(batch) {
    if (!is.character(batch)) {
      return(batch)
    }
    indices <- match(batch, sampleNames)
    if (anyNA(indices)) {
      stop(
        "Unknown sample name(s) in `batchList`: ",
        paste(batch[is.na(indices)], collapse = ", ")
      )
    }
    indices
  })

  calcCytPosGates <- control$calcCytPosGates

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE")
  }
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Install the 'future.apply' package to use `parallel = TRUE`.")
  }

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
