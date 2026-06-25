#' Get example GatingSet
#'
#' Create and save a complete synthetic GatingSet for testing and examples
#' using the simCytExperiment functions, avoiding external data dependencies.
#'
#' @param scenario Character specifying the simulation scenario ("default", "easy", "poorSeparation", "cytPos").
#' @param dirCache Directory to save the GatingSet. If NULL, uses a temporary directory.
#' @param clear Logical indicating whether to force cache clearing.
#' @param nInd Integer. Number of biological samples to simulate.
#' @return A list containing the path to the saved GatingSet, batchList, chnl, and marker names.
#' @export
getExampleData <- function(
  scenario = "default",
  dirCache = NULL,
  clear = FALSE,
  nCell = 1e4,
  nInd = 8
) {
  # Set seed for reproducibility
  set.seed(123)

  if (is.null(dirCache)) {
    dirCache <- tempfile(pattern = "dir_")
    dir.create(dirCache)
  }

  if (dir.exists(dirCache)) {
    if (clear) {
      unlink(dirCache, recursive = TRUE)
    } else {
      pathGs <- file.path(dirCache, "gs")
      pathChnl <- file.path(dirCache, "chnl.rds")
      pathMarker <- file.path(dirCache, "marker.rds")
      pathBatch <- file.path(dirCache, "batchList.rds")

      if (
        dir.exists(pathGs) &&
          file.exists(pathChnl) &&
          file.exists(pathMarker) &&
          file.exists(pathBatch)
      ) {
        return(list(
          pathGs = pathGs,
          chnl = readRDS(pathChnl),
          marker = readRDS(pathMarker),
          batchList = readRDS(pathBatch)
        ))
      } else {
        message("Cache incomplete, regenerating synthetic test data...")
        unlink(dirCache, recursive = TRUE)
      }
    }
  }
  dir.create(dirCache, recursive = TRUE)

  # 1. Define Simulation Matrix Parameters
  nMarker <- 2L
  nCondition <- 2L
  nCluster <- 4L
  nCellByCondition <- nCell # Can be adjusted for speed vs depth
  transformationFunc <- function(x) x
  clusterLabelVec <- c("negNeg", "negPos", "posNeg", "posPos")

  # Base probability of a cell falling into each of the 4 clusters in an Unstimulated state
  probVecUns <- c(0.90, 0.04, 0.04, 0.02)
  # Shift in probabilities under Stimulated conditions
  probResponseVec <- list(c(-0.25, 0.2, 0.04, 0.01))

  # Adjust scenario means and perturbations
  if (scenario %in% c("default", "easy")) {
    meanExprMat <- matrix(
      c(0, 0, 0, 8, 8, 0, 8, 8),
      nrow = nCluster,
      byrow = TRUE
    )
  } else if (scenario == "poorSeparation" || scenario == "cytPos") {
    meanExprMat <- matrix(
      c(0, 0, 0, 4, 4, 0, 4, 4),
      nrow = nCluster,
      byrow = TRUE
    )
  } else {
    stop("Unknown scenario: ", scenario)
  }

  # 2. Run simCytExperiment
  simRes <- simCytExperiment(
    nSample = as.integer(nInd),
    nMarker = nMarker,
    nCondition = nCondition,
    nCluster = nCluster,
    nCellByCondition = nCellByCondition,
    transformationFunc = transformationFunc,
    mixtureType = "gaussianOnly",
    meanExprMat = meanExprMat,
    clusterLabelVec = clusterLabelVec,
    probVecUns = probVecUns,
    probResponseVecByStimCondition = probResponseVec,
    clusterPerturbationSd = 0
  )

  saveRDS(simRes[["labelsList"]], file = file.path(dirCache, "labelsList.rds"))
  saveRDS(
    simRes[["flowFrameList"]],
    file = file.path(dirCache, "flowFrameList.rds")
  )

  # 3. Format generated FlowFrames to mimic historical test channels
  chnlVec <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
  markerVec <- paste0("MarkerF", seq_len(nMarker))

  fsList <- lapply(simRes$flowFrameList, function(ff) {
    # Update matrix col names
    flowCore::colnames(ff) <- chnlVec
    # Update Parameter Annotations
    p <- flowCore::parameters(ff)
    pDataNew <- flowCore::pData(p)
    pDataNew$name <- chnlVec
    pDataNew$desc <- markerVec
    flowCore::pData(p) <- pDataNew
    ff
  })

  # 4. Construct flowSet and GatingSet
  fs <- flowCore::flowSet(fsList)
  gs <- flowWorkspace::GatingSet(fs)

  pathGs <- file.path(dirCache, "gs")
  flowWorkspace::save_gs(gs, path = pathGs)

  # 5. Build matching batchList mapping
  # simCytExperiment returns sequence as: [Unstim_1, Stim_1, Unstim_2, Stim_2...]
  # The historical test structure requires the LAST element in a batch vector to be the Unstim index.
  batchList <- lapply(seq_len(nInd), function(i) {
    idxUnstim <- (i - 1) * nCondition + 1
    idxStim <- (i - 1) * nCondition + 2
    c(idxStim, idxUnstim)
  })

  # Save cache states
  saveRDS(chnlVec, file = file.path(dirCache, "chnl.rds"))
  saveRDS(markerVec, file = file.path(dirCache, "marker.rds"))
  saveRDS(batchList, file = file.path(dirCache, "batchList.rds"))

  list(
    pathGs = pathGs,
    batchList = batchList,
    chnl = chnlVec,
    marker = markerVec,
    pathCache = dirCache
  )
}

#' @rdname getExampleData
#' @export
getTestData <- getExampleData
