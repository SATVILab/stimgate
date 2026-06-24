#' Get example GatingSet
#'
#' Create and save a complete synthetic GatingSet for testing and examples
#' using the simCytExperiment functions, avoiding external data dependencies.
#'
#' @param scenario Character specifying the simulation scenario ("default", "easy", "poor_separation", "cyt_pos").
#' @param dir_cache Directory to save the GatingSet. If NULL, uses a temporary directory.
#' @param clear Logical indicating whether to force cache clearing.
#' @param n_ind Integer. Number of biological samples to simulate.
#' @return A list containing the path to the saved GatingSet, batch_list, chnl, and marker names.
#' @export
get_example_data <- function(
  scenario = "default",
  dir_cache = NULL,
  clear = FALSE,
  n_ind = 8
) {
  # Set seed for reproducibility
  set.seed(123)

  if (is.null(dir_cache)) {
    dir_cache <- testthat::test_path("cache", "test_data", scenario)
  }

  if (dir.exists(dir_cache)) {
    if (clear) {
      unlink(dir_cache, recursive = TRUE)
    } else {
      path_gs <- file.path(dir_cache, "gs")
      path_chnl <- file.path(dir_cache, "chnl.rds")
      path_marker <- file.path(dir_cache, "marker.rds")
      path_batch <- file.path(dir_cache, "batch_list.rds")
      
      if (
        dir.exists(path_gs) &&
        file.exists(path_chnl) &&
        file.exists(path_marker) &&
        file.exists(path_batch)
      ) {
        return(list(
          path_gs = path_gs,
          chnl = readRDS(path_chnl),
          marker = readRDS(path_marker),
          batch_list = readRDS(path_batch)
        ))
      } else {
        message("Cache incomplete, regenerating synthetic test data...")
        unlink(dir_cache, recursive = TRUE)
      }
    }
  }
  dir.create(dir_cache, recursive = TRUE)

  # 1. Define Simulation Matrix Parameters 
  nMarker <- 2L
  nCondition <- 2L
  nCluster <- 4L
  nCellByCondition <- 1000L # Can be adjusted for speed vs depth
  transformationFunc <- function(x) x
  clusterLabelVec <- c("neg_neg", "neg_pos", "pos_neg", "pos_pos")
  
  # Base probability of a cell falling into each of the 4 clusters in an Unstimulated state
  probVecUns <- c(0.90, 0.04, 0.04, 0.02)
  # Shift in probabilities under Stimulated conditions
  probResponseVec <- list(c(-0.10, 0.05, 0.04, 0.01))

  # Adjust scenario means and perturbations
  if (scenario %in% c("default", "easy")) {
    meanExprMat <- matrix(c(0,0, 0,2, 2,0, 2,2), nrow = nCluster, byrow = TRUE)
    clusterPerturbationSd <- 0.1
  } else if (scenario == "poor_separation" || scenario == "cyt_pos") {
    meanExprMat <- matrix(c(0,0, 0,0.5, 0.5,0, 0.5,0.5), nrow = nCluster, byrow = TRUE)
    clusterPerturbationSd <- 0.2
  } else {
    stop("Unknown scenario: ", scenario)
  }

  # 2. Run simCytExperiment 
  sim_res <- simCytExperiment(
    nSample = as.integer(n_ind),
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
    clusterPerturbationSd = clusterPerturbationSd
  )

  # 3. Format generated FlowFrames to mimic historical test channels
  chnl_vec <- c("BC1(La139)Dd", "BC2(Pr141)Dd")
  marker_vec <- paste0("MarkerF", seq_len(nMarker))
  
  fs_list <- lapply(sim_res$flowFrameList, function(ff) {
    # Update matrix col names
    flowCore::colnames(ff) <- chnl_vec
    # Update Parameter Annotations
    p <- flowCore::parameters(ff)
    pData_new <- flowCore::pData(p)
    pData_new$name <- chnl_vec
    pData_new$desc <- marker_vec
    flowCore::pData(p) <- pData_new
    ff
  })

  # 4. Construct flowSet and GatingSet
  fs <- flowCore::flowSet(fs_list)
  gs <- flowWorkspace::GatingSet(fs)
  
  path_gs <- file.path(dir_cache, "gs")
  flowWorkspace::save_gs(gs, path = path_gs)

  # 5. Build matching batch_list mapping
  # simCytExperiment returns sequence as: [Unstim_1, Stim_1, Unstim_2, Stim_2...]
  # The historical test structure requires the LAST element in a batch vector to be the Unstim index.
  batch_list <- lapply(seq_len(n_ind), function(i) {
    idx_unstim <- (i - 1) * nCondition + 1
    idx_stim <- (i - 1) * nCondition + 2
    c(idx_stim, idx_unstim)
  })

  # Save cache states
  saveRDS(chnl_vec, file = file.path(dir_cache, "chnl.rds"))
  saveRDS(marker_vec, file = file.path(dir_cache, "marker.rds"))
  saveRDS(batch_list, file = file.path(dir_cache, "batch_list.rds"))

  list(
    path_gs = path_gs,
    batch_list = batch_list,
    chnl = chnl_vec,
    marker = marker_vec
  )
}

#' @rdname get_example_data
#' @export
get_test_data <- get_example_data
