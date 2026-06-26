.simBandwidth <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  samplePerturbationSd,
  conditionPerturbationSd,
  clusterPerturbationSd,
  backgroundRelativeToResponse,
  ncellUnsRelativeToStim,
  covEvMin = 1,
  covEvMax = 2
) {

  purrr::map_df(seq_len(nIter), function(iterNum) {
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simMiscGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet") # Fixed case
    gs <- flowWorkspace::GatingSet(fs)
    labelsList <- outListExperiment[["labelsList"]] # Fixed case

    pathProject <- file.path(
      tempdir(), "stimgate", "sim-bandwidth",
      paste0("iter-", iterNum, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    )
    on.exit({
      if (dir.exists(pathProject)) {
        unlink(pathProject, recursive = TRUE)
      }
    })
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE)
    
    # FIXED: Corrected parenthesis and variable name 'nSample'
    batchList <- lapply(seq(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    # gate
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      bw = bw
    ))

    stopifnot(file.exists(file.path(pathProject, "gateStats.rds")))

    propBsTblTruth <- purrr::map_df(
      seq_len(nSample),
      function(sampleCurr) {
        labelsVec <- (sampleCurr - 1) * nCondition + seq(1, nCondition)
        purrr::map_df(labelsVec, function(i) {
          labelVec <- labelsList[[i]]
          stimCondition <- ifelse(i %% nCondition == 1, "unstim", "stim")
          f1p <- sum(grepl("^gp$", labelVec)) / length(labelVec)
          tibble::tibble(
            ind = sampleCurr |> as.character(),
            stim = stimCondition,
            chnl = "F1",
            F1 = f1p
          )
        })
      }) |>
        tidyr::pivot_wider(
          names_from = stim,
          values_from = F1
        ) |>
        dplyr::mutate(
          propRespTruth = stim - unstim
        ) |>
        dplyr::rename(
          sample = ind
        ) |>
        dplyr::select(
          all_of(c("chnl", "sample", "propRespTruth"))
        )
      
      # FIXED: Swapped out snake_case to match pathProject and camelCase consistency
      pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
      chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
      chnl <- chnlVec[[1]]
      ind <- 2

    probBsTblEst <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- seq(2, nSample * nCondition, by = 2) |>
        as.character()
      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", ind)
        probSmooth <- file.path(pathInd, "dataMod.rds") |>
          readRDS()
        ncellRespSmooth <- if (nrow(probSmooth) > 0L) {
          sum(probSmooth$probSmooth)
        } else {
          0
        }
        ncellRespPred <- if (nrow(probSmooth) > 0L) {
          sum(probSmooth$pred)
        } else {
          0
        }
        tibble::tibble(
          sample = as.character(as.numeric(ind) / nCondition),
          ind = ind,
          chnl = chnl,
          ncellCondition = nCellStim,
          ncellRespSmooth = ncellRespSmooth,
          ncellRespPred = ncellRespPred,
          propRespSmooth = ncellRespSmooth / ncellCondition,
          propRespPred = ncellRespPred / ncellCondition
        )
      })
    }) |>
      dplyr::select(
        sample,
        chnl,
        propRespSmooth,
        propRespPred
      ) |>
      tidyr::pivot_longer(
        cols = c(propRespSmooth, propRespPred),
        names_to = "method",
        values_to = "propResp"
      )
    
    comparisonTbl <- probBsTblTruth |>
      dplyr::left_join(
        probBsTblEst,
        by = c("sample", "chnl")
      )
    comparisonTbl |>
      dplyr::mutate(
        iter = iterNum
      ) |>
      dplyr::select(iter, chnl, dplyr::everything())
  })
}
