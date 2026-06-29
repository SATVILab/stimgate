.simBandwidthBsFreq <- function(
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
  probExact = FALSE,
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
  covEvMax = 2,
  tolClust = 1e-7
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
      probExact = probExact,
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
      tempdir(),
      "stimgate",
      "sim-bw",
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
      bw = bw,
      tolClust = tolClust
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
      }
    ) |>
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

    pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
    chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
    chnl <- chnlVec[[1]]
    ind <- 2

    propBsTblEstSmooth <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- seq(2, nSample * nCondition, by = 2) |>
        as.character()
      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", ind)
        probSmooth <- file.path(pathInd, "dataMod.rds") |>
          readRDS()
        if (!inherits(probSmooth, "data.frame")) {
          # it's zero if it doesn't inherit anything
          return(tibble::tibble(
            sample = as.character(as.numeric(ind) / nCondition),
            ind = ind,
            chnl = chnl,
            ncellCondition = nCellStim,
            ncellRespSmooth = 0,
            ncellRespPred = 0,
            propRespSmooth = 0,
            propRespPred = 0
          ))
        }
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
        values_to = "propRespEst"
      ) |>
      dplyr::mutate(
        threshold = NA_real_,
        thresholdOrigin = NA_character_,
        locGenerated = NA,
        locGeneratedDirect = NA,
        locSource = NA_character_,
        locReason = NA_character_,
        detailLevel = NA_character_
      )

    propBsTblDetailed <- try(
      getStimGatesDetailed(pathProject),
      silent = TRUE
    )
    propBsTblDetailed <- if (
      inherits(propBsTblDetailed, "try-error") ||
        !is.data.frame(propBsTblDetailed) ||
        nrow(propBsTblDetailed) == 0L
    ) {
      tibble::tibble()
    } else {
      propBsTblDetailed |>
        dplyr::filter(
          .data$detailLevel %in% c("condition", "sample", "cluster_final")
        ) |>
        dplyr::filter(is.finite(.data$propBs)) |>
        dplyr::mutate(
          sample = as.character(as.numeric(.data$ind) / nCondition),
          method = paste0("loc_", .data$detailLevel),
          propRespEst = .data$propBs
        ) |>
        dplyr::select(
          sample,
          chnl,
          method,
          propRespEst,
          threshold,
          thresholdOrigin,
          locGenerated,
          locGeneratedDirect,
          locSource,
          locReason,
          detailLevel
        )
    }

    propBsTblEst <- dplyr::bind_rows(
      propBsTblEstSmooth,
      propBsTblDetailed
    )

    comparisonTbl <- propBsTblTruth |>
      dplyr::left_join(
        propBsTblEst,
        by = c("sample", "chnl")
      )
    comparisonTbl |>
      dplyr::mutate(
        iter = iterNum
      ) |>
      dplyr::select(iter, chnl, dplyr::everything())
  })
}

.simBandwidthEstBw <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  probExact = FALSE,
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
  covEvMax = 2,
  tolClust = 1e-7
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
      probExact = probExact,
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
      tempdir(),
      "stimgate",
      "sim-bw",
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
      bwMtd = bwMtd,
      tolClust = tolClust
    ))

    stopifnot(file.exists(file.path(pathProject, "gateStats.rds")))

    pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
    chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
    chnl <- chnlVec[[1]]
    ind <- 2

    bwTbl <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- seq(2, nSample * nCondition, by = 2) |>
        as.character()
      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", ind)
        bwEst <- file.path(pathInd, "bwCpUnsLoc.rds") |>
          readRDS()
        if (!is.numeric(bwEst)) {
          # it's zero if it doesn't inherit anything
          return(tibble::tibble(
            sample = as.character(as.numeric(ind) / nCondition),
            ind = ind,
            chnl = chnl,
            bw = bwEst
          ))
        }

        tibble::tibble(
          sample = as.character(as.numeric(ind) / nCondition),
          ind = ind,
          chnl = chnl,
          bw = bwEst
        )
      })
    }) |>
      dplyr::mutate(
        iter = iterNum
      ) |>
      dplyr::select(iter, chnl, dplyr::everything())
  })
}

.simPlotGate <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  probExact = FALSE,
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
  covEvMax = 2,
  tolClust = 1e-7,
  markerToPlot = "MarkerF1" # Specify univariate marker here
) {
  allPlots <- lapply(seq_len(nIter), function(iterNum) {
    # 1. Setup Simulation Parameters
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

    # 2. Simulate Cytometry Experiment
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
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

    # 3. Setup Project Directory for gating
    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-plot",
      paste0("iter-", iterNum, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    )

    on.exit(
      {
        if (dir.exists(pathProject)) {
          unlink(pathProject, recursive = TRUE)
        }
      },
      add = TRUE
    )

    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE)

    batchList <- lapply(seq(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    # 4. Execute Gating
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      bw = bw,
      biasUns = biasUns,
      bwMtd = bwMtd,
      bwMin = bwMin,
      bwMax = bwMax,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      bwCluster = bwCluster,
      tolClust = tolClust
    ))

    # 5. Extract plots for each sample
    samplePlots <- lapply(seq_along(batchList), function(sampleIdx) {
      # batchList[[sampleIdx]] grabs both unstim and stim conditions
      indCurr <- batchList[[sampleIdx]]

      # Assuming condition 1 is unstim and condition 2 is stim based on setup
      conditionLabels <- c("Unstimulated", "Stimulated")
      names(conditionLabels) <- indCurr

      # Generate the overlaid plot
      pList <- plotStim(
        ind = indCurr,
        .data = gs,
        pathProject = pathProject,
        marker = markerToPlot,
        indLab = conditionLabels,
        grid = FALSE, # Return raw list of ggplots instead of rendering a cowplot grid
        showGate = TRUE, # Overlay gates
        excMin = TRUE
      )

      # plotStim returns a list (length 1 if univariate). Extract it.
      if (!is.null(pList) && length(pList) > 0) {
        p <- pList[[1]]
        # Update title to track iteration and sample
        p <- p +
          ggplot2::ggtitle(
            sprintf(
              "Iter: %d | Sample: %d | %s",
              iterNum,
              sampleIdx,
              markerToPlot
            )
          )
        return(p)
      } else {
        return(NULL)
      }
    })

    return(samplePlots)
  })

  # Flatten the nested lists so you get a single continuous list of ggplot objects
  do.call(c, allPlots)
}
