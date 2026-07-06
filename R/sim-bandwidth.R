.simBwMtdColVec <- c(
  "hpi0" = "#005a32",
  "hpi1" = "#238b45",
  "hpi2" = "#74c476",
  "hpi3" = "#c7e9c0",
  "nrd0" = "#fe9929",
  "sj" = "#3d8bb3ff",
  "hpi0Norm" = "#54278f",
  "hpi1Norm" = "#756bb1",
  "hpi2Norm" = "#9e9ac8",
  "hpi3Norm" = "#cbc9e2",
  "nrd0Norm" = "#d95f0e",
  "sjNorm" = "#2b8cbe"
)

.simBandwidthReadRdsOrNull <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    readRDS(path),
    error = function(e) NULL
  )
}

.simBandwidthAddMissingColumns <- function(.data, cols) {
  for (nm in names(cols)) {
    if (!nm %in% names(.data)) {
      .data[[nm]] <- cols[[nm]]
    }
  }
  .data
}

.simBandwidthSampleFromInd <- function(ind, nCondition) {
  ind_num <- suppressWarnings(as.numeric(ind))
  as.character(((ind_num - 1) %/% nCondition) + 1)
}

.simBandwidthLocDetailGatePoint <- function(
  detailLevel,
  locGenerated,
  locGeneratedDirect,
  locSource
) {
  dplyr::case_when(
    detailLevel %in%
      "condition" &
      locGenerated %in% TRUE &
      locGeneratedDirect %in% TRUE ~ "condition_direct_local_fdr",
    detailLevel %in%
      "condition" &
      !(locGenerated %in% TRUE) ~ "condition_fallback_high_value",
    detailLevel %in%
      "sample" &
      locGenerated %in% TRUE &
      locSource %in% "combined" ~ "sample_combined_from_other_stim_conditions",
    detailLevel %in%
      "sample" &
      locGenerated %in% TRUE &
      locSource %in% "prejoin" ~ "sample_prejoin_from_joined_stim_conditions",
    detailLevel %in%
      "sample" &
      locGenerated %in% TRUE ~ "sample_final_local_fdr",
    detailLevel %in%
      "sample" &
      !(locGenerated %in% TRUE) ~ "sample_fallback_high_value",
    TRUE ~ NA_character_
  )
}

.simBandwidthReadLocDetails <- function(
  pathProject,
  nSample,
  nCondition
) {
  pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
  if (!dir.exists(pathDirIntInit)) {
    return(tibble::tibble())
  }

  chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
  if (length(chnlVec) == 0L) {
    return(tibble::tibble())
  }

  detailTbl <- purrr::map_df(chnlVec, function(chnl) {
    stageChnl <- file.path(pathDirIntInit, chnl)

    purrr::map_df(seq_len(nSample), function(sampleCurr) {
      indBatch <- seq(
        (sampleCurr - 1L) * nCondition + 1L,
        sampleCurr * nCondition
      )
      indStim <- indBatch[-1L]
      indCombined <- paste0(indStim, collapse = "_")

      conditionTbl <- purrr::map_df(indStim, function(ind) {
        pathInd <- file.path(stageChnl, "ind", as.character(ind))
        out <- .simBandwidthReadRdsOrNull(
          file.path(pathInd, "locDetailCondition.rds")
        )
        if (!is.data.frame(out) || nrow(out) == 0L) {
          return(tibble::tibble())
        }
        out
      })

      sampleTbl <- .simBandwidthReadRdsOrNull(
        file.path(stageChnl, "ind", indCombined, "locDetailSample.rds")
      )
      if (!is.data.frame(sampleTbl) || nrow(sampleTbl) == 0L) {
        sampleTbl <- tibble::tibble()
      }

      dplyr::bind_rows(conditionTbl, sampleTbl)
    })
  })

  if (!is.data.frame(detailTbl) || nrow(detailTbl) == 0L) {
    return(tibble::tibble())
  }

  detailTbl <- .simBandwidthAddMissingColumns(
    detailTbl,
    list(
      detailLevel = NA_character_,
      stage = NA_character_,
      chnl = NA_character_,
      ind = NA_character_,
      threshold = NA_real_,
      thresholdOrigin = NA_character_,
      locGenerated = NA,
      locGeneratedDirect = NA,
      locSource = NA_character_,
      locReason = NA_character_,
      bias = NA_real_,
      propBsEst = NA_real_,
      propBsDiff = NA_real_,
      nCellStim = NA_integer_,
      nCellUns = NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propBs = NA_real_
    )
  )

  conditionThresholdTbl <- detailTbl |>
    dplyr::filter(.data$detailLevel %in% "condition") |>
    dplyr::transmute(
      chnl = .data$chnl,
      ind = as.character(.data$ind),
      thresholdCondition = suppressWarnings(as.numeric(.data$threshold)),
      thresholdOriginCondition = .data$thresholdOrigin,
      locGeneratedCondition = .data$locGenerated %in% TRUE,
      locSourceCondition = .data$locSource,
      locReasonCondition = .data$locReason
    )

  detailTbl |>
    dplyr::mutate(
      ind = as.character(.data$ind),
      sample = .simBandwidthSampleFromInd(.data$ind, nCondition),
      method = paste0("loc_", .data$detailLevel),
      propRespEst = .data$propBs,
      nCellStim = suppressWarnings(as.numeric(.data$nCellStim)),
      nCellUns = suppressWarnings(as.numeric(.data$nCellUns)),
      propStim = suppressWarnings(as.numeric(.data$propStim)),
      propUns = suppressWarnings(as.numeric(.data$propUns)),
      nPosStim = dplyr::if_else(
        is.finite(.data$propStim) & is.finite(.data$nCellStim),
        as.integer(round(.data$propStim * .data$nCellStim)),
        NA_integer_
      ),
      nPosUns = dplyr::if_else(
        is.finite(.data$propUns) & is.finite(.data$nCellUns),
        as.integer(round(.data$propUns * .data$nCellUns)),
        NA_integer_
      ),
      gateReturnPoint = .simBandwidthLocDetailGatePoint(
        detailLevel = .data$detailLevel,
        locGenerated = .data$locGenerated,
        locGeneratedDirect = .data$locGeneratedDirect,
        locSource = .data$locSource
      )
    ) |>
    dplyr::left_join(
      conditionThresholdTbl,
      by = c("chnl", "ind")
    ) |>
    dplyr::mutate(
      thresholdBeforeSampleCombining = dplyr::if_else(
        .data$detailLevel %in% "sample",
        .data$thresholdCondition,
        NA_real_
      ),
      thresholdChangedBySampleCombining = dplyr::case_when(
        .data$detailLevel %in%
          "sample" &
          is.finite(.data$threshold) &
          is.finite(.data$thresholdCondition) ~
          abs(.data$threshold - .data$thresholdCondition) >
            sqrt(.Machine$double.eps),
        .data$detailLevel %in% "sample" ~ NA,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::select(
      sample,
      ind,
      chnl,
      method,
      propRespEst,
      detailLevel,
      stage,
      threshold,
      thresholdOrigin,
      gateReturnPoint,
      locGenerated,
      locGeneratedDirect,
      locSource,
      locReason,
      thresholdBeforeSampleCombining,
      thresholdChangedBySampleCombining,
      thresholdCondition,
      thresholdOriginCondition,
      locGeneratedCondition,
      locSourceCondition,
      locReasonCondition,
      nCellStim,
      nCellUns,
      nPosStim,
      nPosUns,
      propStim,
      propUns,
      propBs,
      propBsEst,
      propBsDiff,
      bias
    )
}

.simBandwidthBsFreq <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwFallback = "auto",
  bwMin = "auto",
  bwMax = "auto",
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 1e2,
  bwNcellMax = 1e5,
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
  tolClust = NULL,
  minCell = 1e2,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locLeftLowRel = 0.25,
  locLeftLowAbs = 0.15,
  locLeftCellFrac = 0.5,
  locLeftLengthFrac = 0.5,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  locTolRefPeak = "highest",
  gateCombn = "min",
  calcCytPosGates = FALSE,
  calcSinglePosGates = FALSE
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

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet")
    gs <- flowWorkspace::GatingSet(fs)
    labelsList <- outListExperiment[["labelsList"]]

    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-bw",
      paste0(
        "iter-",
        iterNum,
        "-",
        Sys.getpid(),
        "-",
        format(Sys.time(), "%Y%m%d%H%M%S"),
        "-",
        sample.int(1e7, 1)
      )
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
    dir.create(pathProject, recursive = TRUE, showWarnings = FALSE)

    batchList <- lapply(seq_len(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      calcCytPosGates = calcCytPosGates,
      calcSinglePosGates = calcSinglePosGates,
      biasUns = biasUns,
      bw = bw,
      bwFallback = bwFallback,
      bwMin = bwMin,
      bwMax = bwMax,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      bwCluster = bwCluster,
      minCell = minCell,
      maxPosProbX = maxPosProbX,
      gateQuant = gateQuant,
      tolClust = tolClust,
      locProbCol = locProbCol,
      locMinPeakProb = locMinPeakProb,
      locDipAlpha = locDipAlpha,
      locAntimodeHeightFrac = locAntimodeHeightFrac,
      locAntimodeLowRel = locAntimodeLowRel,
      locAntimodeLowAbs = locAntimodeLowAbs,
      locFlatDerivFrac = locFlatDerivFrac,
      locFlatHardDerivFrac = locFlatHardDerivFrac,
      locLeftLowRel = locLeftLowRel,
      locLeftLowAbs = locLeftLowAbs,
      locLeftCellFrac = locLeftCellFrac,
      locLeftLengthFrac = locLeftLengthFrac,
      locMarginalPurityRel = locMarginalPurityRel,
      locMarginalCellBinRatio = locMarginalCellBinRatio,
      locMarginalRefQuantile = locMarginalRefQuantile,
      locTolRefPeak = locTolRefPeak,
      gateCombn = gateCombn
    ))

    stopifnot(file.exists(file.path(pathProject, "gateStats.rds")))

    propBsTblTruth <- purrr::map_df(
      seq_len(nSample),
      function(sampleCurr) {
        indUns <- (sampleCurr - 1L) * nCondition + 1L
        indStim <- seq.int(indUns + 1L, sampleCurr * nCondition)
        labelVecUns <- labelsList[[indUns]]
        propUnsTruth <- sum(grepl("^gp$", labelVecUns)) / length(labelVecUns)

        purrr::map_df(indStim, function(ind) {
          labelVecStim <- labelsList[[ind]]
          propStimTruth <- sum(grepl("^gp$", labelVecStim)) /
            length(labelVecStim)
          tibble::tibble(
            sample = as.character(sampleCurr),
            ind = as.character(ind),
            chnl = "F1",
            propStimTruth = propStimTruth,
            propUnsTruth = propUnsTruth,
            propRespTruth = propStimTruth - propUnsTruth
          )
        })
      }
    )

    pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
    chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)

    propBsTblEstSmooth <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- unlist(lapply(seq_len(nSample), function(sampleCurr) {
        indUns <- (sampleCurr - 1L) * nCondition + 1L
        seq.int(indUns + 1L, sampleCurr * nCondition)
      }))

      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", as.character(ind))
        probSmooth <- .simBandwidthReadRdsOrNull(
          file.path(pathInd, "dataMod.rds")
        )
        if (!inherits(probSmooth, "data.frame")) {
          return(tibble::tibble(
            sample = .simBandwidthSampleFromInd(ind, nCondition),
            ind = as.character(ind),
            chnl = chnl,
            method = c("propRespSmooth", "propRespPred"),
            propRespEst = c(0, 0),
            nCellStim = nCellStim,
            nCellUns = nCellUns,
            nPosStim = NA_integer_,
            nPosUns = NA_integer_,
            threshold = NA_real_,
            thresholdOrigin = NA_character_,
            gateReturnPoint = NA_character_,
            locGenerated = NA,
            locGeneratedDirect = NA,
            locSource = NA_character_,
            locReason = NA_character_,
            detailLevel = NA_character_,
            stage = NA_character_,
            thresholdBeforeSampleCombining = NA_real_,
            thresholdChangedBySampleCombining = NA
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
          sample = .simBandwidthSampleFromInd(ind, nCondition),
          ind = as.character(ind),
          chnl = chnl,
          method = c("propRespSmooth", "propRespPred"),
          propRespEst = c(
            ncellRespSmooth / nCellStim,
            ncellRespPred / nCellStim
          ),
          nCellStim = nCellStim,
          nCellUns = nCellUns,
          nPosStim = NA_integer_,
          nPosUns = NA_integer_,
          threshold = NA_real_,
          thresholdOrigin = NA_character_,
          gateReturnPoint = NA_character_,
          locGenerated = NA,
          locGeneratedDirect = NA,
          locSource = NA_character_,
          locReason = NA_character_,
          detailLevel = NA_character_,
          stage = NA_character_,
          thresholdBeforeSampleCombining = NA_real_,
          thresholdChangedBySampleCombining = NA
        )
      })
    })

    propBsTblDetailed <- .simBandwidthReadLocDetails(
      pathProject = pathProject,
      nSample = nSample,
      nCondition = nCondition
    ) |>
      dplyr::filter(.data$detailLevel %in% c("condition", "sample"))

    propBsTblEst <- dplyr::bind_rows(
      propBsTblEstSmooth,
      propBsTblDetailed
    )

    comparisonTbl <- propBsTblTruth |>
      dplyr::left_join(
        propBsTblEst,
        by = c("sample", "ind", "chnl")
      )

    comparisonTbl |>
      dplyr::mutate(
        iter = iterNum,
        nCellStimSim = nCellStim,
        nCellUnsSim = nCellUns,
        biasUns = biasUns,
        bw = if (is.null(bw)) NA_real_ else bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = if (is.null(bwCluster)) NA_real_ else bwCluster,
        samplePerturbationSd = samplePerturbationSd,
        conditionPerturbationSd = conditionPerturbationSd,
        clusterPerturbationSd = clusterPerturbationSd,
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim
      ) |>
      dplyr::select(iter, chnl, sample, ind, dplyr::everything())
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

#' Estimate bandwidths directly from simulated data, without running gateStim
#'
#' This mirrors the bandwidth-estimation part of cp_uns_loc:
#'   1. simulate unstim/stim data
#'   2. optionally exclude the minimum values, as gateStim does by default
#'   3. optionally cap values at the same max-density x used by cp_uns_loc
#'   4. estimate bw separately for stim and unstim
#'   5. return min(bw_stim, bw_uns)
#'
#' @keywords internal
.simBandwidthEstBwDirect <- function(
  nSample = 10L,
  nMarker = 1L,
  nCondition = 2L,
  nCluster = 2L,
  nIter = 10L,
  biasUns = 0.05,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = 1e-10,
  bwMax = 1e10,
  bwFallback = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL, # retained only for signature compatibility
  tolClust = NULL, # retained only for signature compatibility
  probExact = TRUE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  backgroundRelativeToResponse = 0.2,
  ncellUnsRelativeToStim = 1,
  covEvMin = 2,
  covEvMax = 2,
  excMin = TRUE,
  capStimRange = TRUE,
  summarise = TRUE
) {
  if (!identical(as.integer(nMarker), 1L)) {
    stop("This helper currently expects nMarker = 1.")
  }
  if (!identical(as.integer(nCondition), 2L)) {
    stop("This helper currently expects nCondition = 2.")
  }
  if (!identical(as.integer(nCluster), 2L)) {
    stop("This helper currently expects nCluster = 2.")
  }

  nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
  nCellByCondition <- c(nCellUns, nCellStim)

  transformationFunc <- .simBandwidthGetTrans(transformation)

  meanExprMat <- matrix(
    c(0, meanPos),
    byrow = TRUE,
    ncol = 1
  )

  clusterLabelVec <- c("gn", "gp")

  probResponseUns <- probResponse * backgroundRelativeToResponse
  probVecUns <- c(1 - probResponseUns, probResponseUns)

  probResponseVecByStimCondition <- list(
    c(-probResponse, probResponse)
  )

  raw_tbl <- purrr::map_dfr(seq_len(nIter), function(iterNum) {
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
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    flowFrameList <- outListExperiment[["flowFrameList"]]

    purrr::map_dfr(seq_len(nSample), function(sampleCurr) {
      indUns <- (sampleCurr - 1L) * nCondition + 1L
      indStim <- indUns + 1L

      x_uns <- as.numeric(flowCore::exprs(flowFrameList[[indUns]])[, "F1"])
      x_stim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, "F1"])

      # cp_uns_loc applies bias to the unstim expression before density work.
      x_uns <- x_uns + (biasUns %||% 0)

      if (excMin) {
        x_uns <- .simBandwidthExcMin(x_uns)
        x_stim <- .simBandwidthExcMin(x_stim)
      }

      x_cap <- .simBandwidthCapForCpUnsLoc(
        x_stim = x_stim,
        x_uns = x_uns,
        capStimRange = capStimRange
      )

      bw_stim <- .simBandwidthBwOne(
        x = x_cap$x_stim,
        bwMtd = bwMtd,
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwFallback = bwFallback
      )

      bw_stim <- .simBandwidthRemoveFallbackBw(
        bw = bw_stim,
        bwFallback = bwFallback
      )

      bw_uns <- .simBandwidthBwOne(
        x = x_cap$x_uns,
        bwMtd = bwMtd,
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwFallback = bwFallback
      )

      bw_uns <- .simBandwidthRemoveFallbackBw(
        bw = bw_uns,
        bwFallback = bwFallback
      )

      bw_final <- if (!is.null(bw)) {
        bw
      } else {
        bw_vec <- c(bw_stim, bw_uns)
        bw_vec <- bw_vec[!is.na(bw_vec)]
        if (length(bw_vec) == 0) {
          NA_real_
        } else {
          min(bw_vec)
        }
      }

      tibble::tibble(
        transformation = transformation,
        prob_response = probResponse,
        n_cell = nCellStim,
        mean_pos = meanPos,
        bw_mtd = bwMtd,
        iter = iterNum,
        sample = as.character(sampleCurr),
        ind = as.character(indStim),
        chnl = "F1",
        n_cell_uns = nCellUns,
        n_cell_stim = nCellStim,
        n_uns_bw = length(x_cap$x_uns),
        n_stim_bw = length(x_cap$x_stim),
        max_dens_x = x_cap$max_dens_x,
        bw_uns = bw_uns,
        bw_stim = bw_stim,
        bw = bw_final,
        bw_source = dplyr::case_when(
          !is.null(bw) ~ "fixed",
          is.finite(bw_stim) &
            (!is.finite(bw_uns) || bw_stim <= bw_uns) ~ "stim",
          is.finite(bw_uns) ~ "unstim",
          TRUE ~ NA_character_
        )
      )
    })
  })

  if (!summarise) {
    return(raw_tbl)
  }

  .simBandwidthSummariseBw(raw_tbl)
}

#' Estimate bandwidths directly from simulated data, without running gateStim
#'
#' This mirrors the bandwidth-estimation part of cp_uns_loc:
#'   1. simulate unstim/stim data
#'   2. optionally exclude the minimum values, as gateStim does by default
#'   3. optionally cap values at the same max-density x used by cp_uns_loc
#'   4. estimate bw separately for stim and unstim
#'   5. return min(bw_stim, bw_uns)
#'
#' @keywords internal
.simBandwidthEstBwDirectAdaptive <- function(
  nSample = 10L,
  nMarker = 1L,
  nCondition = 2L,
  nCluster = 2L,
  nIter = 10L,
  biasUns = 0.05,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = 1e-10,
  bwMax = 1e10,
  bwFallback = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  tolClust = NULL,
  probExact = TRUE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  backgroundRelativeToResponse = 0.2,
  ncellUnsRelativeToStim = 1,
  covEvMin = 2,
  covEvMax = 2,
  excMin = TRUE,
  capStimRange = TRUE,
  summarise = FALSE
) {
  if (!identical(as.integer(nMarker), 1L)) {
    stop("This helper currently expects nMarker = 1.")
  }
  if (!identical(as.integer(nCondition), 2L)) {
    stop("This helper currently expects nCondition = 2.")
  }
  if (!identical(as.integer(nCluster), 2L)) {
    stop("This helper currently expects nCluster = 2.")
  }

  nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
  nCellByCondition <- c(nCellUns, nCellStim)

  transformationFunc <- .simBandwidthGetTrans(transformation)

  meanExprMat <- matrix(
    c(0, meanPos),
    byrow = TRUE,
    ncol = 1
  )

  clusterLabelVec <- c("gn", "gp")

  probResponseUns <- probResponse * backgroundRelativeToResponse
  probVecUns <- c(1 - probResponseUns, probResponseUns)

  probResponseVecByStimCondition <- list(
    c(-probResponse, probResponse)
  )

  raw_tbl <- purrr::map_dfr(seq_len(nIter), function(iterNum) {
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
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    flowFrameList <- outListExperiment[["flowFrameList"]]

    purrr::map_dfr(seq_len(nSample), function(sampleCurr) {
      indUns <- (sampleCurr - 1L) * nCondition + 1L
      indStim <- indUns + 1L

      x_uns <- as.numeric(flowCore::exprs(flowFrameList[[indUns]])[, "F1"])
      x_stim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, "F1"])

      # cp_uns_loc applies bias to the unstim expression before density work.
      x_uns <- x_uns + (biasUns %||% 0)

      if (excMin) {
        x_uns <- .simBandwidthExcMin(x_uns)
        x_stim <- .simBandwidthExcMin(x_stim)
      }

      x_cap <- .simBandwidthCapForCpUnsLoc(
        x_stim = x_stim,
        x_uns = x_uns,
        capStimRange = capStimRange
      )

      exTblStimThreshold <- tibble::tibble(
        F1 = x_cap$x_stim
      )
      attr(exTblStimThreshold, "chnlCut") <- "F1"
      exTblUnsThreshold <- tibble::tibble(
        F1 = x_cap$x_uns
      )
      attr(exTblUnsThreshold, "chnlCut") <- "F1"

      bwObj <- .getCpUnsLocGetDensRawDensitiesAdaptive(
        exTblStimThreshold = exTblStimThreshold,
        exTblUnsThreshold = exTblUnsThreshold,
        chnlSettings = list(
          bwMtd = bwMtd,
          bwMin = bwMin,
          bwMax = bwMax,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax,
          bwFallback = bwFallback
        )
      )
      bwStimCore <- tryCatch(bwObj$bw$stim$bwCore, error = function(e) NA_real_)
      bwStimExtra <- tryCatch(bwObj$bw$stim$bwExtra, error = function(e) {
        NA_real_
      })
      thresholdStim <- tryCatch(
        bwObj$bw$stim$coreObj$thresholdX,
        error = function(e) NA_real_
      )
      bwUnsCore <- tryCatch(bwObj$bw$uns$bwCore, error = function(e) NA_real_)
      bwUnsExtra <- tryCatch(bwObj$bw$uns$bwExtra, error = function(e) NA_real_)
      thresholdUns <- tryCatch(
        bwObj$bw$uns$coreObj$thresholdX,
        error = function(e) NA_real_
      )
      sharedGrid <- tryCatch(bwObj$bw$sharedGrid, error = function(e) NA_real_)
      densUnsWeight <- tryCatch(bwObj$bw$densUnsWeight, error = function(e) {
        NA_real_
      })
      densStimWeight <- tryCatch(bwObj$bw$densStimWeight, error = function(e) {
        NA_real_
      })
      tibble::tibble(
        transformation = transformation,
        prob_response = probResponse,
        n_cell = nCellStim,
        mean_pos = meanPos,
        bw_mtd = bwMtd,
        iter = iterNum,
        sample = as.character(sampleCurr),
        ind = as.character(indStim),
        chnl = "F1",
        n_cell_uns = nCellUns,
        n_cell_stim = nCellStim,
        n_uns_bw_core = length(x_cap$x_uns),
        n_stim_bw_core = length(x_cap$x_stim),
        max_dens_x = x_cap$max_dens_x,
        bw_uns_core = bwUnsCore,
        bw_stim_core = bwStimCore,
        bw_uns_extra = bwUnsExtra,
        bw_stim_extra = bwStimExtra,
        threshold_uns = thresholdUns,
        threshold_stim = thresholdStim,
        shared_grid = list(sharedGrid),
        dens_uns_weight = list(densUnsWeight),
        dens_stim_weight = list(densStimWeight)
      )
    })
  })

  if (!summarise) {
    return(raw_tbl)
  }

  .simBandwidthSummariseBw(raw_tbl)
}


#' Run the direct bandwidth estimator over a Simulation-Bandwidth sim_grid
#'
#' @keywords internal
.simBandwidthEstBwDirectGrid <- function(
  sim_grid,
  nSample = 10L,
  nIter = 10L,
  biasUns = 0.05,
  bwMin = 1e-10,
  bwMax = 1e10,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  probExact = TRUE,
  backgroundRelativeToResponse = 0.2,
  ncellUnsRelativeToStim = 1,
  covEvMin = 2,
  covEvMax = 2,
  excMin = TRUE,
  capStimRange = TRUE,
  summarise = TRUE
) {
  raw_tbl <- purrr::map_dfr(seq_len(nrow(sim_grid)), function(i) {
    row <- sim_grid[i, , drop = FALSE]

    out <- tryCatch(
      .simBandwidthEstBwDirect(
        nSample = nSample,
        nIter = nIter,
        biasUns = biasUns,
        bwMtd = row$bw_mtd[[1]],
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        probExact = probExact,
        nCellStim = row$n_cell[[1]],
        probResponse = row$prob_response[[1]],
        meanPos = row$mean_pos[[1]],
        transformation = row$transformation[[1]],
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim,
        covEvMin = covEvMin,
        covEvMax = covEvMax,
        excMin = excMin,
        capStimRange = capStimRange,
        summarise = FALSE
      ),
      error = function(e) {
        tibble::tibble(
          iter = NA_integer_,
          sample = NA_character_,
          ind = NA_character_,
          chnl = "F1",
          n_cell_uns = NA_real_,
          n_cell_stim = row$n_cell[[1]],
          n_uns_bw = NA_integer_,
          n_stim_bw = NA_integer_,
          max_dens_x = NA_real_,
          bw_uns = NA_real_,
          bw_stim = NA_real_,
          bw = NA_real_,
          bw_source = NA_character_,
          error = e$message
        )
      }
    )

    # Avoid duplicating scenario columns if the direct helper already added them.
    out <- out |>
      dplyr::select(
        -dplyr::any_of(c(
          "transformation",
          "prob_response",
          "n_cell",
          "mean_pos",
          "bw_mtd"
        ))
      )

    dplyr::bind_cols(
      row[rep(1L, nrow(out)), , drop = FALSE],
      out
    )
  })

  if (!summarise) {
    return(raw_tbl)
  }

  .simBandwidthSummariseBw(raw_tbl)
}


#' Summarise raw bandwidth estimates by simulation scenario
#'
#' @keywords internal
.simBandwidthSummariseBw <- function(.data) {
  group_vars <- intersect(
    c(
      "transformation",
      "prob_response",
      "n_cell",
      "mean_pos",
      "bw_mtd",
      "chnl"
    ),
    names(.data)
  )

  .data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_est = sum(is.finite(.data$bw)),
      n_error = sum(!is.na(.data$error %||% NA_character_)),
      bw_mean = mean(.data$bw, na.rm = TRUE),
      bw_median = stats::median(.data$bw, na.rm = TRUE),
      bw_q05 = stats::quantile(.data$bw, 0.05, na.rm = TRUE),
      bw_q25 = stats::quantile(.data$bw, 0.25, na.rm = TRUE),
      bw_q75 = stats::quantile(.data$bw, 0.75, na.rm = TRUE),
      bw_q95 = stats::quantile(.data$bw, 0.95, na.rm = TRUE),
      bw_min = min(.data$bw, na.rm = TRUE),
      bw_max = max(.data$bw, na.rm = TRUE),
      bw_uns_median = stats::median(.data$bw_uns, na.rm = TRUE),
      bw_stim_median = stats::median(.data$bw_stim, na.rm = TRUE),
      n_source_uns = sum(.data$bw_source == "unstim", na.rm = TRUE),
      n_source_stim = sum(.data$bw_source == "stim", na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$n_cell,
      dplyr::desc(.data$prob_response),
      .data$transformation,
      .data$mean_pos,
      .data$bw_mtd
    )
}


#' Bandwidth estimate for one vector
#'
#' This is the vector-only equivalent of
#' .getCpUnsLocGetDensRawDensitiesBwInit(). It intentionally routes through
#' .bwCalcOne() when available so that the direct bandwidth simulations use the
#' same ordinary and *Norm bandwidth methods as the gating code. In particular,
#' hpi0Norm, hpi1Norm, hpi2Norm, hpi3Norm, sjNorm and nrd0Norm are handled by
#' the shared normalised-bandwidth helper rather than by this wrapper.
#'
#' @keywords internal
.simBandwidthBwOne <- function(
  x,
  bwMtd,
  bwMin,
  bwMax,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwFallback
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(.simBandwidthBwFallbackOrNa(bwFallback))
  }

  bw_calc <- tryCatch(
    {
      if (exists(".bwCalcOne", mode = "function")) {
        .bwCalcOne(
          x = x,
          bwMtd = bwMtd,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax
        )
      } else {
        .simBandwidthBwOneBaseLegacy(
          x = x,
          bwMtd = bwMtd,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax
        )
      }
    },
    error = function(e) NA_real_
  )

  bw_calc <- suppressWarnings(as.numeric(bw_calc)[1])

  if (!is.finite(bw_calc) || bw_calc <= 0) {
    return(.simBandwidthBwFallbackOrNa(bwFallback))
  }

  if (.simBandwidthIsFiniteScalar(bwMin)) {
    bw_calc <- max(as.numeric(bwMin)[1], bw_calc)
  }
  if (.simBandwidthIsFiniteScalar(bwMax)) {
    bw_calc <- min(as.numeric(bwMax)[1], bw_calc)
  }

  bw_calc
}

#' @keywords internal
.simBandwidthRemoveFallbackBw <- function(
  bw,
  bwFallback
) {
  bw <- suppressWarnings(as.numeric(bw)[1])

  if (!is.finite(bw)) {
    return(NA_real_)
  }

  if (!.simBandwidthIsFiniteScalar(bwFallback)) {
    return(bw)
  }

  if (isTRUE(all.equal(bw, as.numeric(bwFallback)[1], tolerance = 0))) {
    return(NA_real_)
  }

  bw
}

#' @keywords internal
.simBandwidthIsFiniteScalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

#' @keywords internal
.simBandwidthBwFallbackOrNa <- function(bwFallback) {
  if (is.null(bwFallback)) {
    return(NA_real_)
  }

  bwFallback <- suppressWarnings(as.numeric(bwFallback)[1])

  if (!is.finite(bwFallback) || bwFallback <= 0) {
    return(NA_real_)
  }

  bwFallback
}

#' @keywords internal
.simBandwidthBwOneBaseLegacy <- function(
  x,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(NA_real_)
  }

  if (.simBandwidthIsFiniteScalar(bwNcellMin) && length(x) < bwNcellMin) {
    sdX <- .simBandwidthRobustSd(x)
    x <- sample(x, replace = TRUE, size = bwNcellMin) +
      stats::rnorm(bwNcellMin, mean = 0, sd = sdX / 10)
  }

  if (.simBandwidthIsFiniteScalar(bwNcellMax) && length(x) > bwNcellMax) {
    x <- sample(x, size = bwNcellMax, replace = FALSE)
  }

  bwMtd <- as.character(bwMtd)[1]
  if (grepl("Norm$", bwMtd)) {
    return(NA_real_)
  }

  bwMtdBase <- bwMtd

  bw_calc <- switch(
    bwMtdBase,
    "nrd0" = try(stats::bw.nrd0(x), silent = TRUE),
    "sj" = try(stats::bw.SJ(x), silent = TRUE),
    {
      derivOrder <- suppressWarnings(as.numeric(gsub("^hpi", "", bwMtdBase)))

      if (!is.finite(derivOrder)) {
        return(NA_real_)
      }

      try(
        suppressWarnings(
          ks::hpi(x = x, deriv.order = derivOrder)
        ),
        silent = TRUE
      )
    }
  )

  if (inherits(bw_calc, "try-error") || !is.finite(bw_calc) || bw_calc <= 0) {
    return(NA_real_)
  }

  as.numeric(bw_calc)[1] * bwAdj
}

#' @keywords internal
.simBandwidthRobustSd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    return(.Machine$double.eps)
  }

  iqrX <- diff(stats::quantile(x, c(0.75, 0.25), na.rm = TRUE))
  sdX <- abs(iqrX) / 1.5

  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(sdX) || sdX <= 0) {
    sdX <- .Machine$double.eps
  }

  sdX
}

#' @keywords internal
.simBandwidthExcMin <- function(x) {
  x <- x[is.finite(x)]
  x[x > min(x, na.rm = TRUE)]
}


#' @keywords internal
.simBandwidthCapForCpUnsLoc <- function(
  x_stim,
  x_uns,
  capStimRange
) {
  if (!capStimRange || length(x_stim) < 2L) {
    return(list(
      x_stim = x_stim,
      x_uns = x_uns,
      max_dens_x = NA_real_
    ))
  }

  range_stim <- range(x_stim, na.rm = TRUE)
  max_dens_x <- max(x_stim, na.rm = TRUE) - 0.05 * diff(range_stim)

  list(
    x_stim = pmin(x_stim, max_dens_x),
    x_uns = pmin(x_uns, max_dens_x),
    max_dens_x = max_dens_x
  )
}


#' @keywords internal
.simBandwidthGetTrans <- function(transformation) {
  if (exists(".simMiscGetTrans", mode = "function")) {
    return(.simMiscGetTrans(transformation))
  }

  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("Transformation not recognised: ", transformation)
  )
}
