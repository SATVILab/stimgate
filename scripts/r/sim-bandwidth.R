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
.simBwMtdLabVec <- c(
  "hpi0" = "HPI0",
  "hpi1" = "HPI1",
  "hpi2" = "HPI2",
  "hpi3" = "HPI3",
  "nrd0" = "NRD0",
  "sj" = "SJ",
  "hpi0Norm" = "HPI0Norm",
  "hpi1Norm" = "HPI1Norm",
  "hpi2Norm" = "HPI2Norm",
  "hpi3Norm" = "HPI3Norm",
  "nrd0Norm" = "NRD0Norm",
  "sjNorm" = "SJNorm"
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


.simBandwidthFiniteNumeric <- function(x) {
  x <- suppressWarnings(as.numeric(x)[1])
  is.finite(x)
}

.simBandwidthHasAdaptiveSetting <- function(
  bwAdaptive = FALSE,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL
) {
  isTRUE(bwAdaptive) ||
    .simBandwidthFiniteNumeric(bwAdaptiveCore) ||
    .simBandwidthFiniteNumeric(bwAdaptiveExtra) ||
    .simBandwidthFiniteNumeric(bwAdaptiveCrossover)
}

.simBandwidthAdaptiveBwMtd <- function(
  bwMtd,
  bwAdaptive = FALSE,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL
) {
  if (is.null(bwMtd) || length(bwMtd) == 0L) {
    bwMtd <- "hpi1Norm"
  }

  bwMtd <- as.character(bwMtd)[1]

  if (is.na(bwMtd) || !nzchar(bwMtd)) {
    bwMtd <- "hpi1Norm"
  }

  if (
    .simBandwidthHasAdaptiveSetting(
      bwAdaptive = bwAdaptive,
      bwAdaptiveCore = bwAdaptiveCore,
      bwAdaptiveExtra = bwAdaptiveExtra,
      bwAdaptiveCrossover = bwAdaptiveCrossover
    ) &&
      !grepl("Norm$", bwMtd)
  ) {
    bwMtd <- paste0(bwMtd, "Norm")
  }

  bwMtd
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


.simBandwidthNegativeShoulderWidth <- function(
  x,
  bw,
  heightFrac = 0.15,
  densityN = 512L
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (
    length(x) < 2L ||
      length(unique(x)) < 2L ||
      length(bw) != 1L ||
      !is.finite(bw) ||
      bw <= 0 ||
      length(heightFrac) != 1L ||
      !is.finite(heightFrac) ||
      heightFrac <= 0 ||
      heightFrac >= 1
  ) {
    return(NA_real_)
  }

  dens <- stats::density(x, bw = bw, n = as.integer(densityN))
  peakIdx <- which.max(dens$y)
  if (peakIdx >= length(dens$y)) {
    return(NA_real_)
  }

  rightIdx <- seq.int(peakIdx + 1L, length(dens$y))
  heightIdx <- rightIdx[dens$y[rightIdx] <= heightFrac * dens$y[[peakIdx]]]
  antimodeCandidates <- rightIdx[rightIdx < length(dens$y)]
  antimodeIdx <- antimodeCandidates[
    dens$y[antimodeCandidates] <= dens$y[antimodeCandidates - 1L] &
      dens$y[antimodeCandidates] < dens$y[antimodeCandidates + 1L]
  ]

  endpointIdx <- min(
    c(
      if (length(heightIdx) > 0L) heightIdx[[1]] else Inf,
      if (length(antimodeIdx) > 0L) antimodeIdx[[1]] else Inf,
      length(dens$x)
    )
  )
  width <- dens$x[[endpointIdx]] - dens$x[[peakIdx]]
  if (is.finite(width) && width > 0) width else NA_real_
}

.simBandwidthBsFreq <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  biasUnsWidthMultiplier = NULL,
  biasUnsWidthHeightFrac = 0.15,
  bw = NULL,
  bwAdaptive = FALSE,
  bwAdaptiveDensityN = NULL,
  bwAdaptivePadFrac = 0.15,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
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
  stimMeanShift = 0,
  stimSdMultiplier = 1,
  stimMeanShiftClusters = NULL,
  stimSdMultiplierClusters = NULL,
  covEvMin = 1,
  covEvMax = 2,
  tolClust = NULL,
  minCell = 1e2,
  maxPosProbX = Inf,
  gateQuant = c(0.25, 0.75),
  locProbCol = "pred",
  locMinPeakProb = 0.25,
  locEnforceShapeThreshold = FALSE,
  locDipAlpha = 0.2,
  locAntimodeHeightFrac = 1 / 6,
  locAntimodeLowRel = 0.25,
  locAntimodeLowAbs = 0.15,
  locFlatDerivFrac = 1 / 2,
  locFlatHardDerivFrac = 1 / 4,
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  gateCombn = "min",
  calcCytPosGates = FALSE
) {
  bwMtdGate <- .simBandwidthAdaptiveBwMtd(
    bwMtd = bwMtd,
    bwAdaptive = bwAdaptive,
    bwAdaptiveCore = bwAdaptiveCore,
    bwAdaptiveExtra = bwAdaptiveExtra,
    bwAdaptiveCrossover = bwAdaptiveCrossover
  )

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

    simArgs <- list(
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
    hasMismatch <- !identical(as.numeric(stimMeanShift), 0) ||
      !identical(as.numeric(stimSdMultiplier), 1) ||
      !is.null(stimMeanShiftClusters) ||
      !is.null(stimSdMultiplierClusters)
    if (hasMismatch) {
      if (!exists(".simCompareSimCytExperiment", mode = "function")) {
        stop(
          "Batch-mismatch simulations require .simCompareSimCytExperiment(). ",
          "Source scripts/r/sim-compare-freq_bs.R before running this scenario."
        )
      }
      simArgs <- c(simArgs, list(
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier,
        stimMeanShiftClusters = stimMeanShiftClusters,
        stimSdMultiplierClusters = stimSdMultiplierClusters
      ))
      outListExperiment <- do.call(.simCompareSimCytExperiment, simArgs)
    } else {
      outListExperiment <- do.call(simcyto::simCytExperiment, simArgs)
    }

    flowFrameList <- outListExperiment[["flowFrameList"]]
    labelsList <- outListExperiment[["labelsList"]]

    biasUnsNegativeWidth <- NA_real_
    biasUnsUse <- biasUns
    if (!is.null(biasUnsWidthMultiplier)) {
      if (is.null(bw) || length(bw) != 1L || !is.finite(bw) || bw <= 0) {
        stop("Width-based biasUns requires one finite positive fixed bandwidth.")
      }
      indUns <- seq.int(1L, length(flowFrameList), by = nCondition)
      xUns <- unlist(lapply(indUns, function(ind) {
        flowCore::exprs(flowFrameList[[ind]])[, 1L]
      }), use.names = FALSE)
      biasUnsNegativeWidth <- .simBandwidthNegativeShoulderWidth(
        x = xUns,
        bw = bw,
        heightFrac = biasUnsWidthHeightFrac
      )
      if (!is.finite(biasUnsNegativeWidth)) {
        stop("Could not calculate a finite negative-population shoulder width.")
      }
      biasUnsUse <- biasUnsWidthMultiplier * biasUnsNegativeWidth
    }

    fs <- as(flowFrameList, "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

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

    oldIntermediate <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA)
    on.exit(
      {
        if (is.na(oldIntermediate)) {
          Sys.unsetenv("STIMGATE_INTERMEDIATE")
        } else {
          Sys.setenv("STIMGATE_INTERMEDIATE" = oldIntermediate)
        }
      },
      add = TRUE
    )
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      calcCytPosGates = calcCytPosGates,
      biasUns = biasUnsUse,
      bw = bw,
      bwAdaptive = bwAdaptive,
      bwAdaptiveDensityN = bwAdaptiveDensityN,
      bwAdaptivePadFrac = bwAdaptivePadFrac,
      bwAdaptiveCore = bwAdaptiveCore,
      bwAdaptiveExtra = bwAdaptiveExtra,
      bwAdaptiveCrossover = bwAdaptiveCrossover,
      bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
      bwFallback = bwFallback,
      bwMin = bwMin,
      bwMax = bwMax,
      bwMtd = bwMtdGate,
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
      locEnforceShapeThreshold = locEnforceShapeThreshold,
      locDipAlpha = locDipAlpha,
      locAntimodeHeightFrac = locAntimodeHeightFrac,
      locAntimodeLowRel = locAntimodeLowRel,
      locAntimodeLowAbs = locAntimodeLowAbs,
      locFlatDerivFrac = locFlatDerivFrac,
      locFlatHardDerivFrac = locFlatHardDerivFrac,
      locMarginalPurityRel = locMarginalPurityRel,
      locMarginalCellBinRatio = locMarginalCellBinRatio,
      locMarginalRefQuantile = locMarginalRefQuantile,
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
        biasUns = biasUnsUse,
        biasUnsNegativeWidth = biasUnsNegativeWidth,
        biasUnsWidthMultiplier = if (is.null(biasUnsWidthMultiplier)) {
          NA_real_
        } else {
          biasUnsWidthMultiplier
        },
        stimMeanShift = stimMeanShift,
        stimSdMultiplier = stimSdMultiplier,
        stimMeanShiftClusters = if (is.null(stimMeanShiftClusters)) {
          NA_character_
        } else {
          paste(stimMeanShiftClusters, collapse = ",")
        },
        stimSdMultiplierClusters = if (is.null(stimSdMultiplierClusters)) {
          NA_character_
        } else {
          paste(stimSdMultiplierClusters, collapse = ",")
        },
        bw = if (is.null(bw)) NA_real_ else bw,
        bwAdaptive = bwAdaptive,
        bwAdaptiveDensityN = if (is.null(bwAdaptiveDensityN)) {
          NA_real_
        } else {
          bwAdaptiveDensityN
        },
        bwAdaptivePadFrac = bwAdaptivePadFrac,
        bwAdaptiveCore = if (is.null(bwAdaptiveCore)) {
          NA_real_
        } else {
          bwAdaptiveCore
        },
        bwAdaptiveExtra = if (is.null(bwAdaptiveExtra)) {
          NA_real_
        } else {
          bwAdaptiveExtra
        },
        bwAdaptiveCrossover = if (is.null(bwAdaptiveCrossover)) {
          NA_real_
        } else {
          bwAdaptiveCrossover
        },
        bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtdGate,
        bwMtdInput = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = if (is.null(bwCluster)) NA_real_ else bwCluster,
        tolClust = if (is.null(tolClust)) NA_real_ else tolClust,
        locEnforceShapeThreshold = locEnforceShapeThreshold,
        calcCytPosGates = calcCytPosGates,
        samplePerturbationSd = samplePerturbationSd,
        conditionPerturbationSd = conditionPerturbationSd,
        clusterPerturbationSd = clusterPerturbationSd,
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim
      ) |>
      dplyr::select(iter, chnl, sample, ind, dplyr::everything())
  })
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
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  normMtd = "moments",
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
    outListExperiment <- simcyto::simCytExperiment(
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
        bwFallback = bwFallback,
        normPeakMinRel = normPeakMinRel,
        normExtraFrac = normExtraFrac,
        normExtraMax = normExtraMax,
        normLambda = normLambda,
        normDensityN = normDensityN,
        normExcessBwMtd = normExcessBwMtd,
        normExcessNcell = normExcessNcell,
        normAdaptiveNcell = normAdaptiveNcell,
        normMtd = normMtd
      )

      bw_stim_norm_fallback <- isTRUE(attr(bw_stim, "normFallback"))
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
        bwFallback = bwFallback,
        normPeakMinRel = normPeakMinRel,
        normExtraFrac = normExtraFrac,
        normExtraMax = normExtraMax,
        normLambda = normLambda,
        normDensityN = normDensityN,
        normExcessBwMtd = normExcessBwMtd,
        normExcessNcell = normExcessNcell,
        normAdaptiveNcell = normAdaptiveNcell,
        normMtd = normMtd
      )

      bw_uns_norm_fallback <- isTRUE(attr(bw_uns, "normFallback"))
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

      bw_source <- dplyr::case_when(
        !is.null(bw) ~ "fixed",
        is.finite(bw_stim) &
          (!is.finite(bw_uns) || bw_stim <= bw_uns) ~ "stim",
        is.finite(bw_uns) ~ "unstim",
        TRUE ~ NA_character_
      )
      bw_norm_fallback <- dplyr::case_when(
        bw_source == "stim" ~ bw_stim_norm_fallback,
        bw_source == "unstim" ~ bw_uns_norm_fallback,
        TRUE ~ NA
      )

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
        bw_source = bw_source,
        bw_norm_fallback_stim = bw_stim_norm_fallback,
        bw_norm_fallback_uns = bw_uns_norm_fallback,
        bw_norm_fallback = bw_norm_fallback
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
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normMtd = "moments",
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

  bwMtdGate <- .simBandwidthAdaptiveBwMtd(
    bwMtd = bwMtd,
    bwAdaptive = TRUE
  )

  raw_tbl <- purrr::map_dfr(seq_len(nIter), function(iterNum) {
    outListExperiment <- simcyto::simCytExperiment(
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
          bwMtd = bwMtdGate,
          bwMin = bwMin,
          bwMax = bwMax,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax,
          bwFallback = bwFallback,
          normPeakMinRel = normPeakMinRel,
          normExtraFrac = normExtraFrac,
          normExtraMax = normExtraMax,
          normLambda = normLambda,
          normDensityN = normDensityN,
          normExcessBwMtd = normExcessBwMtd,
          normExcessNcell = normExcessNcell,
          normAdaptiveNcell = normAdaptiveNcell,
          bwAdaptiveCore = bwAdaptiveCore,
          bwAdaptiveExtra = bwAdaptiveExtra,
          bwAdaptiveCrossover = bwAdaptiveCrossover,
          bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
          normMtd = normMtd
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
        bw_mtd = bwMtdGate,
        bw_mtd_input = bwMtd,
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
      n_error = if ("error" %in% colnames(.data)) {
        sum(!is.na(.data$error %||% NA_character_))
      } else {
        NA_integer_
      },
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
.simBandwidthEnsureCurrentCheckout <- function(pathRoot = NULL) {
  if (is.null(pathRoot)) {
    pathRoot <- normalizePath(".", winslash = "/", mustWork = FALSE)
  }
  pathRoot <- normalizePath(pathRoot, winslash = "/", mustWork = FALSE)

  if (!nzchar(pathRoot)) {
    return(invisible(FALSE))
  }

  if (requireNamespace("stimgate", quietly = TRUE)) {
    ns <- tryCatch(asNamespace("stimgate"), error = function(e) NULL)
    if (!is.null(ns)) {
      ns_path <- normalizePath(
        getNamespaceInfo(ns, "path"),
        winslash = "/",
        mustWork = FALSE
      )
      if (isTRUE(identical(ns_path, pathRoot))) {
        return(invisible(TRUE))
      }
    }
  }

  suppressMessages(devtools::load_all(pathRoot, quiet = TRUE))
  invisible(TRUE)
}

.simBandwidthFiniteMean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x_finite <- x[is.finite(x)]

  if (length(x_finite) == 0L) {
    return(NA_real_)
  }

  mean(x_finite)
}

.simBandwidthResolveCurrentBwCalc <- function(bwMtd, adaptive = FALSE) {
  bwMtd <- as.character(bwMtd)[1]
  ns <- tryCatch(asNamespace("stimgate"), error = function(e) NULL)

  if (!is.null(ns) && exists(".bwCalcOne", mode = "function", envir = ns)) {
    return(get(".bwCalcOne", mode = "function", envir = ns))
  }

  if (isTRUE(grepl("Norm$", bwMtd)) || isTRUE(adaptive)) {
    stop(
      "stimgate::.bwCalcOne is unavailable for the requested bandwidth method '",
      bwMtd,
      "'. Ensure workers are initialised against the current package checkout.",
      call. = FALSE
    )
  }

  .simBandwidthBwOneBaseLegacy
}

.simBandwidthBwOne <- function(
  x,
  bwMtd,
  bwMin,
  bwMax,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwFallback,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normMtd = "moments",
  adaptive = FALSE
) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(.simBandwidthBwFallbackOrNa(bwFallback))
  }

  bwCalcFun <- .simBandwidthResolveCurrentBwCalc(
    bwMtd = bwMtd,
    adaptive = adaptive
  )

  if (identical(bwCalcFun, .simBandwidthBwOneBaseLegacy)) {
    bw_calc <- bwCalcFun(
      x = x,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax
    )
  } else {
    bw_calc <- bwCalcFun(
      x = x,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      normPeakMinRel = normPeakMinRel,
      normExtraFrac = normExtraFrac,
      normExtraMax = normExtraMax,
      normLambda = normLambda,
      normDensityN = normDensityN,
      normExcessBwMtd = normExcessBwMtd,
      normExcessNcell = normExcessNcell,
      normAdaptiveNcell = normAdaptiveNcell,
      bwAdaptiveCore = bwAdaptiveCore,
      bwAdaptiveExtra = bwAdaptiveExtra,
      bwAdaptiveCrossover = bwAdaptiveCrossover,
      bwAdaptiveTransitionWidth = bwAdaptiveTransitionWidth,
      normMtd = normMtd,
      adaptive = adaptive
    )
  }

  norm_fallback <- isTRUE(attr(bw_calc, "normFallback"))
  bw_calc <- suppressWarnings(as.numeric(bw_calc)[1])

  if (!is.finite(bw_calc) || bw_calc <= 0) {
    return(structure(
      .simBandwidthBwFallbackOrNa(bwFallback),
      normFallback = norm_fallback
    ))
  }

  if (.simBandwidthIsFiniteScalar(bwMin)) {
    bw_calc <- max(as.numeric(bwMin)[1], bw_calc)
  }
  if (.simBandwidthIsFiniteScalar(bwMax)) {
    bw_calc <- min(as.numeric(bwMax)[1], bw_calc)
  }

  structure(
    bw_calc,
    normFallback = norm_fallback
  )
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

  bw_calc <- switch(bwMtdBase,
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
  if (is.function(transformation)) {
    return(transformation)
  }

  switch(transformation,
    "gamma" = simcyto::simCytTransformGamma(),
    "gamma_fixed_mean_and_spread" = ,
    "gammaFixed" = simcyto::simCytTransformGammaFixed(),
    "gaussian" = simcyto::simCytTransformGaussian(),
    "identity" = simcyto::simCytTransformIdentity(),
    "skew" = simcyto::simCytTransformSkew(),
    simcyto::simCytGetTransformation(transformation)
  )
}
