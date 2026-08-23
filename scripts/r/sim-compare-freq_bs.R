.simCompareCacheEnv <- new.env(parent = emptyenv())

#' @keywords internal
.simCompareAddMissingColumns <- function(.data, cols) {
  for (nm in names(cols)) {
    if (!nm %in% names(.data)) {
      .data[[nm]] <- cols[[nm]]
    }
  }
  .data
}

#' @keywords internal
.simCompareGetTrans <- function(transformation) {
  if (exists(".simBandwidthGetTrans", mode = "function")) {
    return(.simBandwidthGetTrans(transformation))
  }
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

#' @keywords internal
.simCompareSampleFromInd <- function(ind, nCondition) {
  if (exists(".simBandwidthSampleFromInd", mode = "function")) {
    return(.simBandwidthSampleFromInd(ind, nCondition))
  }
  ind_num <- suppressWarnings(as.numeric(ind))
  as.character(((ind_num - 1) %/% nCondition) + 1)
}

#' @keywords internal
.simCompareReadLocDetails <- function(pathProject, nSample, nCondition) {
  if (!exists(".simBandwidthReadLocDetails", mode = "function")) {
    stop(
      ".simBandwidthReadLocDetails() must be available before running ",
      "the StimGate comparison. Source sim-bandwidth.R first."
    )
  }
  .simBandwidthReadLocDetails(
    pathProject = pathProject,
    nSample = nSample,
    nCondition = nCondition
  )
}

#' Moving mean with leading NA values, matching fbeta.py
#'
#' Locate the fbeta Python script
#'
#' @keywords internal
.simCompareFbetaPath <- function(pathFbeta = NULL) {
  if (!is.null(pathFbeta)) {
    return(pathFbeta)
  }
  if (!requireNamespace("projr", quietly = TRUE)) {
    stop(
      "pathFbeta was not supplied and projr is not available. ",
      "Pass pathFbeta explicitly."
    )
  }
  projr::projr_path_get("project", "scripts", "python", "fbeta.py")
}

#' Patch a temporary copy of fbeta.py only for Python 3 / NumPy compatibility
#'
#' The thresholding itself is still called from the fbeta.py implementation.
#' This wrapper only fixes import/syntax issues that prevent reticulate from
#' loading the historical script in a modern Python session.
#'
#' @keywords internal
.simCompareFbetaCompatPath <- function(pathFbeta, patchPy2Compat = TRUE) {
  if (!file.exists(pathFbeta)) {
    stop("Could not find fbeta.py at: ", pathFbeta)
  }

  if (!isTRUE(patchPy2Compat)) {
    return(pathFbeta)
  }

  pyTxt <- readLines(pathFbeta, warn = FALSE)

  pyTxt <- gsub("normed\\s*=\\s*True", "density=True", pyTxt)
  pyTxt <- gsub("calculate_fscores", "calculate_fscore", pyTxt, fixed = TRUE)
  # Suppress expected divide-by-zero/invalid warnings when both precision
  # and recall are zero. The original code subsequently converts NaN F-scores
  # to zero, so this does not change the numerical result.
  pyTxt <- gsub(
    "    fscores = (1+beta*beta)*(precision*recall)/(beta*beta*precision + recall)",
    paste(
      "    with np.errstate(divide='ignore', invalid='ignore'):",
      "        fscores = (1+beta*beta)*(precision*recall)/(beta*beta*precision + recall)",
      sep = "\n"
    ),
    pyTxt,
    fixed = TRUE
  )

  # The plotting helpers use Python 2 print syntax. They are not used for the
  # threshold, but Python 3 still has to parse them on import.
  pyTxt <- gsub(
    '^([[:space:]]*)print "([^"]*)"$',
    '\\1print("\\2")',
    pyTxt
  )
  pyTxt <- gsub(
    "^([[:space:]]*)print '([^']*)'$",
    '\\1print("\\2")',
    pyTxt
  )
  pyTxt <- gsub(
    "^([[:space:]]*)print '([^']*)', (.*)$",
    '\\1print("\\2", \\3)',
    pyTxt
  )

  # Plotting dependencies are not needed by get_positivity_threshold().
  pyTxt <- gsub(
    "^import matplotlib as mpl$",
    paste(
      "try:",
      "    import matplotlib as mpl",
      "except Exception:",
      "    mpl = None",
      sep = "\n"
    ),
    pyTxt
  )
  pyTxt <- gsub(
    "^from fcm\\.graphics import bilinear_interpolate$",
    paste(
      "try:",
      "    from fcm.graphics import bilinear_interpolate",
      "except Exception:",
      "    bilinear_interpolate = None",
      sep = "\n"
    ),
    pyTxt
  )
  pyTxt <- gsub(
    "^from fcm\\.core\\.transforms import _logicle as logicle$",
    paste(
      "try:",
      "    from fcm.core.transforms import _logicle as logicle",
      "except Exception:",
      "    def logicle(x, *args, **kwargs):",
      "        return x",
      sep = "\n"
    ),
    pyTxt
  )

  pyTxt <- gsub(
    're\\.sub\\("\\\\\\.fcs","",fileName\\)',
    're.sub(r"\\\\.fcs","",fileName)',
    pyTxt
  )

  pathTmp <- file.path(
    tempdir(),
    paste0(
      "fbeta_reticulate_",
      Sys.getpid(),
      "_",
      as.integer(stats::runif(1, 1, 1e9)),
      ".py"
    )
  )
  writeLines(pyTxt, pathTmp)
  pathTmp
}

#' Call the fbeta Python implementation via reticulate
#'
#' @keywords internal
.simCompareFbetaThreshold <- function(
    xUns,
    xStim,
    pathFbeta = NULL,
    patchPy2Compat = TRUE,
    beta = 0.8,
    theta = 2,
    width = 10,
    numBins = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to call fbeta.py.")
  }

  xUns <- as.numeric(xUns)
  xStim <- as.numeric(xStim)
  xUns <- xUns[is.finite(xUns)]
  xStim <- xStim[is.finite(xStim)]

  if (length(xUns) < 2L || length(xStim) < 2L) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_too_few_cells"
    ))
  }

  pathFbeta <- .simCompareFbetaPath(pathFbeta)
  pathInfo <- file.info(pathFbeta)
  cacheName <- make.names(paste(
    "fbeta",
    normalizePath(pathFbeta, winslash = "/", mustWork = FALSE),
    patchPy2Compat,
    pathInfo$mtime,
    sep = "_"
  ))

  if (!exists(cacheName, envir = .simCompareCacheEnv, inherits = FALSE)) {
    pathPyUse <- .simCompareFbetaCompatPath(
      pathFbeta = pathFbeta,
      patchPy2Compat = patchPy2Compat
    )

    pyEnv <- new.env(parent = emptyenv())
    reticulate::source_python(pathPyUse, envir = pyEnv)

    if (!exists("get_positivity_threshold", envir = pyEnv, inherits = FALSE)) {
      stop("fbeta.py did not define get_positivity_threshold().")
    }

    assign(cacheName, pyEnv, envir = .simCompareCacheEnv)
  }

  pyEnv <- get(cacheName, envir = .simCompareCacheEnv, inherits = FALSE)

  negMat <- matrix(xUns, ncol = 1L)
  posMat <- matrix(xStim, ncol = 1L)

  out <- pyEnv$get_positivity_threshold(
    neg = negMat,
    pos = posMat,
    channelIndex = 0L,
    beta = beta,
    theta = theta,
    width = as.integer(width),
    numBins = if (is.null(numBins)) NULL else as.integer(numBins)
  )

  threshold <- suppressWarnings(as.numeric(out[["threshold"]]))[1]
  fscores <- suppressWarnings(as.numeric(out[["fscores"]]))
  metric <- if (length(fscores) > 0L && any(is.finite(fscores))) {
    max(fscores, na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    threshold = threshold,
    thresholdMetric = metric,
    thresholdOrigin = if (is.finite(threshold)) {
      "calculated"
    } else {
      "failed_no_cutpoint"
    },
    fbeta = out
  )
}

#' Call the cytoUtils tailgate implementation
#'
#' If bandwidth is NULL, cytoUtils:::.cytokine_cutpoint() forwards NULL to
#' cytoUtils:::.deriv_density(), whose default behaviour is to estimate the
#' bandwidth with ks::hpi().
#'
#' @keywords internal
.simCompareTailgateThreshold <- function(
    x,
    tailgateSourceFiles = NULL,
    adjust = 1,
    bandwidth = NULL,
    numPeaks = 1,
    refPeak = 1,
    method = c("firstDeriv", "secondDeriv"),
    tol = 1e-2,
    side = "right",
    strict = FALSE,
    autoTol = FALSE) {
  method <- match.arg(method)
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_too_few_unique_cells"
    ))
  }

  if (!requireNamespace("cytoUtils", quietly = TRUE)) {
    stop("Package 'cytoUtils' is required for tailgate comparisons.")
  }

  threshold <- cytoUtils:::.cytokine_cutpoint(
    x = x,
    adjust = adjust,
    num_peaks = numPeaks,
    ref_peak = refPeak,
    method = dplyr::case_when(
      method == "firstDeriv" ~ "first_deriv",
      method == "secondDeriv" ~ "second_deriv",
      .unmatched = "error"
    ),
    tol = tol,
    side = side,
    strict = strict,
    auto_tol = autoTol,
    bandwidth = bandwidth
  )

  list(
    threshold = as.numeric(threshold)[1],
    thresholdMetric = NA_real_,
    thresholdOrigin = if (is.finite(threshold)) {
      "calculated"
    } else {
      "failed_no_cutpoint"
    }
  )
}

.simCompareHighValueGate <- function(xStim, xUns, margin = 0.05) {
  x <- c(xStim, xUns)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(Inf)
  }
  rng <- range(x, na.rm = TRUE)
  rng[2] + max(1, diff(rng)) * margin
}

.simCompareResolveClusterMismatchNames <- function(clusterSpec) {
  if (is.null(clusterSpec) || length(clusterSpec) == 0L) {
    return(character(0))
  }

  spec <- if (length(clusterSpec) == 1L && is.character(clusterSpec)) {
    strsplit(as.character(clusterSpec), ",", fixed = TRUE)[[1]]
  } else {
    as.character(clusterSpec)
  }

  spec <- trimws(spec)
  spec <- spec[nzchar(spec)]
  unique(spec)
}

.simCompareApplyClusterMismatch <- function(
    outListExperiment,
    stimMeanShift = 0,
    stimSdMultiplier = 1,
    stimMeanShiftClusters = NULL,
    stimSdMultiplierClusters = NULL) {
  if (is.null(outListExperiment) || is.null(outListExperiment[["flowFrameList"]])) {
    return(outListExperiment)
  }

  shiftLabels <- .simCompareResolveClusterMismatchNames(stimMeanShiftClusters)
  sdLabels <- .simCompareResolveClusterMismatchNames(stimSdMultiplierClusters)

  if (length(shiftLabels) == 0L && length(sdLabels) == 0L) {
    return(outListExperiment)
  }

  flowFrameList <- outListExperiment[["flowFrameList"]]
  labelsList <- outListExperiment[["labelsList"]]
  nCondition <- if (!is.null(outListExperiment[["nCondition"]])) {
    outListExperiment[["nCondition"]]
  } else {
    2L
  }

  for (idx in seq_along(flowFrameList)) {
    condIndex <- ((idx - 1L) %% nCondition) + 1L
    if (condIndex == 1L) {
      next
    }

    labelVec <- labelsList[[idx]]
    expr <- flowCore::exprs(flowFrameList[[idx]])
    if (ncol(expr) == 0L || length(labelVec) != nrow(expr)) {
      next
    }

    if (length(shiftLabels) > 0L) {
      shiftMask <- labelVec %in% shiftLabels
      if (any(shiftMask)) {
        expr[shiftMask, 1L] <- expr[shiftMask, 1L] + stimMeanShift
        flowCore::exprs(flowFrameList[[idx]]) <- expr
      }
    }

    if (length(sdLabels) > 0L) {
      sdMask <- labelVec %in% sdLabels
      if (any(sdMask)) {
        negVals <- expr[sdMask, 1L]
        center <- mean(negVals)
        expr[sdMask, 1L] <- center + (negVals - center) * stimSdMultiplier
        flowCore::exprs(flowFrameList[[idx]]) <- expr
      }
    }
  }

  outListExperiment[["flowFrameList"]] <- flowFrameList
  outListExperiment
}

.simCompareSimCytExperiment <- function(
    nSample = NULL,
    nMarker = NULL,
    nCondition = NULL,
    nCluster = NULL,
    nCellByCondition = NULL,
    transformationFunc = NULL,
    mixtureType = "gaussianOnly",
    meanExprMat = NA,
    clusterLabelVec = NA,
    probVecUns = NULL,
    probExact = FALSE,
    probResponseVecByStimCondition = NULL,
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    covEvMin = 1,
    covEvMax = 2,
    stimMeanShift = 0,
    stimSdMultiplier = 1,
    stimMeanShiftClusters = NULL,
    stimSdMultiplierClusters = NULL,
    scenario = NULL) {
  callArgs <- list(
    nSample = nSample,
    nMarker = nMarker,
    nCondition = nCondition,
    nCluster = nCluster,
    nCellByCondition = nCellByCondition,
    transformationFunc = transformationFunc,
    mixtureType = mixtureType,
    meanExprMat = meanExprMat,
    clusterLabelVec = clusterLabelVec,
    probVecUns = probVecUns,
    probExact = probExact,
    probResponseVecByStimCondition = probResponseVecByStimCondition,
    samplePerturbationSd = samplePerturbationSd,
    conditionPerturbationSd = conditionPerturbationSd,
    clusterPerturbationSd = clusterPerturbationSd,
    covEvMin = covEvMin,
    covEvMax = covEvMax,
    stimMeanShift = stimMeanShift,
    stimSdMultiplier = stimSdMultiplier,
    scenario = scenario
  )

  simcytoArgs <- names(formals(simcyto::simCytExperiment))
  clusterShiftSupported <- !is.null(simcytoArgs) &&
    "stimMeanShiftClusters" %in% simcytoArgs
  clusterSdSupported <- !is.null(simcytoArgs) &&
    "stimSdMultiplierClusters" %in% simcytoArgs

  if (clusterShiftSupported && !is.null(stimMeanShiftClusters)) {
    callArgs$stimMeanShiftClusters <- stimMeanShiftClusters
  }
  if (clusterSdSupported && !is.null(stimSdMultiplierClusters)) {
    callArgs$stimSdMultiplierClusters <- stimSdMultiplierClusters
  }

  result <- do.call(simcyto::simCytExperiment, callArgs)

  .simCompareApplyClusterMismatch(
    outListExperiment = result,
    stimMeanShift = stimMeanShift,
    stimSdMultiplier = stimSdMultiplier,
    stimMeanShiftClusters = stimMeanShiftClusters,
    stimSdMultiplierClusters = stimSdMultiplierClusters
  )
}

#' @keywords internal
.simCompareEstimateFromThreshold <- function(
    xStim,
    xUns,
    threshold,
    fallbackHighValue = TRUE,
    fallbackMargin = 0.05) {
  xStim <- as.numeric(xStim)
  xUns <- as.numeric(xUns)
  nCellStim <- length(xStim)
  nCellUns <- length(xUns)

  thresholdUsed <- as.numeric(threshold)[1]
  usedFallback <- FALSE
  if (!is.finite(thresholdUsed) && isTRUE(fallbackHighValue)) {
    thresholdUsed <- .simCompareHighValueGate(
      xStim = xStim,
      xUns = xUns,
      margin = fallbackMargin
    )
    usedFallback <- TRUE
  }

  nPosStim <- if (is.finite(thresholdUsed)) {
    sum(xStim > thresholdUsed, na.rm = TRUE)
  } else {
    NA_integer_
  }
  nPosUns <- if (is.finite(thresholdUsed)) {
    sum(xUns > thresholdUsed, na.rm = TRUE)
  } else {
    NA_integer_
  }

  propStim <- nPosStim / nCellStim
  propUns <- nPosUns / nCellUns
  propRespEst <- propStim - propUns

  list(
    threshold = thresholdUsed,
    thresholdFallbackUsed = usedFallback,
    nCellStim = nCellStim,
    nCellUns = nCellUns,
    nPosStim = nPosStim,
    nPosUns = nPosUns,
    propStim = propStim,
    propUns = propUns,
    propRespEst = propRespEst
  )
}

#' @keywords internal
.simCompareTruthTable <- function(
    labelsList,
    nSample,
    nCondition,
    chnl = "F1") {
  purrr::map_df(seq_len(nSample), function(sampleCurr) {
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
        chnl = chnl,
        propStimTruth = propStimTruth,
        propUnsTruth = propUnsTruth,
        propRespTruth = propStimTruth - propUnsTruth
      )
    })
  })
}

#' @keywords internal
.simCompareAlternativeRows <- function(
    flowFrameList,
    labelsList,
    nSample,
    nCondition,
    chnl = "F1",
    biasUns = 0,
    pathFbeta = NULL,
    fbetaPatchPy2Compat = TRUE,
    fbetaBeta = 0.8,
    fbetaTheta = 2,
    fbetaWidth = 10,
    fbetaNumBins = NULL,
    tailgateX = c("stim", "unstim", "combined"),
    tailgateSourceFiles = NULL,
    tailgateAdjust = 1,
    tailgateBandwidth = NULL,
    tailgateNumPeaks = 1,
    tailgateRefPeak = 1,
    tailgateMethod = c("firstDeriv", "secondDeriv"),
    tailgateTol = 1e-2,
    tailgateSide = "right",
    tailgateAutoTol = TRUE,
    fallbackHighValue = TRUE,
    fallbackMargin = 0.05) {
  tailgateX <- match.arg(tailgateX)
  tailgateMethod <- match.arg(tailgateMethod)

  truthTbl <- .simCompareTruthTable(
    labelsList = labelsList,
    nSample = nSample,
    nCondition = nCondition,
    chnl = chnl
  )

  out <- purrr::map_df(seq_len(nSample), function(sampleCurr) {
    indUns <- (sampleCurr - 1L) * nCondition + 1L
    indStimVec <- seq.int(indUns + 1L, sampleCurr * nCondition)

    xUnsRaw <- as.numeric(flowCore::exprs(flowFrameList[[indUns]])[, chnl])
    # Competitor methods receive raw unstimulated data and do not inherit StimGate's biasUns
    xUnsFbeta <- xUnsRaw
    xUnsTailgate <- xUnsRaw

    purrr::map_df(indStimVec, function(indStim) {
      xStim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, chnl])

      fbetaObj <- tryCatch(
        .simCompareFbetaThreshold(
          xUns = xUnsFbeta,
          xStim = xStim,
          pathFbeta = pathFbeta,
          patchPy2Compat = fbetaPatchPy2Compat,
          beta = fbetaBeta,
          theta = fbetaTheta,
          width = fbetaWidth,
          numBins = fbetaNumBins
        ),
        error = function(e) {
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", e$message)
          )
        }
      )

      fbetaEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsFbeta,
        threshold = fbetaObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin
      )

      xTail <- switch(tailgateX,
        "stim" = xStim,
        "unstim" = xUnsTailgate,
        "combined" = c(xUnsTailgate, xStim)
      )

      tailgateObj <- tryCatch(
        .simCompareTailgateThreshold(
          x = xTail,
          tailgateSourceFiles = tailgateSourceFiles,
          adjust = tailgateAdjust,
          bandwidth = tailgateBandwidth,
          numPeaks = tailgateNumPeaks,
          refPeak = tailgateRefPeak,
          method = tailgateMethod,
          tol = tailgateTol,
          side = tailgateSide,
          strict = FALSE,
          autoTol = tailgateAutoTol
        ),
        error = function(e) {
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", e$message)
          )
        }
      )

      tailgateEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsTailgate,
        threshold = tailgateObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin
      )

      tibble::tibble(
        sample = as.character(sampleCurr),
        ind = as.character(indStim),
        chnl = chnl,
        approach = c("fbeta", "tailgate"),
        method = c("fbeta", "tailgate"),
        threshold = c(fbetaEst$threshold, tailgateEst$threshold),
        thresholdOrigin = c(
          fbetaObj$thresholdOrigin,
          tailgateObj$thresholdOrigin
        ),
        gateReturnPoint = c(
          if (isTRUE(fbetaEst$thresholdFallbackUsed)) {
            "fbeta_fallback_high_value"
          } else {
            "fbeta_calculated"
          },
          if (isTRUE(tailgateEst$thresholdFallbackUsed)) {
            "tailgate_fallback_high_value"
          } else {
            "tailgate_calculated"
          }
        ),
        thresholdMetric = c(
          fbetaObj$thresholdMetric %||% NA_real_,
          tailgateObj$thresholdMetric %||% NA_real_
        ),
        thresholdFallbackUsed = c(
          fbetaEst$thresholdFallbackUsed,
          tailgateEst$thresholdFallbackUsed
        ),
        nCellStim = c(fbetaEst$nCellStim, tailgateEst$nCellStim),
        nCellUns = c(fbetaEst$nCellUns, tailgateEst$nCellUns),
        nPosStim = c(fbetaEst$nPosStim, tailgateEst$nPosStim),
        nPosUns = c(fbetaEst$nPosUns, tailgateEst$nPosUns),
        propStim = c(fbetaEst$propStim, tailgateEst$propStim),
        propUns = c(fbetaEst$propUns, tailgateEst$propUns),
        propRespEst = c(fbetaEst$propRespEst, tailgateEst$propRespEst),
        detailLevel = NA_character_,
        locGenerated = NA,
        locGeneratedDirect = NA,
        locSource = NA_character_,
        locReason = NA_character_,
        error = NA_character_
      )
    })
  })

  out |>
    dplyr::left_join(truthTbl, by = c("sample", "ind", "chnl"))
}

#' @keywords internal
.simCompareStimgateFailureRows <- function(truthTbl, errorMessage) {
  truthTbl |>
    dplyr::mutate(
      approach = "stimgate",
      method = "stimgate_error",
      threshold = NA_real_,
      thresholdOrigin = "error",
      gateReturnPoint = "stimgate_error",
      thresholdMetric = NA_real_,
      thresholdFallbackUsed = NA,
      nCellStim = NA_real_,
      nCellUns = NA_real_,
      nPosStim = NA_integer_,
      nPosUns = NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propRespEst = NA_real_,
      detailLevel = NA_character_,
      locGenerated = NA,
      locGeneratedDirect = NA,
      locSource = NA_character_,
      locReason = NA_character_,
      error = errorMessage
    )
}

#' @keywords internal
.simCompareStimgateRows <- function(
    gs,
    labelsList,
    pathProject,
    nSample,
    nCondition,
    nMarker,
    biasUns,
    bw,
    bwFallback = bw,
    bwMin = "none",
    bwMax = "none",
    bwMtd = "hpi1",
    bwAdj = 1,
    bwNcellMin = 1e2,
    bwNcellMax = 1e5,
    bwCluster = NULL,
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
    tolClust = NULL,
    locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE,
    includeLocCondition = FALSE,
    includeLocDetails = includeLocCondition) {
  truthTbl <- .simCompareTruthTable(
    labelsList = labelsList,
    nSample = nSample,
    nCondition = nCondition,
    chnl = "F1"
  )

  out <- tryCatch(
    {
      batchList <- lapply(seq_len(nSample), function(i) {
        seq((i - 1L) * nCondition + 1L, i * nCondition)
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

      invisible(stimgate::gateStim(
        .data = gs,
        pathProject = pathProject,
        popGate = "root",
        batchList = batchList,
        marker = paste0("MarkerF", seq_len(nMarker)),
        calcCytPosGates = calcCytPosGates,
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
        locEnforceShapeThreshold = locEnforceShapeThreshold,
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

      # Extract final cluster-refined StimGate gates and statistics
      gateTblFinal <- tryCatch(
        stimgate::getStimGates(pathProject),
        error = function(e) tibble::tibble()
      )
      statsTblFinal <- tryCatch(
        stimgate::getStimStats(pathProject),
        error = function(e) tibble::tibble()
      )

      stimgatePrimaryTbl <- purrr::map_df(seq_len(nSample), function(sampleCurr) {
        indUns <- (sampleCurr - 1L) * nCondition + 1L
        indStimVec <- seq.int(indUns + 1L, sampleCurr * nCondition)

        purrr::map_df(indStimVec, function(indStim) {
          ind_curr <- as.character(indStim)
          gRow <- if (nrow(gateTblFinal) > 0L && "ind" %in% names(gateTblFinal)) {
            gateTblFinal[as.character(gateTblFinal$ind) == ind_curr & gateTblFinal$chnl == "F1", , drop = FALSE]
          } else {
            tibble::tibble()
          }
          if (nrow(gRow) > 0L) {
            if (any(grepl("Clust$", gRow$gateName))) {
              gRow <- gRow[grepl("Clust$", gRow$gateName), , drop = FALSE]
            } else {
              gRow <- gRow[nrow(gRow), , drop = FALSE]
            }
          }

          sRow <- if (nrow(statsTblFinal) > 0L && "ind" %in% names(statsTblFinal)) {
            statsTblFinal[as.character(statsTblFinal$ind) == ind_curr & grepl("~\\+~", statsTblFinal$cytCombn), , drop = FALSE]
          } else {
            tibble::tibble()
          }
          if (nrow(sRow) > 0L) {
            if (any(grepl("Clust$", sRow$gateName))) {
              sRow <- sRow[grepl("Clust$", sRow$gateName), , drop = FALSE]
            } else {
              sRow <- sRow[1L, , drop = FALSE]
            }
          }

          gateVal <- if (nrow(gRow) > 0L) suppressWarnings(as.numeric(unname(gRow$gate[[1]]))) else NA_real_
          gateNm <- if (nrow(gRow) > 0L) as.character(gRow$gateName[[1]]) else NA_character_

          nCellStimVal <- if (nrow(sRow) > 0L) suppressWarnings(as.numeric(sRow$nCellStim[[1]])) else NA_real_
          nCellUnsVal <- if (nrow(sRow) > 0L) suppressWarnings(as.numeric(sRow$nCellUns[[1]])) else NA_real_
          nPosStimVal <- if (nrow(sRow) > 0L) suppressWarnings(as.integer(sRow$countStim[[1]])) else NA_integer_
          nPosUnsVal <- if (nrow(sRow) > 0L) suppressWarnings(as.integer(sRow$countUns[[1]])) else NA_integer_
          propStimVal <- if (nrow(sRow) > 0L) suppressWarnings(as.numeric(sRow$propStim[[1]])) else NA_real_
          propUnsVal <- if (nrow(sRow) > 0L) suppressWarnings(as.numeric(sRow$propUns[[1]])) else NA_real_
          propBsVal <- if (nrow(sRow) > 0L) suppressWarnings(as.numeric(sRow$propBs[[1]])) else NA_real_

          isClustered <- grepl("Clust$", gateNm %||% "")

          tibble::tibble(
            sample = as.character(sampleCurr),
            ind = as.character(indStim),
            chnl = "F1",
            approach = "stimgate",
            method = "stimgate",
            threshold = gateVal,
            thresholdOrigin = if (is.finite(gateVal)) {
              if (isClustered) "calculated_clustered" else "calculated"
            } else {
              "failed_no_cutpoint"
            },
            gateReturnPoint = if (isClustered) "stimgate_clustered" else "stimgate_calculated",
            thresholdMetric = NA_real_,
            thresholdFallbackUsed = !is.finite(gateVal),
            nCellStim = nCellStimVal,
            nCellUns = nCellUnsVal,
            nPosStim = nPosStimVal,
            nPosUns = nPosUnsVal,
            propStim = propStimVal,
            propUns = propUnsVal,
            propRespEst = propBsVal,
            detailLevel = if (isClustered) "cluster_final" else "sample_final",
            locGenerated = is.finite(gateVal),
            locGeneratedDirect = !isClustered,
            locSource = if (isClustered) "cluster" else "sample",
            locReason = NA_character_,
            error = NA_character_
          )
        })
      })

      detailTbl <- if (isTRUE(includeLocDetails) || isTRUE(includeLocCondition)) {
        tryCatch(
          .simCompareReadLocDetails(
            pathProject = pathProject,
            nSample = nSample,
            nCondition = nCondition
          ),
          error = function(e) tibble::tibble()
        )
      } else {
        tibble::tibble()
      }

      if (nrow(detailTbl) > 0L) {
        if (!isTRUE(includeLocCondition)) {
          detailTbl <- detailTbl |>
            dplyr::filter(.data$detailLevel %in% "sample")
        } else {
          detailTbl <- detailTbl |>
            dplyr::filter(.data$detailLevel %in% c("condition", "sample"))
        }

        detailTbl <- .simCompareAddMissingColumns(
          detailTbl,
          list(
            method = NA_character_,
            propRespEst = NA_real_,
            propBsEst = NA_real_,
            propBs = NA_real_,
            threshold = NA_real_,
            thresholdOrigin = NA_character_,
            gateReturnPoint = NA_character_,
            nCellStim = NA_real_,
            nCellUns = NA_real_,
            nPosStim = NA_integer_,
            nPosUns = NA_integer_,
            propStim = NA_real_,
            propUns = NA_real_,
            detailLevel = NA_character_,
            locGenerated = NA,
            locGeneratedDirect = NA,
            locSource = NA_character_,
            locReason = NA_character_
          )
        )

        detailTbl <- detailTbl |>
          dplyr::mutate(
            approach = "stimgate",
            method = paste0("stimgate_", .data$method),
            propRespEst = dplyr::coalesce(
              suppressWarnings(as.numeric(.data$propRespEst)),
              suppressWarnings(as.numeric(.data$propBsEst)),
              suppressWarnings(as.numeric(.data$propBs))
            ),
            thresholdMetric = NA_real_,
            thresholdFallbackUsed = grepl(
              "fallback_high_value",
              .data$gateReturnPoint %||% ""
            ),
            error = NA_character_
          )
      } else {
        detailTbl <- tibble::tibble()
      }

      dplyr::bind_rows(stimgatePrimaryTbl, detailTbl) |>
        dplyr::left_join(truthTbl, by = c("sample", "ind", "chnl")) |>
        dplyr::select(
          sample,
          ind,
          chnl,
          approach,
          method,
          threshold,
          thresholdOrigin,
          gateReturnPoint,
          thresholdMetric,
          thresholdFallbackUsed,
          nCellStim,
          nCellUns,
          nPosStim,
          nPosUns,
          propStim,
          propUns,
          propRespEst,
          propStimTruth,
          propUnsTruth,
          propRespTruth,
          detailLevel,
          locGenerated,
          locGeneratedDirect,
          locSource,
          locReason,
          error,
          dplyr::everything()
        )
    },
    error = function(e) {
      .simCompareStimgateFailureRows(
        truthTbl = truthTbl,
        errorMessage = e$message
      )
    }
  )

  out
}

#' Compare StimGate, fbeta and tailgate estimates on the same simulated data
#'
#' @keywords internal
.simCompareFreqBs <- function(
    nSample,
    nMarker,
    nCondition,
    nCluster,
    nIter,
    biasUns,
    bw,
    bwFallback = bw,
    bwMin = "none",
    bwMax = "none",
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
    locEnforceShapeThreshold = FALSE,
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
    includeLocCondition = FALSE,
    includeLocDetails = includeLocCondition,
    pathFbeta = NULL,
    fbetaPatchPy2Compat = TRUE,
    fbetaBeta = 0.8,
    fbetaTheta = 2,
    fbetaWidth = 10,
    fbetaNumBins = NULL,
    tailgateX = c("stim", "unstim", "combined"),
    tailgateSourceFiles = NULL,
    tailgateAdjust = 1,
    tailgateBandwidth = NULL,
    tailgateNumPeaks = 1,
    tailgateRefPeak = 1,
    tailgateMethod = c("firstDeriv", "secondDeriv"),
    tailgateTol = 1e-2,
    tailgateSide = "right",
    tailgateAutoTol = FALSE,
    fallbackHighValue = TRUE,
    fallbackMargin = 0.05,
    stimMeanShift = 0,
    stimSdMultiplier = 1,
    stimMeanShiftClusters = NULL,
    stimSdMultiplierClusters = NULL,
    pathProject = NULL) {
  if (!identical(as.integer(nMarker), 1L)) {
    stop("This comparison helper currently expects nMarker = 1.")
  }
  if (!identical(as.integer(nCondition), 2L)) {
    stop("This comparison helper currently expects nCondition = 2.")
  }
  if (!identical(as.integer(nCluster), 2L)) {
    stop("This comparison helper currently expects nCluster = 2.")
  }

  tailgateX <- match.arg(tailgateX)
  tailgateMethod <- match.arg(tailgateMethod)

  purrr::map_df(seq_len(nIter), function(iterNum) {
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simCompareGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    outListExperiment <- .simCompareSimCytExperiment(
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
      covEvMax = covEvMax,
      stimMeanShift = stimMeanShift,
      stimSdMultiplier = stimSdMultiplier,
      stimMeanShiftClusters = stimMeanShiftClusters,
      stimSdMultiplierClusters = stimSdMultiplierClusters
    )

    flowFrameList <- outListExperiment[["flowFrameList"]]
    labelsList <- outListExperiment[["labelsList"]]
    fs <- as(flowFrameList, "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

    pathProjectUse <- if (!is.null(pathProject) && nzchar(pathProject)) {
      if (identical(as.integer(nIter), 1L)) {
        pathProject
      } else {
        file.path(pathProject, paste0("iter-", iterNum))
      }
    } else {
      file.path(
        tempdir(),
        "stimgate-sim-compare",
        paste0(
          "pid-",
          Sys.getpid(),
          "-iter-",
          iterNum,
          "-",
          format(Sys.time(), "%Y%m%d%H%M%OS6"),
          "-",
          sample.int(1e9, 1)
        )
      )
    }
    on.exit(
      {
        if (dir.exists(pathProjectUse)) {
          unlink(pathProjectUse, recursive = TRUE)
        }
      },
      add = TRUE
    )
    if (dir.exists(pathProjectUse)) {
      unlink(pathProjectUse, recursive = TRUE)
    }
    dir.create(pathProjectUse, recursive = TRUE, showWarnings = FALSE)

    stimgateTbl <- .simCompareStimgateRows(
      gs = gs,
      labelsList = labelsList,
      pathProject = pathProjectUse,
      nSample = nSample,
      nCondition = nCondition,
      nMarker = nMarker,
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
      gateCombn = gateCombn,
      tolClust = tolClust,
      locEnforceShapeThreshold = locEnforceShapeThreshold,
      calcCytPosGates = calcCytPosGates,
      includeLocCondition = includeLocCondition,
      includeLocDetails = includeLocDetails
    )

    alternativeTbl <- .simCompareAlternativeRows(
      flowFrameList = flowFrameList,
      labelsList = labelsList,
      nSample = nSample,
      nCondition = nCondition,
      chnl = "F1",
      biasUns = biasUns,
      pathFbeta = pathFbeta,
      fbetaPatchPy2Compat = fbetaPatchPy2Compat,
      fbetaBeta = fbetaBeta,
      fbetaTheta = fbetaTheta,
      fbetaWidth = fbetaWidth,
      fbetaNumBins = fbetaNumBins,
      tailgateX = tailgateX,
      tailgateSourceFiles = tailgateSourceFiles,
      tailgateAdjust = tailgateAdjust,
      tailgateBandwidth = tailgateBandwidth,
      tailgateNumPeaks = tailgateNumPeaks,
      tailgateRefPeak = tailgateRefPeak,
      tailgateMethod = tailgateMethod,
      tailgateTol = tailgateTol,
      tailgateSide = tailgateSide,
      tailgateAutoTol = tailgateAutoTol,
      fallbackHighValue = fallbackHighValue,
      fallbackMargin = fallbackMargin
    )

    dplyr::bind_rows(stimgateTbl, alternativeTbl) |>
      dplyr::mutate(
        iter = iterNum,
        nCellStimSim = nCellStim,
        nCellUnsSim = nCellUns,
        biasUns = biasUns,
        bw = bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = bwCluster %||% NA_real_,
        tolClust = tolClust %||% NA_real_,
        locEnforceShapeThreshold = locEnforceShapeThreshold,
        calcCytPosGates = calcCytPosGates,
        samplePerturbationSd = samplePerturbationSd,
        conditionPerturbationSd = conditionPerturbationSd,
        clusterPerturbationSd = clusterPerturbationSd,
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim,
        fbetaBeta = fbetaBeta,
        fbetaTheta = fbetaTheta,
        fbetaWidth = fbetaWidth,
        fbetaNumBins = fbetaNumBins %||% NA_integer_,
        fbetaPatchPy2Compat = fbetaPatchPy2Compat,
        tailgateX = tailgateX,
        tailgateAdjust = tailgateAdjust,
        tailgateBandwidth = tailgateBandwidth %||% NA_real_,
        tailgateBandwidthEstimated = is.null(tailgateBandwidth),
        tailgateNumPeaks = tailgateNumPeaks,
        tailgateRefPeak = tailgateRefPeak,
        tailgateMethod = tailgateMethod,
        tailgateTol = tailgateTol,
        tailgateSide = tailgateSide,
        tailgateAutoTol = tailgateAutoTol,
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
        }
      ) |>
      dplyr::select(
        iter,
        chnl,
        sample,
        ind,
        approach,
        method,
        dplyr::everything()
      )
  })
}

#' Generate scenario output file path
#'
#' @keywords internal
.simCompareScenarioOutputPath <- function(
    sim_id,
    dirCache,
    sim_grid_chunk_index = NULL,
    sim_grid_n_chunks = NULL) {
  if (is.null(dirCache) || !nzchar(dirCache)) {
    return(character(0))
  }

  if (
    !is.null(sim_grid_chunk_index) && !is.na(sim_grid_chunk_index) &&
      !is.null(sim_grid_n_chunks) && !is.na(sim_grid_n_chunks) &&
      as.integer(sim_grid_n_chunks) > 1L
  ) {
    file.path(
      dirCache,
      sprintf(
        "compare_raw-chunk_%03d-of_%03d-sim_id_%06d.rds",
        as.integer(sim_grid_chunk_index),
        as.integer(sim_grid_n_chunks),
        as.integer(sim_id)
      )
    )
  } else {
    file.path(
      dirCache,
      sprintf(
        "compare_raw-sim_id_%06d.rds",
        as.integer(sim_id)
      )
    )
  }
}

#' Find scenario output files in a directory
#'
#' @keywords internal
.simCompareFindScenarioOutputs <- function(dirCache) {
  if (is.null(dirCache) || !dir.exists(dirCache)) {
    return(character(0))
  }

  files <- list.files(
    dirCache,
    pattern = "^(compare_raw.*|sim_scenario.*|sim_raw.*)sim_id_[0-9]+[.]rds$",
    full.names = TRUE
  )
  sort(unique(files))
}

#' Validate scenario cached output against grid row settings
#'
#' @keywords internal
.simCompareValidateScenarioCache <- function(
    cached,
    row,
    nSample = NULL,
    nIter = NULL,
    retryErrors = FALSE) {
  if (!is.data.frame(cached) || nrow(cached) == 0L) {
    return(FALSE)
  }

  has_error <- "error" %in% names(cached) &&
    any(!is.na(cached$error) & nzchar(as.character(cached$error)))
  if (isTRUE(retryErrors) && has_error) {
    return(FALSE)
  }

  # Ensure selective SD cluster settings are matched even if cached is missing the column
  if ("stim_sd_multiplier_clusters" %in% names(row)) {
    row_val <- row$stim_sd_multiplier_clusters[[1]]
    if (!is.null(row_val) && !is.na(row_val) && nzchar(as.character(row_val))) {
      cached_val <- if ("stim_sd_multiplier_clusters" %in% names(cached)) {
        cached$stim_sd_multiplier_clusters[[1]]
      } else if ("stimSdMultiplierClusters" %in% names(cached)) {
        cached$stimSdMultiplierClusters[[1]]
      } else {
        NA_character_
      }
      if (is.na(cached_val) || as.character(cached_val) != as.character(row_val)) {
        return(FALSE)
      }
    }
  }

  # Ensure selective mean-shift cluster settings are matched even if cached is missing the column
  if ("stim_mean_shift_clusters" %in% names(row)) {
    row_val <- row$stim_mean_shift_clusters[[1]]
    if (!is.null(row_val) && !is.na(row_val) && nzchar(as.character(row_val))) {
      cached_val <- if ("stim_mean_shift_clusters" %in% names(cached)) {
        cached$stim_mean_shift_clusters[[1]]
      } else if ("stimMeanShiftClusters" %in% names(cached)) {
        cached$stimMeanShiftClusters[[1]]
      } else {
        NA_character_
      }
      if (is.na(cached_val) || as.character(cached_val) != as.character(row_val)) {
        return(FALSE)
      }
    }
  }

  for (nm in names(row)) {
    if (nm %in% names(cached)) {
      row_val <- row[[nm]]
      cached_val <- cached[[nm]]

      if (is.null(row_val) || length(row_val) == 0L) {
        next
      }

      row_scalar <- row_val[[1]]
      cached_scalar <- cached_val[[1]]

      if (is.na(row_scalar) != is.na(cached_scalar)) {
        return(FALSE)
      }

      if (!is.na(row_scalar)) {
        if (is.numeric(row_scalar) && is.numeric(cached_scalar)) {
          if (
            !isTRUE(
              all.equal(
                as.numeric(row_scalar),
                as.numeric(cached_scalar),
                tolerance = 1e-7
              )
            )
          ) {
            return(FALSE)
          }
        } else if (nm %in% c("mismatch_type", "mismatchType")) {
          is_all_1 <- as.character(row_scalar) %in%
            c("mean_shift", "mean_shift_all")
          is_all_2 <- as.character(cached_scalar) %in%
            c("mean_shift", "mean_shift_all")
          if (is_all_1 != is_all_2) {
            return(FALSE)
          }
          if (
            !is_all_1 &&
              as.character(row_scalar) != as.character(cached_scalar)
          ) {
            return(FALSE)
          }
        } else if (
          nm %in% c(
            "stim_mean_shift_clusters",
            "stimMeanShiftClusters",
            "stim_sd_multiplier_clusters",
            "stimSdMultiplierClusters"
          )
        ) {
          v1 <- if (is.na(row_scalar) || !nzchar(as.character(row_scalar))) {
            ""
          } else {
            as.character(row_scalar)
          }
          v2 <- if (
            is.na(cached_scalar) || !nzchar(as.character(cached_scalar))
          ) {
            ""
          } else {
            as.character(cached_scalar)
          }
          if (v1 != v2) {
            return(FALSE)
          }
        } else {
          if (as.character(row_scalar) != as.character(cached_scalar)) {
            return(FALSE)
          }
        }
      }
    }
  }

  if (!has_error) {
    if (!is.null(nIter) && "iter" %in% names(cached)) {
      cached_iters <- unique(cached$iter[!is.na(cached$iter)])
      if (length(cached_iters) != as.integer(nIter)) {
        return(FALSE)
      }
    }
    if (!is.null(nSample) && "sample" %in% names(cached)) {
      cached_samples <- unique(cached$sample[!is.na(cached$sample)])
      if (length(cached_samples) != as.integer(nSample)) {
        return(FALSE)
      }
    }
  }

  TRUE
}

#' Format scenario settings string for logging
#'
#' @keywords internal
.simCompareFormatScenarioLog <- function(row, sim_id) {
  parts <- c(
    paste0("sim_id = ", sim_id),
    if ("transformation" %in% names(row)) {
      paste0("trans = ", row$transformation[[1]])
    },
    if ("prob_response" %in% names(row)) {
      paste0("prob = ", row$prob_response[[1]])
    },
    if ("n_cell" %in% names(row)) {
      paste0("n_cell = ", row$n_cell[[1]])
    },
    if ("mean_pos_setting" %in% names(row)) {
      paste0("mean_pos_setting = ", row$mean_pos_setting[[1]])
    },
    if ("mean_pos" %in% names(row)) {
      paste0("mean = ", row$mean_pos[[1]])
    },
    if ("bw_mtd" %in% names(row)) {
      paste0("bw_mtd = ", row$bw_mtd[[1]])
    },
    if ("bias_uns" %in% names(row)) {
      paste0("bias = ", row$bias_uns[[1]])
    },
    if ("mismatch_type" %in% names(row)) {
      paste0("mismatch_type = ", row$mismatch_type[[1]])
    },
    if ("mismatch_val" %in% names(row)) {
      paste0("mismatch_val = ", row$mismatch_val[[1]])
    },
    if ("stim_mean_shift" %in% names(row)) {
      paste0("stim_mean_shift = ", row$stim_mean_shift[[1]])
    },
    if (
      "stim_mean_shift_clusters" %in% names(row) &&
        !is.na(row$stim_mean_shift_clusters[[1]])
    ) {
      paste0("shift_clusters = ", row$stim_mean_shift_clusters[[1]])
    },
    if ("stim_sd_multiplier" %in% names(row)) {
      paste0("stim_sd_multiplier = ", row$stim_sd_multiplier[[1]])
    },
    if (
      "stim_sd_multiplier_clusters" %in% names(row) &&
        !is.na(row$stim_sd_multiplier_clusters[[1]])
    ) {
      paste0("sd_clusters = ", row$stim_sd_multiplier_clusters[[1]])
    }
  )
  paste(parts, collapse = " | ")
}

#' Append a log line safely
#'
#' @keywords internal
.simCompareLogMessage <- function(pathProgress, msg) {
  if (is.null(pathProgress) || !nzchar(pathProgress)) {
    return(invisible(NULL))
  }
  tryCatch(
    {
      dir.create(dirname(pathProgress), recursive = TRUE, showWarnings = FALSE)
      timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S] ")
      cat(paste0(timestamp, msg, "\n"), file = pathProgress, append = TRUE)
    },
    error = function(e) invisible(NULL)
  )
}

#' Ensure the current stimgate checkout is loaded
#'
#' @keywords internal
.simCompareEnsureCurrentCheckout <- function(pathRoot = NULL) {
  if (is.null(pathRoot) || !nzchar(pathRoot)) {
    if (requireNamespace("projr", quietly = TRUE)) {
      pathRoot <- tryCatch(
        projr::projr_path_get("project"),
        error = function(e) normalizePath(".", winslash = "/", mustWork = FALSE)
      )
    } else {
      pathRoot <- normalizePath(".", winslash = "/", mustWork = FALSE)
    }
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

  if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(pathRoot, quiet = TRUE))
  } else if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(pathRoot, quiet = TRUE))
  }
  invisible(TRUE)
}

#' Run a single comparison scenario row
#'
#' @keywords internal
.simCompareRunScenario <- function(
    row,
    nSample = 5,
    nIter = 5,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    covEvMin = 2,
    covEvMax = 2,
    tolClust = NULL,
    locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE,
    includeLocCondition = FALSE,
    includeLocDetails = includeLocCondition,
    dirCache = NULL,
    pathProgress = NULL,
    resume = TRUE,
    sim_grid_chunk_index = NULL,
    sim_grid_n_chunks = NULL,
    retryErrors = FALSE,
    p = NULL,
    ...) {
  .simCompareEnsureCurrentCheckout()

  sim_id <- if ("sim_id" %in% names(row)) {
    row$sim_id[[1]]
  } else {
    1L
  }

  file_output <- if (!is.null(dirCache) && nzchar(dirCache)) {
    .simCompareScenarioOutputPath(
      sim_id = sim_id,
      dirCache = dirCache,
      sim_grid_chunk_index = sim_grid_chunk_index,
      sim_grid_n_chunks = sim_grid_n_chunks
    )
  } else {
    character(0)
  }

  if (isTRUE(resume) && length(file_output) > 0L && file.exists(file_output)) {
    cached <- tryCatch(readRDS(file_output), error = function(e) NULL)
    if (
      .simCompareValidateScenarioCache(
        cached = cached,
        row = row,
        nSample = nSample,
        nIter = nIter,
        retryErrors = retryErrors
      )
    ) {
      if (!is.null(p)) {
        p(sprintf("Skipped existing sim_id: %s", sim_id))
      }
      .simCompareLogMessage(
        pathProgress = pathProgress,
        msg = paste0(
          "Skipped (cached): ",
          .simCompareFormatScenarioLog(row, sim_id)
        )
      )
      return(cached)
    }
  }

  settings_log <- .simCompareFormatScenarioLog(row, sim_id)
  .simCompareLogMessage(pathProgress, paste0("Running: ", settings_log))

  stimMeanShiftVal <- if ("stim_mean_shift" %in% names(row)) {
    row$stim_mean_shift[[1]]
  } else if ("stimMeanShift" %in% names(row)) {
    row$stimMeanShift[[1]]
  } else {
    0
  }

  stimSdMultVal <- if ("stim_sd_multiplier" %in% names(row)) {
    row$stim_sd_multiplier[[1]]
  } else if ("stimSdMultiplier" %in% names(row)) {
    row$stimSdMultiplier[[1]]
  } else if ("sd_multiplier" %in% names(row)) {
    row$sd_multiplier[[1]]
  } else if ("sd_increase" %in% names(row)) {
    1 + row$sd_increase[[1]]
  } else {
    1
  }

  stimMeanShiftClustersVal <- if ("stim_mean_shift_clusters" %in% names(row)) {
    val <- row$stim_mean_shift_clusters[[1]]
    if (is.null(val) || is.na(val) || !nzchar(as.character(val))) {
      NULL
    } else {
      as.character(val)
    }
  } else if ("stimMeanShiftClusters" %in% names(row)) {
    val <- row$stimMeanShiftClusters[[1]]
    if (is.null(val) || is.na(val) || !nzchar(as.character(val))) {
      NULL
    } else {
      as.character(val)
    }
  } else if (identical(row$mismatch_type[[1]], "mean_shift_negative")) {
    "gn"
  } else {
    NULL
  }

  stimSdMultiplierClustersVal <- if (
    "stim_sd_multiplier_clusters" %in% names(row)
  ) {
    val <- row$stim_sd_multiplier_clusters[[1]]
    if (is.null(val) || is.na(val) || !nzchar(as.character(val))) {
      NULL
    } else {
      as.character(val)
    }
  } else if ("stimSdMultiplierClusters" %in% names(row)) {
    val <- row$stimSdMultiplierClusters[[1]]
    if (is.null(val) || is.na(val) || !nzchar(as.character(val))) {
      NULL
    } else {
      as.character(val)
    }
  } else if (identical(row$mismatch_type[[1]], "sd_inflation_negative")) {
    "gn"
  } else {
    NULL
  }

  out <- tryCatch(
    {
      sim_res <- .simCompareFreqBs(
        nSample = nSample,
        nMarker = nMarker,
        nCondition = nCondition,
        nCluster = nCluster,
        nIter = nIter,
        biasUns = if ("bias_uns" %in% names(row)) row$bias_uns[[1]] else 0,
        bw = if ("bw" %in% names(row)) row$bw[[1]] else NULL,
        bwFallback = if ("bw_fallback" %in% names(row)) {
          row$bw_fallback[[1]]
        } else if ("bw" %in% names(row)) {
          row$bw[[1]]
        } else {
          "auto"
        },
        bwMin = if ("bw_min" %in% names(row)) row$bw_min[[1]] else "none",
        bwMax = if ("bw_max" %in% names(row)) row$bw_max[[1]] else "none",
        bwMtd = if ("bw_mtd" %in% names(row)) row$bw_mtd[[1]] else "hpi1",
        bwNcellMax = if ("bw_ncell_max" %in% names(row)) {
          row$bw_ncell_max[[1]]
        } else {
          1e4
        },
        probExact = probExact,
        nCellStim = if ("n_cell" %in% names(row)) row$n_cell[[1]] else 1e4,
        probResponse = if ("prob_response" %in% names(row)) {
          row$prob_response[[1]]
        } else {
          0.05
        },
        meanPos = if ("mean_pos" %in% names(row)) row$mean_pos[[1]] else 5,
        transformation = if ("transformation" %in% names(row)) {
          row$transformation[[1]]
        } else {
          "gaussian"
        },
        samplePerturbationSd = if ("sample_perturbation_sd" %in% names(row)) {
          row$sample_perturbation_sd[[1]]
        } else {
          0
        },
        conditionPerturbationSd = if (
          "condition_perturbation_sd" %in% names(row)
        ) {
          row$condition_perturbation_sd[[1]]
        } else {
          0
        },
        clusterPerturbationSd = if ("cluster_perturbation_sd" %in% names(row)) {
          row$cluster_perturbation_sd[[1]]
        } else {
          0
        },
        backgroundRelativeToResponse = if (
          "background_relative_to_response" %in% names(row)
        ) {
          row$background_relative_to_response[[1]]
        } else {
          0.1
        },
        ncellUnsRelativeToStim = if (
          "n_cell_uns_relative_to_stim" %in% names(row)
        ) {
          row$n_cell_uns_relative_to_stim[[1]]
        } else {
          1
        },
        covEvMin = covEvMin,
        covEvMax = covEvMax,
        tolClust = tolClust,
        locEnforceShapeThreshold = locEnforceShapeThreshold,
        calcCytPosGates = calcCytPosGates,
        includeLocCondition = includeLocCondition,
        includeLocDetails = includeLocDetails,
        stimMeanShift = stimMeanShiftVal,
        stimSdMultiplier = stimSdMultVal,
        stimMeanShiftClusters = stimMeanShiftClustersVal,
        stimSdMultiplierClusters = stimSdMultiplierClustersVal,
        ...
      )

      sim_res <- sim_res |>
        dplyr::select(-dplyr::any_of(names(row)))

      res <- dplyr::bind_cols(
        row[rep(1L, nrow(sim_res)), , drop = FALSE],
        sim_res
      )

      if (length(file_output) > 0L) {
        if (exists(".write_rds_atomic", mode = "function")) {
          .write_rds_atomic(res, file_output)
        } else {
          dir.create(
            dirname(file_output),
            recursive = TRUE,
            showWarnings = FALSE
          )
          saveRDS(res, file_output)
        }
      }

      .simCompareLogMessage(pathProgress, paste0("Completed: ", settings_log))
      if (!is.null(p)) {
        p(sprintf("Completed sim_id: %s", sim_id))
      }

      res
    },
    error = function(e) {
      .simCompareLogMessage(
        pathProgress,
        paste0("Error [", settings_log, "]: ", e$message)
      )
      if (!is.null(p)) {
        p(sprintf("ERROR on sim_id: %s", sim_id))
      }

      err_res <- dplyr::bind_cols(
        row,
        tibble::tibble(
          iter = NA_integer_,
          chnl = "F1",
          sample = NA_character_,
          ind = NA_character_,
          approach = NA_character_,
          method = NA_character_,
          propRespTruth = NA_real_,
          propRespEst = NA_real_,
          threshold = NA_real_,
          thresholdOrigin = NA_character_,
          gateReturnPoint = NA_character_,
          thresholdMetric = NA_real_,
          thresholdFallbackUsed = NA,
          nCellStim = NA_real_,
          nCellUns = NA_real_,
          nPosStim = NA_integer_,
          nPosUns = NA_integer_,
          propStim = NA_real_,
          propUns = NA_real_,
          propStimTruth = NA_real_,
          propUnsTruth = NA_real_,
          detailLevel = NA_character_,
          locGenerated = NA,
          locGeneratedDirect = NA,
          locSource = NA_character_,
          locReason = NA_character_,
          error = e$message
        )
      )

      if (length(file_output) > 0L) {
        if (exists(".write_rds_atomic", mode = "function")) {
          .write_rds_atomic(err_res, file_output)
        } else {
          dir.create(
            dirname(file_output),
            recursive = TRUE,
            showWarnings = FALSE
          )
          saveRDS(err_res, file_output)
        }
      }

      err_res
    }
  )

  out
}

#' Run the comparison helper over a simulation grid
#'
#' @keywords internal
.simCompareFreqBsGrid <- function(
    sim_grid,
    nSample = 5,
    nIter = 5,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    probExact = TRUE,
    covEvMin = 2,
    covEvMax = 2,
    tolClust = NULL,
    locEnforceShapeThreshold = FALSE,
    calcCytPosGates = FALSE,
    includeLocCondition = FALSE,
    includeLocDetails = includeLocCondition,
    parallel = FALSE,
    workers = NULL,
    dirCache = NULL,
    pathProgress = NULL,
    resume = TRUE,
    progress = TRUE,
    sim_grid_chunk_index = NULL,
    sim_grid_n_chunks = NULL,
    retryErrors = FALSE,
    ...) {
  if (nrow(sim_grid) == 0L) {
    return(tibble::tibble())
  }

  if (!"sim_id" %in% names(sim_grid)) {
    sim_grid <- sim_grid |>
      dplyr::mutate(sim_id = dplyr::row_number())
  }

  if (!is.null(dirCache) && nzchar(dirCache)) {
    dir.create(dirCache, recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(pathProgress) && nzchar(pathProgress)) {
    dir.create(dirname(pathProgress), recursive = TRUE, showWarnings = FALSE)
  }

  run_one <- function(i, p = NULL) {
    row <- sim_grid[i, , drop = FALSE]
    .simCompareRunScenario(
      row = row,
      nSample = nSample,
      nIter = nIter,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      probExact = probExact,
      covEvMin = covEvMin,
      covEvMax = covEvMax,
      tolClust = tolClust,
      locEnforceShapeThreshold = locEnforceShapeThreshold,
      calcCytPosGates = calcCytPosGates,
      includeLocCondition = includeLocCondition,
      includeLocDetails = includeLocDetails,
      dirCache = dirCache,
      pathProgress = pathProgress,
      resume = resume,
      sim_grid_chunk_index = sim_grid_chunk_index,
      sim_grid_n_chunks = sim_grid_n_chunks,
      retryErrors = retryErrors,
      p = p,
      ...
    )
  }

  if (
    !requireNamespace("furrr", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE)
  ) {
    out_list <- lapply(seq_len(nrow(sim_grid)), function(i) run_one(i, p = NULL))
    res_tbl <- purrr::list_rbind(out_list)
    if ("sim_id" %in% names(res_tbl)) {
      res_tbl <- res_tbl |>
        dplyr::arrange(
          .data$sim_id,
          dplyr::across(dplyr::any_of(c("iter", "sample", "ind", "method")))
        )
    }
    return(res_tbl)
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (isTRUE(parallel)) {
    workers_use <- if (!is.null(workers)) {
      max(1L, as.integer(workers))
    } else {
      max(1L, .simGetCores())
    }
    if (workers_use > 1L) {
      future::plan(future::multisession, workers = workers_use)
    } else {
      future::plan(future::sequential)
    }
  } else {
    future::plan(future::sequential)
  }

  furrr_opts <- furrr::furrr_options(seed = TRUE)

  out_list <- if (
    isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)
  ) {
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(sim_grid))
      furrr::future_map(
        seq_len(nrow(sim_grid)),
        function(i) run_one(i, p = p),
        .options = furrr_opts
      )
    })
  } else {
    furrr::future_map(
      seq_len(nrow(sim_grid)),
      function(i) run_one(i, p = NULL),
      .options = furrr_opts
    )
  }

  res_tbl <- purrr::list_rbind(out_list)

  if ("sim_id" %in% names(res_tbl)) {
    res_tbl <- res_tbl |>
      dplyr::arrange(
        .data$sim_id,
        dplyr::across(dplyr::any_of(c("iter", "sample", "ind", "method")))
      )
  }

  res_tbl
}

#' Collate scenario output files
#'
#' @keywords internal
.simCompareCollateScenarioOutputs <- function(
    dirCache = NULL,
    pathList = NULL,
    sim_grid = NULL) {
  if (is.null(pathList) || length(pathList) == 0L) {
    if (is.null(dirCache) || !dir.exists(dirCache)) {
      return(tibble::tibble())
    }
    pathList <- .simCompareFindScenarioOutputs(dirCache)
  }

  if (length(pathList) == 0L) {
    return(tibble::tibble())
  }

  out_list <- purrr::map(
    pathList,
    function(path) {
      tryCatch(
        readRDS(path),
        error = function(e) {
          warning(
            "Could not read scenario output file: ",
            path,
            " (",
            e$message,
            ")"
          )
          NULL
        }
      )
    }
  ) |>
    purrr::compact()

  if (length(out_list) == 0L) {
    return(tibble::tibble())
  }

  collated <- purrr::list_rbind(out_list)

  if (
    !is.null(sim_grid) &&
      "sim_id" %in% names(sim_grid) &&
      "sim_id" %in% names(collated)
  ) {
    collated <- collated |>
      dplyr::filter(.data$sim_id %in% sim_grid$sim_id)
  }

  if ("sim_id" %in% names(collated)) {
    collated <- collated |>
      dplyr::arrange(
        .data$sim_id,
        dplyr::across(dplyr::any_of(c("iter", "sample", "ind", "method")))
      )
  }

  collated
}

#' Summarise comparison runs by scenario and method
#'
#' @keywords internal
.simCompareSummariseFreqBs <- function(
    .data,
    scenarioCols = NULL,
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  if (!"error" %in% names(.data)) {
    .data$error <- NA_character_
  }

  if (is.null(scenarioCols)) {
    scenarioCols <- intersect(
      c(
        "mismatch_type",
        "mismatch_val",
        "stim_mean_shift",
        "stimMeanShift",
        "stim_mean_shift_clusters",
        "stimMeanShiftClusters",
        "stim_sd_multiplier",
        "stimSdMultiplier",
        "stim_sd_multiplier_clusters",
        "stimSdMultiplierClusters",
        "sd_increase",
        "transformation",
        "prob_response",
        "n_cell",
        "mean_pos",
        "mean_pos_setting",
        "scenario_desc",
        "bw",
        "bias_uns",
        "sample_perturbation_sd",
        "condition_perturbation_sd",
        "cluster_perturbation_sd",
        "background_relative_to_response",
        "n_cell_uns_relative_to_stim",
        "approach",
        "method"
      ),
      names(.data)
    )
  }

  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::mutate(
      run_error = .data$error,
      freq_error = .data$propRespEst - .data$propRespTruth,
      abs_error = abs(.data$freq_error),
      sq_error = .data$freq_error^2,
      rel_error = dplyr::if_else(
        .data$propRespTruth != 0,
        .data$freq_error / .data$propRespTruth,
        NA_real_
      )
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_est = sum(is.finite(.data$propRespEst)),
      n_run_error = sum(!is.na(.data$run_error) & .data$run_error != ""),
      n_threshold_fallback = sum(
        .data$thresholdFallbackUsed %in% TRUE,
        na.rm = TRUE
      ),
      propRespTruth_mean = mean(.data$propRespTruth, na.rm = TRUE),
      propRespEst_mean = mean(.data$propRespEst, na.rm = TRUE),
      propStim_mean = mean(.data$propStim, na.rm = TRUE),
      propStim_median = stats::median(.data$propStim, na.rm = TRUE),
      propUns_mean = mean(.data$propUns, na.rm = TRUE),
      propUns_median = stats::median(.data$propUns, na.rm = TRUE),
      bias = mean(.data$freq_error, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      median_abs_error = stats::median(.data$abs_error, na.rm = TRUE),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      med_abs_rel_error = stats::median(abs(.data$rel_error), na.rm = TRUE),
      q90_abs_rel_error = stats::quantile(
        abs(.data$rel_error),
        probs = 0.9,
        na.rm = TRUE
      ),
      q95_abs_rel_error = stats::quantile(
        abs(.data$rel_error),
        probs = 0.95,
        na.rm = TRUE
      ),
      max_abs_rel_error = max(abs(.data$rel_error), na.rm = TRUE),
      threshold_mean = mean(.data$threshold, na.rm = TRUE),
      threshold_median = stats::median(.data$threshold, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::across(
        dplyr::any_of(
          c(
            "n_cell",
            "prob_response",
            "transformation",
            "mean_pos",
            "mismatch_type",
            "mismatch_val",
            "method"
          )
        )
      )
    )
}
