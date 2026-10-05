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

  switch(
    transformation,
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

#' Locate the fbeta Python script
#'
#' @keywords internal
.simCompareFbetaPath <- function(pathFbeta = NULL) {
  if (!is.null(pathFbeta)) {
    return(pathFbeta)
  }
  if (requireNamespace("projr", quietly = TRUE)) {
    pathCandidate <- projr::projr_path_get("project", "scripts", "python", "fbeta.py")
    if (file.exists(pathCandidate)) {
      return(pathCandidate)
    }
  }
  # If projr did not locate the script, look above the working directory
  # (e.g. analysis tests run from analysis/tests/testthat).
  pathDir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  repeat {
    pathCandidate <- file.path(pathDir, "scripts", "python", "fbeta.py")
    if (file.exists(pathCandidate)) {
      return(pathCandidate)
    }
    pathParent <- dirname(pathDir)
    if (identical(pathParent, pathDir)) {
      break
    }
    pathDir <- pathParent
  }
  stop(
    "pathFbeta was not supplied and ",
    "scripts/python/fbeta.py was not found above the working directory. ",
    "Pass pathFbeta explicitly."
  )
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

#' Load the fbeta Python implementation for one R run
#'
#' The returned environment contains Reticulate-backed Python objects and must
#' remain local to the R process that created it. In particular, it must not be
#' stored in a global cache that can be serialised to multisession workers.
#'
#' @keywords internal
.simCompareFbetaEnvironment <- function(
  pathFbeta = NULL,
  patchPy2Compat = TRUE
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to call fbeta.py.")
  }

  pathFbeta <- .simCompareFbetaPath(pathFbeta)
  pathPyUse <- .simCompareFbetaCompatPath(
    pathFbeta = pathFbeta,
    patchPy2Compat = patchPy2Compat
  )

  pyEnv <- new.env(parent = emptyenv())
  reticulate::source_python(pathPyUse, envir = pyEnv)

  if (!exists("get_positivity_threshold", envir = pyEnv, inherits = FALSE)) {
    stop("fbeta.py did not define get_positivity_threshold().")
  }

  pyEnv
}

#' Call the fbeta Python implementation via reticulate
#'
#' @keywords internal
.simCompareFbetaThreshold <- function(
  xUns,
  xStim,
  pathFbeta = NULL,
  patchPy2Compat = TRUE,
  fbetaEnv = NULL,
  beta = 0.8,
  theta = 2,
  width = 10,
  numBins = NULL
) {
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

  if (is.null(fbetaEnv)) {
    fbetaEnv <- .simCompareFbetaEnvironment(
      pathFbeta = pathFbeta,
      patchPy2Compat = patchPy2Compat
    )
  }

  negMat <- matrix(xUns, ncol = 1L)
  posMat <- matrix(xStim, ncol = 1L)

  out <- fbetaEnv$get_positivity_threshold(
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
  autoTol = FALSE
) {
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

  bandwidthUse <- if (is.null(bandwidth)) {
    suppressWarnings(
      tryCatch(
        ks::hpi(x, deriv.order = 1L),
        error = function(e) NA_real_
      )
    )
  } else {
    bandwidth
  }

  if (
    length(bandwidthUse) != 1L ||
      !is.finite(bandwidthUse) ||
      bandwidthUse <= 0
  ) {
    return(list(
      threshold = NA_real_,
      thresholdMetric = NA_real_,
      thresholdOrigin = "failed_bandwidth_nonfinite"
    ))
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
    bandwidth = bandwidthUse
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
  stimSdMultiplierClusters = NULL
) {
  if (
    is.null(outListExperiment) || is.null(outListExperiment[["flowFrameList"]])
  ) {
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
  scenario = NULL
) {
  simcytoArgs <- names(formals(simcyto::simCytExperiment))
  clusterShiftSupported <- !is.null(simcytoArgs) &&
    "stimMeanShiftClusters" %in% simcytoArgs
  clusterSdSupported <- !is.null(simcytoArgs) &&
    "stimSdMultiplierClusters" %in% simcytoArgs

  selectiveShift <- !is.null(stimMeanShiftClusters)
  selectiveSd <- !is.null(stimSdMultiplierClusters)

  # If the installed simcyto supports selective mismatch, let simcyto apply it
  # exactly once. For older simcyto versions, neutralise the global mismatch
  # during simulation and apply the selective mismatch locally afterwards.
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
    stimMeanShift = if (selectiveShift && !clusterShiftSupported) {
      0
    } else {
      stimMeanShift
    },
    stimSdMultiplier = if (selectiveSd && !clusterSdSupported) {
      1
    } else {
      stimSdMultiplier
    },
    scenario = scenario
  )

  if (clusterShiftSupported && selectiveShift) {
    callArgs$stimMeanShiftClusters <- stimMeanShiftClusters
  }
  if (clusterSdSupported && selectiveSd) {
    callArgs$stimSdMultiplierClusters <- stimSdMultiplierClusters
  }

  result <- do.call(simcyto::simCytExperiment, callArgs)

  if (
    (selectiveShift && !clusterShiftSupported) ||
      (selectiveSd && !clusterSdSupported)
  ) {
    result <- .simCompareApplyClusterMismatch(
      outListExperiment = result,
      stimMeanShift = if (selectiveShift && !clusterShiftSupported) {
        stimMeanShift
      } else {
        0
      },
      stimSdMultiplier = if (selectiveSd && !clusterSdSupported) {
        stimSdMultiplier
      } else {
        1
      },
      stimMeanShiftClusters = if (selectiveShift && !clusterShiftSupported) {
        stimMeanShiftClusters
      } else {
        NULL
      },
      stimSdMultiplierClusters = if (selectiveSd && !clusterSdSupported) {
        stimSdMultiplierClusters
      } else {
        NULL
      }
    )
  }

  result
}

#' @keywords internal
.simCompareEstimateFromThreshold <- function(
  xStim,
  xUns,
  threshold,
  fallbackHighValue = TRUE,
  fallbackMargin = 0.05,
  labelsStim = NULL
) {
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

  c(
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
    ),
    if (!is.null(labelsStim)) {
      .simCompareConfusionCounts(xStim, labelsStim, thresholdUsed)
    }
  )
}

# Confusion-matrix counts of one stimulated tube against its simulation labels.
# A cell is classified positive when its expression is strictly above the gate
# (`x > gate`), as in the package (`R/pos_ind.R` and the statistics helpers),
# so a cell exactly at the gate is negative. Cells labelled `positiveLabel`
# ("gp") are genuine positives, whether induced by stimulation or part of the
# positive background; all other cells ("gn") are genuine negatives. Counts
# are NA when the gate is not finite.
.simCompareConfusionCounts <- function(
    x,
    labels,
    threshold,
    positiveLabel = "gp") {
  x <- as.numeric(x)
  labels <- as.character(labels)
  if (length(x) != length(labels)) {
    stop("Expression and label vectors must have the same length.")
  }
  threshold <- as.numeric(threshold)[1]
  if (length(threshold) == 0L || !is.finite(threshold)) {
    return(list(
      nTruePos = NA_integer_, nFalsePos = NA_integer_,
      nFalseNeg = NA_integer_, nTrueNeg = NA_integer_
    ))
  }
  pos <- x > threshold
  truth <- labels %in% positiveLabel
  list(
    nTruePos = sum(pos & truth),
    nFalsePos = sum(pos & !truth),
    nFalseNeg = sum(!pos & truth),
    nTrueNeg = sum(!pos & !truth)
  )
}

# Replicate-level classification outcomes from the confusion-matrix counts.
# Proportions with a zero denominator are NA, not zero: FDP is undefined when
# the gate selects no stimulated cells, sensitivity when the tube has no
# genuine positives and the false-positive rate when it has no genuine
# negatives. `gate_status` separates failed runs, fallback gates and
# calculated gates, each split by whether any stimulated cell was selected.
.simCompareClassificationMetrics <- function(.data) {
  ratio <- function(num, den) {
    dplyr::if_else(!is.na(den) & den > 0, num / den, NA_real_)
  }
  has_error <- if ("error" %in% names(.data)) {
    !is.na(.data$error) & nzchar(as.character(.data$error))
  } else {
    rep(FALSE, nrow(.data))
  }
  fallback <- if ("thresholdFallbackUsed" %in% names(.data)) {
    .data$thresholdFallbackUsed %in% TRUE
  } else {
    rep(FALSE, nrow(.data))
  }
  .data |>
    dplyr::mutate(
      n_selected = .data$nTruePos + .data$nFalsePos,
      n_genuine_pos = .data$nTruePos + .data$nFalseNeg,
      n_genuine_neg = .data$nFalsePos + .data$nTrueNeg,
      n_classified = .data$n_selected + .data$nFalseNeg + .data$nTrueNeg,
      fdp = ratio(.data$nFalsePos, .data$n_selected),
      sensitivity = ratio(.data$nTruePos, .data$n_genuine_pos),
      false_positive_rate = ratio(.data$nFalsePos, .data$n_genuine_neg),
      selected_fraction = ratio(.data$n_selected, .data$n_classified),
      gate_empty = .data$n_selected == 0L,
      gate_status = dplyr::case_when(
        has_error | is.na(.data$n_selected) ~ "failed",
        fallback & .data$gate_empty ~ "fallback_empty",
        fallback ~ "fallback_selected",
        .data$gate_empty ~ "calculated_empty",
        TRUE ~ "calculated_selected"
      )
    )
}

#' @keywords internal
.simCompareTruthTable <- function(
  labelsList,
  nSample,
  nCondition,
  chnl = "F1"
) {
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
  fallbackMargin = 0.05
) {
  tailgateX <- match.arg(tailgateX)
  tailgateMethod <- match.arg(tailgateMethod)

  fbetaEnv <- .simCompareFbetaEnvironment(
    pathFbeta = pathFbeta,
    patchPy2Compat = fbetaPatchPy2Compat
  )

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
      labelsStim <- labelsList[[indStim]]

      fbetaError <- NA_character_
      fbetaObj <- tryCatch(
        .simCompareFbetaThreshold(
          xUns = xUnsFbeta,
          xStim = xStim,
          pathFbeta = pathFbeta,
          patchPy2Compat = fbetaPatchPy2Compat,
          fbetaEnv = fbetaEnv,
          beta = fbetaBeta,
          theta = fbetaTheta,
          width = fbetaWidth,
          numBins = fbetaNumBins
        ),
        error = function(e) {
          fbetaError <<- conditionMessage(e)
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", fbetaError)
          )
        }
      )

      fbetaEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsFbeta,
        threshold = fbetaObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin,
        labelsStim = labelsStim
      )

      xTail <- switch(
        tailgateX,
        "stim" = xStim,
        "unstim" = xUnsTailgate,
        "combined" = c(xUnsTailgate, xStim)
      )

      tailgateError <- NA_character_
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
          tailgateError <<- conditionMessage(e)
          list(
            threshold = NA_real_,
            thresholdMetric = NA_real_,
            thresholdOrigin = paste0("error: ", tailgateError)
          )
        }
      )

      tailgateEst <- .simCompareEstimateFromThreshold(
        xStim = xStim,
        xUns = xUnsTailgate,
        threshold = tailgateObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin,
        labelsStim = labelsStim
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
          if (!is.na(fbetaError)) {
            "fbeta_error_fallback_high_value"
          } else if (isTRUE(fbetaEst$thresholdFallbackUsed)) {
            "fbeta_fallback_high_value"
          } else {
            "fbeta_calculated"
          },
          if (!is.na(tailgateError)) {
            "tailgate_error_fallback_high_value"
          } else if (isTRUE(tailgateEst$thresholdFallbackUsed)) {
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
        nTruePos = c(fbetaEst$nTruePos, tailgateEst$nTruePos),
        nFalsePos = c(fbetaEst$nFalsePos, tailgateEst$nFalsePos),
        nFalseNeg = c(fbetaEst$nFalseNeg, tailgateEst$nFalseNeg),
        nTrueNeg = c(fbetaEst$nTrueNeg, tailgateEst$nTrueNeg),
        detailLevel = NA_character_,
        locGenerated = NA,
        locGeneratedDirect = NA,
        locSource = NA_character_,
        locReason = NA_character_,
        error = c(fbetaError, tailgateError)
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
      nTruePos = NA_integer_,
      nFalsePos = NA_integer_,
      nFalseNeg = NA_integer_,
      nTrueNeg = NA_integer_,
      detailLevel = NA_character_,
      locGenerated = NA,
      locGeneratedDirect = NA,
      locSource = NA_character_,
      locReason = NA_character_,
      error = errorMessage
    )
}

#' Extract final StimGate threshold provenance for comparison outputs
#'
#' @keywords internal
.simCompareStimgateGateProvenance <- function(gRow, gateVal, isClustered) {
  has_row <- is.data.frame(gRow) && nrow(gRow) > 0L

  locGenerated <- if (
    has_row &&
      "locGenerated" %in% names(gRow) &&
      !is.na(gRow$locGenerated[[1]])
  ) {
    isTRUE(gRow$locGenerated[[1]])
  } else {
    is.finite(gateVal)
  }

  locGeneratedDirect <- if (
    has_row &&
      "locGeneratedDirect" %in% names(gRow) &&
      !is.na(gRow$locGeneratedDirect[[1]])
  ) {
    isTRUE(gRow$locGeneratedDirect[[1]])
  } else {
    isTRUE(locGenerated) && !isTRUE(isClustered)
  }

  locSource <- if (
    has_row &&
      "locSource" %in% names(gRow) &&
      !is.na(gRow$locSource[[1]])
  ) {
    as.character(gRow$locSource[[1]])
  } else if (isTRUE(isClustered)) {
    "cluster"
  } else if (isTRUE(locGenerated)) {
    "sample"
  } else {
    "not_calculated"
  }

  locReason <- if (
    has_row &&
      "locReason" %in% names(gRow) &&
      !is.na(gRow$locReason[[1]])
  ) {
    as.character(gRow$locReason[[1]])
  } else {
    NA_character_
  }

  thresholdFallbackUsed <- !isTRUE(locGenerated)

  list(
    thresholdOrigin = if (thresholdFallbackUsed) {
      "fallback_high_value"
    } else if (isTRUE(isClustered)) {
      "calculated_clustered"
    } else {
      "calculated"
    },
    gateReturnPoint = if (thresholdFallbackUsed) {
      "stimgate_fallback_high_value"
    } else if (isTRUE(isClustered)) {
      "stimgate_clustered"
    } else {
      "stimgate_calculated"
    },
    thresholdFallbackUsed = thresholdFallbackUsed,
    locGenerated = locGenerated,
    locGeneratedDirect = locGeneratedDirect,
    locSource = locSource,
    locReason = locReason
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
  biasUnsFactor = 1,
  bwFallback = bw,
  bwMin = "none",
  bwMax = "none",
  bwMtd = "hpi1",
  bwAdj = 1,
  bwNcellMin = 1e2,
  bwNcellMax = 1e5,
  bwCluster = NULL,
  minCell = 1e2,
  # Retained for callers/manifests; gateStim no longer uses these settings.
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
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
  clusterGates = FALSE,
  locEnforceShapeThreshold = FALSE,
  calcCytPosGates = FALSE,
  includeLocCondition = FALSE,
  includeLocDetails = includeLocCondition
) {
  if (!is.logical(clusterGates) || length(clusterGates) != 1L || is.na(clusterGates)) {
    stop("clusterGates must be a single TRUE or FALSE.")
  }

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
        biasUns = biasUns,
        bw = bw,
        control = stimgate::stimControl(
          biasUnsFactor = biasUnsFactor,
          bwFallback = bwFallback,
          bwMin = bwMin,
          bwMax = bwMax,
          bwMtd = bwMtd,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax,
          bwCluster = bwCluster,
          clusterGates = clusterGates,
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
          calcCytPosGates = calcCytPosGates,
          minCell = minCell
        )
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

      stimgatePrimaryTbl <- purrr::map_df(
        seq_len(nSample),
        function(sampleCurr) {
          indUns <- (sampleCurr - 1L) * nCondition + 1L
          indStimVec <- seq.int(indUns + 1L, sampleCurr * nCondition)

          purrr::map_df(indStimVec, function(indStim) {
            ind_curr <- as.character(indStim)
            gRow <- if (
              nrow(gateTblFinal) > 0L && "ind" %in% names(gateTblFinal)
            ) {
              gateTblFinal[
                as.character(gateTblFinal$ind) == ind_curr &
                  gateTblFinal$chnl == "F1",
                ,
                drop = FALSE
              ]
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

            sRow <- if (
              nrow(statsTblFinal) > 0L && "ind" %in% names(statsTblFinal)
            ) {
              statsTblFinal[
                as.character(statsTblFinal$ind) == ind_curr &
                  grepl("~\\+~", statsTblFinal$cytCombn),
                ,
                drop = FALSE
              ]
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

            gateVal <- if (nrow(gRow) > 0L) {
              suppressWarnings(as.numeric(unname(gRow$gate[[1]])))
            } else {
              NA_real_
            }
            gateNm <- if (nrow(gRow) > 0L) {
              as.character(gRow$gateName[[1]])
            } else {
              NA_character_
            }

            nCellStimVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.numeric(sRow$nCellStim[[1]]))
            } else {
              NA_real_
            }
            nCellUnsVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.numeric(sRow$nCellUns[[1]]))
            } else {
              NA_real_
            }
            nPosStimVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.integer(sRow$countStim[[1]]))
            } else {
              NA_integer_
            }
            nPosUnsVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.integer(sRow$countUns[[1]]))
            } else {
              NA_integer_
            }
            propStimVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.numeric(sRow$propStim[[1]]))
            } else {
              NA_real_
            }
            propUnsVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.numeric(sRow$propUns[[1]]))
            } else {
              NA_real_
            }
            propBsVal <- if (nrow(sRow) > 0L) {
              suppressWarnings(as.numeric(sRow$propBs[[1]]))
            } else {
              NA_real_
            }

            # Classify the stimulated cells with the final (cluster-refined
            # where applicable) gate, exactly as StimGate's statistics do.
            # Read the GatingSet copy that StimGate gated: it stores
            # expression in single precision, so a cell next to the gate can
            # fall on the other side of it in the double-precision simulation.
            xStim <- flowCore::exprs(
              flowWorkspace::gh_pop_get_data(gs[[indStim]], "root")
            )[, "F1"]
            counts <- .simCompareConfusionCounts(
              xStim, labelsList[[indStim]], gateVal
            )

            isClustered <- grepl("Clust$", gateNm %||% "")
            provenance <- .simCompareStimgateGateProvenance(
              gRow = gRow,
              gateVal = gateVal,
              isClustered = isClustered
            )

            tibble::tibble(
              sample = as.character(sampleCurr),
              ind = as.character(indStim),
              chnl = "F1",
              approach = "stimgate",
              method = "stimgate",
              threshold = gateVal,
              thresholdOrigin = provenance$thresholdOrigin,
              gateReturnPoint = provenance$gateReturnPoint,
              thresholdMetric = NA_real_,
              thresholdFallbackUsed = provenance$thresholdFallbackUsed,
              nCellStim = nCellStimVal,
              nCellUns = nCellUnsVal,
              nPosStim = nPosStimVal,
              nPosUns = nPosUnsVal,
              propStim = propStimVal,
              propUns = propUnsVal,
              propRespEst = propBsVal,
              nTruePos = counts$nTruePos,
              nFalsePos = counts$nFalsePos,
              nFalseNeg = counts$nFalseNeg,
              nTrueNeg = counts$nTrueNeg,
              detailLevel = if (isClustered) {
                "cluster_final"
              } else {
                "sample_final"
              },
              locGenerated = provenance$locGenerated,
              locGeneratedDirect = provenance$locGeneratedDirect,
              locSource = provenance$locSource,
              locReason = provenance$locReason,
              error = NA_character_
            )
          })
        }
      )

      detailTbl <- if (
        isTRUE(includeLocDetails) || isTRUE(includeLocCondition)
      ) {
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
          dplyr::any_of(c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg")),
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
  biasUnsFactor = 1,
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
  clusterGates = FALSE,
  locEnforceShapeThreshold = FALSE,
  minCell = 1e2,
  # Retained for callers/manifests; gateStim no longer uses these settings.
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
  locMarginalPurityRel = 0.5,
  locMarginalCellBinRatio = 2,
  locMarginalRefQuantile = 0.75,
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
  pathProject = NULL,
  keepCells = FALSE
) {
  if (!is.logical(clusterGates) || length(clusterGates) != 1L || is.na(clusterGates)) {
    stop("clusterGates must be a single TRUE or FALSE.")
  }

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

  # One seed per iteration, drawn before any method runs. The methods consume
  # different amounts of randomness in different mismatch settings, so without
  # this the data of later iterations would not be paired across settings.
  # Drawing with replacement makes the first seeds independent of `nIter`.
  iterSeeds <- sample.int(.Machine$integer.max, nIter, replace = TRUE)
  cellsList <- list()

  out <- purrr::map_df(seq_len(nIter), function(iterNum) {
    set.seed(iterSeeds[[iterNum]])
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
    # Fingerprint of each sample's unstimulated tube, which no mismatch
    # changes: equal values across mismatch settings show the data are paired.
    unsTbl <- tibble::tibble(
      sample = as.character(seq_len(nSample)),
      unsExprSum = vapply(seq_len(nSample), function(sampleCurr) {
        sum(flowCore::exprs(
          flowFrameList[[(sampleCurr - 1L) * nCondition + 1L]]
        )[, "F1"])
      }, numeric(1))
    )
    if (isTRUE(keepCells)) {
      cellsList[[iterNum]] <<- .simCompareCellTable(
        flowFrameList, labelsList, nSample, nCondition, iterNum
      )
    }
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
      biasUnsFactor = biasUnsFactor,
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
      locMarginalPurityRel = locMarginalPurityRel,
      locMarginalCellBinRatio = locMarginalCellBinRatio,
      locMarginalRefQuantile = locMarginalRefQuantile,
      clusterGates = clusterGates,
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
      dplyr::left_join(unsTbl, by = "sample") |>
      dplyr::mutate(
        iter = iterNum,
        nCellStimSim = nCellStim,
        nCellUnsSim = nCellUns,
        # NULL biasUns: StimGate sets it from the initial bandwidth estimate.
        biasUns = biasUns %||% NA_real_,
        biasUnsFactor = biasUnsFactor,
        bw = bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwCluster = bwCluster %||% NA_real_,
        clusterGates = clusterGates,
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

  if (isTRUE(keepCells)) {
    attr(out, "cells") <- dplyr::bind_rows(cellsList)
  }
  out
}

# Cell-level expression and simulation labels of every tube, for diagnostics
# that need to show the distributions behind a gate.
.simCompareCellTable <- function(
    flowFrameList,
    labelsList,
    nSample,
    nCondition,
    iter = 1L,
    chnl = "F1") {
  purrr::map_df(seq_len(nSample), function(sampleCurr) {
    indUns <- (sampleCurr - 1L) * nCondition + 1L
    purrr::map_df(seq.int(indUns, sampleCurr * nCondition), function(ind) {
      # Evaluate before `tibble()`, whose `ind` column would mask the index.
      expr <- as.numeric(flowCore::exprs(flowFrameList[[ind]])[, chnl])
      label <- as.character(labelsList[[ind]])
      tibble::tibble(
        iter = as.integer(iter),
        sample = as.character(sampleCurr),
        ind = as.character(ind),
        condition = if (ind == indUns) "unstim" else "stim",
        expr = expr,
        label = label
      )
    })
  })
}

#' Generate scenario output file path
#'
#' @keywords internal
.simCompareScenarioOutputPath <- function(
  sim_id,
  dirCache,
  sim_grid_chunk_index = NULL,
  sim_grid_n_chunks = NULL
) {
  if (is.null(dirCache) || !nzchar(dirCache)) {
    return(character(0))
  }

  if (
    !is.null(sim_grid_chunk_index) &&
      !is.na(sim_grid_chunk_index) &&
      !is.null(sim_grid_n_chunks) &&
      !is.na(sim_grid_n_chunks) &&
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

#' Check that a scenario has one complete primary result per replicate and method
#'
#' @keywords internal
.simComparePrimaryOutputComplete <- function(
  .data,
  nSample,
  nIter,
  methods = c("stimgate", "fbeta", "tailgate")
) {
  # The label-based confusion-matrix counts (and the unstimulated-data
  # fingerprint used to check pairing) are required, so outputs saved before
  # they were recorded cannot satisfy the classification analysis.
  required_cols <- c(
    "iter",
    "sample",
    "method",
    "propRespTruth",
    "propRespEst",
    "nCellStim",
    "nPosStim",
    .simCompareCountCols,
    "unsExprSum"
  )
  if (
    !is.data.frame(.data) ||
      nrow(.data) == 0L ||
      !all(required_cols %in% names(.data))
  ) {
    return(FALSE)
  }

  if (
    "error" %in% names(.data) &&
      any(!is.na(.data$error) & nzchar(as.character(.data$error)))
  ) {
    return(FALSE)
  }

  primary <- .data |>
    dplyr::filter(.data$method %in% methods)

  expected_n_per_method <- as.integer(nSample) * as.integer(nIter)
  if (nrow(primary) != expected_n_per_method * length(methods)) {
    return(FALSE)
  }

  key_counts <- primary |>
    dplyr::count(.data$iter, .data$sample, .data$method, name = "n")

  if (any(key_counts$n != 1L)) {
    return(FALSE)
  }

  method_counts <- primary |>
    dplyr::count(.data$method, name = "n")

  if (
    !setequal(as.character(method_counts$method), methods) ||
      any(method_counts$n != expected_n_per_method)
  ) {
    return(FALSE)
  }

  all(is.finite(primary$propRespTruth)) &&
    all(is.finite(primary$propRespEst)) &&
    .simCompareCountsConsistent(primary) &&
    all(is.finite(primary$unsExprSum))
}

.simCompareCountCols <- c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg")

# TRUE when every row has complete confusion-matrix counts that reproduce the
# method's own gated stimulated count (`nPosStim`) and stimulated cell count.
.simCompareCountsConsistent <- function(.data) {
  if (!all(c(.simCompareCountCols, "nPosStim", "nCellStim") %in% names(.data))) {
    return(FALSE)
  }
  counts <- as.matrix(.data[, .simCompareCountCols])
  if (anyNA(counts) || any(counts < 0)) {
    return(FALSE)
  }
  n_selected <- .data$nTruePos + .data$nFalsePos
  isTRUE(all(n_selected == .data$nPosStim)) &&
    isTRUE(all(rowSums(counts) == .data$nCellStim))
}

#' Validate scenario cached output against grid row settings
#'
#' @keywords internal
.simCompareValidateScenarioCache <- function(
  cached,
  row,
  nSample = NULL,
  nIter = NULL,
  retryErrors = FALSE
) {
  if (!is.data.frame(cached) || nrow(cached) == 0L) {
    return(FALSE)
  }

  has_error <- "error" %in%
    names(cached) &&
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
      if (
        is.na(cached_val) || as.character(cached_val) != as.character(row_val)
      ) {
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
      if (
        is.na(cached_val) || as.character(cached_val) != as.character(row_val)
      ) {
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
          nm %in%
            c(
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
    if (
      !is.null(nIter) &&
        !is.null(nSample) &&
        "method" %in% names(cached) &&
        !.simComparePrimaryOutputComplete(
          cached,
          nSample = nSample,
          nIter = nIter
        )
    ) {
      return(FALSE)
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
    if ("sim_seed" %in% names(row)) {
      paste0("sim_seed = ", row$sim_seed[[1]])
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
      "stim_mean_shift_clusters" %in%
        names(row) &&
        !is.na(row$stim_mean_shift_clusters[[1]])
    ) {
      paste0("shift_clusters = ", row$stim_mean_shift_clusters[[1]])
    },
    if ("stim_sd_multiplier" %in% names(row)) {
      paste0("stim_sd_multiplier = ", row$stim_sd_multiplier[[1]])
    },
    if (
      "stim_sd_multiplier_clusters" %in%
        names(row) &&
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
  clusterGates = FALSE,
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
  dirJobs = NULL,
  totalSims = NULL,
  progressHeading = "COMPARISON SIMULATION PROGRESS",
  ...
) {
  if (!is.logical(clusterGates) || length(clusterGates) != 1L || is.na(clusterGates)) {
    stop("clusterGates must be a single TRUE or FALSE.")
  }

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

  # With `dirJobs`, keep per-row marker files and rewrite `pathProgress` as the
  # same summary the bandwidth analyses use; otherwise append log lines.
  use_dashboard <- !is.null(dirJobs) && nzchar(dirJobs)
  report <- function(status, msg) {
    if (!use_dashboard) {
      .simCompareLogMessage(pathProgress, msg)
      return(invisible(NULL))
    }
    markers <- file.path(
      dirJobs, paste0(c("running-", "completed-", "error-"), sim_id)
    )
    names(markers) <- c("running", "completed", "error")
    dir.create(dirJobs, recursive = TRUE, showWarnings = FALSE)
    unlink(markers)
    file.create(markers[[status]])
    if (!is.null(pathProgress) && nzchar(pathProgress)) {
      .update_progress_summary(
        path_progress_file = pathProgress,
        dir_jobs_chunk = dirJobs,
        total_sims = totalSims,
        sim_grid_chunk_index = sim_grid_chunk_index,
        sim_grid_n_chunks = sim_grid_n_chunks,
        dir_output = dirCache,
        heading = progressHeading
      )
    }
    invisible(NULL)
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
      report(
        "completed",
        paste0(
          "Skipped (cached): ",
          .simCompareFormatScenarioLog(row, sim_id)
        )
      )
      return(cached)
    }
  }

  if (
    "sim_seed" %in% names(row) &&
      length(row$sim_seed) > 0L &&
      is.finite(as.numeric(row$sim_seed[[1]]))
  ) {
    rng_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rng_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      do.call(RNGkind, as.list(rng_kind))
      if (had_seed) {
        assign(".Random.seed", rng_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    sim_seed <- as.integer(row$sim_seed[[1]])
    set.seed(
      sim_seed, kind = "Mersenne-Twister",
      normal.kind = "Inversion", sample.kind = "Rejection"
    )
  }

  settings_log <- .simCompareFormatScenarioLog(row, sim_id)
  report("running", paste0("Running: ", settings_log))

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
  } else if ("mismatch_type" %in% names(row) &&
    identical(row$mismatch_type[[1]], "mean_shift_negative")) {
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
  } else if ("mismatch_type" %in% names(row) &&
    identical(row$mismatch_type[[1]], "sd_inflation_negative")) {
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
        # A missing bias_uns lets StimGate derive it from the bandwidth.
        biasUns = if (!"bias_uns" %in% names(row)) {
          0
        } else if (is.na(row$bias_uns[[1]])) {
          NULL
        } else {
          row$bias_uns[[1]]
        },
        biasUnsFactor = if ("bias_uns_factor" %in% names(row)) {
          row$bias_uns_factor[[1]]
        } else {
          1
        },
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
        clusterGates = clusterGates,
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

      cells <- attr(sim_res, "cells")
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

      report("completed", paste0("Completed: ", settings_log))
      if (!is.null(p)) {
        p(sprintf("Completed sim_id: %s", sim_id))
      }

      # Cell-level data (only with `keepCells = TRUE`) are returned, not cached.
      if (!is.null(cells)) {
        attr(res, "cells") <- cells
      }
      res
    },
    error = function(e) {
      report("error", paste0("Error [", settings_log, "]: ", e$message))
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
          nTruePos = NA_integer_,
          nFalsePos = NA_integer_,
          nFalseNeg = NA_integer_,
          nTrueNeg = NA_integer_,
          unsExprSum = NA_real_,
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
  clusterGates = FALSE,
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
  dirJobs = NULL,
  progressHeading = "COMPARISON SIMULATION PROGRESS",
  ...
) {
  if (!is.logical(clusterGates) || length(clusterGates) != 1L || is.na(clusterGates)) {
    stop("clusterGates must be a single TRUE or FALSE.")
  }

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
      clusterGates = clusterGates,
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
      dirJobs = dirJobs,
      totalSims = nrow(sim_grid),
      progressHeading = progressHeading,
      ...
    )
  }

  if (
    !requireNamespace("furrr", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE)
  ) {
    out_list <- lapply(seq_len(nrow(sim_grid)), function(i) {
      run_one(i, p = NULL)
    })
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
  sim_grid = NULL
) {
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

#' Validate primary comparison output coverage for a simulation grid
#'
#' @keywords internal
.simCompareGridOutputStatus <- function(
  .data,
  sim_grid,
  nSample,
  nIter
) {
  expected_ids <- if (
    is.data.frame(sim_grid) &&
      "sim_id" %in% names(sim_grid)
  ) {
    sort(unique(as.integer(sim_grid$sim_id)))
  } else {
    integer()
  }

  observed_ids <- if (
    is.data.frame(.data) &&
      nrow(.data) > 0L &&
      "sim_id" %in% names(.data)
  ) {
    sort(unique(as.integer(.data$sim_id[!is.na(.data$sim_id)])))
  } else {
    integer()
  }

  extra_ids <- setdiff(observed_ids, expected_ids)
  missing_ids <- setdiff(expected_ids, observed_ids)

  if (length(expected_ids) == 0L) {
    return(list(
      expected_ids = expected_ids,
      observed_ids = observed_ids,
      completed_ids = integer(),
      failed_ids = integer(),
      missing_ids = integer(),
      extra_ids = extra_ids,
      collate_ok = length(extra_ids) == 0L,
      validation_ok = length(extra_ids) == 0L
    ))
  }

  completed_ids <- integer()
  failed_ids <- integer()

  for (sim_id in intersect(expected_ids, observed_ids)) {
    sim_data <- .data[as.integer(.data$sim_id) == sim_id, , drop = FALSE]
    has_error <- "error" %in% names(sim_data) &&
      any(!is.na(sim_data$error) & nzchar(as.character(sim_data$error)))
    complete <- !has_error &&
      .simComparePrimaryOutputComplete(
        sim_data,
        nSample = nSample,
        nIter = nIter
      )

    if (isTRUE(complete)) {
      completed_ids <- c(completed_ids, sim_id)
    } else {
      failed_ids <- c(failed_ids, sim_id)
    }
  }

  collate_ok <- length(missing_ids) == 0L &&
    length(extra_ids) == 0L
  validation_ok <- collate_ok &&
    length(failed_ids) == 0L &&
    identical(sort(completed_ids), expected_ids)

  list(
    expected_ids = expected_ids,
    observed_ids = observed_ids,
    completed_ids = sort(completed_ids),
    failed_ids = sort(failed_ids),
    missing_ids = sort(missing_ids),
    extra_ids = sort(extra_ids),
    collate_ok = collate_ok,
    validation_ok = validation_ok
  )
}

#' Summarise comparison runs by scenario and method
#'
#' @keywords internal
.simCompareSummariseFreqBs <- function(
  .data,
  scenarioCols = NULL,
  keepMethods = c("stimgate", "fbeta", "tailgate")
) {
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

# Paired comparisons use the simulated dataset (iter), not its dependent tubes,
# as the independent unit. Scenario columns must include every varying setting.
# Method summaries use finite tubes; dataset-pair coverage stays explicit.
.simCompareDatasetDifferences <- function(
    .data, scenarioCols,
    outcomes = c("abs_error", "abs_rel_error"),
    competitors = c("fbeta", "tailgate")) {
  scenarioCols <- setdiff(
    scenarioCols, c("method", "approach", "sim_id", "sim_seed", "iter", "sample", "ind")
  )
  allowed <- c("abs_error", "abs_rel_error", "fdp", "sensitivity")
  if (!length(outcomes) || any(!outcomes %in% allowed)) {
    stop("Unknown dataset comparison outcome")
  }
  required <- c(scenarioCols, "iter", "method", "sample")
  if (!all(required %in% names(.data)) || !is.numeric(.data$iter) || any(!is.finite(.data$iter))) {
    stop("Dataset comparisons require scenario columns, method, sample and finite iter IDs")
  }
  primary <- .data[.data$method %in% c("stimgate", competitors), , drop = FALSE]
  keys <- c(scenarioCols, "iter", "method", "sample")
  if (anyDuplicated(primary[keys])) {
    stop("Dataset comparisons require one primary row per tube and method")
  }
  if (any(outcomes %in% c("fdp", "sensitivity"))) {
    primary <- .simCompareClassificationMetrics(primary)
  }
  if (any(outcomes %in% c("abs_error", "abs_rel_error"))) {
    primary <- primary |>
      dplyr::mutate(
        abs_error = abs(.data$propRespEst - .data$propRespTruth),
        abs_rel_error = dplyr::if_else(
          .data$propRespTruth != 0,
          .data$abs_error / abs(.data$propRespTruth), NA_real_
        )
      )
  }
  dataset <- purrr::map_dfr(outcomes, function(outcome) {
    primary |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(scenarioCols, "iter", "method")))) |>
      dplyr::summarise(
        n_tube = dplyr::n(),
        n_finite = sum(is.finite(.data[[outcome]])),
        value = {
          x <- .data[[outcome]]
          x <- x[is.finite(x)]
          if (!length(x)) NA_real_ else if (outcome %in% c("abs_error", "abs_rel_error")) {
            mean(x)
          } else {
            stats::median(x)
          }
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(outcome = outcome)
  })
  join_cols <- c(scenarioCols, "iter", "outcome")
  reference <- dataset |>
    dplyr::filter(.data$method == "stimgate") |>
    dplyr::select(-"method")
  purrr::map_dfr(competitors, function(competitor) {
    other <- dataset |>
      dplyr::filter(.data$method == competitor) |>
      dplyr::select(-"method")
    dplyr::full_join(reference, other, by = join_cols, suffix = c("_stim", "_other")) |>
      dplyr::mutate(difference = .data$value_stim - .data$value_other) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(scenarioCols, "outcome")))) |>
      dplyr::summarise(
        n_dataset = dplyr::n(),
        n_pair = sum(is.finite(.data$difference)),
        tube_coverage_stimgate = if (sum(.data$n_tube_stim, na.rm = TRUE) > 0) {
          sum(.data$n_finite_stim, na.rm = TRUE) / sum(.data$n_tube_stim, na.rm = TRUE)
        } else NA_real_,
        tube_coverage_competitor = if (sum(.data$n_tube_other, na.rm = TRUE) > 0) {
          sum(.data$n_finite_other, na.rm = TRUE) / sum(.data$n_tube_other, na.rm = TRUE)
        } else NA_real_,
        mean_difference = {
          x <- .data$difference[is.finite(.data$difference)]
          if (length(x)) mean(x) else NA_real_
        },
        half_width = {
          x <- .data$difference[is.finite(.data$difference)]
          if (length(x) >= 5L) 1.96 * stats::sd(x) / sqrt(length(x)) else NA_real_
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(
        method = competitor,
        pair_coverage = .data$n_pair / .data$n_dataset,
        lower = .data$mean_difference - .data$half_width,
        upper = .data$mean_difference + .data$half_width
      ) |>
      dplyr::select(-"half_width")
  })
}

# Promote only a complete cross-chunk comparison grid.
.simComparePromoteIfReady <- function(
    run_ctx,
    sim_grid_all,
    total_sims,
    completed_sims,
    failed_sims,
    nSample,
    nIter) {
  if (isTRUE(run_ctx$read_only) || !.analysis_can_promote(run_ctx)) {
    return(invisible(FALSE))
  }

  scenario_paths <- list.files(
    run_ctx$staging_run_dir,
    pattern = "^(compare_raw.*|sim_scenario.*|sim_raw.*)sim_id_[0-9]+[.]rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  compare_raw_full <- .simCompareCollateScenarioOutputs(
    pathList = scenario_paths,
    sim_grid = sim_grid_all
  )
  expected_sim_ids <- sort(unique(as.integer(sim_grid_all$sim_id)))
  full_check <- .simCompareGridOutputStatus(
    compare_raw_full,
    sim_grid = sim_grid_all,
    nSample = nSample,
    nIter = nIter
  )
  full_collate_ok <-
    length(scenario_paths) == length(expected_sim_ids) &&
    isTRUE(full_check$collate_ok)
  full_validation_ok <-
    full_collate_ok &&
    isTRUE(full_check$validation_ok)

  if (!isTRUE(full_validation_ok)) {
    error_message <- paste0(
      "Refusing to promote comparison: canonical collation did not contain ",
      "exactly the complete error-free simulation grid."
    )
    .analysis_mark_chunk(
      run_ctx = run_ctx,
      total_sims = total_sims,
      completed_sims = completed_sims,
      failed_sims = failed_sims,
      collate_ok = full_collate_ok,
      validation_ok = FALSE,
      error_message = error_message
    )
    stop(error_message)
  }

  path_rds_full <- file.path(run_ctx$staging_collated_dir, "compare_raw.rds")
  .write_rds_atomic(compare_raw_full, path_rds_full)
  invisible(isTRUE(.analysis_promote_run(run_ctx)))
}

# ---------------------------------------------------------------------------
# Figure helpers for analyses 7 and 8. These need `analysis-plot-style.R` (and,
# for the signed-error plots, `sim-bandwidth-analysis-plot.R`) sourced first.
# ---------------------------------------------------------------------------

# Method sets shown for every method-comparison figure. Tailgate performs very
# poorly without more tuning and squashes the other methods' scale, so each
# figure is also drawn without it.
.simCompareMethodSets <- function() {
  list(
    all_methods = list(
      heading = "All methods",
      methods = c("stimgate", "tailgate", "fbeta")
    ),
    no_tailgate = list(
      heading = "Without Tailgate",
      methods = c("stimgate", "fbeta")
    )
  )
}

.simCompareMismatchLabels <- c(
  mean_shift_all = "Shift all stimulated cells",
  mean_shift_negative = "Shift stimulated negatives only",
  sd_inflation = "Inflate stimulated background SD"
)

# Horizontal axis label for each mismatch mechanism.
.simCompareMismatchAxisLabels <- c(
  mean_shift_all =
    "Additive shift of all stimulated cells (transformed expression scale)",
  mean_shift_negative =
    "Additive shift of stimulated negatives (transformed expression scale)",
  sd_inflation = "Fractional increase in stimulated background SD"
)

# Print, save and describe one figure per method set, mean position setting
# and (optionally) a third setting, with a heading for each loop level. Use in
# a `results: asis` chunk. `level` is the heading level of the method set;
# mean position headings are one deeper and the optional third loop deeper
# still. Figures go to `dir/<method set>/<file_fn(pos, extra)>`.
.simCompareFigureLoop <- function(
    data,
    make_plot,
    dir,
    file_fn,
    height,
    level,
    extra_col = NULL,
    extra_heading = function(x) as.character(x),
    method_col = "method",
    pos_col = "mean_pos_setting",
    allow_tall = FALSE) {
  if (level + 1L + as.integer(!is.null(extra_col)) > 6L) {
    stop("Figure loop headings would be deeper than level 6.")
  }
  sets <- .simCompareMethodSets()
  for (set_name in names(sets)) {
    set <- sets[[set_name]]
    .analysis_heading(set$heading, level)
    set_data <- data[data[[method_col]] %in% set$methods, , drop = FALSE]
    for (pos in unique(as.character(data[[pos_col]]))) {
      .analysis_heading(paste0("Mean position: ", pos), level + 1L)
      in_pos <- as.character(data[[pos_col]]) == pos
      pos_data <- set_data[as.character(set_data[[pos_col]]) == pos, ,
        drop = FALSE
      ]
      # One heading for all `extra_col` columns together: a single column is
      # passed on as a value, several as a list (one element per column), so
      # that the deepest heading can name them all within level 6.
      extras <- if (is.null(extra_col)) {
        list(NA)
      } else {
        .simCompareLoopExtras(data[in_pos, extra_col, drop = FALSE])
      }
      for (extra in extras) {
        curr <- pos_data
        if (!is.null(extra_col)) {
          .analysis_heading(extra_heading(extra), level + 2L)
          keep <- rep(TRUE, nrow(curr))
          for (col in extra_col) {
            value <- if (length(extra_col) == 1L) extra else extra[[col]]
            keep <- keep & curr[[col]] == value
          }
          curr <- curr[keep, , drop = FALSE]
        }
        if (nrow(curr) == 0L) next
        p <- make_plot(curr)
        .analysis_save_fig(
          p, file.path(dir, set_name, file_fn(pos, extra)),
          height = height, allow_tall = allow_tall
        )
        .analysis_print_fig(p)
      }
    }
  }
  invisible(NULL)
}

# Distinct, sorted combinations of the loop columns, as a list.
.simCompareLoopExtras <- function(extra_data) {
  extra_data <- unique(extra_data)
  extra_data <- extra_data[do.call(order, unname(as.list(extra_data))), ,
    drop = FALSE
  ]
  if (ncol(extra_data) == 1L) {
    return(as.list(extra_data[[1]]))
  }
  lapply(seq_len(nrow(extra_data)), function(i) as.list(extra_data[i, ]))
}

# Error statistic against stimulated cell count, one panel per response
# frequency and transformation.
.simComparePlotByCell <- function(
    data, y, y_label, zero_line = FALSE, free_y = FALSE, percent_y = FALSE) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = n_cell, y = .data[[y]], colour = method)
  )
  if (zero_line) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0, colour = "gray25", linetype = "dashed"
    )
  }
  p +
    ggplot2::geom_line(alpha = 0.75) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_x_log10(
      breaks = sort(unique(data$n_cell)),
      labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    ) +
    ggplot2::scale_y_continuous(
      labels = if (percent_y) .analysis_label_percent else .analysis_label_number
    ) +
    .analysis_scale_method() +
    ggplot2::facet_grid(
      prob_response ~ transformation,
      scales = if (free_y) "free_y" else "fixed",
      labeller = ggplot2::labeller(
        prob_response = .analysis_labeller_percent()
      )
    ) +
    ggplot2::labs(
      x = "Number of stimulated cells", y = y_label, colour = "Method"
    ) +
    .analysis_theme()
}

# Estimated against true background-subtracted frequency.
.simComparePlotEstVsTruth <- function(data) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  top <- max(data$propRespTruth, na.rm = TRUE) * 2
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = propRespTruth, y = propRespEst, colour = method)
  ) +
    ggplot2::geom_point(alpha = 0.35, size = 1) +
    ggplot2::expand_limits(x = c(0, top), y = c(0, top)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::scale_x_continuous(labels = .analysis_label_percent) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
    .analysis_scale_method() +
    ggplot2::facet_wrap(~transformation) +
    ggplot2::labs(
      x = "True background-subtracted response frequency",
      y = "Estimated background-subtracted response frequency",
      colour = "Method"
    ) +
    .analysis_theme()
}

# Share of estimates that used a threshold fallback.
.simComparePlotFallback <- function(data) {
  .simComparePlotByCell(
    data, "fallback_rate", "Share of estimates using a threshold fallback",
    percent_y = TRUE
  )
}

# Histograms of finite thresholds by method; `data` has `approach`.
.simComparePlotThresholdDensity <- function(data) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = threshold, y = ggplot2::after_stat(density),
      fill = approach, colour = approach
    )
  ) +
    ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    ggplot2::facet_grid(
      prob_response ~ transformation, scales = "free",
      labeller = ggplot2::labeller(
        prob_response = .analysis_labeller_percent()
      )
    ) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number) +
    .analysis_scale_method(aesthetics = c("colour", "fill")) +
    ggplot2::labs(
      x = "Threshold", y = "Density", colour = "Method", fill = "Method"
    ) +
    .analysis_theme()
}

# Shared pieces of the analysis 8 mismatch plots.
.simCompareMismatchScales <- function(x_label) {
  list(
    ggplot2::scale_x_continuous(labels = .analysis_label_number),
    .analysis_scale_method(),
    ggplot2::labs(x = x_label, colour = "Method", subtitle = "Tube-level summaries"),
    .analysis_theme()
  )
}

.simCompareStripWrap <- function() ggplot2::label_wrap_gen(width = 22)

# Share of unstimulated cells above the gate against the mean shift.
.simComparePlotPlacement <- function(data) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = mismatch_val, y = propUns_mean, colour = method,
      linetype = mismatch_type,
      group = interaction(method, mismatch_type)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 2, alpha = 0.75) +
    ggplot2::facet_wrap(
      ~scenario_desc, scales = "free_y", labeller = .simCompareStripWrap()
    ) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
    ggplot2::scale_linetype_discrete(labels = .simCompareMismatchLabels) +
    .simCompareMismatchScales("Additive mean shift (transformed expression scale)") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::labs(
      y = "Unstimulated cells above the gate",
      linetype = "Mismatch variant"
    )
}

# 90th percentile error for every mismatch mechanism.
.simComparePlotUpperTail <- function(data) {
  data$mismatch_type <- factor(
    data$mismatch_type, levels = names(.simCompareMismatchLabels)
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = mismatch_val, y = q90_abs_rel_error, colour = method, group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 2, alpha = 0.75) +
    ggplot2::facet_grid(
      mismatch_type ~ scenario_desc, scales = "free",
      labeller = ggplot2::labeller(
        mismatch_type = .simCompareMismatchLabels,
        scenario_desc = .simCompareStripWrap()
      )
    ) +
    ggplot2::scale_y_continuous(
      transform = scales::asinh_trans(), labels = .analysis_label_number
    ) +
    .simCompareMismatchScales("Mismatch size") +
    ggplot2::labs(y = "Tube-level 90th percentile absolute relative error (asinh scale)")
}

# Over- and under-estimate summary of signed relative error per scenario and
# method, from raw comparison rows, with the same rows and estimand as
# `.simCompareSummariseFreqBs()`.
.simCompareSignedErrorSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::mutate(
      rel_error = dplyr::if_else(
        .data$propRespTruth != 0,
        (.data$propRespEst - .data$propRespTruth) / .data$propRespTruth,
        NA_real_
      )
    ) |>
    .simBandwidthSignedErrorSummary(scenarioCols)
}

# Size of the relative error, |estimate - truth| / truth, per scenario and
# method: median, 95th percentile and maximum over the sample-level errors.
# Uses the same rows and estimand as `.simCompareSummariseFreqBs()`; samples
# with a true frequency of zero have no relative error and are left out.
.simCompareUnsignedErrorSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::mutate(
      abs_rel_error = dplyr::if_else(
        .data$propRespTruth != 0,
        abs((.data$propRespEst - .data$propRespTruth) / .data$propRespTruth),
        NA_real_
      )
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      median = stats::median(.data$abs_rel_error, na.rm = TRUE),
      q95 = unname(stats::quantile(
        .data$abs_rel_error, probs = 0.95, na.rm = TRUE
      )),
      max = suppressWarnings(max(.data$abs_rel_error, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(
      c("median", "q95", "max"),
      ~ dplyr::if_else(is.finite(.x), .x, NA_real_)
    ))
}

# Average the scenario statistics equally over the scenarios in each group
# (each scenario counts once, however many samples it has). Scenarios are
# averaged only over the columns left out of `group_cols`.
.simCompareErrorAverage <- function(tbl, group_cols) {
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(c("median", "q95", "max")),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

# Panels for the mismatch plots: rows are statistics, columns are
# transformations (and response probabilities when `by_prob`).
.simCompareMismatchFacet <- function(by_prob) {
  ggplot2::facet_grid(
    if (by_prob) statistic ~ transformation + prob_response else
      statistic ~ transformation,
    scales = "free",
    labeller = ggplot2::labeller(
      prob_response = .analysis_labeller_percent("Response probability: ")
    )
  )
}

# Size of the relative error against mismatch size, coloured by method.
# `tbl` comes from `.simCompareUnsignedErrorSummary()` (optionally averaged).
.simComparePlotMismatchError <- function(
    tbl,
    x_label = "Mismatch size",
    by_prob = FALSE,
    stat_cols = c(median = "Median", q95 = "95th percentile", max = "Maximum")) {
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols)
  # Errors above 1500% (16 times the truth) are drawn at the cap.
  capped <- .simBandwidthSignedErrorIsCapped(tbl$value)
  tbl$value_shown <- .simBandwidthSignedErrorSquish(tbl$value)
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = mismatch_val, y = value_shown, colour = method, group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 1.5, alpha = 0.75) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    .simBandwidthAbsErrorLayers(capped = capped) +
    .analysis_scale_method() +
    .simCompareMismatchFacet(by_prob) +
    ggplot2::labs(x = x_label, colour = "Method", subtitle = "Tube-level error distributions") +
    .analysis_theme()
}

# Signed relative error against `x`: rows of panels are statistics, columns
# are transformations (and response probabilities when `by_prob`);
# over-estimates sit above zero and under-estimates below, and line weight is
# the share of estimates in that direction.
.simComparePlotSignedError <- function(
    tbl,
    x = "n_cell",
    x_label = "Number of stimulated cells",
    x_log = TRUE,
    by_prob = FALSE,
    stat_cols = c(median = "Median", q95 = "95th percentile", max = "Maximum")) {
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols)
  x_scale <- if (x_log) {
    ggplot2::scale_x_log10(
      breaks = sort(unique(tbl[[x]])),
      labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    )
  } else {
    ggplot2::scale_x_continuous(labels = .analysis_label_number)
  }
  line_cols <- c(
    "statistic", "transformation", if (by_prob) "prob_response",
    "method", "direction"
  )
  # Over-estimates above +1500% (16 times the truth) are drawn at the cap.
  capped <- .simBandwidthSignedErrorIsCapped(tbl$value)
  tbl$value_shown <- .simBandwidthSignedErrorSquish(tbl$value)
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = .data[[x]], y = value_shown, colour = method,
      group = interaction(method, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers("Tube-level relative error", capped = capped) +
    .simBandwidthSignedErrorSegmentLayer(
      tbl, x, "value_shown", line_cols, alpha = 0.75
    ) +
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    x_scale +
    .analysis_scale_method() +
    .simCompareMismatchFacet(by_prob) +
    ggplot2::labs(x = x_label, colour = "Method", subtitle = "Tube-level error distributions") +
    .analysis_theme()
}

# ---------------------------------------------------------------------------
# Gate purity and detection (analysis 8): tube-level classification of
# stimulated cells against the simulation labels.
# ---------------------------------------------------------------------------

# Quantile of the finite values, NA when there are none.
.simCompareQuantileFinite <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  unname(stats::quantile(x, probs = prob, names = FALSE))
}

# Tube-level distribution per scenario and method. Each
# stimulated tube in each dataset gives one FDP, sensitivity and
# false-positive rate; these are summarised directly, never pooled over cells.
# FDP summaries use only tubes whose gate selected at least one cell
# (`n_fdp_defined`), whereas an empty gate contributes zero sensitivity.
.simCompareClassificationSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  q <- .simCompareQuantileFinite
  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    .simCompareClassificationMetrics() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_valid = sum(.data$gate_status != "failed"),
      n_failed = sum(.data$gate_status == "failed"),
      n_fdp_defined = sum(is.finite(.data$fdp)),
      n_empty = sum(.data$gate_empty %in% TRUE),
      n_fallback = sum(
        .data$gate_status %in% c("fallback_empty", "fallback_selected")
      ),
      n_fallback_empty = sum(.data$gate_status == "fallback_empty"),
      fdp_median = q(.data$fdp, 0.5),
      fdp_q90 = q(.data$fdp, 0.9),
      sensitivity_median = q(.data$sensitivity, 0.5),
      sensitivity_q10 = q(.data$sensitivity, 0.1),
      fpr_median = q(.data$false_positive_rate, 0.5),
      fpr_q90 = q(.data$false_positive_rate, 0.9),
      selected_fraction_median = q(.data$selected_fraction, 0.5),
      prevalence_median = q(.data$n_genuine_pos / .data$n_classified, 0.5),
      .groups = "drop"
    )
}

# Pairing check: within each baseline scenario, replicate and sample, every
# mismatch setting must use the same simulated draws. The unstimulated tube is
# never changed by a mismatch and the labels are not changed by it, so the
# unstimulated-expression fingerprint, the number of genuine positives and the
# stimulated cell count must each take one value. Returns one row per group.
.simComparePairingCheck <- function(
    .data,
    pairCols = "base_scenario_id",
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(pairCols, "iter", "sample")))) |>
    dplyr::summarise(
      n_settings = dplyr::n_distinct(.data$sim_id),
      n_uns_values = dplyr::n_distinct(.data$unsExprSum),
      n_genuine_pos_values = dplyr::n_distinct(.data$nTruePos + .data$nFalseNeg),
      n_cell_values = dplyr::n_distinct(.data$nCellStim),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      paired = .data$n_uns_values == 1L &
        .data$n_genuine_pos_values == 1L &
        .data$n_cell_values == 1L
    )
}

# Agreement of the zero-mismatch rows of each mechanism with the zero-shift
# "shift all stimulated cells" rows on the same replicate, sample and method.
# The shift variants add exactly zero, so they must agree exactly; zero SD
# inflation rescales by one, which can change the last bit of a value.
.simCompareZeroMismatchAgreement <- function(
    .data,
    pairCols = "base_scenario_id",
    reference = "mean_shift_all",
    keepMethods = c("stimgate", "fbeta", "tailgate")) {
  keys <- c(pairCols, "iter", "sample", "method")
  zero <- .data |>
    dplyr::filter(.data$method %in% keepMethods, .data$mismatch_val == 0) |>
    dplyr::select(dplyr::all_of(c(
      keys, "mismatch_type", "threshold", .simCompareCountCols
    )))
  ref <- zero |>
    dplyr::filter(.data$mismatch_type == reference) |>
    dplyr::select(-"mismatch_type") |>
    dplyr::rename_with(~ paste0(.x, "_ref"), -dplyr::all_of(keys))
  zero |>
    dplyr::filter(.data$mismatch_type != reference) |>
    dplyr::inner_join(ref, by = keys) |>
    dplyr::mutate(
      same_threshold = .data$threshold == .data$threshold_ref,
      same_counts = .data$nTruePos == .data$nTruePos_ref &
        .data$nFalsePos == .data$nFalsePos_ref &
        .data$nFalseNeg == .data$nFalseNeg_ref &
        .data$nTrueNeg == .data$nTrueNeg_ref
    ) |>
    dplyr::group_by(.data$mismatch_type, .data$method) |>
    dplyr::summarise(
      n_compared = dplyr::n(),
      n_same_threshold = sum(.data$same_threshold %in% TRUE),
      n_same_counts = sum(.data$same_counts %in% TRUE),
      max_abs_threshold_diff = max(abs(.data$threshold - .data$threshold_ref)),
      .groups = "drop"
    )
}

# Facet label for each baseline scenario: transformation and description,
# ordered Gaussian, Skew, Gamma and then by baseline scenario.
.simCompareScenarioLabel <- function(data) {
  trans <- .analysis_trans_factor(data$transformation)
  lab <- paste0(as.character(trans), ": ", data$scenario_desc)
  ord <- order(as.integer(trans), data$base_scenario_id)
  factor(lab, levels = unique(lab[ord]))
}

.simCompareClassificationOutcomes <- list(
  fdp = c(
    label = "False discovery proportion",
    median = "fdp_median", tail = "fdp_q90"
  ),
  sensitivity = c(
    label = "Sensitivity",
    median = "sensitivity_median", tail = "sensitivity_q10"
  ),
  fpr = c(
    label = "False-positive rate",
    median = "fpr_median", tail = "fpr_q90"
  )
)

# Classification outcomes against mismatch size from
# `.simCompareClassificationSummary()`. Rows of panels are outcomes, columns
# are baseline scenarios, each with its own horizontal range; colour is the
# method and line type the statistic (median, or the worse tail: 90th
# percentile for FDP and false-positive rate, 10th for sensitivity). With
# `unit_scale`, every vertical scale is fixed at 0-100%.
.simComparePlotClassification <- function(
    tbl,
    outcomes = c("fdp", "sensitivity"),
    x_label = "Mismatch size",
    unit_scale = TRUE) {
  spec <- .simCompareClassificationOutcomes[outcomes]
  long <- purrr::map_df(names(spec), function(nm) {
    s <- spec[[nm]]
    dplyr::bind_rows(
      dplyr::mutate(tbl, outcome = s[["label"]], statistic = "median",
        value = .data[[s[["median"]]]]),
      dplyr::mutate(tbl, outcome = s[["label"]], statistic = "tail",
        value = .data[[s[["tail"]]]])
    )
  })
  long$outcome <- factor(
    long$outcome,
    levels = vapply(spec, function(s) s[["label"]], character(1))
  )
  long$scenario <- .simCompareScenarioLabel(long)
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = mismatch_val, y = value, colour = method, linetype = statistic,
      group = interaction(method, statistic)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::facet_grid(
      outcome ~ scenario,
      scales = if (unit_scale) "free_x" else "free",
      labeller = ggplot2::labeller(
        scenario = ggplot2::label_wrap_gen(width = 16),
        outcome = ggplot2::label_wrap_gen(width = 14)
      )
    ) +
    ggplot2::scale_x_continuous(
      transform = "sqrt", labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    ) +
    .analysis_scale_method() +
    ggplot2::scale_linetype_manual(
      values = c(median = "solid", tail = "22"),
      labels = c(
        median = "Tube-level median",
        tail = "Tube-level 90th percentile (FDP, FPR) or 10th (sensitivity)"
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2, ncol = 1)
    ) +
    ggplot2::labs(
      x = x_label, y = NULL, colour = "Method", linetype = "Tube-level distribution"
    ) +
    .analysis_theme()
  if (unit_scale) {
    p + ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = .analysis_label_percent
    )
  } else {
    p + ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
      ggplot2::expand_limits(y = 0)
  }
}

# Counts behind each plotted setting, one column per method:
# "valid / FDP defined / empty gates / fallbacks".
.simCompareClassificationCountTable <- function(summary_tbl) {
  method_labels <- .analysis_method_labels[
    intersect(names(.analysis_method_labels), unique(summary_tbl$method))
  ]
  summary_tbl |>
    dplyr::mutate(
      scenario = .simCompareScenarioLabel(summary_tbl),
      counts = paste(
        .data$n_valid, .data$n_fdp_defined, .data$n_empty, .data$n_fallback,
        sep = " / "
      ),
      method = factor(.data$method, levels = names(method_labels))
    ) |>
    dplyr::arrange(.data$scenario, .data$mismatch_val, .data$method) |>
    dplyr::mutate(
      method = unname(method_labels[as.character(.data$method)]),
      mismatch_val = .analysis_label_number(.data$mismatch_val)
    ) |>
    dplyr::select("scenario", "mismatch_val", "method", "counts") |>
    tidyr::pivot_wider(names_from = "method", values_from = "counts") |>
    dplyr::rename(Scenario = "scenario", `Mismatch size` = "mismatch_val")
}

# ---------------------------------------------------------------------------
# Gate diagnostic for one baseline scenario: the cells, their labels and the
# final gates, on matched replicates across mismatch settings.
# ---------------------------------------------------------------------------

# Rerun `rows` (from the full grid, so with their production `sim_seed`) for
# one replicate, keeping the cells and StimGate's sample-level details. With
# the same `nSample`, the replicate equals iteration 1 of the production run
# of each row. Returns `list(results, cells)`.
.simCompareGateDiagnosticRun <- function(rows, nSample, ...) {
  out <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    res <- .simCompareRunScenario(
      row = row,
      nSample = nSample,
      nIter = 1L,
      resume = FALSE,
      includeLocDetails = TRUE,
      keepCells = TRUE,
      ...
    )
    cells <- attr(res, "cells")
    if (is.null(cells)) {
      stop("Gate diagnostic for sim_id ", row$sim_id[[1]], " failed: ",
        paste(unique(stats::na.omit(res$error)), collapse = "; "))
    }
    attr(res, "cells") <- NULL
    list(
      results = res,
      cells = dplyr::bind_cols(
        row[rep(1L, nrow(cells)), c("sim_id", "mismatch_type", "mismatch_val")],
        cells
      )
    )
  })
  list(
    results = purrr::list_rbind(purrr::map(out, "results")),
    cells = purrr::list_rbind(purrr::map(out, "cells"))
  )
}

# Stimulated cells stacked by true label, the unstimulated tube as an outline
# and each method's final gate, for one sample. Rows of panels are mismatch
# sizes and columns the mismatch mechanisms. `gates` has `mismatch_type`,
# `mismatch_val`, `method` and `threshold`.
.simComparePlotGateDiagnostic <- function(
    cells,
    gates,
    x_label = "Expression (transformed scale)",
    bins = 80) {
  binwidth <- diff(range(cells$expr)) / bins
  boundary <- min(cells$expr)
  type_levels <- intersect(names(.simCompareMismatchLabels), cells$mismatch_type)
  prep <- function(d) {
    d$mismatch_type <- factor(d$mismatch_type, levels = type_levels)
    d
  }
  stim <- prep(cells[cells$condition == "stim", , drop = FALSE])
  stim$label <- factor(
    ifelse(stim$label == "gp", "Genuine positive", "Genuine negative"),
    levels = c("Genuine negative", "Genuine positive")
  )
  uns <- prep(cells[cells$condition == "unstim", , drop = FALSE])
  gates <- prep(gates)
  # A distinct group per line keeps coincident gates visible.
  gates$line_id <- seq_len(nrow(gates))
  ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = stim,
      ggplot2::aes(x = expr, fill = label),
      binwidth = binwidth, boundary = boundary, position = "stack"
    ) +
    ggplot2::geom_freqpoly(
      data = uns,
      ggplot2::aes(x = expr, linetype = "Unstimulated tube"),
      binwidth = binwidth, boundary = boundary, colour = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      data = gates,
      ggplot2::aes(xintercept = threshold, colour = method, group = line_id),
      linewidth = 0.7
    ) +
    ggplot2::facet_grid(
      mismatch_val ~ mismatch_type,
      labeller = ggplot2::labeller(
        mismatch_val = function(x) paste0("Shift: ", .analysis_label_number(as.numeric(x))),
        mismatch_type = .simCompareMismatchLabels
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c("Genuine negative" = "grey78", "Genuine positive" = "grey35")
    ) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    # A square-root scale keeps the small positive component visible.
    ggplot2::scale_y_sqrt(labels = .analysis_label_number) +
    .analysis_scale_method() +
    ggplot2::labs(
      x = x_label, y = "Number of cells (square-root scale)",
      fill = "Stimulated cells",
      colour = "Final gate", linetype = NULL
    ) +
    .analysis_theme()
}
