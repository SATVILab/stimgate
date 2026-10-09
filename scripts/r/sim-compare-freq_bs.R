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
#' With `removeZero = TRUE`, cells with an expression of exactly zero are
#' dropped from both tubes before the histograms are built (in CyTOF data they
#' are clearly negative), and each pdf is scaled by its tube's retained
#' fraction so that it integrates to that fraction rather than to one.
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
  numBins = NULL,
  removeZero = FALSE
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to call fbeta.py.")
  }

  xUns <- as.numeric(xUns)
  xStim <- as.numeric(xStim)
  xUns <- xUns[is.finite(xUns)]
  xStim <- xStim[is.finite(xStim)]

  negScale <- 1
  posScale <- 1
  if (isTRUE(removeZero)) {
    nUns <- length(xUns)
    nStim <- length(xStim)
    xUns <- xUns[xUns != 0]
    xStim <- xStim[xStim != 0]
    negScale <- length(xUns) / nUns
    posScale <- length(xStim) / nStim
  }

  if (length(xUns) < 2L || length(xStim) < 2L) {
    stop("F-beta estimation requires at least two finite cells per tube.")
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
    numBins = if (is.null(numBins)) NULL else as.integer(numBins),
    negScale = negScale,
    posScale = posScale
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
#' `tol` is an absolute bound on the first derivative of the density, so it
#' depends on the scale of `x`; `autoTol = TRUE` replaces it with 1% of the
#' largest absolute derivative. A finite cutpoint is moved up by `bias`, as
#' StimGate's gates are. With `removeZero = TRUE`, cells with an expression of
#' exactly zero are dropped before the density is estimated.
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
  autoTol = FALSE,
  bias = 0,
  removeZero = FALSE
) {
  method <- match.arg(method)
  if (length(bias) != 1L || !is.finite(bias)) {
    stop("Tailgate bias must be a finite scalar.")
  }
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (isTRUE(removeZero)) {
    x <- x[x != 0]
  }

  if (length(x) < 2L || length(unique(x)) < 2L) {
    stop("Tailgate estimation requires at least two distinct finite cells.")
  }

  if (!requireNamespace("cytoUtils", quietly = TRUE)) {
    stop("Package 'cytoUtils' is required for tailgate comparisons.")
  }

  bandwidthUse <- if (is.null(bandwidth)) {
    suppressWarnings(ks::hpi(x, deriv.order = 1L))
  } else {
    bandwidth
  }

  if (
    length(bandwidthUse) != 1L ||
      !is.finite(bandwidthUse) ||
      bandwidthUse <= 0
  ) {
    stop("Tailgate bandwidth must be a finite positive scalar.")
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

  threshold <- as.numeric(threshold)[1] + bias
  list(
    threshold = threshold,
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
# negatives. The F1 score, 2TP / (2TP + FP + FN), is undefined only when the
# tube has no genuine positives and the gate selects no cells. `gate_status` separates failed runs, fallback gates and
# calculated gates, each split by whether any stimulated cell was selected.
# Runtime failures belong to the primary method's failure cohort. Keep their
# error/provenance fields while recognizing the historical diagnostic label.
.simComparePrimaryMethodRows <- function(.data) {
  if (!"method" %in% names(.data)) return(.data)
  .data$method <- as.character(.data$method)
  .data$method[.data$method %in% "stimgate_error"] <- "stimgate"
  .data
}

.simCompareClassificationMetrics <- function(.data) {
  .data <- .simComparePrimaryMethodRows(.data)
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
      dplyr::across(
        dplyr::all_of(.simCompareCountCols),
        ~ dplyr::if_else(has_error, NA_integer_, .x)
      ),
      n_selected = .data$nTruePos + .data$nFalsePos,
      n_genuine_pos = .data$nTruePos + .data$nFalseNeg,
      n_genuine_neg = .data$nFalsePos + .data$nTrueNeg,
      n_classified = .data$n_selected + .data$nFalseNeg + .data$nTrueNeg,
      fdp = ratio(.data$nFalsePos, .data$n_selected),
      sensitivity = ratio(.data$nTruePos, .data$n_genuine_pos),
      false_positive_rate = ratio(.data$nFalsePos, .data$n_genuine_neg),
      f1 = ratio(2 * .data$nTruePos, 2 * .data$nTruePos + .data$nFalsePos + .data$nFalseNeg),
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
  tailgateBias = 0,
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
        fallbackHighValue = is.na(fbetaError) && isTRUE(fallbackHighValue),
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
          autoTol = tailgateAutoTol,
          bias = tailgateBias
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
        fallbackHighValue = is.na(tailgateError) && isTRUE(fallbackHighValue),
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
            "fbeta_error"
          } else if (isTRUE(fbetaEst$thresholdFallbackUsed)) {
            "fbeta_fallback_high_value"
          } else {
            "fbeta_calculated"
          },
          if (!is.na(tailgateError)) {
            "tailgate_error"
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
.simCompareStimgateFailureRows <- function(
  truthTbl,
  errorMessage,
  locThresholdMethod = NA_character_
) {
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
      locThresholdMethod = as.character(.env$locThresholdMethod),
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

#' Read the local-FDR threshold method StimGate saved for one marker
#'
#' Errors when the saved channel settings lack a valid method, so output rows
#' are never labelled with a method StimGate did not record.
#'
#' @keywords internal
.simCompareStimgateLocThresholdMethod <- function(
  pathProject,
  marker,
  chnl = "F1"
) {
  # Saved settings are keyed by marker label; fall back to the channel.
  settings <- stimgate::stimgateMetaReadSettingsChnls(pathProject)
  entry <- settings[[marker]] %||% purrr::detect(
    settings,
    function(x) identical(x[["chnlCut"]], chnl)
  )
  method <- entry[["locThresholdMethod"]]
  if (
    !is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("region", "match", "cap")
  ) {
    stop(
      "StimGate saved no valid locThresholdMethod for marker ", marker, "."
    )
  }
  method
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
  bwScope = "cytokine",
  bwAdj = 1,
  bwNcellMin = bwNcellMax,
  bwNcellMax = 1e4,
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
  locThresholdMethod = "region",
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
          bwScope = bwScope,
          bwAdj = bwAdj,
          bwNcellMin = bwNcellMin,
          bwNcellMax = bwNcellMax,
          bwCluster = bwCluster,
          clusterGates = clusterGates,
          locProbCol = locProbCol,
          locMinPeakProb = locMinPeakProb,
          locEnforceShapeThreshold = locEnforceShapeThreshold,
          locThresholdMethod = locThresholdMethod,
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

      # Record the threshold method StimGate saved for the marker, so output
      # rows show the method actually used rather than only the request.
      locThresholdMethodUsed <- .simCompareStimgateLocThresholdMethod(
        pathProject = pathProject,
        marker = "MarkerF1"
      )

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
              locThresholdMethod = locThresholdMethodUsed,
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
            locReason = NA_character_,
            locThresholdMethod = NA_character_
          )
        )

        detailTbl <- detailTbl |>
          dplyr::mutate(
            approach = "stimgate",
            locThresholdMethod = dplyr::coalesce(
              as.character(.data$locThresholdMethod),
              locThresholdMethodUsed
            ),
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
          locThresholdMethod,
          error,
          dplyr::everything()
        )
    },
    error = function(e) {
      .simCompareStimgateFailureRows(
        truthTbl = truthTbl,
        errorMessage = e$message,
        locThresholdMethod = locThresholdMethod
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
  bwScope = "cytokine",
  bwAdj = 1,
  bwNcellMin = bwNcellMax,
  bwNcellMax = 1e4,
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
  locThresholdMethod = "region",
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
  tailgateBias = 0,
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
      bwScope = bwScope,
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
      locThresholdMethod = locThresholdMethod,
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
      tailgateBias = tailgateBias,
      fallbackHighValue = fallbackHighValue,
      fallbackMargin = fallbackMargin
    )

    dplyr::bind_rows(stimgateTbl, alternativeTbl) |>
      dplyr::left_join(unsTbl, by = "sample") |>
      dplyr::mutate(
        iter = iterNum,
        nCellStimSim = nCellStim,
        nCellUnsSim = nCellUns,
        # NULL biasUns: StimGate sets it from the shared bandwidth.
        biasUns = biasUns %||% NA_real_,
        biasUnsFactor = biasUnsFactor,
        bw = bw,
        bwFallback = bwFallback,
        bwMin = bwMin,
        bwMax = bwMax,
        bwMtd = bwMtd,
        bwScope = bwScope,
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
        tailgateBias = tailgateBias,
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

# Explicit comparator exceptions are completed observations with missing scientific
# outputs, not empty gates or whole-scenario infrastructure failures.
.simCompareRecordedComparatorErrors <- function(data) {
  required <- c("method", "error", "threshold", "thresholdOrigin", "gateReturnPoint",
    "thresholdFallbackUsed", "propRespEst", "nPosStim", .simCompareCountCols)
  if (!is.data.frame(data) || !all(required %in% names(data))) return(rep(FALSE, nrow(data)))
  missing_cols <- intersect(c("threshold", "thresholdMetric", "propRespEst", "propStim", "propUns",
    "nPosStim", "nPosUns", .simCompareCountCols), names(data))
  missing_outputs <- rowSums(!is.na(as.matrix(data[, missing_cols, drop = FALSE]))) == 0L
  recorded <- data$method %in% c("fbeta", "tailgate") &
    !is.na(data$error) & nzchar(trimws(as.character(data$error))) &
    !is.na(data$thresholdOrigin) & startsWith(data$thresholdOrigin, "error: ") &
    !is.na(data$gateReturnPoint) & data$gateReturnPoint == paste0(data$method, "_error") &
    !is.na(data$thresholdFallbackUsed) & !data$thresholdFallbackUsed & missing_outputs
  recorded %in% TRUE
}

# Give one concise reason for invalid structure; method failures remain separate
# from completion so promotion and resume use exactly the same contract.
.simComparePrimaryOutputStatus <- function(
    data, nSample, nIter, methods = c("stimgate", "fbeta", "tailgate")) {
  invalid <- function(reason) list(complete = FALSE, reason = reason)
  required <- c("iter", "sample", "method", "propRespTruth", "propRespEst", "nCellStim",
    "nPosStim", .simCompareCountCols, "unsExprSum")
  if (!is.data.frame(data) || !nrow(data) || !all(required %in% names(data))) {
    return(invalid("missing primary outcome columns or rows"))
  }
  recorded <- .simCompareRecordedComparatorErrors(data)
  error <- if ("error" %in% names(data)) !is.na(data$error) & nzchar(as.character(data$error)) else rep(FALSE, nrow(data))
  if (any(error & !recorded)) return(invalid("unrecognised, malformed or whole-scenario runtime error"))
  primary <- data[data$method %in% methods, , drop = FALSE]
  expected <- as.integer(nSample) * as.integer(nIter)
  if (nrow(primary) != expected * length(methods)) return(invalid("missing or extra primary method/sample rows"))
  if (anyNA(primary$iter) || anyNA(primary$sample) ||
      !all(as.character(primary$iter) %in% as.character(seq_len(nIter))) ||
      !all(as.character(primary$sample) %in% as.character(seq_len(nSample)))) {
    return(invalid("sample or iteration IDs differ from the intended grid"))
  }
  if (anyDuplicated(primary[c("iter", "sample", "method")])) return(invalid("duplicate primary method/sample rows"))
  method_counts <- table(factor(primary$method, levels = methods))
  if (any(method_counts != expected)) return(invalid("incomplete primary method coverage"))
  if (!all(is.finite(primary$propRespTruth)) || !all(is.finite(primary$unsExprSum)) ||
      !all(is.finite(primary$nCellStim) & primary$nCellStim > 0)) {
    return(invalid("missing simulated truth, cell counts or pairing fingerprint"))
  }
  shared_truth <- primary |>
    dplyr::group_by(.data$iter, .data$sample) |>
    dplyr::summarise(dplyr::across(c("propRespTruth", "nCellStim", "unsExprSum"),
      dplyr::n_distinct),
      n_genuine = dplyr::n_distinct(.data$nTruePos + .data$nFalseNeg, na.rm = TRUE),
      .groups = "drop")
  if (any(shared_truth$n_genuine > 1L) || any(as.matrix(shared_truth[c("propRespTruth", "nCellStim", "unsExprSum")]) != 1L)) {
    return(invalid("simulated truth, biological counts or pairing fingerprints differ across methods"))
  }
  valid <- primary[!.simCompareRecordedComparatorErrors(primary), , drop = FALSE]
  if (!all(is.finite(valid$propRespEst)) || !.simCompareCountsConsistent(valid)) {
    return(invalid("unlabelled missing estimates or inconsistent gate counts"))
  }
  list(complete = TRUE, reason = NA_character_)
}

.simComparePrimaryOutputComplete <- function(
    .data, nSample, nIter, methods = c("stimgate", "fbeta", "tailgate")) {
  .simComparePrimaryOutputStatus(.data, nSample, nIter, methods)$complete
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

# TRUE when every StimGate row records `locThresholdMethod` equal to the
# requested method.
.simCompareCacheLocThresholdMethodOk <- function(cached, locThresholdMethod) {
  if (!all(c("approach", "locThresholdMethod") %in% names(cached))) {
    return(FALSE)
  }
  isStimgate <- cached$approach %in% "stimgate"
  if (!any(isStimgate)) {
    return(FALSE)
  }
  used <- as.character(cached$locThresholdMethod[isStimgate])
  !anyNA(used) && all(used == locThresholdMethod)
}

#' Validate scenario cached output against grid row settings
#'
#' @keywords internal
.simCompareValidateScenarioCache <- function(
  cached,
  row,
  nSample = NULL,
  nIter = NULL,
  retryErrors = FALSE,
  locThresholdMethod = NULL
) {
  if (!is.data.frame(cached) || nrow(cached) == 0L) {
    return(FALSE)
  }

  # Outputs from before the threshold method was recorded, or made with
  # another method, are not reused.
  if (!is.null(locThresholdMethod) && !.simCompareCacheLocThresholdMethodOk(
    cached, locThresholdMethod
  )) {
    return(FALSE)
  }

  has_error <- "error" %in%
    names(cached) &&
    any(!is.na(cached$error) & nzchar(as.character(cached$error)))
  # retryErrors retries scenario/infrastructure failures. A fully recorded
  # comparator failure is valid missing scientific data and is reused.
  if (has_error && any(
    !is.na(cached$error) & nzchar(as.character(cached$error)) &
      !.simCompareRecordedComparatorErrors(cached)
  )) return(FALSE)

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
      !.simComparePrimaryOutputComplete(
        cached,
        nSample = nSample,
        nIter = nIter
      )
  ) {
    return(FALSE)
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
  locThresholdMethod = "region",
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
        retryErrors = retryErrors,
        locThresholdMethod = locThresholdMethod
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
    identical(row$mismatch_type[[1]], "sd_inflation")) {
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
        bwScope = if ("bw_scope" %in% names(row)) row$bw_scope[[1]] else "cytokine",
        bwNcellMax = if ("bw_ncell_max" %in% names(row)) {
          row$bw_ncell_max[[1]]
        } else {
          1e4
        },
        bwNcellMin = if ("bw_ncell_min" %in% names(row)) {
          row$bw_ncell_min[[1]]
        } else if ("bw_ncell_max" %in% names(row)) {
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
        locThresholdMethod = locThresholdMethod,
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
          locThresholdMethod = NA_character_,
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
  locThresholdMethod = "region",
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
      locThresholdMethod = locThresholdMethod,
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
  failure_reasons <- character()

  for (sim_id in intersect(expected_ids, observed_ids)) {
    sim_data <- .data[as.integer(.data$sim_id) == sim_id, , drop = FALSE]
    status <- .simComparePrimaryOutputStatus(sim_data, nSample, nIter)
    complete <- status$complete

    if (isTRUE(complete)) {
      completed_ids <- c(completed_ids, sim_id)
    } else {
      failed_ids <- c(failed_ids, sim_id)
      failure_reasons[as.character(sim_id)] <- status$reason
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
    failure_reasons = failure_reasons,
    missing_ids = sort(missing_ids),
    extra_ids = sort(extra_ids),
    collate_ok = collate_ok,
    validation_ok = validation_ok
  )
}

# Counts aggregate across biological scenarios without pretending they share a
# bootstrap family. This validation table contains no performance intervals.
.simCompareMethodOutcomeCounts <- function(raw) {
  raw <- .simComparePrimaryMethodRows(raw)
  if (!"error" %in% names(raw)) raw$error <- NA_character_
  if (!"thresholdOrigin" %in% names(raw)) raw$thresholdOrigin <- NA_character_
  raw |>
    dplyr::filter(.data$method %in% c("stimgate", "fbeta", "tailgate")) |>
    .simCompareClassificationMetrics() |>
    dplyr::group_by(.data$method) |>
    dplyr::summarise(n = dplyr::n(),
      n_valid = sum(.data$gate_status != "failed"),
      n_run_error = sum(!is.na(.data$error) & nzchar(.data$error)),
      n_no_cutpoint = sum(.data$thresholdOrigin %in% "failed_no_cutpoint"),
      n_fallback = sum(.data$gate_status %in% c("fallback_empty", "fallback_selected")),
      .groups = "drop")
}

# sim_seed identifies biological draws shared across methods and deterministic
# mismatches. Fallback biological keys support small fixtures without grid seeds.
.simCompareBootstrapContext <- function(data, unit = "iter") {
  .simCompareRequireUnit(data, unit)
  if ("sim_seed" %in% names(data)) {
    if (any(!is.finite(data$sim_seed))) stop("Bootstrap requires finite biological sim_seed values.")
    family <- paste0("sim_seed:", data$sim_seed)
  } else {
    keys <- intersect(c("base_scenario_id", "transformation", "mean_pos_setting", "mean_pos",
      "prob_response", "n_cell", "sample_perturbation_sd", "condition_perturbation_sd",
      "cluster_perturbation_sd", "background_relative_to_response", "n_cell_uns_relative_to_stim"), names(data))
    family <- if (length(keys)) do.call(paste, c(lapply(data[keys], as.character), sep = "|")) else rep("fixture", nrow(data))
  }
  units <- lapply(split(data[[unit]], family), function(ids) {
    ids <- ids[!is.na(ids)]
    if (is.numeric(ids) && length(ids) && all(ids == as.integer(ids) & ids > 0)) {
      seq_len(max(ids))
    } else sort(unique(ids))
  })
  data$.bootstrap_family <- family
  data$.bootstrap_units <- unname(units[family])
  data
}

.simComparePooledStats <- function(data, scenarioCols, spec, unit, mcse) {
  data <- .simCompareBootstrapContext(data, unit)
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::group_modify(function(rows, key) {
      families <- unique(rows$.bootstrap_family)
      if (length(families) != 1L) stop("A pooled scenario must have one biological bootstrap family.")
      intended <- rows$.bootstrap_units[[1]]
      missing <- setdiff(as.character(intended), as.character(rows[[unit]]))
      units <- c(as.character(rows[[unit]]), missing)
      out <- dplyr::bind_cols(purrr::imap(spec, function(item, name) {
        x <- c(rows[[item$column]], rep(NA_real_, length(missing)))
        .analysis_mcse_pooled_cols(x, units, item$stat, name, families, mcse)
      }))
      out$n_dataset_total <- length(intended)
      out
    }) |>
    dplyr::ungroup()
}

#' Summarise comparison runs by scenario and method
#'
#' @keywords internal
# With `mcse`, the plotted statistics also get `<stat>_mcse`, `<stat>_lower`
# and `<stat>_upper` (`analysis-mcse.R`): samples within one simulated
# dataset are not independent (bandwidth and bias are estimated across the
# dataset's samples and gates can be clustered), so the MCSE comes from the
# dataset-block bootstrap of the plotted pooled statistic. At least five
# contributing datasets and 95% finite bootstrap statistics are required.
# Point estimates are identical with intervals on or off.
.simCompareSummariseFreqBs <- function(
  .data,
  scenarioCols = NULL,
  keepMethods = c("stimgate", "fbeta", "tailgate"),
  mcse = FALSE,
  unit = "iter"
) {
  .data <- .simComparePrimaryMethodRows(.data)
  if (!"error" %in% names(.data)) {
    .data$error <- NA_character_
  }

  if (!"thresholdOrigin" %in% names(.data)) {
    .data$thresholdOrigin <- NA_character_
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

  scored <- .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    dplyr::mutate(
      run_error = .data$error,
      dplyr::across(
        dplyr::all_of(c("propRespEst", "propStim", "propUns", "threshold")),
        ~ dplyr::if_else(
          !is.na(.data$run_error) & nzchar(.data$run_error), NA_real_, .x
        )
      ),
      freq_error = .data$propRespEst - .data$propRespTruth,
      abs_error = abs(.data$freq_error),
      sq_error = .data$freq_error^2,
      rel_error = dplyr::if_else(
        .data$propRespTruth != 0,
        .data$freq_error / .data$propRespTruth,
        NA_real_
      )
    )
  out <- scored |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_est = sum(is.finite(.data$propRespEst)),
      n_run_error = sum(!is.na(.data$run_error) & .data$run_error != ""),
      n_no_cutpoint = sum(.data$thresholdOrigin %in% "failed_no_cutpoint"),
      n_threshold_fallback = sum(
        .data$thresholdFallbackUsed %in% TRUE &
          (is.na(.data$run_error) | !nzchar(.data$run_error)),
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
  q <- .simCompareQuantileFinite
  fmean <- function(v) { v <- v[is.finite(v)]; if (length(v)) mean(v) else NA_real_ }
  item <- function(column, stat) list(column = column, stat = stat)
  scored$fallback <- scored$thresholdFallbackUsed %in% TRUE &
    (is.na(scored$run_error) | !nzchar(scored$run_error))
  spec <- list(
    mean_abs_error = item("abs_error", fmean),
    med_abs_rel_error = item("rel_error", function(v) q(abs(v), 0.5)),
    q90_abs_rel_error = item("rel_error", function(v) q(abs(v), 0.9)),
    q95_abs_rel_error = item("rel_error", function(v) q(abs(v), 0.95)),
    fallback_rate = item("fallback", fmean),
    propUns_mean = item("propUns", fmean)
  )
  bootstrap <- .simComparePooledStats(scored, scenarioCols, spec, unit, mcse)
  out <- dplyr::select(out, -dplyr::any_of(c(names(spec), "max_abs_rel_error")))
  dplyr::left_join(out, bootstrap, by = scenarioCols)
}

# Stop unless `.data` has the dataset column used as the Monte Carlo unit.
.simCompareRequireUnit <- function(.data, unit) {
  if (!is.character(unit) || length(unit) != 1L || !unit %in% names(.data)) {
    stop(
      "Monte Carlo errors need the dataset column `", paste(unit, collapse = ""),
      "` in the comparison results."
    )
  }
  invisible(TRUE)
}

# Validate cross-setting invariants on the complete table before promotion.
.simCompareValidateMismatch <- function(compare_raw) {
  counts_ok <- .simCompareCountsConsistent(
    compare_raw[compare_raw$method %in% c("stimgate", "fbeta", "tailgate") &
      !.simCompareRecordedComparatorErrors(compare_raw), ]
  )
  if (!isTRUE(counts_ok)) {
    stop("Confusion-matrix counts do not reproduce the gated stimulated counts.")
  }

  pairing_check <- .simComparePairingCheck(compare_raw)
  if (!all(pairing_check$paired)) {
    stop(
      "Simulated data differ across mismatch settings within a baseline ",
      "scenario, replicate and sample, so the settings are not paired."
    )
  }

  zero_agreement <- .simCompareZeroMismatchAgreement(compare_raw)
  zero_shift <- zero_agreement[
    zero_agreement$mismatch_type == "mean_shift_negative",
  ]
  if (
    nrow(zero_shift) == 0L ||
      any(zero_shift$n_same_threshold != zero_shift$n_compared) ||
      any(zero_shift$n_same_counts != zero_shift$n_compared)
  ) {
    stop(
      "With no mismatch, shifting only the stimulated negatives did not ",
      "reproduce shifting all stimulated cells."
    )
  }
  invisible(list(pairing_check = pairing_check, zero_agreement = zero_agreement))
}

# Paired comparisons use the simulated dataset (iter), not its dependent tubes,
# as the independent unit. Scenario columns must include every varying setting.
# Method summaries use finite tubes; dataset-pair coverage stays explicit.
.simCompareDatasetDifferences <- function(
    .data, scenarioCols,
    outcomes = c("abs_error", "abs_rel_error"),
    competitors = c("fbeta", "tailgate")) {
  .data <- .simComparePrimaryMethodRows(.data)
  scenarioCols <- setdiff(
    scenarioCols, c("method", "approach", "sim_id", "sim_seed", "iter", "sample", "ind")
  )
  allowed <- c("abs_error", "abs_rel_error", "fdp", "sensitivity", "f1")
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
  if (any(outcomes %in% c("fdp", "sensitivity", "f1"))) {
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
    nIter,
    validate_full = NULL) {
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
      "exactly the complete simulation grid. ",
      paste(names(full_check$failure_reasons), full_check$failure_reasons, collapse = "; ")
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

  if (!is.null(validate_full)) {
    tryCatch(validate_full(compare_raw_full), error = function(e) {
      .analysis_mark_chunk(
        run_ctx = run_ctx, total_sims = total_sims,
        completed_sims = completed_sims, failed_sims = failed_sims,
        collate_ok = TRUE, validation_ok = FALSE,
        error_message = conditionMessage(e)
      )
      stop(e)
    })
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
    allow_tall = FALSE,
    ratio_twins = FALSE,
    mcse_mode = NULL,
    after_plot = NULL,
    fit_panels = FALSE) {
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
        plot_height <- if (isTRUE(fit_panels)) .analysis_facet_height(p) else height
        .analysis_print_save_fig(
          p, file.path(dir, set_name, file_fn(pos, extra)),
          height = plot_height, allow_tall = allow_tall, mcse_mode = mcse_mode,
          fit_panels = fit_panels
        )
        if (!is.null(after_plot)) after_plot(curr, p)
        if (isTRUE(ratio_twins)) {
          .simBandwidthPrintRatioTwin(p, file.path(dir, set_name, file_fn(pos, extra)),
            height = plot_height, allow_tall = allow_tall || fit_panels, mcse_mode = mcse_mode)
        }
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
# `mcse`: draw the `<y>_lower`/`<y>_upper` Monte Carlo intervals.
.simComparePlotByCell <- function(
    data, y, y_label, zero_line = FALSE, free_y = FALSE, percent_y = FALSE,
    mcse = FALSE) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = n_cell, y = .data[[y]], colour = method, shape = method, linetype = method)
  )
  if (zero_line) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0, colour = "gray25", linetype = "dashed"
    )
  }
  bounds <- paste0(y, c("_lower", "_upper"))
  p +
    ggplot2::geom_line(alpha = 0.75) +
    ggplot2::geom_point(alpha = 0.75) +
    (if (isTRUE(mcse) && all(bounds %in% names(data))) {
      .analysis_mcse_errorbar(data, bounds[[1]], bounds[[2]])
    }) +
    ggplot2::scale_x_log10(
      breaks = sort(unique(data$n_cell)),
      labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    ) +
    ggplot2::scale_y_continuous(
      labels = if (percent_y) .analysis_label_percent else .analysis_label_number
    ) +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    ggplot2::facet_wrap(
      ggplot2::vars(prob_response, transformation),
      ncol = length(unique(data$transformation)),
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
.simComparePlotEstVsTruth <- function(data, maxwidth = 0.2, lower_limit = NULL) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  response_levels <- sort(unique(data$prob_response))
  methods <- .simCompareMethodSets()$all_methods$methods
  methods <- methods[methods %in% data$method]
  offsets <- if (length(methods) == 1L) 0 else seq(-0.25, 0.25, length.out = length(methods))
  spacing <- if (length(methods) == 1L) 0.5 else min(diff(offsets))
  if (length(maxwidth) != 1L || !is.finite(maxwidth) || maxwidth <= 0 || maxwidth >= spacing) {
    stop("maxwidth must be positive and smaller than the method spacing (", spacing, ").")
  }
  if (is.null(lower_limit)) {
    positive <- c(data$propRespEst, response_levels)
    lower_limit <- min(positive[is.finite(positive) & positive > 0]) / 10
  }
  if (length(lower_limit) != 1L || !is.finite(lower_limit) || lower_limit <= 0) {
    stop("lower_limit must be one finite positive frequency.")
  }
  data$response_position <- match(data$prob_response, response_levels)
  data$method_position <- data$response_position + offsets[match(data$method, methods)]
  data$estimate_shown <- pmax(data$propRespEst, lower_limit)
  data$plot_floor <- lower_limit
  truth <- data |>
    dplyr::distinct(.data$transformation, .data$prob_response, .data$response_position)
  ggplot2::ggplot(data,
    ggplot2::aes(x = method_position, y = estimate_shown, colour = method, shape = method,
      group = interaction(prob_response, method))) +
    ggforce::geom_sina(maxwidth = maxwidth, scale = "width", orientation = "x",
      position = "identity", seed = 271L, jitter_y = FALSE, alpha = 0.35, size = 1) +
    # Truth segments go on top so dense groups cannot hide them.
    ggplot2::geom_segment(data = truth,
      ggplot2::aes(x = response_position - 0.4, xend = response_position + 0.4,
        y = prob_response, yend = prob_response), inherit.aes = FALSE,
      colour = "black", linetype = "dashed", linewidth = 0.6) +
    ggplot2::scale_x_continuous(breaks = seq_along(response_levels),
      labels = .analysis_label_percent(response_levels)) +
    ggplot2::scale_y_log10(labels = function(x) {
        lab <- .analysis_label_percent(x)
        floor_at <- !is.na(x) & abs(x / lower_limit - 1) < 1e-8
        lab[floor_at] <- paste0("\u2264 ", lab[floor_at])
        lab
      },
      breaks = function(limits) sort(unique(c(lower_limit, scales::log_breaks()(limits)))),
      limits = c(lower_limit, NA), oob = scales::squish) +
    .analysis_scale_method(c("colour", "shape")) +
    ggplot2::facet_wrap(~transformation) +
    ggplot2::labs(
      x = "Simulated response frequency (evenly spaced levels)",
      y = "Estimated background-subtracted response frequency",
      colour = "Method",
      caption = paste0("Estimates at or below ", .analysis_label_percent(lower_limit),
        " (including zero and negative estimates) are squished to that lower limit; ",
        sum(is.finite(data$propRespEst) & data$propRespEst <= lower_limit),
        " method/sample estimates shown at the limit. Dashed segments mark true response frequencies.")
    ) + .analysis_theme()
}

# Share of estimates that used a threshold fallback.
.simComparePlotFallback <- function(data, mcse = FALSE) {
  .simComparePlotByCell(
    data, "fallback_rate", "Share of estimates using a threshold fallback",
    percent_y = TRUE, mcse = mcse
  )
}

# Histograms of finite thresholds by method; `data` has `approach`.
# Reference expression densities are rescaled within each panel so their peak
# matches the tallest threshold-histogram bar: they show shape only, not
# density values. A dotted line marks the stimulated positive-component mean,
# which is often too rare to see in the density.
.simComparePlotThresholdDensity <- function(data, densities = NULL) {
  data$transformation <- .analysis_trans_factor(data$transformation)
  facet_vars <- c("prob_response", "transformation")
  build_plot <- function(densities) {
    p <- ggplot2::ggplot(data,
      ggplot2::aes(x = threshold, y = ggplot2::after_stat(density),
        fill = approach, colour = approach))
    if (!is.null(densities)) {
      p <- p + ggplot2::geom_line(data = densities,
        ggplot2::aes(x = expression, y = reference_y,
          linetype = condition, group = condition), inherit.aes = FALSE,
        colour = "gray25", linewidth = 0.5) +
        ggplot2::scale_linetype_manual(values = c(unstimulated = "dashed", stimulated = "solid"),
          labels = c(unstimulated = "Unstimulated (all cells)", stimulated = "Stimulated (all cells)"))
      if ("positive_mean" %in% names(densities)) {
        positive_means <- dplyr::distinct(densities,
          dplyr::across(dplyr::all_of(facet_vars)), .data$positive_mean)
        p <- p + ggplot2::geom_vline(data = positive_means,
          ggplot2::aes(xintercept = positive_mean), inherit.aes = FALSE,
          colour = "gray10", linetype = "dotted", linewidth = 0.6) +
          ggplot2::labs(caption = paste("Dotted line: mean of the stimulated reference",
            "tube's positive-component (responding) cells."))
      }
    }
    p + ggplot2::geom_histogram(alpha = 0.2, position = "identity", bins = 30) +
      ggplot2::facet_wrap(~ prob_response + transformation, scales = "free", ncol = 3,
        labeller = ggplot2::labeller(prob_response = .analysis_labeller_percent())) +
      ggplot2::scale_x_continuous(labels = .analysis_label_number) +
      ggplot2::scale_y_continuous(labels = .analysis_label_number) +
      .analysis_scale_method(aesthetics = c("colour", "fill")) +
      ggplot2::labs(x = "Threshold / marker expression", y = "Threshold density",
        colour = "Method", fill = "Method", linetype = "Reference tube (density shape, peak matched)") +
      .analysis_theme()
  }
  if (is.null(densities)) {
    return(build_plot(NULL))
  }
  keys <- intersect(c("transformation", "prob_response", "mean_pos",
    "sample_perturbation_sd", "condition_perturbation_sd",
    "cluster_perturbation_sd", "background_relative_to_response"), names(data))
  densities$transformation <- .analysis_trans_factor(densities$transformation)
  densities <- dplyr::semi_join(densities, dplyr::distinct(data, dplyr::across(dplyr::all_of(keys))), by = keys)
  densities$reference_y <- pmax(densities$density, 0)
  # Histogram bins depend only on each panel's x range, which the reference
  # lines also train, so one build gives the bar heights of the final plot.
  plot <- build_plot(densities)
  built <- ggplot2::ggplot_build(plot)
  hist_data <- built$data[[length(plot$layers)]]
  panel_key <- function(x) do.call(paste, c(lapply(x[facet_vars], as.character), sep = "\r"))
  layout <- built$layout$layout
  hist_max <- tapply(hist_data$y, hist_data$PANEL, max, na.rm = TRUE)
  layout$hist_max <- as.numeric(hist_max[as.character(layout$PANEL)])
  ref_max <- tapply(densities$reference_y, panel_key(densities), max, na.rm = TRUE)
  key <- panel_key(densities)
  scale <- layout$hist_max[match(key, panel_key(layout))] / as.numeric(ref_max[key])
  scale[!is.finite(scale) | scale <= 0] <- 1
  densities$reference_y <- densities$reference_y * scale
  build_plot(densities)
}

# Shared pieces of the analysis 8 mismatch plots.
.simCompareMismatchScales <- function(x_label, aesthetics = "colour") {
  list(
    ggplot2::scale_x_continuous(labels = .analysis_label_number),
    .analysis_scale_method(aesthetics),
    ggplot2::labs(x = x_label, colour = "Method"),
    .analysis_theme()
  )
}

.simCompareStripWrap <- function() ggplot2::label_wrap_gen(width = 22)

# Share of unstimulated cells above the gate against the mean shift.
.simComparePlotPlacement <- function(data, mcse = FALSE) {
  bars <- if (isTRUE(mcse) && "propUns_mean_lower" %in% names(data)) {
    .analysis_mcse_errorbar(data, "propUns_mean_lower", "propUns_mean_upper")
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = mismatch_val, y = propUns_mean, colour = method, shape = method,
      linetype = method,
      group = interaction(method, mismatch_type)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 2, alpha = 0.75) +
    bars +
    ggplot2::facet_wrap(
      ~scenario_desc, scales = "free_y", ncol = 2, labeller = .simCompareStripWrap()
    ) +
    ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
    .analysis_y_floor() +
    .simCompareMismatchScales("Additive mean shift (transformed expression scale)", c("colour", "shape", "linetype")) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::labs(
      y = "Unstimulated cells above the gate",
      linetype = "Method"
    )
}

# 90th percentile error for every mismatch mechanism.
.simComparePlotUpperTail <- function(data, mcse = FALSE) {
  data$mismatch_type <- factor(
    data$mismatch_type, levels = names(.simCompareMismatchLabels)
  )
  bars <- if (isTRUE(mcse) && "q90_abs_rel_error_lower" %in% names(data)) {
    .analysis_mcse_errorbar(
      data, "q90_abs_rel_error_lower", "q90_abs_rel_error_upper"
    )
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = mismatch_val, y = q90_abs_rel_error, colour = method, shape = method, linetype = method, group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 2, alpha = 0.75) +
    bars +
    ggplot2::facet_wrap(~scenario_desc, scales = "free", ncol = 2,
      labeller = .simCompareStripWrap()) +
    ggplot2::scale_y_continuous(
      transform = scales::asinh_trans(), labels = .analysis_label_percent,
      breaks = scales::breaks_pretty(n = 4)
    ) +
    .analysis_y_floor() +
    .simCompareMismatchScales("Mismatch size", c("colour", "shape", "linetype")) +
    ggplot2::labs(y = "Tube-level 90th percentile absolute relative error (asinh scale)")
}

# Over- and under-estimate summary of signed relative error per scenario and
# method, from raw comparison rows, with the same rows and estimand as
# `.simCompareSummariseFreqBs()`. With `mcse`, the Monte Carlo errors come
# from whole-dataset (`unit`) resampling of those pooled statistics.
.simCompareSignedErrorSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate"),
    mcse = FALSE,
    unit = "iter") {
  .data <- .simComparePrimaryMethodRows(.data)
  rows <- .simCompareBootstrapContext(.data, unit) |>
    dplyr::filter(.data$method %in% keepMethods)
  if (!"error" %in% names(rows)) rows$error <- NA_character_
  rows <- rows |> dplyr::mutate(rel_error = dplyr::if_else(
    .data$propRespTruth != 0 & (is.na(.data$error) | !nzchar(.data$error)),
    (.data$propRespEst - .data$propRespTruth) / .data$propRespTruth, NA_real_))
  rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::group_modify(function(data, key) {
      family <- unique(data$.bootstrap_family)
      if (length(family) != 1L) stop("A signed scenario must have one bootstrap family.")
      missing <- setdiff(as.character(data$.bootstrap_units[[1]]), as.character(data[[unit]]))
      .simBandwidthSignedErrorSides(
        c(data$rel_error, rep(NA_real_, length(missing))), mcse = mcse,
        unit = c(as.character(data[[unit]]), missing), bootstrap_family = family
      )
    }) |>
    dplyr::ungroup()
}

# Size of the relative error, |estimate - truth| / truth, per scenario and
# method: median and 95th percentile over pooled sample-level errors.
# Uses the same rows and estimand as `.simCompareSummariseFreqBs()`; samples
# with a true frequency of zero have no relative error and are left out.
# Bounds use whole-dataset resampling of the same pooled percentiles. Maxima
# belong to `.simCompareDatasetMaxSummary()`, not these main figures.
.simCompareUnsignedErrorSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate"),
    mcse = FALSE,
    unit = "iter") {
  .data <- .simComparePrimaryMethodRows(.data)
  rows <- .data |> dplyr::filter(.data$method %in% keepMethods)
  if (!"error" %in% names(rows)) rows$error <- NA_character_
  rows <- rows |> dplyr::mutate(abs_rel_error = dplyr::if_else(
    .data$propRespTruth != 0 & (is.na(.data$error) | !nzchar(.data$error)),
    abs((.data$propRespEst - .data$propRespTruth) / .data$propRespTruth), NA_real_))
  q <- .simCompareQuantileFinite
  spec <- list(
    median = list(column = "abs_rel_error", stat = function(v) q(v, 0.5)),
    q95 = list(column = "abs_rel_error", stat = function(v) q(v, 0.95))
  )
  .simComparePooledStats(rows, scenarioCols, spec, unit, mcse)
}

# Bootstrap the complete plotted scenario average, retaining CRN covariance.
.simComparePerformanceAverage <- function(
    raw, scenarioCols, group_cols, signed = FALSE, mcse = FALSE) {
  summary <- if (isTRUE(signed)) .simCompareSignedErrorSummary(raw, scenarioCols, mcse = mcse) else
    .simCompareUnsignedErrorSummary(raw, scenarioCols, mcse = mcse)
  stats <- if (isTRUE(signed)) c("prop", "median", "q90", "q95") else c("median", "q95")
  groups <- c(group_cols, if (isTRUE(signed)) "direction")
  out <- .analysis_mcse_bootstrap_average(summary, groups, stats, mcse)
  if (isTRUE(signed)) .simBandwidthSignedErrorClipSide(out) else out
}

# Fixed-size maxima use complete datasets only. All intended dataset IDs remain
# in the joint bootstrap, including incomplete and unaffected datasets.
.simCompareDatasetMaxSummary <- function(
    raw, scenarioCols, expected_samples = 20L, expected_datasets = NULL,
    mcse = FALSE) {
  raw <- .simComparePrimaryMethodRows(raw)
  if (length(expected_samples) != 1L || !is.finite(expected_samples) ||
      expected_samples < 1L || expected_samples != as.integer(expected_samples)) {
    stop("Dataset maxima require a positive integer expected sample count.")
  }
  required <- c(scenarioCols, "iter", "sample", "method", "propRespTruth", "propRespEst")
  if (!all(required %in% names(raw))) stop("Dataset maxima need scenario keys, iter, sample and frequency outcomes.")
  rows <- .simCompareBootstrapContext(raw) |> dplyr::filter(.data$method %in% c("stimgate", "fbeta", "tailgate"))
  if (!"error" %in% names(rows)) rows$error <- NA_character_
  rows <- rows |> dplyr::mutate(rel_error = dplyr::if_else(
    .data$propRespTruth != 0 & (is.na(.data$error) | !nzchar(.data$error)),
    (.data$propRespEst - .data$propRespTruth) / .data$propRespTruth, NA_real_))
  out <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::group_modify(function(data, key) {
      families <- unique(data$.bootstrap_family)
      if (length(families) != 1L) stop("A maximum scenario must have one biological bootstrap family.")
      intended <- expected_datasets
      if (is.null(intended)) intended <- data$.bootstrap_units[[1]]
      if (length(intended) == 1L && is.numeric(intended) && intended > 1L) intended <- seq_len(intended)
      intended <- sort(unique(as.character(intended)))
      if (!length(intended) || anyNA(intended)) stop("Dataset maxima need complete intended dataset IDs.")
      if (anyNA(data$iter) || any(!as.character(data$iter) %in% intended)) stop("Dataset IDs fall outside the intended maximum cohort.")
      datasets <- purrr::map_dfr(intended, function(id) {
        block <- data[as.character(data$iter) == id, , drop = FALSE]
        eligible <- nrow(block) == expected_samples && !anyNA(block$sample) &&
          !anyDuplicated(block$sample) &&
          setequal(as.character(block$sample), as.character(seq_len(expected_samples))) &&
          all(is.finite(block$rel_error))
        tibble::tibble(id = id, eligible = eligible,
          over = if (eligible) max(c(0, block$rel_error)) else NA_real_,
          under = if (eligible) max(c(0, -block$rel_error)) else NA_real_)
      })
      purrr::map_dfr(c("over", "under"), function(direction) {
        magnitude <- datasets[[direction]]
        affected <- datasets$eligible & is.finite(magnitude) & magnitude > 0
        n_eligible <- sum(datasets$eligible)
        occurrence_fn <- function(index) {
          eligible <- datasets$eligible[index]
          if (any(eligible)) mean(affected[index][eligible]) else NA_real_
        }
        severity_fn <- function(index) {
          selected <- affected[index]
          if (any(selected)) mean(magnitude[index][selected]) else NA_real_
        }
        index <- seq_len(nrow(datasets))
        occurrence_draws <- if (isTRUE(mcse)) .analysis_mcse_block_draws(index, datasets$id, occurrence_fn, families) else numeric()
        severity_draws <- if (isTRUE(mcse)) .analysis_mcse_block_draws(index, datasets$id, severity_fn, families) else numeric()
        result <- dplyr::bind_cols(
          tibble::tibble(direction = direction, n_dataset_total = nrow(datasets),
            n_eligible = n_eligible, n_incomplete = nrow(datasets) - n_eligible,
            n_affected = sum(affected), eligible_dataset_ids = paste(datasets$id[datasets$eligible], collapse = ",")),
          .analysis_mcse_bootstrap_cols(occurrence_fn(index), occurrence_draws,
            "occurrence", nrow(datasets), n_eligible, mcse),
          .analysis_mcse_bootstrap_cols(severity_fn(index), severity_draws,
            "severity", nrow(datasets), sum(affected), mcse)
        )
        result
      })
    }) |>
    dplyr::ungroup()
  if ("method" %in% names(out)) {
    keys <- setdiff(c(scenarioCols, "direction"), c("method", "approach"))
    reference <- out |> dplyr::filter(.data$method == "stimgate") |>
      dplyr::select(dplyr::all_of(keys), reference_ids = eligible_dataset_ids)
    out <- dplyr::left_join(out, reference, by = keys) |>
      dplyr::mutate(cohort_matches_stimgate = .data$eligible_dataset_ids == .data$reference_ids) |>
      dplyr::select(-"reference_ids")
  }
  out
}

# Average the scenario statistics equally over the scenarios in each group
# (each scenario counts once, however many samples it has). Scenarios are
# averaged only over the columns left out of `group_cols`.
# Dataset-bootstrap summaries average aligned bootstrap draws, preserving
# covariance between scenarios with common biological seeds. Legacy independent
# sample summaries retain their independent-scenario MCSE calculation.
.simCompareErrorAverage <- function(tbl, group_cols) {
  if (".boot_median" %in% names(tbl)) {
    return(.analysis_mcse_bootstrap_average(tbl, group_cols, c("median", "q95"),
      mcse = any(lengths(tbl$.boot_median) > 0L)))
  }
  if (any(c("median_mcse", "q95_mcse", "max_mcse") %in% names(tbl))) {
    out <- .analysis_mcse_average_cols(tbl, group_cols, c("median", "q95", "max"))
    for (s in c("median", "q95", "max")) {
      out[[paste0(s, "_lower")]] <- pmax(out[[paste0(s, "_lower")]], 0)
    }
    return(out)
  }
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_scenario = dplyr::n(),
      dplyr::across(
        dplyr::all_of(c("median", "q95", "max")),
        ~ sum(is.finite(.x)), .names = "n_scenario_{.col}"
      ),
      dplyr::across(
        dplyr::all_of(c("median", "q95", "max")),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

# Panels for mismatch plots, with independent ranges for each statistic,
# transformation and (when `by_prob`) response probability.
.simCompareMismatchFacet <- function(by_prob, ncol = 3L) {
  ggplot2::facet_wrap(
    if (by_prob) ggplot2::vars(statistic, transformation, prob_response) else
      ggplot2::vars(statistic, transformation),
    ncol = ncol,
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
    stat_cols = c(median = "Median", q95 = "95th percentile"),
    mcse = FALSE) {
  tbl$transformation <- .analysis_trans_factor(tbl$transformation)
  tbl <- .simBandwidthErrorStatLong(tbl, stat_cols)
  # Errors above 1500% (16 times the truth) are drawn at the cap.
  capped <- .simBandwidthSignedErrorIsCapped(tbl$value)
  tbl$value_shown <- .simBandwidthSignedErrorSquish(tbl$value)
  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = mismatch_val, y = value_shown, colour = method, shape = method, linetype = method, group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.75) +
    ggplot2::geom_point(size = 1.5, alpha = 0.75) +
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl)) +
    ggplot2::scale_x_continuous(labels = .analysis_label_number) +
    .simBandwidthAbsErrorLayers(capped = capped) +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    .simCompareMismatchFacet(by_prob, length(unique(tbl$transformation)) *
      if (by_prob) length(unique(tbl$prob_response)) else 1L) +
    ggplot2::labs(x = x_label, colour = "Method") +
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
    stat_cols = c(median = "Median", q95 = "95th percentile"),
    mcse = FALSE) {
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
      x = .data[[x]], y = value_shown, colour = method, shape = method, linetype = method,
      group = interaction(method, direction)
    )
  ) +
    .simBandwidthSignedErrorLayers("Tube-level relative error", capped = capped) +
    .simBandwidthSignedErrorSegmentLayer(
      tbl, x, "value_shown", line_cols, alpha = 0.75
    ) +
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    (if (isTRUE(mcse)) .simBandwidthSignedErrorBars(tbl)) +
    x_scale +
    .analysis_scale_method(c("colour", "shape", "linetype")) +
    .simCompareMismatchFacet(by_prob, length(unique(tbl$transformation)) *
      if (by_prob) length(unique(tbl$prob_response)) else 1L) +
    ggplot2::labs(x = x_label, colour = "Method") +
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
# With `mcse`, each plotted percentile gets `<stat>_mcse`, `<stat>_lower` and
# `<stat>_upper`: the percentile is also computed within each dataset
# (`unit`) is resampled as a block and the pooled percentile recomputed.
# At least five contributing datasets and 95% finite draws are required.
# Percentile-bootstrap bounds stay within 0-100%; the pooled points are unchanged.
.simCompareClassificationSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate"),
    mcse = FALSE,
    unit = "iter") {
  .data <- .simComparePrimaryMethodRows(.data)
  if (isTRUE(mcse)) {
    .simCompareRequireUnit(.data, unit)
  }
  q <- .simCompareQuantileFinite
  if (!"error" %in% names(.data)) {
    .data$error <- NA_character_
  }
  if (!"thresholdOrigin" %in% names(.data)) {
    .data$thresholdOrigin <- NA_character_
  }
  out <- .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    .simCompareClassificationMetrics() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_valid = sum(.data$gate_status != "failed"),
      n_failed = sum(.data$gate_status == "failed"),
      n_run_error = sum(!is.na(.data$error) & nzchar(.data$error)),
      n_no_cutpoint = sum(.data$thresholdOrigin %in% "failed_no_cutpoint"),
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
      f1_median = q(.data$f1, 0.5),
      f1_q10 = q(.data$f1, 0.1),
      fpr_median = q(.data$false_positive_rate, 0.5),
      fpr_q90 = q(.data$false_positive_rate, 0.9),
      selected_fraction_median = q(.data$selected_fraction, 0.5),
      prevalence_median = q(.data$n_genuine_pos / .data$n_classified, 0.5),
      .groups = "drop"
    )
  metrics <- .data |>
    dplyr::filter(.data$method %in% keepMethods) |>
    .simCompareClassificationMetrics()
  item <- function(column, p) list(column = column, stat = function(v) q(v, p))
  spec <- list(
    fdp_median = item("fdp", 0.5), fdp_q90 = item("fdp", 0.9),
    sensitivity_median = item("sensitivity", 0.5), sensitivity_q10 = item("sensitivity", 0.1),
    f1_median = item("f1", 0.5), f1_q10 = item("f1", 0.1),
    fpr_median = item("false_positive_rate", 0.5), fpr_q90 = item("false_positive_rate", 0.9)
  )
  out <- dplyr::select(out, -dplyr::any_of(names(spec)))
  dplyr::left_join(out, .simComparePooledStats(metrics, scenarioCols, spec, unit, mcse), by = scenarioCols)
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
      n_genuine_pos_values = dplyr::n_distinct(.data$nTruePos + .data$nFalseNeg, na.rm = TRUE),
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
# The shift variants add exactly zero and the counts must agree exactly.
# Thresholds may differ in the last bits when the settings ran on different
# compute nodes, so they are compared with a small relative tolerance.
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
      same_threshold = .data$threshold == .data$threshold_ref |
        abs(.data$threshold - .data$threshold_ref) <=
          1e-10 * pmax(1, abs(.data$threshold_ref)),
      same_counts = .data$nTruePos == .data$nTruePos_ref &
        .data$nFalsePos == .data$nFalsePos_ref &
        .data$nFalseNeg == .data$nFalseNeg_ref &
        .data$nTrueNeg == .data$nTrueNeg_ref
    ) |>
    dplyr::group_by(.data$mismatch_type, .data$method) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      n_failed_pairs = sum(is.na(.data$same_threshold) | is.na(.data$same_counts)),
      n_compared = sum(!is.na(.data$same_threshold) & !is.na(.data$same_counts)),
      n_same_threshold = sum(.data$same_threshold %in% TRUE),
      n_same_counts = sum(.data$same_counts %in% TRUE),
      max_abs_threshold_diff = {
        diff <- abs(.data$threshold - .data$threshold_ref)
        if (any(is.finite(diff))) max(diff[is.finite(diff)]) else NA_real_
      },
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
  f1 = c(
    label = "F1 score",
    median = "f1_median", tail = "f1_q10"
  ),
  fpr = c(
    label = "False-positive rate",
    median = "fpr_median", tail = "fpr_q90"
  )
)

# Classification outcomes against mismatch size from
# `.simCompareClassificationSummary()`. Each scenario has adjacent outcome
# panels, each with its own horizontal range; colour and shape identify the
# method and line type the statistic (median, or the worse tail: 90th
# percentile for FDP and false-positive rate, 10th for sensitivity and F1). With
# `unit_scale`, every vertical scale is fixed at 0-100%.
.simComparePlotClassification <- function(
    tbl,
    outcomes = c("fdp", "sensitivity"),
    x_label = "Mismatch size",
    unit_scale = TRUE,
    mcse = FALSE) {
  spec <- .simCompareClassificationOutcomes[outcomes]
  tail_label <- if (length(outcomes) == 1L) {
    if (outcomes %in% c("sensitivity", "f1")) "10th percentile" else "90th percentile"
  } else {
    paste(c(fdp = "90th FDP", sensitivity = "10th sensitivity", f1 = "10th F1", fpr = "90th FPR")[outcomes],
      collapse = " / ")
  }
  # One statistic's rows, with its Monte Carlo bounds when present.
  stat_rows <- function(label, statistic, col) {
    bound <- function(suffix) {
      b <- paste0(col, suffix)
      if (b %in% names(tbl)) tbl[[b]] else rep(NA_real_, nrow(tbl))
    }
    dplyr::mutate(
      tbl,
      outcome = label, statistic = statistic, value = tbl[[col]],
      lower = bound("_lower"), upper = bound("_upper")
    )
  }
  long <- purrr::map_df(names(spec), function(nm) {
    s <- spec[[nm]]
    dplyr::bind_rows(
      stat_rows(s[["label"]], "median", s[["median"]]),
      stat_rows(s[["label"]], "tail", s[["tail"]])
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
      x = mismatch_val, y = value, colour = method, shape = method, linetype = statistic,
      group = interaction(method, statistic)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8, na.rm = TRUE) +
    (if (isTRUE(mcse)) .analysis_mcse_errorbar(long)) +
    ggplot2::facet_wrap(
      ggplot2::vars(scenario, outcome), ncol = 2,
      scales = if (unit_scale) "free_x" else "free",
      labeller = ggplot2::labeller(
        scenario = ggplot2::label_wrap_gen(width = 28),
        outcome = ggplot2::label_wrap_gen(width = 24)
      )
    ) +
    ggplot2::scale_x_continuous(
      transform = "sqrt", labels = .analysis_label_number,
      guide = ggplot2::guide_axis(angle = 45)
    ) +
    .analysis_scale_method(c("colour", "shape")) +
    ggplot2::scale_linetype_manual(
      values = c(median = "solid", tail = "22"),
      labels = c(
        median = "Tube-level median",
        tail = tail_label
      )
    ) +
    # Same order for colour and shape, so they merge into one Method legend.
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2, ncol = 1)
    ) +
    ggplot2::labs(
      x = paste0(x_label, " (square-root spacing)"), y = NULL, colour = "Method", linetype = "Tube-level distribution"
    ) +
    .analysis_theme()
  if (unit_scale) {
    p + ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = .analysis_label_percent
    ) + .analysis_y_floor(c(0, 1))
  } else {
    p + ggplot2::scale_y_continuous(labels = .analysis_label_percent) +
      .analysis_y_floor()
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
      ggplot2::aes(x = expr, alpha = "Unstimulated tube"),
      binwidth = binwidth, boundary = boundary, colour = "black",
      linewidth = 0.4, linetype = "dotted"
    ) +
    ggplot2::geom_vline(
      data = gates,
      ggplot2::aes(xintercept = threshold, colour = method, linetype = method, group = line_id),
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
    .analysis_scale_method(c("colour", "linetype")) +
    ggplot2::scale_alpha_manual(values = c("Unstimulated tube" = 1), name = "Reference tube") +
    ggplot2::labs(
      x = x_label, y = "Number of cells (square-root scale)",
      fill = "Stimulated cells",
      colour = "Method", linetype = "Method"
    ) +
    .analysis_theme()
}

# Unconditional pooled tube percentiles with the same dataset context and exclusions.
.simCompareSignedPercentileSummary <- function(
    .data,
    scenarioCols,
    keepMethods = c("stimgate", "fbeta", "tailgate"),
    mcse = FALSE,
    unit = "iter") {
  .data <- .simComparePrimaryMethodRows(.data)
  rows <- .simCompareBootstrapContext(.data, unit) |>
    dplyr::filter(.data$method %in% keepMethods)
  if (!"error" %in% names(rows)) rows$error <- NA_character_
  rows <- rows |> dplyr::mutate(rel_error = dplyr::if_else(
    .data$propRespTruth != 0 & (is.na(.data$error) | !nzchar(.data$error)),
    (.data$propRespEst - .data$propRespTruth) / .data$propRespTruth, NA_real_))
  rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scenarioCols))) |>
    dplyr::group_modify(function(data, key) {
      family <- unique(data$.bootstrap_family)
      if (length(family) != 1L) stop("A signed scenario must have one bootstrap family.")
      missing <- setdiff(as.character(data$.bootstrap_units[[1]]), as.character(data[[unit]]))
      .simBandwidthSignedErrorPercentiles(
        c(data$rel_error, rep(NA_real_, length(missing))), mcse = mcse,
        unit = c(as.character(data[[unit]]), missing), bootstrap_family = family
      )
    }) |>
    dplyr::ungroup()
}

.simCompareSignedPercentileAverage <- function(raw, scenarioCols, group_cols, mcse = FALSE) {
  .analysis_mcse_bootstrap_average(
    .simCompareSignedPercentileSummary(raw, scenarioCols, mcse = mcse),
    group_cols, names(.simBandwidthSignedErrorProbs), mcse = TRUE)
}
