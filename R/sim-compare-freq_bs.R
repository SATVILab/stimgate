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

  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("Transformation not recognised: ", transformation)
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

  # fcm is only needed by plotting/tick helpers, not by get_positivity_threshold.
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

#' Locate the tailgate helper scripts if they have not already been sourced
#'
#' @keywords internal
.simCompareTailgateFindSourceFiles <- function(tailgateSourceFiles = NULL) {
  if (!is.null(tailgateSourceFiles)) {
    return(tailgateSourceFiles[file.exists(tailgateSourceFiles)])
  }

  root <- tryCatch(
    {
      if (requireNamespace("projr", quietly = TRUE)) {
        projr::projr_path_get("project")
      } else {
        getwd()
      }
    },
    error = function(e) getwd()
  )

  candidate <- list.files(
    root,
    pattern = "^(openCyto-find_peaks_and_valleys|cytoUtils-cytokine_cutpoint).*\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(candidate) == 0L) {
    return(character())
  }

  candidate[order(
    !grepl("openCyto-find_peaks_and_valleys", basename(candidate))
  )]
}

#' Source and return the real tailgate helper environment
#'
#' @keywords internal
.simCompareTailgateEnv <- function(tailgateSourceFiles = NULL) {
  if (
    exists(".cytokineCutpoint", mode = "function") &&
      exists(".findPeaks", mode = "function") &&
      exists(".findValleys", mode = "function")
  ) {
    return(environment(get(".cytokineCutpoint", mode = "function")))
  }

  tailgateSourceFiles <- .simCompareTailgateFindSourceFiles(tailgateSourceFiles)

  if (length(tailgateSourceFiles) == 0L) {
    stop(
      "Could not find the tailgate helper scripts. Either source the ",
      "openCyto/cytoUtils scripts before calling this helper, or pass ",
      "tailgateSourceFiles."
    )
  }

  cacheName <- make.names(paste(
    "tailgate",
    paste(
      normalizePath(tailgateSourceFiles, winslash = "/", mustWork = FALSE),
      collapse = "_"
    ),
    sep = "_"
  ))
  if (exists(cacheName, envir = .simCompareCacheEnv, inherits = FALSE)) {
    return(get(cacheName, envir = .simCompareCacheEnv, inherits = FALSE))
  }

  tailgateEnv <- new.env(parent = parent.frame())
  for (pathCurr in tailgateSourceFiles) {
    source(pathCurr, local = tailgateEnv)
  }

  required <- c(
    ".cytokineCutpoint",
    ".findPeaks",
    ".findValleys",
    ".derivDensity"
  )
  missing <- required[
    !vapply(
      required,
      exists,
      logical(1),
      envir = tailgateEnv,
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(missing) > 0L) {
    stop(
      "Tailgate helper scripts were sourced, but the following functions ",
      "were not found: ",
      paste(missing, collapse = ", ")
    )
  }

  assign(cacheName, tailgateEnv, envir = .simCompareCacheEnv)

  tailgateEnv
}

#' Call the cytoUtils/openCyto tailgate implementation
#'
#' If bandwidth is NULL, .cytokineCutpoint() forwards NULL to .derivDensity(),
#' whose default behaviour is to estimate the bandwidth with ks::hpi().
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

  tailgateEnv <- .simCompareTailgateEnv(tailgateSourceFiles)
  cytokineCutpoint <- get(
    ".cytokineCutpoint",
    envir = tailgateEnv,
    inherits = TRUE
  )

  threshold <- cytokineCutpoint(
    x = x,
    adjust = adjust,
    numPeaks = numPeaks,
    refPeak = refPeak,
    method = method,
    tol = tol,
    side = side,
    strict = strict,
    autoTol = autoTol,
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

#' @keywords internal
.simCompareEstimateFromThreshold <- function(
  xStim,
  xUns,
  threshold,
  fallbackHighValue = TRUE,
  fallbackMargin = 0.05
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
  tailgateX = c("unstim", "combined", "stim"),
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
    xUnsGate <- xUnsRaw + (biasUns %||% 0)

    purrr::map_df(indStimVec, function(indStim) {
      xStim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, chnl])

      fbetaObj <- tryCatch(
        .simCompareFbetaThreshold(
          xUns = xUnsGate,
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
        xUns = xUnsGate,
        threshold = fbetaObj$threshold,
        fallbackHighValue = fallbackHighValue,
        fallbackMargin = fallbackMargin
      )

      xTail <- switch(
        tailgateX,
        "unstim" = xUnsGate,
        "combined" = c(xUnsGate, xStim),
        "stim" = xStim
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
        xUns = xUnsGate,
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
  calcCytPosGates = FALSE,
  calcSinglePosGates = FALSE,
  includeLocCondition = TRUE
) {
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

      detailTbl <- .simCompareReadLocDetails(
        pathProject = pathProject,
        nSample = nSample,
        nCondition = nCondition
      )

      if (!includeLocCondition) {
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

      detailTbl |>
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
        ) |>
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
  calcSinglePosGates = FALSE,
  includeLocCondition = TRUE,
  pathFbeta = NULL,
  fbetaPatchPy2Compat = TRUE,
  fbetaBeta = 0.8,
  fbetaTheta = 2,
  fbetaWidth = 10,
  fbetaNumBins = NULL,
  tailgateX = c("unstim", "combined", "stim"),
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
  fallbackMargin = 0.05
) {
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

    flowFrameList <- outListExperiment[["flowFrameList"]]
    labelsList <- outListExperiment[["labelsList"]]
    fs <- as(flowFrameList, "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-compare-freq-bs",
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

    stimgateTbl <- .simCompareStimgateRows(
      gs = gs,
      labelsList = labelsList,
      pathProject = pathProject,
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
      calcCytPosGates = calcCytPosGates,
      calcSinglePosGates = calcSinglePosGates,
      includeLocCondition = includeLocCondition
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
        tailgateAutoTol = tailgateAutoTol
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
  includeLocCondition = TRUE,
  ...
) {
  purrr::map_df(seq_len(nrow(sim_grid)), function(i) {
    row <- sim_grid[i, , drop = FALSE]

    out <- .simCompareFreqBs(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nIter = nIter,
      biasUns = row$bias_uns[[1]],
      bw = row$bw[[1]],
      bwFallback = row$bw[[1]],
      bwMin = "none",
      bwMax = "none",
      probExact = probExact,
      nCellStim = row$n_cell[[1]],
      probResponse = row$prob_response[[1]],
      meanPos = row$mean_pos[[1]],
      transformation = row$transformation[[1]],
      samplePerturbationSd = row$sample_perturbation_sd[[1]],
      conditionPerturbationSd = row$condition_perturbation_sd[[1]],
      clusterPerturbationSd = row$cluster_perturbation_sd[[1]],
      backgroundRelativeToResponse = row$background_relative_to_response[[1]],
      ncellUnsRelativeToStim = row$n_cell_uns_relative_to_stim[[1]],
      covEvMin = covEvMin,
      covEvMax = covEvMax,
      tolClust = tolClust,
      includeLocCondition = includeLocCondition,
      ...
    )

    out <- out |>
      dplyr::select(-dplyr::any_of(names(row)))

    dplyr::bind_cols(
      row[rep(1L, nrow(out)), , drop = FALSE],
      out
    )
  })
}

#' Summarise comparison runs by scenario and method
#'
#' @keywords internal
.simCompareSummariseFreqBs <- function(
  .data,
  scenarioCols = NULL,
  keepMethods = c("stimgate_loc_sample", "fbeta", "tailgate")
) {
  if (!"error" %in% names(.data)) {
    .data$error <- NA_character_
  }

  if (is.null(scenarioCols)) {
    scenarioCols <- intersect(
      c(
        "transformation",
        "prob_response",
        "n_cell",
        "mean_pos",
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
      bias = mean(.data$freq_error, na.rm = TRUE),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      median_abs_error = stats::median(.data$abs_error, na.rm = TRUE),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      mean_rel_error = mean(.data$rel_error, na.rm = TRUE),
      threshold_mean = mean(.data$threshold, na.rm = TRUE),
      threshold_median = stats::median(.data$threshold, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$n_cell,
      dplyr::desc(.data$prob_response),
      .data$transformation,
      .data$mean_pos,
      .data$method
    )
}
