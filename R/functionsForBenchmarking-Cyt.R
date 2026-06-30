#' Detect whether a transformation should receive upper-population ratio correction
#'
#' This is currently applied to the raw gamma transformation and to the skew
#' transformation.  The lower component is never rescaled by this correction;
#' only upper components are moved and rescaled so that their observed-scale
#' separation from the pre-perturbation lower reference matches the corresponding
#' raw-scale separation after sample, condition, and cluster perturbations.
#'
#' @keywords internal
.simCytUsesUpperRatioCorrection <- function(transformationFunc) {
  transName <- attr(transformationFunc, "sim_transformation")
  isTRUE(transName %in% c("gamma", "skew")) ||
    (exists("calc_gamma", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_gamma))) ||
    (exists("calc_skew", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_skew)))
}

#' Backwards-compatible alias for older notebooks/scripts
#'
#' @keywords internal
.simCytUsesGammaRatioCorrection <- function(transformationFunc) {
  .simCytUsesUpperRatioCorrection(transformationFunc)
}

#' Lower reference mean after transformation
#'
#' For gamma this is the standard gamma transform of the raw lower reference.
#' For raw skew, calc_skew() is stochastic and depends on the full vector, so
#' the transformed lower reference should be estimated from the realised lower
#' population for that condition and marker.
#'
#' @keywords internal
.simCytLowerMeanTransReference <- function(
  transformationFunc,
  lowerMeanRawReference,
  yLower = NULL
) {
  transName <- attr(transformationFunc, "sim_transformation")
  isSkew <- isTRUE(identical(transName, "skew")) ||
    (exists("calc_skew", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_skew)))

  if (isSkew) {
    yLower <- as.numeric(yLower)
    yLower <- yLower[is.finite(yLower)]
    if (length(yLower) > 0L) {
      return(mean(yLower))
    }
    return(NA_real_)
  }

  as.numeric(transformationFunc(lowerMeanRawReference))[1]
}

#' Solve the correction factor with combined lower and upper variance
#'
#' We adjust only the upper population.  If c is the mean-distance multiplier,
#' the upper SD multiplier is 1 / c.  c is chosen so that
#'
#'   c * delta / sqrt(sd_lower^2 + (sd_upper / c)^2) == target_ratio.
#'
#' @keywords internal
.simCytCombinedRatioFactor <- function(
  delta,
  sdLower,
  sdUpper,
  targetRatio,
  eps = sqrt(.Machine$double.eps)
) {
  delta <- as.numeric(delta)[1]
  sdLower <- as.numeric(sdLower)[1]
  sdUpper <- as.numeric(sdUpper)[1]
  targetRatio <- as.numeric(targetRatio)[1]

  if (
    !is.finite(delta) ||
      delta <= eps ||
      !is.finite(sdLower) ||
      sdLower < 0 ||
      !is.finite(sdUpper) ||
      sdUpper <= eps ||
      !is.finite(targetRatio) ||
      targetRatio <= eps
  ) {
    return(NA_real_)
  }

  a <- delta^2
  b <- -(targetRatio^2) * sdLower^2
  cc <- -(targetRatio^2) * sdUpper^2

  disc <- b^2 - 4 * a * cc
  if (!is.finite(disc) || disc < 0) {
    return(NA_real_)
  }

  z <- (-b + sqrt(disc)) / (2 * a)
  if (!is.finite(z) || z <= eps) {
    return(NA_real_)
  }

  sqrt(z)
}

#' Apply the upper-population ratio correction
#'
#' The target ratio is measured on the pre-transformation scale as the distance from the
#' upper mean to the lower reference mean, divided by the combined spread of the
#' lower and upper populations.  The achieved ratio is measured after the
#' transformation in the same way.  Only the upper population is adjusted: its mean
#' distance is multiplied by c and its SD is multiplied by 1 / c.
#'
#' @keywords internal
.simCytRatioAdjustUpper <- function(
  yUpper,
  xUpperRaw,
  xLowerRaw,
  yLower,
  lowerMeanRawReference,
  lowerMeanTransReference,
  eps = sqrt(.Machine$double.eps)
) {
  yUpper <- as.numeric(yUpper)
  xUpperRaw <- as.numeric(xUpperRaw)
  xLowerRaw <- as.numeric(xLowerRaw)
  yLower <- as.numeric(yLower)

  okUpper <- is.finite(yUpper) & is.finite(xUpperRaw)
  okLower <- is.finite(xLowerRaw) & is.finite(yLower)

  if (sum(okUpper) < 2L || sum(okLower) < 2L) {
    return(yUpper)
  }

  yUpperMean <- mean(yUpper[okUpper])
  yUpperSd <- stats::sd(yUpper[okUpper])
  xUpperMean <- mean(xUpperRaw[okUpper])
  xUpperSd <- stats::sd(xUpperRaw[okUpper])
  xLowerSd <- stats::sd(xLowerRaw[okLower])
  yLowerSd <- stats::sd(yLower[okLower])

  if (
    !is.finite(yUpperSd) ||
      yUpperSd <= eps ||
      !is.finite(xUpperSd) ||
      xUpperSd <= eps ||
      !is.finite(xLowerSd) ||
      xLowerSd < 0 ||
      !is.finite(yLowerSd) ||
      yLowerSd < 0 ||
      !is.finite(lowerMeanRawReference) ||
      !is.finite(lowerMeanTransReference)
  ) {
    return(yUpper)
  }

  targetRatio <- (xUpperMean - lowerMeanRawReference) /
    sqrt(xLowerSd^2 + xUpperSd^2)

  deltaTrans <- yUpperMean - lowerMeanTransReference

  cFactor <- .simCytCombinedRatioFactor(
    delta = deltaTrans,
    sdLower = yLowerSd,
    sdUpper = yUpperSd,
    targetRatio = targetRatio,
    eps = eps
  )

  if (!is.finite(cFactor) || cFactor <= eps) {
    return(yUpper)
  }

  yNewMean <- lowerMeanTransReference + cFactor * deltaTrans
  yNewSd <- yUpperSd / cFactor

  out <- yUpper
  out[okUpper] <- (yUpper[okUpper] - yUpperMean) / yUpperSd * yNewSd + yNewMean
  out
}

#' Apply a transformation to a simulated cluster before ratio correction
#'
#' Ratio correction is done at the condition level, because it needs the lower
#' component spread as well as the upper component spread.
#'
#' @keywords internal
.simCytTransformMatrix <- function(simData, transformationFunc) {
  simDataMat <- as.matrix(simData)
  if (is.null(dim(simDataMat))) {
    simDataMat <- matrix(as.numeric(simData), ncol = 1L)
  }
  out <- apply(simDataMat, 2, transformationFunc)
  out <- as.matrix(out)
  if (ncol(out) != ncol(simDataMat)) {
    out <- matrix(out, ncol = ncol(simDataMat))
  }
  out
}

#' @title Simulate a set of stimulation conditions for multiple biological samples
#
#' @description Simulate stimulation conditions (e.g., stimulated and unstimulated) for
#' multiple biological samples, where each sample has one unstimulated condition
#' and one or more stimulated conditions. The outputs are returned as lists of `flowFrame`
#' objects representing each sample-condition combination, and matching lists of cellular
#' cluster labels.
#
#' @param nSample Integer. Number of biological samples to simulate.
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCondition Integer. Number of conditions per sample. The first
#'   condition is unstimulated and the rest are stimulated.
#' @param nCluster Integer. Number of clusters. Must equal `2^nMarker`.
#' @param nCellByCondition Numeric or integer vector. Number of cells per
#'   condition. If length is 1, the value is recycled across all conditions.
#' @param transformationFunc Function. Transformation applied marker-wise to
#'   simulated expression values.
#' @param mixtureType Character. Mixture distribution used for simulation.
#'   Supported values are "gaussianOnly", "tOnly", and "tPlusGauss".
#' @param meanExprMat Numeric matrix. Baseline cluster means with dimensions
#'   `nCluster x nMarker`.
#' @param clusterLabelVec Character vector. Cluster labels of length `nCluster`.
#' @param probVecUns Numeric vector. Baseline cluster probabilities for the
#'   unstimulated condition. Must have length `nCluster` and sum to 1.
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution. Default is `FALSE`.
#' @param probResponseVecByStimCondition NULL or list. If provided, must be a list
#'   of length `nCondition - 1`, where each element is a numeric vector of
#'   length `nCluster`. Each vector is added to `probVecUns` to construct the
#'   stimulated-condition cluster probabilities.
#' @param samplePerturbationSd Numeric. Standard deviation of sample-level
#'   perturbations added to cluster means. Default is 0.
#' @param conditionPerturbationSd Numeric. Standard deviation of condition-level
#'   perturbations added to cluster means within each sample. Default is 0.
#' @param clusterPerturbationSd Numeric. Standard deviation of cluster-level
#'   perturbations applied during cell-level simulation. Default is 0.
#
#' @return A list with two elements:
#'   - `flowFrameList`: A named list of `flowCore::flowFrame` objects.
#'   - `labelsList`: A named list of character vectors of per-cell cluster labels.
#'   Names of list elements are formatted as `sampleXXxUnstim`, `sampleXXX_stim1`, etc.
#
#' @export
simCytExperiment <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns,
  probExact = FALSE,
  probResponseVecByStimCondition = NULL,
  samplePerturbationSd = 0,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2
) {
  # Coerce inputs
  stopifnot(is.numeric(nSample))
  nSample <- as.integer(nSample)
  nMarker <- as.integer(nMarker)
  nCondition <- as.integer(nCondition)
  nCluster <- as.integer(nCluster)

  # Validate inputs using helper
  validateExperimentInputs(
    nSample,
    nMarker,
    nCondition,
    nCluster,
    nCellByCondition,
    transformationFunc
  )

  meanExprMatReference <- meanExprMat

  # Begin Simulation Logic
  nSampleXCondition <- nSample * nCondition
  sampleConditionLabelVec <- lapply(seq_len(nSample), function(currentSample) {
    if (currentSample < 10) {
      sampleName <- paste0("sample00", currentSample)
    } else if (currentSample < 100) {
      sampleName <- paste0("sample0", currentSample)
    } else {
      sampleName <- paste0("sample", currentSample)
    }
    sampleName |>
      paste0("_", c("unstim", paste0("stim", seq_len(nCondition - 1L))))
  }) |>
    unlist()

  flowFrameList <- lapply(seq_len(nSampleXCondition), function(i) NULL) |>
    stats::setNames(sampleConditionLabelVec)
  labelsList <- lapply(seq_len(nSampleXCondition), function(i) NULL) |>
    stats::setNames(sampleConditionLabelVec)

  lapply(seq_len(nSample), function(sampleInd) {
    idxLower <- (sampleInd - 1) * nCondition + 1
    meanExprMatCurrent <- if (samplePerturbationSd == 0L) {
      meanExprMat
    } else {
      meanExprMat +
        matrix(
          rep(
            rnorm(nMarker, mean = 0, sd = samplePerturbationSd),
            each = nCluster
          ),
          byrow = FALSE,
          nrow = nCluster,
          ncol = nMarker
        )
    }

    outListSample <- simCytSample(
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = mixtureType,
      meanExprMat = meanExprMatCurrent,
      meanExprMatReference = meanExprMatReference,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    for (condInd in seq_len(nCondition)) {
      flowFrameList[[idxLower + condInd - 1]] <<- outListSample$flowFrameList[[
        condInd
      ]]
      labelsList[[
        idxLower + condInd - 1
      ]] <<- outListSample$conditionLabelsList[[condInd]]
    }
    NULL
  })

  list(
    flowFrameList = flowFrameList,
    labelsList = labelsList
  )
}

#' @title Simulate all cytokine combinations for a set of stimulation conditions for a single biological sample
#'
#' @description Simulate a set of stimulation conditions (e.g., stimulated and unstimulated) for
#' a single biological sample, where the first condition is always the unstimulated condition.
#'
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCondition Integer. Number of stimulation conditions (must be >= 2).
#' @param nCluster Integer. Number of clusters (must be a power of 2 between 2 and 1024).
#' @param nCellByCondition Integer or numeric vector. Number of cells per condition. If a single
#'   value, it is recycled for all conditions.
#' @param transformationFunc Function. Transformation to apply to simulated data (e.g., identity or
#'   gamma transformation).
#' @param mixtureType Character. Type of mixture distribution: "gaussianOnly", "tOnly", or
#'   "tPlusGauss".
#' @param meanExprMat Numeric matrix. Cluster mean vectors (nCluster x nMarker).
#' @param clusterLabelVec Character vector. Labels for each cluster (length nCluster).
#' @param probVecUns Numeric vector. Probability distribution for unstimulated condition
#'   (length nCluster, sums to 1).
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution.
#' Default is `FALSE`.
#' @param probResponseVecByStimCondition List. Probability response vectors for each stimulated
#'   condition (length nCondition - 1, each of length nCluster).
#' @param conditionPerturbationSd Numeric. Standard deviation of condition-level perturbations
#'   to cluster means.
#' @param clusterPerturbationSd Numeric. Standard deviation of cluster-level perturbations
#'   within each condition.
#'
#' @return A list of length nCondition with named elements "unstim", "stim1", etc. Each element
#'   contains:
#'   - `conditionMatrix`: Numeric matrix of simulated data (nCell x nMarker).
#'   - `conditionLabels`: Character vector of cluster labels for each cell.
#'
#' @keywords internal
simCytSample <- function(
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns,
  probExact,
  probResponseVecByStimCondition = NULL,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
) {
  # Validate inputs using helper
  validateSampleInputs(
    nCondition,
    nMarker,
    nCluster,
    nCellByCondition,
    transformationFunc,
    probResponseVecByStimCondition,
    probVecUns,
    probExact,
    conditionPerturbationSd,
    clusterPerturbationSd,
    meanExprMat,
    clusterLabelVec,
    covEvMin,
    covEvMax
  )

  # Begin Simulation Logic
  nCellByCondition <- if (length(nCellByCondition) == 1L) {
    rep(nCellByCondition, nCondition)
  } else {
    nCellByCondition
  }

  probVecByCondition <- list(probVecUns)
  if (!is.null(probResponseVecByStimCondition)) {
    probVecByCondition <- probVecByCondition |>
      append(lapply(probResponseVecByStimCondition, function(probResponseVec) {
        probVecUns + probResponseVec
      }))
  }

  conditionLabelVec <- c("unstim", paste0("stim", seq_len(nCondition - 1L)))
  flowList <- lapply(seq_len(nCondition), function(i) NULL) |>
    stats::setNames(conditionLabelVec)
  labelsList <- lapply(seq_len(nCondition), function(i) NULL) |>
    stats::setNames(conditionLabelVec)

  lapply(seq_len(nCondition), function(i) {
    meanExprMat <- if (conditionPerturbationSd == 0L) {
      meanExprMat
    } else {
      meanExprMat +
        matrix(
          rep(
            rnorm(nMarker, mean = 0, sd = conditionPerturbationSd),
            each = nCluster
          ),
          byrow = FALSE,
          nrow = nCluster,
          ncol = nMarker
        )
    }

    outListCondition <- simCytCondition(
      nMarker = nMarker,
      nCell = nCellByCondition[[i]],
      transformationFunc = transformationFunc,
      mixtureType = mixtureType,
      meanExprMat = meanExprMat,
      meanExprMatReference = meanExprMatReference,
      clusterLabelVec = clusterLabelVec,
      probVec = probVecByCondition[[i]],
      probExact = probExact,
      clusterPerturbationSd = clusterPerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    # Create annotated data frame
    paramMeta <- data.frame(
      name = paste0("F", seq_len(nMarker)),
      desc = paste0("MarkerF", seq_len(nMarker)),
      range = apply(outListCondition$conditionMatrix, 2, max),
      minRange = apply(outListCondition$conditionMatrix, 2, min),
      maxRange = apply(outListCondition$conditionMatrix, 2, max)
    )
    rownames(paramMeta) <- paste0("$P", seq_len(nMarker))

    paramAnnotated <- Biobase::AnnotatedDataFrame(
      data = paramMeta,
      varMetadata = data.frame(
        labelDescription = c(
          "Name of instrument channel",
          "Actual marker description",
          "Range of values",
          "Minimum binary value",
          "Maximum binary value"
        ),
        row.names = c("name", "desc", "range", "minRange", "maxRange")
      )
    )

    ff <- flowCore::flowFrame(
      exprs = outListCondition$conditionMatrix,
      parameters = paramAnnotated
    )
    flowList[[i]] <<- ff
    labelsList[[i]] <<- outListCondition$conditionLabels
    NULL
  })

  list(
    "flowFrameList" = flowList,
    "conditionLabelsList" = labelsList
  )
}

#' @title Simulate cytometric data for a single stimulation condition
#'
#' @description Simulate flow cytometric data for a single stimulation condition.
#'
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCell Integer. Total number of cells to simulate.
#' @param transformationFunc Function. Transformation to apply to simulated data.
#' @param mixtureType Character. Type of mixture distribution.
#' @param meanExprMat Numeric matrix. Cluster mean vectors.
#' @param probVec Character vector. Cluster labels.
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution.
#' Default is `FALSE`.
#' @param clusterPerturbationSd Numeric. Cluster-level perturbation SD.
#'
#' @return A list with `conditionMatrix` and `conditionLabels`.
#'
#' @keywords internal
simCytCondition <- function(
  nMarker,
  nCell,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVec,
  probExact = FALSE,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
) {
  numClusters <- nrow(meanExprMat)

  if (is.null(meanExprMatReference)) {
    meanExprMatReference <- meanExprMat
  }
  stopifnot(is.matrix(meanExprMatReference))
  stopifnot(all(dim(meanExprMatReference) == dim(meanExprMat)))

  ratioCorrection <- .simCytUsesUpperRatioCorrection(transformationFunc)
  lowerMeanExprVecReference <- apply(meanExprMatReference, 2, min)

  nCellVec <- if (probExact) {
    initAllocVec <- round(nCell * probVec)
    nAlloc <- sum(initAllocVec)
    if (nAlloc != nCell) {
      diffAlloc <- nCell - nAlloc
      if (diffAlloc > 0) {
        probOrder <- order(probVec, decreasing = TRUE)
        for (i in seq_len(diffAlloc)) {
          initAllocVec[probOrder[i]] <- initAllocVec[probOrder[i]] + 1
        }
      } else {
        probOrder <- order(probVec, decreasing = FALSE)
        for (i in seq_len(-diffAlloc)) {
          initAllocVec[probOrder[i]] <- initAllocVec[probOrder[i]] - 1
        }
      }
    }
    initAllocVec
  } else {
    as.vector(t(stats::rmultinom(1, nCell, probVec)))
  }
  nCellVecCum <- cumsum(nCellVec)

  nCellVecObserved <- nCellVec[nCellVec > 0L]
  clusterLabelVecObserved <- clusterLabelVec[nCellVec > 0L]
  cellLabelVec <- lapply(seq_along(nCellVecObserved), function(i) {
    rep(clusterLabelVecObserved[i], nCellVecObserved[i])
  }) |>
    unlist()

  rawDataList <- vector("list", numClusters)
  transformedDataList <- vector("list", numClusters)
  outDataIndClusterList <- vector("list", numClusters)

  for (clusterNumber in seq_len(numClusters)) {
    nCellCluster <- nCellVec[[clusterNumber]]
    if (nCellCluster == 0L) {
      next
    }
    outDataIndClusterLower <- if (clusterNumber == 1L) {
      1L
    } else {
      nCellVecCum[clusterNumber - 1] + 1
    }
    outDataIndClusterUpper <- nCellVecCum[[clusterNumber]]
    outDataIndClusterVec <- seq.int(
      outDataIndClusterLower,
      outDataIndClusterUpper
    )
    meanExprVec <- as.numeric(meanExprMat[clusterNumber, , drop = TRUE])
    simData <- simCytCluster(
      nMarker = nMarker,
      nCell = nCellCluster,
      meanExprVec = meanExprVec,
      perturbationSd = clusterPerturbationSd,
      mixtureType = mixtureType,
      clusterNumber = clusterNumber,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )
    simData <- as.matrix(simData)
    if (ncol(simData) != nMarker) {
      simData <- matrix(simData, ncol = nMarker)
    }

    rawDataList[[clusterNumber]] <- simData
    transformedDataList[[clusterNumber]] <- .simCytTransformMatrix(
      simData = simData,
      transformationFunc = transformationFunc
    )
    outDataIndClusterList[[clusterNumber]] <- outDataIndClusterVec
  }

  if (isTRUE(ratioCorrection)) {
    for (markerInd in seq_len(nMarker)) {
      lowerMeanRawReference <- lowerMeanExprVecReference[[markerInd]]

      lowerClusterInd <- which(
        meanExprMatReference[, markerInd] <=
          lowerMeanRawReference + sqrt(.Machine$double.eps)
      )
      lowerClusterInd <- lowerClusterInd[
        vapply(
          lowerClusterInd,
          function(i) {
            !is.null(rawDataList[[i]]) && nrow(rawDataList[[i]]) > 0L
          },
          logical(1)
        )
      ]

      if (length(lowerClusterInd) == 0L) {
        next
      }

      xLowerRaw <- unlist(lapply(
        lowerClusterInd,
        function(i) rawDataList[[i]][, markerInd, drop = TRUE]
      ))
      yLower <- unlist(lapply(
        lowerClusterInd,
        function(i) transformedDataList[[i]][, markerInd, drop = TRUE]
      ))

      lowerMeanTransReference <- .simCytLowerMeanTransReference(
        transformationFunc = transformationFunc,
        lowerMeanRawReference = lowerMeanRawReference,
        yLower = yLower
      )
      if (!is.finite(lowerMeanTransReference)) {
        next
      }

      for (clusterNumber in seq_len(numClusters)) {
        if (is.null(rawDataList[[clusterNumber]])) {
          next
        }
        isUpperPopulation <- meanExprMatReference[clusterNumber, markerInd] >
          lowerMeanRawReference + sqrt(.Machine$double.eps)
        if (!isTRUE(isUpperPopulation)) {
          next
        }

        transformedDataList[[clusterNumber]][, markerInd] <-
          .simCytRatioAdjustUpper(
            yUpper = transformedDataList[[clusterNumber]][, markerInd],
            xUpperRaw = rawDataList[[clusterNumber]][, markerInd],
            xLowerRaw = xLowerRaw,
            yLower = yLower,
            lowerMeanRawReference = lowerMeanRawReference,
            lowerMeanTransReference = lowerMeanTransReference
          )
      }
    }
  }

  outData <- if (nMarker == 1L) {
    rep(NA_real_, nCell)
  } else {
    matrix(NA_real_, nrow = nCell, ncol = nMarker)
  }

  for (clusterNumber in seq_len(numClusters)) {
    if (is.null(transformedDataList[[clusterNumber]])) {
      next
    }
    outDataIndClusterVec <- outDataIndClusterList[[clusterNumber]]
    if (nMarker == 1L) {
      outData[outDataIndClusterVec] <- transformedDataList[[clusterNumber]][,
        1L
      ]
    } else {
      outData[outDataIndClusterVec, ] <- transformedDataList[[clusterNumber]]
    }
  }

  reorderVec <- sample.int(nCell)
  if (nMarker == 1L) {
    outData <- outData[reorderVec]
    outData <- matrix(outData, ncol = 1)
  } else {
    outData <- outData[reorderVec, ]
  }
  colnames(outData) <- paste0("F", seq_len(nMarker))
  cellLabelVec <- cellLabelVec[reorderVec]
  list(
    conditionMatrix = outData,
    conditionLabels = cellLabelVec
  )
}

simCytCluster <- function(
  nMarker,
  nCell,
  meanExprVec,
  perturbationSd = 0,
  mixtureType,
  clusterNumber,
  covEvMin = 1,
  covEvMax = 2
) {
  conditionPerturbationVec <- if (perturbationSd == 0L) {
    meanExprVec
  } else {
    meanExprVec + rnorm(nMarker, mean = 0, sd = perturbationSd)
  }
  currentSigma <- posDef(nMarker, covEvMin, covEvMax)
  simCytClusterData(
    mixtureType = mixtureType,
    clusterNumber = clusterNumber,
    nCell = nCell,
    muVec = conditionPerturbationVec,
    sigmaMat = currentSigma
  )
}

#' @title Generate simulated data from mixture component
#'
#' @description Sample from multivariate normal or t distribution.
#'
#' @param mixtureType Character. "gaussianOnly", "tOnly", or "tPlusGauss".
#' @param clusterNumber Integer. Used for alternating distributions in "tPlusGauss".
#' @param nCell Integer. Number of samples.
#' @param muVec Numeric vector. Mean vector.
#' @param sigmaMat Numeric matrix. Covariance matrix.
#'
#' @return Numeric matrix of sampled data (nCell x length(muVec)).
#'
#' @keywords internal
simCytClusterData <- function(
  mixtureType,
  clusterNumber,
  nCell,
  muVec,
  sigmaMat
) {
  if (mixtureType == "tPlusGauss") {
    if ((clusterNumber %% 2) == 0) {
      MASS::mvrnorm(
        nCell,
        mu = muVec,
        Sigma = sigmaMat
      )
    } else {
      mvtnorm::rmvt(
        nCell,
        delta = muVec,
        sigma = sigmaMat,
        df = 2
      )
    }
  } else if (mixtureType == "tOnly") {
    mvtnorm::rmvt(
      nCell,
      delta = muVec,
      sigma = sigmaMat,
      df = 2
    )
  } else {
    MASS::mvrnorm(
      nCell,
      mu = muVec,
      Sigma = sigmaMat
    )
  }
}

#' @title Validate inputs for simCytExperiment
#' @keywords internal
validateExperimentInputs <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc
) {
  stopifnot(nSample > 0L)
  stopifnot(nMarker > 0L)
  stopifnot(nCondition > 1L)
  stopifnot(nCluster > 0L)
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, nCondition))
  stopifnot(all(nCellByCondition > 0))
}

#' @title Validate inputs for simCytSample
#' @keywords internal
validateSampleInputs <- function(
  nCondition,
  nMarker,
  nCluster,
  nCellByCondition,
  transformationFunc,
  probResponseVecByStimCondition,
  probVecUns,
  probExact,
  conditionPerturbationSd,
  clusterPerturbationSd,
  meanExprMat,
  clusterLabelVec,
  covEvMin,
  covEvMax
) {
  stopifnot(is.logical(probExact))
  stopifnot(is.integer(nCondition))
  stopifnot(is.integer(nMarker))
  stopifnot(nCondition > 1L)
  stopifnot(is.integer(nCluster))
  stopifnot(nCluster > 0L)
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, nCondition))
  stopifnot(all(nCellByCondition > 0))

  if (!is.null(probResponseVecByStimCondition)) {
    stopifnot(is.list(probResponseVecByStimCondition))
    stopifnot(all(sapply(probResponseVecByStimCondition, is.numeric)))
    stopifnot(length(probResponseVecByStimCondition) == (nCondition - 1L))
    stopifnot(all(sapply(probResponseVecByStimCondition, length) == nCluster))
  }

  stopifnot(is.numeric(probVecUns))
  stopifnot(length(probVecUns) == nCluster)
  stopifnot(all(probVecUns >= 0))
  stopifnot(all(probVecUns <= 1))
  stopifnot(abs(sum(probVecUns) - 1) < 1e-6)

  stopifnot(is.numeric(conditionPerturbationSd))
  stopifnot(length(conditionPerturbationSd) == 1L)
  stopifnot(conditionPerturbationSd >= 0)

  stopifnot(is.numeric(clusterPerturbationSd))
  stopifnot(length(clusterPerturbationSd) == 1L)
  stopifnot(clusterPerturbationSd >= 0)

  stopifnot(is.matrix(meanExprMat))
  stopifnot(nrow(meanExprMat) == nCluster)
  stopifnot(ncol(meanExprMat) == nMarker)
  stopifnot(!any(is.na(meanExprMat)))
  stopifnot(is.numeric(meanExprMat))

  stopifnot(is.character(clusterLabelVec))
  stopifnot(length(clusterLabelVec) == nCluster)
}
