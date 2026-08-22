# Local-FDR probability smoothing
#
# Fits the monotone response-probability curve, supplies fallback smoothers, and
# stores the finite-difference derivative evaluated from the fitted curve.

.getCpUnsLocGetProbSmooth <- function(
    dataMod,
    stage,
    pathProject,
    chnl,
    chnlSettings = list()) {
  stageChnl <- file.path(stage, chnl)
  retainedAttrs <- c(
    "locDensityBw",
    "locStimDensity",
    "locDensityComparison",
    "locPeakX",
    "locWindowWidth"
  )
  retainedValues <- stats::setNames(
    lapply(retainedAttrs, function(name) attr(dataMod, name)),
    retainedAttrs
  )

  if (!.getCpUnsLocGetProbSmoothCheckNCell(dataMod)) {
    .intSaveNm(
      "not_enough_cells_to_smooth",
      NULL,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    dataModOut <- .getCpUnsLocGetProbSmoothCheckNCellOut(dataMod)
    for (name in retainedAttrs) {
      if (!is.null(retainedValues[[name]])) {
        attr(dataModOut, name) <- retainedValues[[name]]
      }
    }
    .intSaveNm(
      "probSmoothOut",
      dataModOut,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    return(dataModOut)
  }

  smoothObj <- .getCpUnsLocGetProbSmoothActual(
    dataMod = dataMod,
    stage = stage,
    chnlSettings = chnlSettings
  )

  predVec <- .getCpUnsLocGetProbSmoothObjPred(smoothObj)

  dataModOut <- dataMod
  dataModOut$pred <- predVec

  dataModOut <- .getCpUnsLocGetProbSmoothAttachDeriv(
    dataMod = dataModOut,
    smoothObj = smoothObj
  )

  for (name in retainedAttrs) {
    if (!is.null(retainedValues[[name]])) {
      attr(dataModOut, name) <- retainedValues[[name]]
    }
  }

  .intSaveNm(
    "probSmoothOut",
    dataModOut,
    .getInd(dataMod),
    stageChnl,
    pathProject
  )

  dataModOut
}


#' @keywords internal
.getCpUnsLocGetProbSmoothCheckNCell <- function(dataMod) {
  is.data.frame(dataMod) && nrow(dataMod) > 10
}


#' @keywords internal
.getCpUnsLocGetProbSmoothCheckNCellOut <- function(dataMod) {
  if (is.data.frame(dataMod)) {
    dataMod$pred <- dataMod$probSmooth - 1e-4
  }
  dataMod
}


#' @keywords internal
.getCpUnsLocGetProbSmoothObjPred <- function(smoothObj) {
  if (is.list(smoothObj) && "pred" %in% names(smoothObj)) {
    return(smoothObj$pred)
  }
  smoothObj
}


#' @keywords internal
.getCpUnsLocGetProbSmoothAttachDeriv <- function(dataMod, smoothObj) {
  if (!is.list(smoothObj)) {
    return(dataMod)
  }

  if (!is.null(smoothObj$derivTbl)) {
    attr(dataMod, "locProbDerivTbl") <- smoothObj$derivTbl
  }

  if (!is.null(smoothObj$method)) {
    attr(dataMod, "locProbSmoothMethod") <- smoothObj$method
  }

  dataMod
}


#' @keywords internal
.getCpUnsLocGetProbSmoothActual <- function(
    dataMod,
    stage,
    chnlSettings = list()) {
  fit1 <- .getCpUnsLocGetProbSmoothActualFirst(
    dataMod,
    stage
  )

  .getCpUnsLocGetProbSmoothActualFirstResponse(
    fit = fit1,
    dataMod = dataMod,
    stage = stage,
    chnlSettings = chnlSettings
  )
}


#' Construct the minimal prediction data required by the smoother
#'
#' The fitted SCAM contains only the expression channel as a predictor, so
#' prediction does not require copying every column of dataMod.
#' @keywords internal
.getCpUnsLocGetProbSmoothNewData <- function(chnl, x) {
  out <- data.frame(
    value = x,
    check.names = FALSE
  )
  names(out) <- chnl
  out
}


#' Predict one fitted smoother over the full model data
#'
#' This evaluates the expensive full-data prediction once and calculates the
#' fit diagnostic used to decide whether the smoother is acceptable. The
#' derivative is deliberately not calculated here because a rejected fit does
#' not need one.
#' @keywords internal
.getCpUnsLocGetProbSmoothFitEval <- function(fit, dataMod) {
  if (inherits(fit, "try-error") || is.null(fit)) {
    return(NULL)
  }

  chnl <- .getCpUnsLocGetChnl(dataMod)
  x <- suppressWarnings(
    as.numeric(.getCut(dataMod))
  )

  newData <- .getCpUnsLocGetProbSmoothNewData(
    chnl = chnl,
    x = x
  )

  predVec <- try(
    stats::predict(
      fit,
      newdata = newData,
      type = "response"
    ),
    silent = TRUE
  )

  if (inherits(predVec, "try-error")) {
    return(NULL)
  }

  predVec <- as.numeric(predVec)
  meanAbsError <- mean(
    abs(predVec - dataMod$probSmooth)
  )

  list(
    pred = predVec,
    meanAbsError = meanAbsError
  )
}


#' Test whether an already evaluated smoother is acceptable
#' @keywords internal
.getCpUnsLocGetProbSmoothFitEvalCheck <- function(fitEval) {
  if (is.null(fitEval)) {
    return(FALSE)
  }

  !(all(fitEval$pred > 0.99) ||
    fitEval$meanAbsError > 0.3)
}


.getCpUnsLocGetProbSmoothActualFirst <- function(dataMod, stage) {
  .debug("Smoothing I")

  idxMod <- attr(dataMod, "idxMod") %||%
    seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, , drop = FALSE]
  dataMod$probSmooth <- pmin(
    dataMod$probSmooth,
    0.999
  )
  dataMod$probSmooth <- pmax(
    dataMod$probSmooth,
    0.001
  )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(try(stop(), silent = TRUE)) # nolint
  }

  k <- min(n - 1L, 20L)

  fml <- stats::as.formula(
    paste0(
      "probSmooth ~ s(`",
      chnl,
      "`, bs = 'mpi', k = ",
      k,
      ", m = c(2, 1))"
    )
  )

  try(
    scam::scam(
      fml,
      family = "quasibinomial",
      data = dataMod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(
          steptol.bfgs = 1e-1
        ),
        maxit = 1e1
      )
    ),
    silent = TRUE
  )
}


#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponse <- function(
    fit,
    dataMod,
    stage,
    chnlSettings = list()) {
  # Evaluate the full-data prediction once. Previously the check performed a
  # full prediction and derivative calculation, then the success path repeated
  # both operations.
  fitEval <- .getCpUnsLocGetProbSmoothFitEval(
    fit = fit,
    dataMod = dataMod
  )

  if (.getCpUnsLocGetProbSmoothFitEvalCheck(fitEval)) {
    .debug("Smoothed") # nolint

    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_mpi",
        predVec = fitEval$pred,
        meanAbsError = fitEval$meanAbsError
      )
    )
  }

  .getCpUnsLocGetProbSmoothActualFirstResponseFailure(
    stage = stage,
    dataMod = dataMod,
    chnlSettings = chnlSettings
  )
}


#' Check whether a fitted smoother is acceptable
#'
#' Retained as a separate helper for existing internal callers and tests. Unlike
#' the previous implementation, this does not calculate a derivative table just
#' to decide whether the fit should be retained.
#' @keywords internal
.getCpUnsLocGetProbSmoothActualCheck <- function(fit, dataMod) {
  fitEval <- .getCpUnsLocGetProbSmoothFitEval(
    fit = fit,
    dataMod = dataMod
  )

  .getCpUnsLocGetProbSmoothFitEvalCheck(fitEval)
}


#' Return the full successful smoothing result
#'
#' predVec and meanAbsError may be supplied when they were already calculated
#' during fit validation. This avoids repeating the full-data prediction.
#' @keywords internal
.getCpUnsLocGetProbSmoothActualResponseSuccess <- function(
    fit,
    dataMod,
    chnlSettings = list(),
    method = NA_character_,
    predVec = NULL,
    meanAbsError = NULL) {
  if (is.null(predVec)) {
    fitEval <- .getCpUnsLocGetProbSmoothFitEval(
      fit = fit,
      dataMod = dataMod
    )

    if (is.null(fitEval)) {
      stop("Failed to predict from local-FDR probability smoother.")
    }

    predVec <- fitEval$pred

    if (is.null(meanAbsError)) {
      meanAbsError <- fitEval$meanAbsError
    }
  }

  if (is.null(meanAbsError)) {
    meanAbsError <- mean(
      abs(predVec - dataMod$probSmooth)
    )
  }

  derivTbl <- .getCpUnsLocGetProbSmoothDerivativeTbl(
    fit = fit,
    dataMod = dataMod,
    chnlSettings = chnlSettings
  )

  list(
    "pred" = predVec,
    "meanAbsError" = meanAbsError,
    "derivTbl" = derivTbl,
    "method" = method
  )
}


#' @keywords internal
.getCpUnsLocGetProbSmoothDerivativeTbl <- function(
    fit,
    dataMod,
    chnlSettings = list()) {
  chnl <- .getCpUnsLocGetChnl(dataMod)

  x <- suppressWarnings(
    as.numeric(.getCut(dataMod))
  )
  x <- x[is.finite(x)]

  if (
    length(x) < 4L ||
      diff(range(x)) <= 0
  ) {
    return(NULL)
  }

  nGrid <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivGridN",
    512L
  )
  nGrid <- suppressWarnings(
    as.integer(nGrid[1])
  )

  if (!is.finite(nGrid) || nGrid < 25L) {
    nGrid <- 512L
  }

  epsFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivEpsFrac",
    1e-5
  )
  epsFrac <- suppressWarnings(
    as.numeric(epsFrac[1])
  )

  if (!is.finite(epsFrac) || epsFrac <= 0) {
    epsFrac <- 1e-5
  }

  xRange <- range(x, na.rm = TRUE)
  xWidth <- diff(xRange)

  xGrid <- seq(
    xRange[1],
    xRange[2],
    length.out = nGrid
  )

  eps <- max(
    xWidth * epsFrac,
    sqrt(.Machine$double.eps) *
      max(abs(xRange), 1, na.rm = TRUE)
  )

  if (!is.finite(eps) || eps <= 0) {
    return(NULL)
  }

  xLeft <- pmax(
    xRange[1],
    xGrid - eps
  )
  xRight <- pmin(
    xRange[2],
    xGrid + eps
  )

  denom <- xRight - xLeft

  if (any(!is.finite(denom) | denom <= 0)) {
    return(NULL)
  }

  # Predict all three derivative-grid locations in one call rather than making
  # three separate calls to predict.scam().
  #
  # Ordering:
  #   1:nGrid                   -> xGrid
  #   (nGrid + 1):(2 * nGrid)  -> xLeft
  #   (2 * nGrid + 1):(3*nGrid)-> xRight
  predX <- c(
    xGrid,
    xLeft,
    xRight
  )

  newData <- .getCpUnsLocGetProbSmoothNewData(
    chnl = chnl,
    x = predX
  )

  predAll <- try(
    stats::predict(
      fit,
      newdata = newData,
      type = "response"
    ),
    silent = TRUE
  )

  if (inherits(predAll, "try-error")) {
    return(NULL)
  }

  predAll <- as.numeric(predAll)

  if (length(predAll) != 3L * nGrid) {
    return(NULL)
  }

  idx <- seq_len(nGrid)

  predGrid <- predAll[idx]
  predLeft <- predAll[nGrid + idx]
  predRight <- predAll[2L * nGrid + idx]

  deriv <- (predRight - predLeft) / denom

  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)

  tibble::tibble(
    x = xGrid,
    pred = predGrid,
    deriv = deriv
  )
}


#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponseFailure <- function(
    stage,
    dataMod,
    chnlSettings = list()) {
  fit2 <- .getCpUnsLocGetProbSmoothActualSecond(
    dataMod,
    stage
  )

  fitEval <- .getCpUnsLocGetProbSmoothFitEval(
    fit = fit2,
    dataMod = dataMod
  )

  if (.getCpUnsLocGetProbSmoothFitEvalCheck(fitEval)) {
    .debug("Smoothed") # nolint

    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit2,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_micv",
        predVec = fitEval$pred,
        meanAbsError = fitEval$meanAbsError
      )
    )
  }

  .getCpUnsLocGetProbSmoothActualThird(
    dataMod,
    stage
  )
}


#' @keywords internal
.getCpUnsLocGetProbSmoothActualSecond <- function(dataMod, stage) {
  .debug("Smoothing II") # nolint

  idxMod <- attr(dataMod, "idxMod") %||%
    seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, , drop = FALSE]
  dataMod$probSmooth <- pmin(
    dataMod$probSmooth,
    0.999
  )
  dataMod$probSmooth <- pmax(
    dataMod$probSmooth,
    0.001
  )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(NULL)
  }

  k <- min(n - 1L, 20L)

  fml <- stats::as.formula(
    paste0(
      "probSmooth ~ s(`",
      chnl,
      "`, bs = 'micv', k = ",
      k,
      ", m = c(2, 1))"
    )
  )

  try(
    suppressWarnings(
      scam::scam(
        fml,
        family = "binomial",
        data = dataMod,
        control = scam::scam.control(
          print.warn = FALSE,
          trace = FALSE,
          devtol.fit = 0.5,
          steptol.fit = 1e-1,
          maxHalf = 5,
          bfgs = list(
            steptol.bfgs = 1e-1
          ),
          maxit = 1e1
        )
      )
    ),
    silent = TRUE
  )
}


#' @keywords internal
.getCpUnsLocGetProbSmoothActualThird <- function(dataMod, stage) {
  .debug("Failed to smooth") # nolint

  list(
    "pred" = dataMod$probSmooth - 0.0001,
    "meanAbsError" = NA_real_,
    "derivTbl" = NULL,
    "method" = "probSmooth_fallback"
  )
}


.getCpUnsLocGetCpTrimSetting <- function(
    chnlSettings,
    nm,
    default) {
  .getCpUnsLocSetting(
    chnlSettings,
    nm,
    default
  )
}
