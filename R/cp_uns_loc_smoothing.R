# Local-FDR probability smoothing
#
# Fits the monotone response-probability curve, supplies fallback smoothers, and
# stores the finite-difference derivative evaluated from the fitted curve.

.getCpUnsLocGetProbSmooth <- function(
  dataMod,
  stage,
  pathProject,
  chnl,
  chnlSettings = list()
) {
  stageChnl <- file.path(stage, chnl)
  densityBw <- attr(dataMod, "locDensityBw")

  if (!.getCpUnsLocGetProbSmoothCheckNCell(dataMod)) {
    .intSaveNm(
      "not_enough_cells_to_smooth",
      NULL,
      .getInd(dataMod),
      stageChnl,
      pathProject
    )
    dataModOut <- .getCpUnsLocGetProbSmoothCheckNCellOut(dataMod)
    attr(dataModOut, "locDensityBw") <- densityBw
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
  dataModOut <- dataMod |> dplyr::mutate(pred = predVec)
  dataModOut <- .getCpUnsLocGetProbSmoothAttachDeriv(
    dataMod = dataModOut,
    smoothObj = smoothObj
  )
  attr(dataModOut, "locDensityBw") <- densityBw
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
    dataMod |> dplyr::mutate(pred = probSmooth - 1e-4) # nolint
  } else {
    dataMod
  }
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
  chnlSettings = list()
) {
  fit1 <- .getCpUnsLocGetProbSmoothActualFirst(dataMod, stage)
  .getCpUnsLocGetProbSmoothActualFirstResponse(
    fit = fit1,
    dataMod = dataMod,
    stage = stage,
    chnlSettings = chnlSettings
  )
}

.getCpUnsLocGetProbSmoothActualFirst <- function(dataMod, stage) {
  .debug("Smoothing I")

  idxMod <- attr(dataMod, "idxMod") %||% seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, ] |>
    dplyr::mutate(
      probSmooth = pmin(probSmooth, 0.999),
      probSmooth = pmax(probSmooth, 0.001)
    )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(try(stop(), silent = TRUE)) # nolint
  }

  k <- min(n - 1L, 20L)

  fml <- as.formula(paste0(
    "probSmooth ~ s(`",
    chnl,
    "`, bs = 'mpi', k = ",
    k,
    ", m = c(2, 1))"
  ))

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
        bfgs = list(steptol.bfgs = 1e-1),
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
  chnlSettings = list()
) {
  if (.getCpUnsLocGetProbSmoothActualCheck(fit, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_mpi"
      )
    )
  }
  .getCpUnsLocGetProbSmoothActualFirstResponseFailure(
    stage = stage,
    dataMod = dataMod,
    chnlSettings = chnlSettings
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualCheck <- function(fit, dataMod) {
  if (inherits(fit, "try-error") || is.null(fit)) {
    return(FALSE)
  }
  outList <- .getCpUnsLocGetProbSmoothActualResponseSuccess(
    fit = fit,
    dataMod = dataMod
  )
  !(all(outList$pred > 0.99) || outList$meanAbsError > 0.3) # nolint
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualResponseSuccess <- function(
  fit,
  dataMod,
  chnlSettings = list(),
  method = NA_character_
) {
  # predict on full dataset
  predVec <- predict(fit, newdata = dataMod, type = "response")
  meanAbsError <- mean(abs(predVec - dataMod$probSmooth))
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
  chnlSettings = list()
) {
  chnl <- .getCpUnsLocGetChnl(dataMod)
  x <- suppressWarnings(as.numeric(.getCut(dataMod)))
  x <- x[is.finite(x)]
  if (length(x) < 4L || diff(range(x)) <= 0) {
    return(NULL)
  }

  nGrid <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivGridN",
    512L
  )
  nGrid <- suppressWarnings(as.integer(nGrid[1]))
  if (!is.finite(nGrid) || nGrid < 25L) {
    nGrid <- 512L
  }

  epsFrac <- .getCpUnsLocGetCpTrimSetting(
    chnlSettings,
    "locFlatDerivEpsFrac",
    1e-5
  )
  epsFrac <- suppressWarnings(as.numeric(epsFrac[1]))
  if (!is.finite(epsFrac) || epsFrac <= 0) {
    epsFrac <- 1e-5
  }

  xRange <- range(x, na.rm = TRUE)
  xWidth <- diff(xRange)
  xGrid <- seq(xRange[1], xRange[2], length.out = nGrid)
  eps <- max(
    xWidth * epsFrac,
    sqrt(.Machine$double.eps) * max(abs(xRange), 1, na.rm = TRUE)
  )
  if (!is.finite(eps) || eps <= 0) {
    return(NULL)
  }

  xLeft <- pmax(xRange[1], xGrid - eps)
  xRight <- pmin(xRange[2], xGrid + eps)
  denom <- xRight - xLeft
  if (any(!is.finite(denom) | denom <= 0)) {
    return(NULL)
  }

  newData <- dataMod[rep(1L, length(xGrid)), , drop = FALSE]
  newData[[chnl]] <- xGrid
  predGrid <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  newData[[chnl]] <- xLeft
  predLeft <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  newData[[chnl]] <- xRight
  predRight <- try(
    predict(fit, newdata = newData, type = "response"),
    silent = TRUE
  )
  if (
    inherits(predGrid, "try-error") ||
      inherits(predLeft, "try-error") ||
      inherits(predRight, "try-error")
  ) {
    return(NULL)
  }

  deriv <- (as.numeric(predRight) - as.numeric(predLeft)) / denom
  deriv[!is.finite(deriv)] <- 0
  deriv <- pmax(0, deriv)
  tibble::tibble(
    x = xGrid,
    pred = as.numeric(predGrid),
    deriv = deriv
  )
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualFirstResponseFailure <- function(
  stage,
  dataMod,
  chnlSettings = list()
) {
  fit2 <- .getCpUnsLocGetProbSmoothActualSecond(dataMod, stage)
  if (.getCpUnsLocGetProbSmoothActualCheck(fit2, dataMod)) {
    .debug("Smoothed") # nolint
    return(
      .getCpUnsLocGetProbSmoothActualResponseSuccess(
        fit = fit2,
        dataMod = dataMod,
        chnlSettings = chnlSettings,
        method = "scam_micv"
      )
    )
  }
  .getCpUnsLocGetProbSmoothActualThird(dataMod, stage)
}

#' @keywords internal
.getCpUnsLocGetProbSmoothActualSecond <- function(dataMod, stage) {
  .debug("Smoothing II") # nolint

  idxMod <- attr(dataMod, "idxMod") %||% seq_len(nrow(dataMod))
  chnl <- .getCpUnsLocGetChnl(dataMod)

  dataMod <- dataMod[idxMod, ] |>
    dplyr::mutate(
      probSmooth = pmin(probSmooth, 0.999),
      probSmooth = pmax(probSmooth, 0.001)
    )

  n <- nrow(dataMod)

  if (n <= 4L) {
    return(NULL)
  }

  k <- min(n - 1L, 20L)

  fml <- as.formula(paste0(
    "probSmooth ~ s(`",
    chnl,
    "`, bs = 'micv', k = ",
    k,
    ", m = c(2, 1))"
  ))

  try(
    suppressWarnings(scam::scam(
      fml,
      family = "binomial",
      data = dataMod,
      control = scam::scam.control(
        print.warn = FALSE,
        trace = FALSE,
        devtol.fit = 0.5,
        steptol.fit = 1e-1,
        maxHalf = 5,
        bfgs = list(steptol.bfgs = 1e-1),
        maxit = 1e1
      )
    )),
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


.getCpUnsLocGetCpTrimSetting <- function(chnlSettings, nm, default) {
  .getCpUnsLocSetting(chnlSettings, nm, default)
}
