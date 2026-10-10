# Local-FDR probability smoothing
#
# Fits the monotone response-probability curve, supplies fallback smoothers, and
# stores the finite-difference derivative evaluated from the fitted curve.

#' Smooth response probabilities while retaining negative-width diagnostics
#' @keywords internal
.getCpUnsLocGetProbSmooth <- function(
  dataMod,
  stage,
  pathProject,
  chnl,
  chnlSettings = list()
) {
  stageChnl <- file.path(stage, chnl)
  retainedAttrs <- c(
    "locDensityBw",
    "locStimDensity",
    "locDensityComparison",
    "locPeakX",
    "locWindowWidth",
    "locWindowWidthInfo"
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
  } else {
    smoothObj <- .getCpUnsLocGetProbSmoothFit(
      dataMod = dataMod,
      chnlSettings = chnlSettings
    )
    dataModOut <- dataMod
    dataModOut$pred <- smoothObj$pred
    if (!is.null(smoothObj$derivTbl)) {
      attr(dataModOut, "locProbDerivTbl") <- smoothObj$derivTbl
    }
    attr(dataModOut, "locProbSmoothMethod") <- smoothObj$method
  }

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


#' Fit the probability smoother, trying each SCAM specification in turn
#'
#' Returns the first acceptable fit's predictions and derivative table, or the
#' probSmooth fallback when no fit is acceptable.
#' @keywords internal
.getCpUnsLocGetProbSmoothFit <- function(dataMod, chnlSettings = list()) {
  specList <- list(
    list(
      msg = "Smoothing I", bs = "mpi", family = "quasibinomial",
      quiet = FALSE, method = "scam_mpi"
    ),
    list(
      msg = "Smoothing II", bs = "micv", family = "binomial",
      quiet = TRUE, method = "scam_micv"
    )
  )
  for (spec in specList) {
    .debug(spec$msg) # nolint
    fit <- .fitScam(
      dataMod = dataMod,
      bs = spec$bs,
      family = spec$family,
      quiet = spec$quiet
    )
    # Evaluate the full-data prediction once; only an accepted fit needs the
    # derivative table.
    fitEval <- .getCpUnsLocGetProbSmoothFitEval(
      fit = fit,
      dataMod = dataMod
    )
    if (.getCpUnsLocGetProbSmoothFitEvalCheck(fitEval)) {
      .debug("Smoothed") # nolint
      return(list(
        "pred" = fitEval$pred,
        "meanAbsError" = fitEval$meanAbsError,
        "derivTbl" = .getCpUnsLocGetProbSmoothDerivativeTbl(
          fit = fit,
          dataMod = dataMod,
          chnlSettings = chnlSettings
        ),
        "method" = spec$method
      ))
    }
  }
  .getCpUnsLocGetProbSmoothFallback(dataMod)
}


#' Fit a monotone increasing SCAM to the modelled probabilities
#'
#' Returns NULL when there are too few points and a try-error when the fit
#' fails. quiet suppresses warnings raised while fitting. The expression
#' values are fitted under the fixed column name `x`, because channel names
#' such as `PE-A` are not valid in a model formula.
#' @keywords internal
.fitScam <- function(dataMod, bs, family, quiet) {
  idxMod <- attr(dataMod, "idxMod") %||%
    seq_len(nrow(dataMod))

  dataMod <- dataMod[idxMod, , drop = FALSE]
  dataMod <- data.frame(
    x = as.numeric(.getCut(dataMod)),
    probSmooth = dataMod$probSmooth
  )
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
      "probSmooth ~ s(x, bs = '",
      bs,
      "', k = ",
      k,
      ", m = c(2, 1))"
    )
  )

  fitScam <- function() {
    scam::scam(
      fml,
      family = family,
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
  }

  try(
    if (quiet) suppressWarnings(fitScam()) else fitScam(),
    silent = TRUE
  )
}


#' Construct the minimal prediction data required by the smoother
#'
#' The fitted SCAM has only the expression values, named `x`, as a predictor
#' (see `.fitScam()`), so prediction does not require copying dataMod.
#' @keywords internal
.getCpUnsLocGetProbSmoothNewData <- function(x) {
  data.frame(x = x)
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

  x <- suppressWarnings(
    as.numeric(.getCut(dataMod))
  )

  newData <- .getCpUnsLocGetProbSmoothNewData(x = x)

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


#' @keywords internal
.getCpUnsLocGetProbSmoothDerivativeTbl <- function(
  fit,
  dataMod,
  chnlSettings = list()
) {

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

  nGrid <- .getCpUnsLocSetting(
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

  epsFrac <- .getCpUnsLocSetting(
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

  xRange <- range(x)
  xWidth <- diff(xRange)

  xGrid <- seq(
    xRange[1],
    xRange[2],
    length.out = nGrid
  )

  eps <- max(
    xWidth * epsFrac,
    sqrt(.Machine$double.eps) *
      max(abs(xRange), 1)
  )

  xLeft <- pmax(
    xRange[1],
    xGrid - eps
  )
  xRight <- pmin(
    xRange[2],
    xGrid + eps
  )

  denom <- xRight - xLeft

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

  newData <- .getCpUnsLocGetProbSmoothNewData(x = predX)

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
.getCpUnsLocGetProbSmoothFallback <- function(dataMod) {
  .debug("Failed to smooth") # nolint

  list(
    "pred" = dataMod$probSmooth - 0.0001,
    "meanAbsError" = NA_real_,
    "derivTbl" = NULL,
    "method" = "probSmooth_fallback"
  )
}
