# Get cutpoints using clustering approach
# Clusters thresholds from similar distributions and imputes local-FDR-missing
# thresholds from equivalent tolerance values.
#' @keywords internal
.getCpCluster <- function(
  .data,
  gateTbl,
  gateStatsTbl,
  gateTblCtrl,
  chnlSettings,
  stage,
  pathProject,
  control = list(),
  filterOtherCytPos,
  calcCytPosGates,
  indBatchList
) {
  .browse(seq_along(.data)) # nolint
  stageChnl <- file.path(stage, chnlSettings$chnlCut)
  control <- .getCpClusterControlUpdate(control) # nolint

  gateTbl <- .getCpClusterLocGateTblPrepare(gateTbl)
  cpMin <- .getCpClusterCpGetMin(gateTbl, gateTblCtrl)

  exLookup <- .getCpClusterLocExLookup(
    .data = .data,
    indBatchList = indBatchList,
    chnlSettings = chnlSettings,
    filterOtherCytPos = filterOtherCytPos,
    calcCytPosGates = calcCytPosGates,
    gateTbl = gateTbl,
    pathProject = pathProject
  )

  gateTblStim <- .getCpClusterLocGateTblStim(
    gateTbl = gateTbl,
    exLookup = exLookup
  )

  if (nrow(gateTblStim) == 0L) {
    cpTbl <- tibble::tibble()
    cpTbl <- .getCpClusterLocAddFinalDetail(cpTbl, exLookup)
    .intSaveNm("locDetailClusterFinal", cpTbl, "all", stageChnl, pathProject)
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  if (all(gateTblStim$locGenerated %in% TRUE)) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "all_local_fdr_thresholds_available_skip_cluster"
    )
    cpTbl <- .getCpClusterLocAddFinalDetail(cpTbl, exLookup)
    .intSaveNm("locDetailClusterFinal", cpTbl, "all", stageChnl, pathProject)
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  if (!any(gateTblStim$locGenerated %in% TRUE)) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "no_generated_local_fdr_thresholds_for_tol_imputation"
    )
    cpTbl <- .getCpClusterLocAddFinalDetail(cpTbl, exLookup)
    .intSaveNm("locDetailClusterFinal", cpTbl, "all", stageChnl, pathProject)
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  commonBw <- .getCpClusterLocCommonBw(
    gateTblStim = gateTblStim,
    exLookup = exLookup,
    chnlSettings = chnlSettings
  )
  .intSaveNm("locClusterCommonBw", commonBw, "all", stageChnl, pathProject)

  exprRange <- .getCpClusterLocExprRange(exLookup)
  minThreshold <- .getCpClusterDensTblGetMinThreshold(
    gateTbl = gateTblStim,
    control = control
  )

  densTblUns <- .getCpClusterLocDensTbl(
    exLookup = exLookup,
    type = "uns",
    minThreshold = minThreshold,
    exprMin = exprRange[["min"]],
    exprMax = exprRange[["max"]],
    bw = commonBw,
    chnlCut = chnlSettings$chnlCut
  ) |>
    .getCpClusterLocDensTblAddGrp()

  densTblStim <- .getCpClusterLocDensTbl(
    exLookup = exLookup,
    type = "stim",
    minThreshold = minThreshold,
    exprMin = exprRange[["min"]],
    exprMax = exprRange[["max"]],
    bw = commonBw,
    chnlCut = chnlSettings$chnlCut
  ) |>
    .getCpClusterLocDensTblAddGrp()

  .intSave("all", stageChnl, pathProject, densTblUns, densTblStim)

  locTbl <- .getCpClusterLocTbl(
    gateTblStim = gateTblStim,
    densTblUns = densTblUns,
    densTblStim = densTblStim,
    exLookup = exLookup,
    commonBw = commonBw,
    cpMin = cpMin,
    chnlSettings = chnlSettings
  )
  .intSaveNm("locClusterTolTbl", locTbl, "all", stageChnl, pathProject)

  cpTbl <- purrr::map_df(seq_len(nrow(locTbl)), function(i) {
    .getCpClusterLocImputeRow(
      row = locTbl[i, , drop = FALSE],
      locTbl = locTbl,
      exLookup = exLookup,
      commonBw = commonBw,
      cpMin = cpMin,
      chnlSettings = chnlSettings
    )
  }) |>
    dplyr::arrange(ind)

  cpTbl <- .getCpClusterLocAddFinalDetail(cpTbl, exLookup)
  .intSaveNm("locDetailClusterFinal", cpTbl, "all", stageChnl, pathProject)
  .intSave("all", stageChnl, pathProject, cpTbl)
  cpTbl
}

#' @keywords internal
.getCpClusterLocGateTblPrepare <- function(gateTbl) {
  if (!"locGenerated" %in% names(gateTbl)) {
    gateTbl$locGenerated <- !is.na(gateTbl$gate)
  }
  if (!"locGeneratedDirect" %in% names(gateTbl)) {
    gateTbl$locGeneratedDirect <- gateTbl$locGenerated
  }
  if (!"locSource" %in% names(gateTbl)) {
    gateTbl$locSource <- NA_character_
  }
  if (!"locReason" %in% names(gateTbl)) {
    gateTbl$locReason <- NA_character_
  }
  gateTbl |>
    dplyr::mutate(
      ind = as.character(.data$ind),
      locGenerated = .data$locGenerated %in% TRUE,
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE
    )
}

#' @keywords internal
.getCpClusterLocGateTblStim <- function(gateTbl, exLookup) {
  stimInd <- names(exLookup)
  gateTbl |>
    dplyr::filter(.data$ind %in% .env$stimInd) |>
    dplyr::filter(!(.data$locSource %in% "unstim_summary"))
}

#' @keywords internal
.getCpClusterLocExLookup <- function(
  .data,
  indBatchList,
  chnlSettings,
  filterOtherCytPos,
  calcCytPosGates,
  gateTbl,
  pathProject
) {
  exPairs <- purrr::map(seq_along(indBatchList), function(i) {
    batch <- names(indBatchList)[i]
    exList <- .getExList(
      .data = .data,
      indBatch = indBatchList[[i]],
      pop = chnlSettings$popGate,
      chnlCut = chnlSettings$chnlCut,
      batch = batch,
      pathProject = pathProject
    )

    exListStim <- if (filterOtherCytPos) {
      .getCpClusterDensTblGetBatchPrepExListFilter(
        exList = exList,
        chnlCut = chnlSettings$chnlCut,
        gateTbl = gateTbl,
        calcCytPosGates = calcCytPosGates
      )
    } else {
      exList[-1]
    }

    purrr::map(names(exListStim), function(indCurr) {
      list(
        ind = as.character(indCurr),
        batch = batch,
        stim = exListStim[[indCurr]],
        uns = exList[[1]]
      )
    })
  }) |>
    purrr::flatten()
  stats::setNames(exPairs, purrr::map_chr(exPairs, "ind"))
}

#' @keywords internal
.getCpClusterLocCommonBw <- function(gateTblStim, exLookup, chnlSettings) {
  indGenerated <- gateTblStim |>
    dplyr::filter(.data$locGenerated %in% TRUE) |>
    dplyr::pull("ind") |>
    as.character()

  bwVec <- purrr::map_dbl(indGenerated, function(indCurr) {
    exPair <- exLookup[[indCurr]]
    if (is.null(exPair)) {
      return(NA_real_)
    }
    bwStim <- .getCpClusterLocBwOne(.getCut(exPair$stim), chnlSettings)
    bwUns <- .getCpClusterLocBwOne(.getCut(exPair$uns), chnlSettings)
    suppressWarnings(median(c(bwStim, bwUns), na.rm = TRUE))
  })
  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
  if (length(bwVec) > 0L) {
    return(median(bwVec, na.rm = TRUE))
  }

  bwFallback <- suppressWarnings(as.numeric(chnlSettings$bwCluster))[1]
  if (!is.finite(bwFallback) || bwFallback <= 0) {
    bwFallback <- suppressWarnings(as.numeric(chnlSettings$bw))[1]
  }
  if (is.finite(bwFallback) && bwFallback > 0) {
    return(bwFallback)
  }

  allExpr <- unlist(purrr::map(exLookup, function(x) {
    c(.getCut(x$stim), .getCut(x$uns))
  }))
  .getCpClusterLocBwOne(allExpr, chnlSettings)
}

#' @keywords internal
.getCpClusterLocBwOne <- function(x, chnlSettings) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 5L || length(unique(x)) < 3L) {
    return(NA_real_)
  }
  if (
    !is.null(chnlSettings$bwNcellMax) && length(x) > chnlSettings$bwNcellMax
  ) {
    x <- sample(x, size = chnlSettings$bwNcellMax, replace = FALSE)
  }
  if (
    !is.null(chnlSettings$bwNcellMin) && length(x) < chnlSettings$bwNcellMin
  ) {
    iqrX <- diff(stats::quantile(x, c(0.75, 0.25), na.rm = TRUE))
    sdX <- abs(iqrX) / 1.5
    x <- sample(x, replace = TRUE, size = chnlSettings$bwNcellMin) +
      stats::rnorm(chnlSettings$bwNcellMin, mean = 0, sd = sdX / 10)
  }

  bwMtd <- chnlSettings$bwMtd %||% "hpi1"
  bwAdj <- chnlSettings$bwAdj %||% 1
  bwMin <- suppressWarnings(as.numeric(chnlSettings$bwMin))[1]
  bwMax <- suppressWarnings(as.numeric(chnlSettings$bwMax))[1]
  if (!is.finite(bwMin) || bwMin <= 0) {
    bwMin <- .Machine$double.eps
  }
  if (!is.finite(bwMax) || bwMax <= 0) {
    bwMax <- Inf
  }

  bwCalc <- switch(
    bwMtd,
    "nrd0" = try(stats::bw.nrd0(x), silent = TRUE),
    "sj" = try(stats::bw.SJ(x), silent = TRUE),
    try(
      suppressWarnings(
        ks::hpi(x, deriv.order = as.numeric(gsub("hpi", "", bwMtd)))
      ),
      silent = TRUE
    )
  )
  if (inherits(bwCalc, "try-error") || !is.finite(bwCalc) || bwCalc <= 0) {
    bwCalc <- try(stats::bw.nrd0(x), silent = TRUE)
  }
  if (inherits(bwCalc, "try-error") || !is.finite(bwCalc) || bwCalc <= 0) {
    return(NA_real_)
  }
  max(bwMin, min(bwCalc * bwAdj, bwMax))
}

#' @keywords internal
.getCpClusterLocExprRange <- function(exLookup) {
  rangeTbl <- purrr::map_df(exLookup, function(exPair) {
    x <- c(.getCut(exPair$stim), .getCut(exPair$uns))
    x <- x[is.finite(x)]
    if (length(x) <= 5L) {
      return(NULL)
    }
    q <- stats::quantile(x, c(0.0025, 0.999), na.rm = TRUE)
    tibble::tibble(lb = q[[1]], ub = 3 * q[[2]])
  })
  if (nrow(rangeTbl) == 0L) {
    return(c(min = 0, max = 1))
  }
  c(
    min = stats::quantile(rangeTbl$lb, 0.0025, na.rm = TRUE)[[1]],
    max = max(rangeTbl$ub, na.rm = TRUE)
  )
}

#' @keywords internal
.getCpClusterLocDensTbl <- function(
  exLookup,
  type,
  minThreshold,
  exprMin,
  exprMax,
  bw,
  chnlCut
) {
  densTbl <- purrr::map_df(exLookup, function(exPair) {
    ex <- switch(type, stim = exPair$stim, uns = exPair$uns)
    .getCpClusterDensTblGetActualInd(
      exprVec = .getCut(ex),
      batch = exPair$batch,
      ind = exPair$ind,
      minThreshold = minThreshold,
      chnlCut = chnlCut,
      exprMin = exprMin,
      exprMax = exprMax,
      bw = bw
    )
  }) |>
    dplyr::filter(!is.na(x1)) # nolint

  if (nrow(densTbl) == 0L) {
    return(densTbl)
  }
  notAllNaVecInd <- purrr::map_lgl(
    seq_len(ncol(densTbl)),
    function(i) !all(is.na(densTbl[[i]]))
  )
  densTbl[, notAllNaVecInd, drop = FALSE]
}

#' @keywords internal
.getCpClusterLocDensTblAddGrp <- function(densTbl) {
  if (is.null(densTbl) || nrow(densTbl) == 0L) {
    return(tibble::tibble(ind = character(), grp = character()))
  }
  xCols <- grepl("^x\\d+", colnames(densTbl))
  if (!any(xCols) || nrow(densTbl) < 3L) {
    return(densTbl |> dplyr::mutate(grp = "1"))
  }
  densMat <- as.matrix(densTbl[, xCols, drop = FALSE])
  if (nrow(unique(densMat)) <= 1L) {
    return(densTbl |> dplyr::mutate(grp = "1"))
  }
  nClus <- .getCpClusterNClus(densTbl)
  clusVec <- .getCpClusterClus(densTbl, nClus)
  densTbl |> dplyr::mutate(grp = as.character(clusVec))
}

#' @keywords internal
.getCpClusterLocTbl <- function(
  gateTblStim,
  densTblUns,
  densTblStim,
  exLookup,
  commonBw,
  cpMin,
  chnlSettings
) {
  locTbl <- gateTblStim |>
    dplyr::left_join(
      densTblUns |> dplyr::select(ind, grpUns = grp) |> dplyr::distinct(),
      by = "ind"
    ) |>
    dplyr::left_join(
      densTblStim |> dplyr::select(ind, grpStim = grp) |> dplyr::distinct(),
      by = "ind"
    )

  tolTbl <- purrr::map_df(seq_len(nrow(locTbl)), function(i) {
    row <- locTbl[i, , drop = FALSE]
    indCurr <- as.character(row$ind[1])
    exPair <- exLookup[[indCurr]]
    if (is.null(exPair) || !(row$locGenerated[1] %in% TRUE)) {
      return(tibble::tibble(
        ind = indCurr,
        locTolSignedUns = NA_real_,
        locTolAbsUns = NA_real_,
        locLog10TolAbsUns = NA_real_,
        locDerivUns = NA_real_,
        locDerivSignUns = NA_real_,
        locTolSignedStim = NA_real_,
        locTolAbsStim = NA_real_,
        locLog10TolAbsStim = NA_real_,
        locDerivStim = NA_real_,
        locDerivSignStim = NA_real_
      ))
    }
    tolUns <- .getCpClusterLocTolAtCp(
      x = .getCut(exPair$uns),
      cp = row$gate[1],
      bw = commonBw,
      cpMin = cpMin,
      refPeak = chnlSettings$locTolRefPeak %||% "highest"
    )
    tolStim <- .getCpClusterLocTolAtCp(
      x = .getCut(exPair$stim),
      cp = row$gate[1],
      bw = commonBw,
      cpMin = cpMin,
      refPeak = chnlSettings$locTolRefPeak %||% "highest"
    )
    tibble::tibble(
      ind = indCurr,
      locTolSignedUns = tolUns$signedTol,
      locTolAbsUns = tolUns$tolAbs,
      locLog10TolAbsUns = tolUns$log10TolAbs,
      locDerivUns = tolUns$deriv,
      locDerivSignUns = tolUns$derivSign,
      locTolSignedStim = tolStim$signedTol,
      locTolAbsStim = tolStim$tolAbs,
      locLog10TolAbsStim = tolStim$log10TolAbs,
      locDerivStim = tolStim$deriv,
      locDerivSignStim = tolStim$derivSign
    )
  })

  locTbl |>
    dplyr::left_join(tolTbl, by = "ind")
}

#' @keywords internal
.getCpClusterLocTolAtCp <- function(
  x,
  cp,
  bw,
  cpMin = -Inf,
  n = 512,
  refPeak = c("highest", "first")
) {
  refPeak <- match.arg(refPeak)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  cp <- suppressWarnings(as.numeric(cp))[1]
  bw <- suppressWarnings(as.numeric(bw))[1]
  cpMin <- suppressWarnings(as.numeric(cpMin))[1]
  if (
    length(x) < 5L ||
      length(unique(x)) < 3L ||
      !is.finite(cp) ||
      !is.finite(bw) ||
      bw <= 0
  ) {
    return(.getCpClusterLocTolEmpty(cp = cp, bw = bw))
  }

  dens <- try(
    suppressWarnings(stats::density(x, bw = bw, n = n)),
    silent = TRUE
  )
  if (inherits(dens, "try-error")) {
    return(.getCpClusterLocTolEmpty(cp = cp, bw = bw))
  }

  derObj <- .getCpClusterLocDerivative(dens)
  peakInd <- .getCpClusterLocPeakInd(dens, refPeak = refPeak)
  rightInd <- derObj$x >= max(dens$x[peakInd], cpMin, na.rm = TRUE)
  if (!any(rightInd)) {
    return(.getCpClusterLocTolEmpty(cp = cp, bw = bw))
  }

  derRight <- derObj$deriv[rightInd]
  refDeriv <- derRight[which.max(abs(derRight))]
  if (!is.finite(refDeriv) || refDeriv == 0) {
    return(.getCpClusterLocTolEmpty(cp = cp, bw = bw))
  }

  derivAtCp <- stats::approx(
    x = derObj$x,
    y = derObj$deriv,
    xout = cp,
    rule = 2
  )$y
  tolAbs <- abs(derivAtCp) / abs(refDeriv)
  derivSign <- sign(derivAtCp)
  if (!is.finite(tolAbs) || tolAbs <= 0 || !is.finite(derivSign)) {
    return(.getCpClusterLocTolEmpty(cp = cp, bw = bw))
  }
  signedTol <- derivSign * tolAbs
  list(
    signedTol = signedTol,
    tolAbs = tolAbs,
    log10TolAbs = log10(tolAbs),
    deriv = derivAtCp,
    derivSign = derivSign,
    refDeriv = refDeriv,
    bw = bw,
    cp = cp
  )
}

#' @keywords internal
.getCpClusterLocTolEmpty <- function(cp, bw) {
  list(
    signedTol = NA_real_,
    tolAbs = NA_real_,
    log10TolAbs = NA_real_,
    deriv = NA_real_,
    derivSign = NA_real_,
    refDeriv = NA_real_,
    bw = bw,
    cp = cp
  )
}

#' @keywords internal
.getCpClusterLocDerivative <- function(dens) {
  list(
    x = dens$x[-1] - diff(dens$x) / 2,
    deriv = diff(dens$y) / diff(dens$x)
  )
}

#' @keywords internal
.getCpClusterLocPeakInd <- function(dens, refPeak = "highest") {
  if (refPeak == "highest") {
    return(which.max(dens$y))
  }
  peakInd <- which(
    dens$y[-c(1, length(dens$y))] >
      dens$y[-c(length(dens$y) - 1L, length(dens$y))] &
      dens$y[-c(1, length(dens$y))] >= dens$y[-c(1, 2)]
  ) +
    1L
  if (length(peakInd) == 0L) {
    return(which.max(dens$y))
  }
  peakInd[1]
}

#' @keywords internal
.getCpClusterLocImputeRow <- function(
  row,
  locTbl,
  exLookup,
  commonBw,
  cpMin,
  chnlSettings
) {
  indCurr <- as.character(row$ind[1])
  cpOrig <- as.numeric(row$gate[1])
  cpMedianUns <- .getCpClusterLocMedianCp(
    locTbl = locTbl,
    grpCol = "grpUns",
    grpVal = row$grpUns[1]
  )
  cpMedianStim <- .getCpClusterLocMedianCp(
    locTbl = locTbl,
    grpCol = "grpStim",
    grpVal = row$grpStim[1]
  )

  if (row$locGenerated[1] %in% TRUE) {
    return(.getCpClusterLocRowOut(
      row = row,
      cp = cpOrig,
      cpUns = NA_real_,
      cpStim = NA_real_,
      cpMedianUns = cpMedianUns,
      cpMedianStim = cpMedianStim,
      commonBw = commonBw,
      reason = "local_fdr_available"
    ))
  }

  exPair <- exLookup[[indCurr]]
  if (is.null(exPair)) {
    return(.getCpClusterLocRowOut(
      row = row,
      cp = cpOrig,
      cpUns = NA_real_,
      cpStim = NA_real_,
      cpMedianUns = cpMedianUns,
      cpMedianStim = cpMedianStim,
      commonBw = commonBw,
      reason = "expression_pair_missing"
    ))
  }

  tolUns <- .getCpClusterLocMedianSignedTol(
    locTbl = locTbl,
    grpCol = "grpUns",
    grpVal = row$grpUns[1],
    tolCol = "locTolSignedUns"
  )
  tolStim <- .getCpClusterLocMedianSignedTol(
    locTbl = locTbl,
    grpCol = "grpStim",
    grpVal = row$grpStim[1],
    tolCol = "locTolSignedStim"
  )

  cpUns <- .getCpClusterLocCpFromSignedTol(
    x = .getCut(exPair$uns),
    signedTol = tolUns,
    bw = commonBw,
    cpMin = cpMin,
    refPeak = chnlSettings$locTolRefPeak %||% "highest"
  )
  cpStim <- .getCpClusterLocCpFromSignedTol(
    x = .getCut(exPair$stim),
    signedTol = tolStim,
    bw = commonBw,
    cpMin = cpMin,
    refPeak = chnlSettings$locTolRefPeak %||% "highest"
  )

  candUns <- .getCpClusterLocCandidateBound(cpUns, cpMedianUns)
  candStim <- .getCpClusterLocCandidateBound(cpStim, cpMedianStim)
  cp <- .getCpClusterLocChooseCp(
    cpUns = candUns,
    cpStim = candStim,
    cpMedianUns = cpMedianUns,
    cpMedianStim = cpMedianStim,
    cpOrig = cpOrig
  )

  .getCpClusterLocRowOut(
    row = row,
    cp = cp,
    cpUns = candUns,
    cpStim = candStim,
    cpMedianUns = cpMedianUns,
    cpMedianStim = cpMedianStim,
    commonBw = commonBw,
    reason = "imputed_from_cluster_median_signed_tol"
  )
}

#' @keywords internal
.getCpClusterLocMedianCp <- function(locTbl, grpCol, grpVal) {
  if (is.na(grpVal)) {
    return(NA_real_)
  }
  x <- locTbl |>
    dplyr::filter(.data[[grpCol]] == .env$grpVal) |>
    dplyr::filter(.data$locGenerated %in% TRUE) |>
    dplyr::pull("gate")
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

#' @keywords internal
.getCpClusterLocMedianSignedTol <- function(locTbl, grpCol, grpVal, tolCol) {
  if (is.na(grpVal)) {
    return(NA_real_)
  }
  tolTbl <- locTbl |>
    dplyr::filter(.data[[grpCol]] == .env$grpVal) |>
    dplyr::filter(.data$locGenerated %in% TRUE)
  x <- tolTbl[[tolCol]]
  x <- x[is.finite(x) & x != 0]
  if (length(x) == 0L) {
    return(NA_real_)
  }

  signVec <- sign(x)
  tab <- table(signVec)
  signUse <- as.numeric(names(tab)[which.max(tab)])
  xUse <- x[signVec == signUse]
  if (length(xUse) == 0L) {
    signUse <- sign(stats::median(x, na.rm = TRUE))
    xUse <- x
  }
  signUse * 10^stats::median(log10(abs(xUse)), na.rm = TRUE)
}

#' @keywords internal
.getCpClusterLocCpFromSignedTol <- function(
  x,
  signedTol,
  bw,
  cpMin = -Inf,
  n = 512,
  refPeak = c("highest", "first")
) {
  refPeak <- match.arg(refPeak)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  signedTol <- suppressWarnings(as.numeric(signedTol))[1]
  bw <- suppressWarnings(as.numeric(bw))[1]
  cpMin <- suppressWarnings(as.numeric(cpMin))[1]
  if (
    length(x) < 5L ||
      length(unique(x)) < 3L ||
      !is.finite(signedTol) ||
      signedTol == 0 ||
      !is.finite(bw) ||
      bw <= 0
  ) {
    return(NA_real_)
  }

  dens <- try(
    suppressWarnings(stats::density(x, bw = bw, n = n)),
    silent = TRUE
  )
  if (inherits(dens, "try-error")) {
    return(NA_real_)
  }

  derObj <- .getCpClusterLocDerivative(dens)
  peakInd <- .getCpClusterLocPeakInd(dens, refPeak = refPeak)
  rightInd <- derObj$x >= max(dens$x[peakInd], cpMin, na.rm = TRUE)
  if (!any(rightInd)) {
    return(NA_real_)
  }

  xRight <- derObj$x[rightInd]
  derRight <- derObj$deriv[rightInd]
  refInd <- which.max(abs(derRight))
  refDeriv <- derRight[refInd]
  if (!is.finite(refDeriv) || refDeriv == 0) {
    return(NA_real_)
  }

  tolAbs <- abs(signedTol)
  targetSign <- sign(signedTol)
  ratio <- abs(derRight) / abs(refDeriv)
  afterRef <- seq_along(xRight) >= refInd
  sameSign <- sign(derRight) == targetSign

  if (targetSign < 0) {
    cand <- which(afterRef & sameSign & ratio <= tolAbs)
    if (length(cand) > 0L) {
      return(max(min(xRight[cand[1]], na.rm = TRUE), cpMin, na.rm = TRUE))
    }
  }

  cand <- which(afterRef & sameSign)
  if (length(cand) == 0L) {
    cand <- which(afterRef)
  }
  if (length(cand) == 0L) {
    return(NA_real_)
  }
  best <- cand[which.min(abs(ratio[cand] - tolAbs))]
  max(xRight[best], cpMin, na.rm = TRUE)
}

#' @keywords internal
.getCpClusterLocCandidateBound <- function(cp, cpMedian) {
  if (!is.finite(cp)) {
    return(NA_real_)
  }
  if (!is.finite(cpMedian)) {
    return(cp)
  }
  max(cp, cpMedian)
}

#' @keywords internal
.getCpClusterLocChooseCp <- function(
  cpUns,
  cpStim,
  cpMedianUns,
  cpMedianStim,
  cpOrig
) {
  cand <- tibble::tibble(
    cp = c(cpUns, cpStim),
    cpMedian = c(cpMedianUns, cpMedianStim),
    source = c("uns", "stim")
  ) |>
    dplyr::filter(is.finite(.data$cp)) |>
    dplyr::mutate(
      dist = ifelse(is.finite(.data$cpMedian), .data$cp - .data$cpMedian, Inf),
      dist = ifelse(.data$dist < 0, Inf, .data$dist)
    )
  if (nrow(cand) == 0L) {
    return(cpOrig)
  }
  cand |>
    dplyr::arrange(.data$dist) |>
    dplyr::slice(1) |>
    dplyr::pull("cp")
}


#' @keywords internal
.getCpClusterLocThresholdOrigin <- function(row) {
  reason <- as.character(row$locClusterReason[1] %||% NA_character_)
  locSource <- as.character(row$locSource[1] %||% NA_character_)
  locGenerated <- row$locGenerated[1] %in% TRUE
  locGeneratedDirect <- row$locGeneratedDirect[1] %in% TRUE

  if (reason %in% "imputed_from_cluster_median_signed_tol") {
    return("cluster_imputed_from_similar_conditions")
  }
  if (locGenerated && locGeneratedDirect && locSource %in% "direct") {
    return("condition_detected_response")
  }
  if (locGenerated && locSource %in% "combined") {
    return("sample_imputed_from_other_stim_conditions")
  }
  if (locGenerated && locSource %in% "prejoin") {
    return("prejoin_generated_from_joined_stim_conditions")
  }
  if (reason %in% "local_fdr_available") {
    return("local_fdr_available")
  }
  if (!locGenerated) {
    return("not_generated_fallback")
  }
  "generated_unknown_source"
}

#' @keywords internal
.getCpClusterLocPropBsAtCp <- function(cp, exPair) {
  cp <- suppressWarnings(as.numeric(cp))[1]
  if (
    !is.finite(cp) ||
      is.null(exPair) ||
      is.null(exPair$stim) || is.null(exPair$uns) ||
      !is.data.frame(exPair$stim) || !is.data.frame(exPair$uns) ||
      nrow(exPair$stim) == 0L || nrow(exPair$uns) == 0L
  ) {
    return(tibble::tibble(
      locFinalNCellStim = if (!is.null(exPair$stim) && is.data.frame(exPair$stim)) nrow(exPair$stim) else NA_integer_,
      locFinalNCellUns = if (!is.null(exPair$uns) && is.data.frame(exPair$uns)) nrow(exPair$uns) else NA_integer_,
      locFinalPropStim = NA_real_,
      locFinalPropUns = NA_real_,
      locFinalPropBs = NA_real_
    ))
  }
  propStim <- sum(.getCut(exPair$stim) >= cp, na.rm = TRUE) / nrow(exPair$stim)
  propUns <- sum(.getCut(exPair$uns) >= cp, na.rm = TRUE) / nrow(exPair$uns)
  tibble::tibble(
    locFinalNCellStim = nrow(exPair$stim),
    locFinalNCellUns = nrow(exPair$uns),
    locFinalPropStim = propStim,
    locFinalPropUns = propUns,
    locFinalPropBs = propStim - propUns
  )
}

#' @keywords internal
.getCpClusterLocAddFinalDetail <- function(cpTbl, exLookup) {
  if (!is.data.frame(cpTbl) || nrow(cpTbl) == 0L) {
    return(cpTbl)
  }
  detailTbl <- purrr::map_df(seq_len(nrow(cpTbl)), function(i) {
    row <- cpTbl[i, , drop = FALSE]
    cp <- suppressWarnings(as.numeric(row$cpJoinLseOrigMeanTg[1]))
    if (!is.finite(cp)) {
      cp <- suppressWarnings(as.numeric(row$cpJoinTgOrig[1]))
    }
    exPair <- exLookup[[as.character(row$ind[1])]]
    .getCpClusterLocPropBsAtCp(cp = cp, exPair = exPair) |>
      dplyr::mutate(
        locFinalThreshold = cp,
        locFinalThresholdOrigin = .getCpClusterLocThresholdOrigin(row),
        locFinalThresholdGenerated = is.finite(cp)
      )
  })
  dplyr::bind_cols(cpTbl, detailTbl)
}

#' @keywords internal
.getCpClusterLocSkipOut <- function(gateTblStim, reason) {
  purrr::map_df(seq_len(nrow(gateTblStim)), function(i) {
    .getCpClusterLocRowOut(
      row = gateTblStim[i, , drop = FALSE],
      cp = gateTblStim$gate[i],
      cpUns = NA_real_,
      cpStim = NA_real_,
      cpMedianUns = NA_real_,
      cpMedianStim = NA_real_,
      commonBw = NA_real_,
      reason = reason
    )
  })
}

#' @keywords internal
.getCpClusterLocRowOut <- function(
  row,
  cp,
  cpUns,
  cpStim,
  cpMedianUns,
  cpMedianStim,
  commonBw,
  reason
) {
  grpUns <- suppressWarnings(row$grpUns[1]) %||% NA_character_
  grpStim <- suppressWarnings(row$grpStim[1]) %||% NA_character_
  grp <- ifelse(!is.na(grpStim), grpStim, grpUns)
  cpJoin <- ifelse(is.finite(cpMedianStim), cpMedianStim, cpMedianUns)
  tibble::tibble(
    grp = grp,
    grpUns = grpUns,
    grpStim = grpStim,
    ind = as.character(row$ind[1]),
    cpOrigQuantMin = row$gate[1],
    cpJoin = cpJoin,
    cpJoinLse = NA_real_,
    cpJoinLseOrig = NA_real_,
    cpJoinLseOrigMean = NA_real_,
    cpJoinTgOrig = cp,
    cpJoinTgOrigMean = cp,
    cpJoinLseOrigMeanTg = cp,
    cpTolUns = cpUns,
    cpTolStim = cpStim,
    cpMedianUns = cpMedianUns,
    cpMedianStim = cpMedianStim,
    locGenerated = suppressWarnings(row$locGenerated[1] %in% TRUE),
    locGeneratedDirect = suppressWarnings(row$locGeneratedDirect[1] %in% TRUE),
    locSource = suppressWarnings(row$locSource[1]) %||% NA_character_,
    locReason = suppressWarnings(row$locReason[1]) %||% NA_character_,
    locClusterReason = reason,
    locClusterBw = commonBw,
    locTolSignedUns = suppressWarnings(row$locTolSignedUns[1]) %||% NA_real_,
    locTolSignedStim = suppressWarnings(row$locTolSignedStim[1]) %||% NA_real_,
    locDerivSignUns = suppressWarnings(row$locDerivSignUns[1]) %||% NA_real_,
    locDerivSignStim = suppressWarnings(row$locDerivSignStim[1]) %||% NA_real_,
    propBsOrig = NA_real_,
    propBsCpDiff = NA_real_,
    propBsCpDiffSd = NA_real_,
    propBsCp = NA_real_
  )
}
