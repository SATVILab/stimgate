# Cluster local-FDR thresholds using the joint stimulated/unstimulated
# left-region density shapes. Directly generated thresholds are the only donors.
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
  stageChnl <- file.path(stage, chnlSettings$chnlCut)
  control <- .getCpClusterControlUpdate(control)

  gateTbl <- .getCpClusterLocGateTblPrepare(gateTbl)
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
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  directAvailable <- gateTblStim$locGeneratedDirect %in% TRUE &
    is.finite(gateTblStim$gate)
  if (!any(directAvailable)) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "no_direct_local_fdr_threshold_donors"
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  commonBw <- .getCpClusterLocCommonBw(
    gateTblStim = gateTblStim,
    exLookup = exLookup,
    chnlSettings = chnlSettings
  )
  .intSaveNm("locClusterCommonBw", commonBw, "all", stageChnl, pathProject)

  if (!is.finite(commonBw) || commonBw <= 0) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "common_density_bandwidth_unavailable",
      commonBw = commonBw
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  featureRange <- .getCpClusterLocFeatureRange(
    gateTblStim = gateTblStim,
    exLookup = exLookup,
    control = control
  )
  featureObj <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    indKeep = gateTblStim$ind,
    exprMin = featureRange[["min"]],
    exprMax = featureRange[["max"]],
    bw = commonBw,
    n = control$locClusterDensityN
  )
  featureTbl <- featureObj$featureTbl

  .intSaveNm(
    "locClusterDensityGrid",
    featureObj$x,
    "all",
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "locClusterJointDensityTbl",
    featureTbl,
    "all",
    stageChnl,
    pathProject
  )

  if (nrow(featureTbl) == 0L) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "no_complete_joint_density_features",
      commonBw = commonBw,
      featureMinX = featureRange[["min"]],
      featureMaxX = featureRange[["max"]]
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  locTblFeature <- gateTblStim |>
    dplyr::inner_join(featureTbl, by = "ind")
  nDirectFeature <- sum(
    locTblFeature$locGeneratedDirect %in% TRUE &
      is.finite(locTblFeature$gate)
  )

  if (nDirectFeature < control$locClusterMinDirect) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "fewer_than_minimum_direct_donors_with_density_features",
      commonBw = commonBw,
      featureMinX = featureRange[["min"]],
      featureMaxX = featureRange[["max"]]
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  grpObj <- .getCpClusterLocJointClusters(
    locTblFeature = locTblFeature,
    control = control
  )
  grpTbl <- tibble::tibble(
    ind = locTblFeature$ind,
    grpInitial = as.character(grpObj$initial),
    grp = as.character(grpObj$final)
  )

  locTbl <- gateTblStim |>
    dplyr::left_join(grpTbl, by = "ind")
  clusterSummary <- .getCpClusterLocDirectSummary(
    locTbl = locTbl,
    control = control
  )
  locTbl <- locTbl |>
    dplyr::left_join(clusterSummary, by = "grp")

  .intSaveNm("locClusterThresholdTbl", locTbl, "all", stageChnl, pathProject)

  cpTbl <- purrr::map_df(seq_len(nrow(locTbl)), function(i) {
    .getCpClusterLocFinaliseRow(
      row = locTbl[i, , drop = FALSE],
      commonBw = commonBw,
      featureMinX = featureRange[["min"]],
      featureMaxX = featureRange[["max"]]
    )
  }) |>
    dplyr::arrange(ind)

  .intSave("all", stageChnl, pathProject, cpTbl)
  cpTbl
}

#' @keywords internal
.getCpClusterLocGateTblPrepare <- function(gateTbl) {
  if (!"locGenerated" %in% names(gateTbl)) {
    gateTbl$locGenerated <- is.finite(gateTbl$gate)
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
      locGenerated = .data$locGenerated %in% TRUE &
        is.finite(.data$gate),
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE &
        is.finite(.data$gate)
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
      .getCpClusterLocExListFilter(
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
.getCpClusterLocExListFilter <- function(
  exList,
  chnlCut,
  gateTbl,
  calcCytPosGates
) {
  if (length(exList) <= 1L) {
    return(exList[0])
  }

  purrr::map(exList[-1], function(exTbl) {
    indCurr <- as.character(attr(exTbl, "ind"))
    gateTblInd <- gateTbl[gateTbl$ind == indCurr, , drop = FALSE]
    posInd <- .get_pos_ind_but_single_pos_for_one_cyt(
      ex = exTbl,
      gateTbl = gateTblInd,
      chnlSingleExc = chnlCut,
      chnl = NULL,
      gateTypeCytPos = if (calcCytPosGates) "cyt" else "base",
      gateTypeSinglePos = "base"
    )
    exTbl[!posInd, , drop = FALSE]
  })
}

#' @keywords internal
.getCpClusterLocCommonBw <- function(gateTblStim, exLookup, chnlSettings) {
  indDirect <- gateTblStim |>
    dplyr::filter(.data$locGeneratedDirect %in% TRUE) |>
    dplyr::pull("ind") |>
    as.character()

  bwVec <- purrr::map_dbl(indDirect, function(indCurr) {
    exPair <- exLookup[[indCurr]]
    if (is.null(exPair)) {
      return(NA_real_)
    }
    bwStim <- .getCpClusterLocBwOne(.getCut(exPair$stim), chnlSettings)
    bwUns <- .getCpClusterLocBwOne(.getCut(exPair$uns), chnlSettings)
    suppressWarnings(stats::median(c(bwStim, bwUns), na.rm = TRUE))
  })
  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
  if (length(bwVec) > 0L) {
    return(stats::median(bwVec, na.rm = TRUE))
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

  bwMtd <- chnlSettings$bwMtd %||% "hpi1"
  bwAdj <- chnlSettings$bwAdj %||% 1
  bwMin <- suppressWarnings(as.numeric(chnlSettings$bwMin))[1]
  bwMax <- suppressWarnings(as.numeric(chnlSettings$bwMax))[1]
  bwFallback <- suppressWarnings(as.numeric(chnlSettings$bwFallback))[1]

  if (!is.finite(bwMin)) {
    bwMin <- .Machine$double.eps
  }
  if (!is.finite(bwMax) || bwMax <= 0) {
    bwMax <- Inf
  }

  bwCalc <- .bwCalcOne(
    x = x,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax
  )

  if (!is.finite(bwCalc) || bwCalc <= 0) {
    if (is.finite(bwFallback) && bwFallback > 0) {
      return(max(bwMin, min(bwFallback, bwMax)))
    }
    return(NA_real_)
  }

  max(bwMin, min(bwCalc, bwMax))
}

#' @keywords internal
.getCpClusterLocFeatureRange <- function(gateTblStim, exLookup, control) {
  directGate <- gateTblStim$gate[
    gateTblStim$locGeneratedDirect %in% TRUE &
      is.finite(gateTblStim$gate)
  ]
  featureMax <- stats::quantile(
    directGate,
    probs = control$locClusterLeftQuantile,
    na.rm = TRUE,
    names = FALSE
  )
  featureMax <- control$locClusterLeftFraction * featureMax

  lowExpr <- purrr::map_dbl(exLookup, function(exPair) {
    x <- c(.getCut(exPair$stim), .getCut(exPair$uns))
    x <- x[is.finite(x)]
    if (length(x) < 5L) {
      return(NA_real_)
    }
    stats::quantile(x, 0.0025, na.rm = TRUE, names = FALSE)
  })
  lowExpr <- lowExpr[is.finite(lowExpr)]
  featureMin <- if (length(lowExpr) > 0L) {
    stats::quantile(lowExpr, 0.0025, na.rm = TRUE, names = FALSE)
  } else {
    NA_real_
  }

  if (!is.finite(featureMin) || !is.finite(featureMax) ||
      featureMin >= featureMax) {
    allExpr <- unlist(purrr::map(exLookup, function(exPair) {
      c(.getCut(exPair$stim), .getCut(exPair$uns))
    }))
    allExpr <- allExpr[is.finite(allExpr)]
    if (length(allExpr) >= 5L) {
      featureMin <- stats::quantile(
        allExpr,
        0.0025,
        na.rm = TRUE,
        names = FALSE
      )
      upperFallback <- stats::quantile(
        allExpr,
        0.75,
        na.rm = TRUE,
        names = FALSE
      )
      if (!is.finite(featureMax) || featureMax <= featureMin) {
        featureMax <- upperFallback
      }
    }
  }

  c(min = featureMin, max = featureMax)
}

#' @keywords internal
.getCpClusterLocJointFeatureTbl <- function(
  exLookup,
  indKeep,
  exprMin,
  exprMax,
  bw,
  n = 128L
) {
  n <- max(16L, as.integer(n)[1])
  if (!is.finite(exprMin) || !is.finite(exprMax) || exprMin >= exprMax) {
    return(list(
      featureTbl = tibble::tibble(ind = character()),
      x = numeric()
    ))
  }

  xGrid <- seq(exprMin, exprMax, length.out = n)
  featureRows <- purrr::map(indKeep, function(indCurr) {
    exPair <- exLookup[[as.character(indCurr)]]
    if (is.null(exPair)) {
      return(NULL)
    }
    densUns <- .getCpClusterLocDensityFeature(
      x = .getCut(exPair$uns),
      xGrid = xGrid,
      bw = bw
    )
    densStim <- .getCpClusterLocDensityFeature(
      x = .getCut(exPair$stim),
      xGrid = xGrid,
      bw = bw
    )
    if (is.null(densUns) || is.null(densStim)) {
      return(NULL)
    }
    c(
      list(ind = as.character(indCurr)),
      stats::setNames(as.list(densUns), sprintf("uns%03d", seq_len(n))),
      stats::setNames(as.list(densStim), sprintf("stim%03d", seq_len(n)))
    )
  }) |>
    purrr::compact()

  if (length(featureRows) == 0L) {
    return(list(
      featureTbl = tibble::tibble(ind = character()),
      x = xGrid
    ))
  }

  list(
    featureTbl = dplyr::bind_rows(featureRows),
    x = xGrid
  )
}

#' @keywords internal
.getCpClusterLocDensityFeature <- function(x, xGrid, bw) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 5L || length(unique(x)) < 3L) {
    return(NULL)
  }

  xAboveMin <- x[x > min(x)]
  if (length(xAboveMin) >= 5L && length(unique(xAboveMin)) >= 3L) {
    x <- xAboveMin
  }

  dens <- try(
    suppressWarnings(stats::density(
      x,
      bw = bw,
      n = length(xGrid),
      from = min(xGrid),
      to = max(xGrid)
    )),
    silent = TRUE
  )
  if (inherits(dens, "try-error")) {
    return(NULL)
  }

  y <- as.numeric(dens$y)
  y[!is.finite(y) | y < 0] <- 0
  ySum <- sum(y)
  if (!is.finite(ySum) || ySum <= 0) {
    return(NULL)
  }
  y / ySum
}

#' @keywords internal
.getCpClusterLocFeatureCols <- function(tbl) {
  grep("^(uns|stim)\\d{3}$", names(tbl), value = TRUE)
}

#' @keywords internal
.getCpClusterLocJointClusters <- function(locTblFeature, control) {
  featureCols <- .getCpClusterLocFeatureCols(locTblFeature)
  featureMat <- as.matrix(locTblFeature[, featureCols, drop = FALSE])
  storage.mode(featureMat) <- "double"

  direct <- locTblFeature$locGeneratedDirect %in% TRUE &
    is.finite(locTblFeature$gate)
  nDirect <- sum(direct)
  nDistinct <- nrow(unique(as.data.frame(featureMat)))
  kMax <- min(
    control$locClusterMaxClusters,
    floor(nDirect / control$locClusterMinDirect),
    nrow(featureMat),
    nDistinct
  )
  kMax <- max(1L, as.integer(kMax))

  if (kMax == 1L) {
    grp <- rep(1L, nrow(featureMat))
    return(list(initial = grp, final = grp, selectedK = 1L))
  }

  selectedK <- .getCpClusterLocSelectK(
    featureMat = featureMat,
    kMax = kMax,
    control = control
  )
  initial <- .getCpClusterLocKmeans(
    featureMat = featureMat,
    k = selectedK,
    control = control
  )
  final <- .getCpClusterLocMergeSparseDirect(
    featureMat = featureMat,
    grp = initial,
    direct = direct,
    minDirect = control$locClusterMinDirect
  )

  list(initial = initial, final = final, selectedK = selectedK)
}

#' @keywords internal
.getCpClusterLocSelectK <- function(featureMat, kMax, control) {
  if (kMax <= 1L) {
    return(1L)
  }

  clusterFun <- function(x, k) {
    stats::kmeans(
      x,
      centers = k,
      nstart = control$locClusterNstart,
      iter.max = 100
    )
  }
  gapObj <- .getCpClusterLocWithSeed(
    seed = control$locClusterSeed,
    code = try(
      cluster::clusGap(
        featureMat,
        FUNcluster = clusterFun,
        K.max = kMax,
        B = control$locClusterGapB
      ),
      silent = TRUE
    )
  )
  if (inherits(gapObj, "try-error")) {
    return(1L)
  }

  selected <- try(
    cluster::maxSE(
      gapObj$Tab[, "gap"],
      gapObj$Tab[, "SE.sim"],
      method = "firstSEmax"
    ),
    silent = TRUE
  )
  if (inherits(selected, "try-error") || !is.finite(selected)) {
    return(1L)
  }
  max(1L, min(as.integer(selected), kMax))
}

#' @keywords internal
.getCpClusterLocKmeans <- function(featureMat, k, control) {
  if (k <= 1L) {
    return(rep(1L, nrow(featureMat)))
  }
  fit <- .getCpClusterLocWithSeed(
    seed = control$locClusterSeed,
    code = stats::kmeans(
      featureMat,
      centers = k,
      nstart = control$locClusterNstart,
      iter.max = 100
    )
  )
  as.integer(fit$cluster)
}

#' @keywords internal
.getCpClusterLocWithSeed <- function(seed, code) {
  hadSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (hadSeed) {
    oldSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (hadSeed) {
      assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(as.integer(seed)[1])
  force(code)
}

#' @keywords internal
.getCpClusterLocMergeSparseDirect <- function(
  featureMat,
  grp,
  direct,
  minDirect
) {
  grp <- as.integer(factor(grp))
  repeat {
    grpLevels <- sort(unique(grp))
    if (length(grpLevels) <= 1L) {
      break
    }
    directCount <- vapply(grpLevels, function(g) {
      sum(direct[grp == g])
    }, integer(1))
    if (all(directCount >= minDirect)) {
      break
    }

    sourceGrp <- grpLevels[which.min(directCount)]
    sourceCentre <- colMeans(
      featureMat[grp == sourceGrp, , drop = FALSE]
    )
    targetGrp <- setdiff(grpLevels, sourceGrp)
    targetDist <- vapply(targetGrp, function(g) {
      targetCentre <- colMeans(featureMat[grp == g, , drop = FALSE])
      sqrt(sum((sourceCentre - targetCentre)^2))
    }, numeric(1))
    mergeInto <- targetGrp[which.min(targetDist)]
    grp[grp == sourceGrp] <- mergeInto
    grp <- as.integer(factor(grp))
  }
  grp
}

#' @keywords internal
.getCpClusterLocDirectSummary <- function(locTbl, control) {
  locTbl |>
    dplyr::filter(!is.na(.data$grp)) |>
    dplyr::group_by(.data$grp) |>
    dplyr::summarise(
      locClusterN = dplyr::n(),
      locClusterNDirect = sum(
        .data$locGeneratedDirect %in% TRUE &
          is.finite(.data$gate)
      ),
      locClusterQ15 = .getCpClusterLocQuantile(
        .data$gate[
          .data$locGeneratedDirect %in% TRUE &
            is.finite(.data$gate)
        ],
        control$locClusterWinsorLower
      ),
      locClusterQ60 = .getCpClusterLocQuantile(
        .data$gate[
          .data$locGeneratedDirect %in% TRUE &
            is.finite(.data$gate)
        ],
        control$locClusterImputeQuantile
      ),
      locClusterQ85 = .getCpClusterLocQuantile(
        .data$gate[
          .data$locGeneratedDirect %in% TRUE &
            is.finite(.data$gate)
        ],
        control$locClusterWinsorUpper
      ),
      .groups = "drop"
    )
}

#' @keywords internal
.getCpClusterLocQuantile <- function(x, prob) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  stats::quantile(
    x,
    probs = prob,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
}

#' @keywords internal
.getCpClusterLocFinaliseRow <- function(
  row,
  commonBw,
  featureMinX,
  featureMaxX
) {
  cpOrig <- suppressWarnings(as.numeric(row$gate[1]))
  q15 <- suppressWarnings(as.numeric(row$locClusterQ15[1]))
  q60 <- suppressWarnings(as.numeric(row$locClusterQ60[1]))
  q85 <- suppressWarnings(as.numeric(row$locClusterQ85[1]))
  direct <- row$locGeneratedDirect[1] %in% TRUE && is.finite(cpOrig)

  if (!is.finite(q15) || !is.finite(q60) || !is.finite(q85)) {
    return(.getCpClusterLocRowOut(
      row = row,
      cp = cpOrig,
      commonBw = commonBw,
      reason = if (is.na(row$grp[1])) {
        "joint_density_features_unavailable"
      } else {
        "cluster_direct_donor_summary_unavailable"
      },
      featureMinX = featureMinX,
      featureMaxX = featureMaxX
    ))
  }

  if (!direct) {
    return(.getCpClusterLocRowOut(
      row = row,
      cp = q60,
      commonBw = commonBw,
      reason = "assigned_cluster_direct_q60",
      featureMinX = featureMinX,
      featureMaxX = featureMaxX
    ))
  }

  cp <- min(max(cpOrig, q15), q85)
  reason <- if (cpOrig < q15) {
    "direct_winsorised_to_cluster_q15"
  } else if (cpOrig > q85) {
    "direct_winsorised_to_cluster_q85"
  } else {
    "direct_within_cluster_winsor_limits"
  }
  .getCpClusterLocRowOut(
    row = row,
    cp = cp,
    commonBw = commonBw,
    reason = reason,
    featureMinX = featureMinX,
    featureMaxX = featureMaxX
  )
}

#' @keywords internal
.getCpClusterLocSkipOut <- function(
  gateTblStim,
  reason,
  commonBw = NA_real_,
  featureMinX = NA_real_,
  featureMaxX = NA_real_
) {
  purrr::map_df(seq_len(nrow(gateTblStim)), function(i) {
    .getCpClusterLocRowOut(
      row = gateTblStim[i, , drop = FALSE],
      cp = gateTblStim$gate[i],
      commonBw = commonBw,
      reason = reason,
      featureMinX = featureMinX,
      featureMaxX = featureMaxX
    )
  })
}

#' @keywords internal
.getCpClusterLocRowOut <- function(
  row,
  cp,
  commonBw,
  reason,
  featureMinX = NA_real_,
  featureMaxX = NA_real_
) {
  grp <- .getCpClusterLocRowValue(row, "grp", NA_character_)
  grpInitial <- .getCpClusterLocRowValue(row, "grpInitial", NA_character_)
  q15 <- .getCpClusterLocRowValue(row, "locClusterQ15", NA_real_)
  q60 <- .getCpClusterLocRowValue(row, "locClusterQ60", NA_real_)
  q85 <- .getCpClusterLocRowValue(row, "locClusterQ85", NA_real_)
  nCluster <- .getCpClusterLocRowValue(row, "locClusterN", NA_integer_)
  nDirect <- .getCpClusterLocRowValue(
    row,
    "locClusterNDirect",
    NA_integer_
  )

  tibble::tibble(
    grp = as.character(grp),
    grpUns = as.character(grp),
    grpStim = as.character(grp),
    ind = as.character(row$ind[1]),
    cpOrigQuantMin = suppressWarnings(as.numeric(row$gate[1])),
    cpJoin = as.numeric(q60),
    cpJoinLse = NA_real_,
    cpJoinLseOrig = NA_real_,
    cpJoinLseOrigMean = NA_real_,
    cpJoinTgOrig = cp,
    cpJoinTgOrigMean = cp,
    cpJoinLseOrigMeanTg = cp,
    cpTolUns = NA_real_,
    cpTolStim = NA_real_,
    cpMedianUns = as.numeric(q60),
    cpMedianStim = as.numeric(q60),
    locGenerated = row$locGenerated[1] %in% TRUE,
    locGeneratedDirect = row$locGeneratedDirect[1] %in% TRUE,
    locSource = .getCpClusterLocRowValue(row, "locSource", NA_character_),
    locReason = .getCpClusterLocRowValue(row, "locReason", NA_character_),
    locClusterReason = reason,
    locClusterBw = commonBw,
    locClusterInitialGrp = as.character(grpInitial),
    locClusterN = as.integer(nCluster),
    locClusterNDirect = as.integer(nDirect),
    locClusterQ15 = as.numeric(q15),
    locClusterQ60 = as.numeric(q60),
    locClusterQ85 = as.numeric(q85),
    locClusterFeatureMinX = featureMinX,
    locClusterFeatureMaxX = featureMaxX,
    locTolSignedUns = NA_real_,
    locTolSignedStim = NA_real_,
    locDerivSignUns = NA_real_,
    locDerivSignStim = NA_real_,
    propBsOrig = NA_real_,
    propBsCpDiff = NA_real_,
    propBsCpDiffSd = NA_real_,
    propBsCp = NA_real_
  )
}

#' @keywords internal
.getCpClusterLocRowValue <- function(row, name, default) {
  if (!name %in% names(row) || length(row[[name]]) == 0L) {
    return(default)
  }
  value <- row[[name]][1]
  if (length(value) == 0L) {
    return(default)
  }
  value
}
