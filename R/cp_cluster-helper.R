# Helpers for joint density clustering and within-cluster threshold transfer

#' @keywords internal
.getCpClusterControlUpdate <- function(control) {
  defaults <- list(
    leftThresholdFrac = 0.8,
    leftThresholdQuantile = 0.1,
    nGrid = 128L,
    winsorLower = 0.15,
    imputeQuantile = 0.60,
    winsorUpper = 0.85,
    minDirectForWinsor = 3L,
    maxClusters = 6L,
    gapBootstraps = 50L,
    kmeansNstart = 25L,
    seed = 1L
  )
  for (name in names(defaults)) {
    if (is.null(control[[name]])) {
      control[[name]] <- defaults[[name]]
    }
  }

  stopifnot(
    is.numeric(control$leftThresholdFrac),
    length(control$leftThresholdFrac) == 1L,
    is.finite(control$leftThresholdFrac),
    control$leftThresholdFrac > 0,
    is.numeric(control$leftThresholdQuantile),
    length(control$leftThresholdQuantile) == 1L,
    control$leftThresholdQuantile >= 0,
    control$leftThresholdQuantile <= 1,
    is.numeric(control$nGrid),
    length(control$nGrid) == 1L,
    control$nGrid >= 16,
    is.numeric(control$winsorLower),
    is.numeric(control$imputeQuantile),
    is.numeric(control$winsorUpper),
    0 <= control$winsorLower,
    control$winsorLower < control$imputeQuantile,
    control$imputeQuantile < control$winsorUpper,
    control$winsorUpper <= 1,
    is.numeric(control$minDirectForWinsor),
    length(control$minDirectForWinsor) == 1L,
    control$minDirectForWinsor >= 3,
    is.numeric(control$maxClusters),
    control$maxClusters >= 1,
    is.numeric(control$gapBootstraps),
    control$gapBootstraps >= 1,
    is.numeric(control$kmeansNstart),
    control$kmeansNstart >= 1
  )

  control$nGrid <- as.integer(control$nGrid)
  control$minDirectForWinsor <- as.integer(control$minDirectForWinsor)
  control$maxClusters <- as.integer(control$maxClusters)
  control$gapBootstraps <- as.integer(control$gapBootstraps)
  control$kmeansNstart <- as.integer(control$kmeansNstart)
  control$seed <- as.integer(control$seed)
  control
}

#' @keywords internal
.getCpClusterLocCommonBw <- function(gateTblStim, exLookup, chnlSettings) {
  indDirect <- gateTblStim |>
    dplyr::filter(
      .data$locGeneratedDirect %in% TRUE,
      is.finite(suppressWarnings(as.numeric(.data$gate)))
    ) |>
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
.getCpClusterLocExprRange <- function(exLookup) {
  rangeTbl <- purrr::map_df(exLookup, function(exPair) {
    x <- c(.getCut(exPair$stim), .getCut(exPair$uns))
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) <= 5L) {
      return(NULL)
    }
    q <- stats::quantile(x, c(0.0025, 0.999), na.rm = TRUE)
    tibble::tibble(lb = q[[1]], ub = q[[2]])
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
.getCpClusterLocLeftUpperX <- function(gateTblStim, control) {
  directGates <- gateTblStim |>
    dplyr::filter(.data$locGeneratedDirect %in% TRUE) |>
    dplyr::pull("gate")
  directGates <- suppressWarnings(as.numeric(directGates))
  directGates <- directGates[is.finite(directGates)]
  if (length(directGates) == 0L) {
    return(NA_real_)
  }
  control$leftThresholdFrac * stats::quantile(
    directGates,
    probs = control$leftThresholdQuantile,
    na.rm = TRUE
  )[[1]]
}

#' @keywords internal
.getCpClusterLocDensityGrid <- function(exprMin, leftUpperX, nGrid) {
  exprMin <- suppressWarnings(as.numeric(exprMin))[1]
  leftUpperX <- suppressWarnings(as.numeric(leftUpperX))[1]
  if (!is.finite(exprMin) || !is.finite(leftUpperX)) {
    return(numeric())
  }
  if (leftUpperX <= exprMin) {
    pad <- max(abs(exprMin), 1) * sqrt(.Machine$double.eps)
    leftUpperX <- exprMin + pad
  }
  seq(exprMin, leftUpperX, length.out = as.integer(nGrid))
}

#' @keywords internal
.getCpClusterLocDensityFeature <- function(x, densityGrid, bw) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) >= 1L) {
    x <- x[x > min(x)]
  }
  if (
    length(x) < 3L ||
      length(unique(x)) < 3L ||
      length(densityGrid) < 2L ||
      !is.finite(bw) ||
      bw <= 0
  ) {
    return(rep(NA_real_, length(densityGrid)))
  }

  dens <- try(
    suppressWarnings(stats::density(
      x,
      bw = bw,
      from = densityGrid[1],
      to = densityGrid[length(densityGrid)],
      n = length(densityGrid)
    )),
    silent = TRUE
  )
  if (inherits(dens, "try-error")) {
    return(rep(NA_real_, length(densityGrid)))
  }

  y <- suppressWarnings(as.numeric(dens$y))
  total <- sum(y, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(densityGrid)))
  }
  y / total
}

#' @keywords internal
.getCpClusterLocJointFeatureTbl <- function(
  exLookup,
  densityGrid,
  bw
) {
  if (length(densityGrid) < 2L) {
    return(tibble::tibble(ind = character(), batch = character()))
  }

  featureTbl <- purrr::map_df(exLookup, function(exPair) {
    uns <- .getCpClusterLocDensityFeature(
      x = .getCut(exPair$uns),
      densityGrid = densityGrid,
      bw = bw
    )
    stim <- .getCpClusterLocDensityFeature(
      x = .getCut(exPair$stim),
      densityGrid = densityGrid,
      bw = bw
    )
    values <- c(uns, stim)
    names(values) <- c(
      sprintf("uns_x%03d", seq_along(uns)),
      sprintf("stim_x%03d", seq_along(stim))
    )
    tibble::as_tibble_row(c(
      list(
        ind = as.character(exPair$ind),
        batch = as.character(exPair$batch)
      ),
      as.list(values)
    ))
  })

  featureCols <- .getCpClusterLocFeatureCols(featureTbl)
  if (length(featureCols) == 0L) {
    return(featureTbl[0, , drop = FALSE])
  }
  complete <- stats::complete.cases(featureTbl[, featureCols, drop = FALSE])
  featureTbl[complete, , drop = FALSE]
}

#' @keywords internal
.getCpClusterLocFeatureCols <- function(featureTbl) {
  grep("^(uns|stim)_x\\d+$", names(featureTbl), value = TRUE)
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
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

#' @keywords internal
.getCpClusterLocInitialNClusters <- function(
  featureTbl,
  control
) {
  featureCols <- .getCpClusterLocFeatureCols(featureTbl)
  x <- as.matrix(featureTbl[, featureCols, drop = FALSE])
  nUnique <- nrow(unique(as.data.frame(x)))
  nForClustering <- nrow(x)
  kMaxByN <- max(1L, nForClustering - 1L)
  kMax <- min(
    control$maxClusters,
    kMaxByN,
    nUnique
  )
  kMax <- max(1L, as.integer(kMax))
  if (kMax == 1L) {
    return(1L)
  }

  gap <- .getCpClusterLocWithSeed(control$seed, {
    cluster::clusGap(
      x,
      FUNcluster = function(x, k) {
        stats::kmeans(
          x,
          centers = k,
          nstart = control$kmeansNstart
        )
      },
      K.max = kMax,
      B = control$gapBootstraps
    )
  })
  as.integer(cluster::maxSE(
    gap$Tab[, "gap"],
    gap$Tab[, "SE.sim"]
  ))
}

#' @keywords internal
.getCpClusterLocClusters <- function(featureTbl, control) {
  featureCols <- .getCpClusterLocFeatureCols(featureTbl)
  x <- as.matrix(featureTbl[, featureCols, drop = FALSE])
  nClusters <- .getCpClusterLocInitialNClusters(
    featureTbl = featureTbl,
    control = control
  )

  grp <- if (nClusters == 1L) {
    rep(1L, nrow(featureTbl))
  } else {
    .getCpClusterLocWithSeed(control$seed, {
      stats::kmeans(
        x,
        centers = nClusters,
        nstart = control$kmeansNstart
      )$cluster
    })
  }
  list(
    clusterTbl = tibble::tibble(
      ind = as.character(featureTbl$ind),
      grp = as.character(grp)
    ),
    nInitialClusters = nClusters,
    nFinalClusters = length(unique(grp))
  )
}

# Empirical quantile used by an intercept-only quantile regression. This is the
# inverse empirical distribution function, and therefore minimises the same
# check loss as `quantreg::rq(x ~ 1, tau = tau)` without requiring quantreg.
# When the minimiser is an interval, the lower endpoint is returned.
#' @keywords internal
.getCpClusterLocRqQuantile <- function(x, tau) {
  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x)])
  tau <- suppressWarnings(as.numeric(tau))[1]

  if (length(x) == 0L || !is.finite(tau) || tau < 0 || tau > 1) {
    return(NA_real_)
  }
  if (tau <= 0) {
    return(x[[1]])
  }

  index <- min(length(x), max(1L, ceiling(length(x) * tau)))
  x[[index]]
}

#' @keywords internal
.getCpClusterLocRowUnchanged <- function(row, reason) {
  cp <- suppressWarnings(as.numeric(row$gate[1]))
  tibble::tibble(
    grp = NA_character_,
    grpUns = NA_character_,
    grpStim = NA_character_,
    ind = as.character(row$ind[1]),
    cpOrigQuantMin = cp,
    cpJoin = NA_real_,
    cpJoinLse = cp,
    cpJoinLseOrig = cp,
    cpJoinLseOrigMean = cp,
    cpJoinTgOrig = cp,
    cpJoinTgOrigMean = cp,
    cpJoinLseOrigMeanTg = cp,
    cpTolUns = NA_real_,
    cpTolStim = NA_real_,
    cpMedianUns = NA_real_,
    cpMedianStim = NA_real_,
    locGenerated = suppressWarnings(row$locGenerated[1] %in% TRUE),
    locGeneratedDirect = suppressWarnings(
      row$locGeneratedDirect[1] %in% TRUE
    ),
    locSource = as.character(row$locSource[1] %||% NA_character_),
    locReason = as.character(row$locReason[1] %||% NA_character_),
    locClusterReason = reason,
    locClusterAction = "unchanged",
    locClusterAdjusted = FALSE,
    locClusterBw = NA_real_,
    locClusterNDirect = NA_integer_,
    locClusterQ15 = NA_real_,
    locClusterQ60 = NA_real_,
    locClusterQ85 = NA_real_,
    locClusterNInitial = NA_integer_,
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
.getCpClusterDensTblGetBatchPrepExListFilter <- function(
  exList,
  chnlCut,
  gateTbl,
  calcCytPosGates
) {
  .debug("Filtering other cytokine positive cells")
  exListFilter <- purrr::map(seq_along(exList), function(i) {
    if (i == 1L) {
      return(exList[[i]])
    }
    .getCpClusterDensTblGetBatchPrepExListFilterInd(
      exTbl = exList[[i]],
      gateTbl = gateTbl,
      chnlCut = chnlCut,
      calcCytPosGates = calcCytPosGates
    )
  }) |>
    stats::setNames(names(exList))
  exListFilter[-1]
}

#' @keywords internal
.getCpClusterDensTblGetBatchPrepExListFilterInd <- function(
  exTbl,
  gateTbl,
  chnlCut,
  calcCytPosGates
) {
  posInd <- .get_pos_ind_but_single_pos_for_one_cyt(
    ex = exTbl,
    gateTbl = gateTbl[gateTbl[["ind"]] == attr(exTbl, "ind"), ],
    chnlSingleExc = chnlCut,
    chnl = NULL,
    gateTypeCytPos = if (calcCytPosGates) "cyt" else "base",
    gateTypeSinglePos = "base"
  )
  exTbl[!posInd, , drop = FALSE]
}
