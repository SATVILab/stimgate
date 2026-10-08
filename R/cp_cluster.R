# Get cutpoints using joint density clustering
#
# Samples are clustered from their paired unstimulated and stimulated densities
# on one common absolute-expression grid over the left-hand region. Only
# responders (.getCpShareResponder(), fixed at the per-sample step) are donors,
# using their gates after the batch step. Within each cluster, donor
# thresholds are winsorised to the 15th and 85th percentiles when at least
# three donors are available. Every other threshold is replaced by the 60th
# percentile whenever its cluster has at least one donor. Clusters without a
# donor retain their original thresholds. Lowered thresholds are then limited
# as in the batch step (.getCpClusterLocLimit()).
#' @keywords internal
.getCpCluster <- function(
  .data,
  gateTbl,
  chnlSettings,
  stage,
  pathProject,
  control = list(),
  filterOtherCytPos,
  calcCytPosGates,
  indBatchList,
  exLookup = NULL
) {
  stageChnl <- file.path(stage, chnlSettings$chnlCut)
  control <- .getCpClusterControlUpdate(control)
  gateTbl <- .getCpClusterLocGateTblPrepare(gateTbl)

  if (is.null(exLookup)) {
    exLookup <- .getCpClusterLocExLookup(
      .data = .data,
      indBatchList = indBatchList,
      chnlSettings = chnlSettings,
      filterOtherCytPos = filterOtherCytPos,
      calcCytPosGates = calcCytPosGates,
      gateTbl = gateTbl,
      pathProject = pathProject
    )
  }
  gateTblStim <- gateTbl |>
    dplyr::filter(
      .data$ind %in% names(.env$exLookup),
      !(.data$locSource %in% "unstim_summary")
    )

  if (nrow(gateTblStim) == 0L) {
    cpTbl <- .getCpClusterLocSkipOut(gateTblStim, "no_stimulated_samples")
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  direct <- gateTblStim$locResponder %in% TRUE &
    is.finite(suppressWarnings(as.numeric(gateTblStim$gate)))
  if (!any(direct)) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "no_direct_threshold_donors"
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  commonBw <- .getCpClusterLocCommonBw(
    indDirect = as.character(gateTblStim$ind[direct]),
    exLookup = exLookup,
    chnlSettings = chnlSettings
  )
  if (!is.finite(commonBw) || commonBw <= 0) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "common_density_bandwidth_unavailable"
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }
  .intSaveNm("locClusterCommonBw", commonBw, "all", stageChnl, pathProject)

  exprRange <- .getCpClusterLocExprRange(exLookup)
  leftUpperX <- control$leftThresholdFrac * stats::quantile(
    suppressWarnings(as.numeric(gateTblStim$gate[direct])),
    probs = control$leftThresholdQuantile,
    na.rm = TRUE
  )[[1]]
  densityGrid <- .getCpClusterLocDensityGrid(
    exprMin = exprRange[["min"]],
    leftUpperX = leftUpperX,
    nGrid = control$nGrid
  )

  featureTbl <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    densityGrid = densityGrid,
    bw = commonBw
  )
  .intSaveNm(
    "locClusterJointDensityFeatures",
    featureTbl,
    "all",
    stageChnl,
    pathProject
  )

  featureCols <- .getCpClusterLocFeatureCols(featureTbl)
  clusterable <- gateTblStim$ind %in% featureTbl$ind
  nDirectClusterable <- sum(direct & clusterable)
  if (length(featureCols) == 0L || nDirectClusterable < 1L) {
    cpTbl <- .getCpClusterLocSkipOut(
      gateTblStim = gateTblStim,
      reason = "no_clusterable_direct_threshold_donors"
    )
    .intSave("all", stageChnl, pathProject, cpTbl)
    return(cpTbl)
  }

  clusterObj <- .getCpClusterLocClusters(
    featureTbl = featureTbl,
    control = control
  )
  clusterTbl <- clusterObj$clusterTbl
  .intSaveNm(
    "locClusterAssignments",
    clusterTbl,
    "all",
    stageChnl,
    pathProject
  )

  locTbl <- gateTblStim |>
    dplyr::left_join(clusterTbl, by = "ind")

  cpTbl <- .getCpClusterLocApplyQuantiles(
    locTbl = locTbl,
    commonBw = commonBw,
    control = control,
    nInitialClusters = clusterObj$nInitialClusters
  ) |>
    .getCpClusterLocLimit(
      exLookup = exLookup,
      shareCap = chnlSettings$locShareCap %||% 1.5,
      cellCap = chnlSettings$locShareCellCap %||% 0.5
    ) |>
    dplyr::arrange(.data$ind)

  .intSaveNm("locClusterQuantileTbl", cpTbl, "all", stageChnl, pathProject)
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
  if (!"locResponder" %in% names(gateTbl)) {
    gateTbl$locResponder <- gateTbl$locGeneratedDirect
  }
  if (!"propBsEst" %in% names(gateTbl)) {
    gateTbl$propBsEst <- NA_real_
  }
  gateTbl |>
    dplyr::mutate(
      ind = as.character(.data$ind),
      locGenerated = .data$locGenerated %in% TRUE,
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE,
      locResponder = .data$locResponder %in% TRUE,
      propBsEst = suppressWarnings(as.numeric(.data$propBsEst))
    )
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
.getCpClusterLocApplyQuantiles <- function(
  locTbl,
  commonBw,
  control,
  nInitialClusters
) {
  clusterSummary <- locTbl |>
    dplyr::mutate(
      gateNumeric = suppressWarnings(as.numeric(.data$gate))
    ) |>
    dplyr::filter(
      !is.na(.data$grp),
      .data$locResponder %in% TRUE,
      is.finite(.data$gateNumeric)
    ) |>
    dplyr::group_by(.data$grp) |>
    dplyr::summarise(
      locClusterNDirect = dplyr::n(),
      {
        q <- .getCpClusterLocRqQuantile(
          .data$gateNumeric,
          tau = c(
            control$winsorLower[1], control$imputeQuantile[1],
            control$winsorUpper[1]
          )
        )
        if (dplyr::n() < control$minDirectForWinsor) {
          q[c(1L, 3L)] <- NA_real_
        }
        tibble::tibble(
          locClusterQ15 = q[[1]],
          locClusterQ60 = q[[2]],
          locClusterQ85 = q[[3]]
        )
      },
      .groups = "drop"
    )

  locTbl |>
    dplyr::left_join(clusterSummary, by = "grp") |>
    dplyr::mutate(
      cpOrig = suppressWarnings(as.numeric(.data$gate)),
      isDirectDonor = .data$locResponder %in% TRUE &
        is.finite(.data$cpOrig),
      clusterHasDirect = !is.na(.data$grp) &
        .data$locClusterNDirect >= 1L &
        is.finite(.data$locClusterQ60),
      clusterWinsorises = .data$clusterHasDirect &
        .data$locClusterNDirect >= control$minDirectForWinsor &
        is.finite(.data$locClusterQ15) &
        is.finite(.data$locClusterQ85),
      cpFinal = dplyr::case_when(
        .data$clusterWinsorises & .data$isDirectDonor ~
          pmax(
            .data$locClusterQ15,
            pmin(.data$cpOrig, .data$locClusterQ85)
          ),
        .data$isDirectDonor ~ .data$cpOrig,
        .data$clusterHasDirect ~ .data$locClusterQ60,
        TRUE ~ .data$cpOrig
      ),
      locClusterAction = dplyr::case_when(
        is.na(.data$grp) ~ "unchanged_no_cluster",
        !.data$clusterHasDirect ~
          "unchanged_no_direct_threshold_in_cluster",
        .data$clusterWinsorises & .data$isDirectDonor &
          .data$cpOrig < .data$locClusterQ15 ~
          "direct_winsorised_to_q15",
        .data$clusterWinsorises & .data$isDirectDonor &
          .data$cpOrig > .data$locClusterQ85 ~
          "direct_winsorised_to_q85",
        .data$clusterWinsorises & .data$isDirectDonor ~
          "direct_retained_within_winsor_limits",
        .data$isDirectDonor ~
          "direct_retained_fewer_than_three_direct_thresholds",
        TRUE ~ "non_direct_replaced_by_q60"
      ),
      locClusterAdjusted = .data$clusterHasDirect &
        (
          !.data$isDirectDonor |
            !is.finite(.data$cpOrig) |
            .data$cpFinal != .data$cpOrig
        ),
      locGenerated = dplyr::if_else(
        .data$clusterHasDirect,
        TRUE,
        .data$locGenerated %in% TRUE
      ),
      locSource = dplyr::if_else(
        .data$clusterHasDirect & !.data$isDirectDonor,
        "cluster_q60",
        as.character(.data$locSource)
      ),
      locReason = dplyr::if_else(
        .data$clusterHasDirect & !.data$isDirectDonor,
        "replaced_by_cluster_direct_threshold_q60",
        as.character(.data$locReason)
      )
    ) |>
    dplyr::transmute(
      grp = as.character(.data$grp),
      grpUns = as.character(.data$grp),
      grpStim = as.character(.data$grp),
      ind = as.character(.data$ind),
      cpOrigQuantMin = .data$cpOrig,
      cpJoin = .data$locClusterQ60,
      cpJoinLse = .data$cpFinal,
      cpJoinLseOrig = .data$cpFinal,
      cpJoinLseOrigMean = .data$cpFinal,
      cpJoinTgOrig = .data$cpFinal,
      cpJoinTgOrigMean = .data$cpFinal,
      cpJoinLseOrigMeanTg = .data$cpFinal,
      cpTolUns = NA_real_,
      cpTolStim = NA_real_,
      cpMedianUns = .data$locClusterQ60,
      cpMedianStim = .data$locClusterQ60,
      locGenerated = .data$locGenerated %in% TRUE,
      locGeneratedDirect = .data$isDirectDonor &
        .data$locGeneratedDirect %in% TRUE,
      locResponder = .data$locResponder %in% TRUE,
      propBsEst = .data$propBsEst,
      locShareLimit = "none",
      locShareProposed = .data$cpFinal,
      locSource = as.character(.data$locSource),
      locReason = as.character(.data$locReason),
      locClusterReason = dplyr::if_else(
        .data$clusterHasDirect,
        "cluster_direct_threshold_quantile_transfer",
        dplyr::if_else(
          is.na(.data$grp),
          "cluster_unavailable",
          "cluster_has_no_direct_threshold_original_retained"
        )
      ),
      locClusterAction = .data$locClusterAction,
      locClusterAdjusted = .data$locClusterAdjusted,
      locClusterBw = commonBw,
      locClusterNDirect = .data$locClusterNDirect,
      locClusterQ15 = .data$locClusterQ15,
      locClusterQ60 = .data$locClusterQ60,
      locClusterQ85 = .data$locClusterQ85,
      locClusterNInitial = as.integer(nInitialClusters),
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
.getCpClusterLocSkipOut <- function(gateTblStim, reason) {
  cp <- suppressWarnings(as.numeric(gateTblStim$gate))
  tibble::tibble(
    .rows = nrow(gateTblStim),
    grp = NA_character_,
    grpUns = NA_character_,
    grpStim = NA_character_,
    ind = as.character(gateTblStim$ind),
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
    locGenerated = suppressWarnings(gateTblStim$locGenerated %in% TRUE),
    locGeneratedDirect = suppressWarnings(
      gateTblStim$locGeneratedDirect %in% TRUE
    ),
    locResponder = gateTblStim$locResponder %in% TRUE,
    propBsEst = suppressWarnings(as.numeric(
      gateTblStim$propBsEst %||% rep(NA_real_, nrow(gateTblStim))
    )),
    locShareLimit = "none",
    locShareProposed = cp,
    locSource = as.character(
      gateTblStim$locSource %||% rep(NA_character_, nrow(gateTblStim))
    ),
    locReason = as.character(
      gateTblStim$locReason %||% rep(NA_character_, nrow(gateTblStim))
    ),
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

# Limit lowered cluster gates as in the batch step (.getCpShareApply()):
# donors use the rule for responders, other tubes in clusters with a donor the
# rule for non-responders, with the median frequency of the cluster's donors
# at their current gates. Frequencies use the stimulated and unstimulated
# expression in `exLookup`.
#' @keywords internal
.getCpClusterLocLimit <- function(cpTbl, exLookup, shareCap, cellCap) {
  if (nrow(cpTbl) == 0L) {
    return(cpTbl)
  }
  getX <- function(ind, type) .getCut(exLookup[[ind]][[type]])
  donor <- cpTbl$locResponder %in% TRUE & is.finite(cpTbl$cpOrigQuantMin)
  hasDonor <- cpTbl$locClusterReason %in%
    "cluster_direct_threshold_quantile_transfer"
  freqCurr <- vapply(seq_len(nrow(cpTbl)), function(i) {
    if (!donor[[i]]) {
      return(NA_real_)
    }
    .getCpShareFreq(
      cpTbl$cpOrigQuantMin[[i]],
      getX(cpTbl$ind[[i]], "stim"),
      getX(cpTbl$ind[[i]], "uns")
    )
  }, numeric(1))
  freqDonor <- stats::ave(
    freqCurr,
    dplyr::coalesce(cpTbl$grp, ""),
    FUN = function(x) stats::median(x, na.rm = TRUE)
  )

  gate <- cpTbl$cpJoinTgOrig
  limit <- cpTbl$locShareLimit
  for (i in which(hasDonor)) {
    res <- .getCpShareApply(
      gs = cpTbl$cpJoinTgOrig[[i]],
      gc = cpTbl$cpOrigQuantMin[[i]],
      responder = donor[[i]],
      propBsEst = cpTbl$propBsEst[[i]],
      freqDonor = freqDonor[[i]],
      xStim = getX(cpTbl$ind[[i]], "stim"),
      xUns = getX(cpTbl$ind[[i]], "uns"),
      shareCap = shareCap,
      cellCap = cellCap
    )
    gate[[i]] <- res$gate
    limit[[i]] <- res$limit
  }

  limited <- !(limit %in% "none")
  cpTbl$locShareLimit <- limit
  cpTbl$locReason[limited] <- paste0(
    cpTbl$locReason[limited], "_limited_by_", limit[limited]
  )
  cpTbl$locClusterAction[limited] <- paste0(
    cpTbl$locClusterAction[limited], "_limited_by_", limit[limited]
  )
  cpTbl$locClusterAdjusted <- hasDonor &
    (!donor | !is.finite(cpTbl$cpOrigQuantMin) | gate != cpTbl$cpOrigQuantMin)
  for (col in c(
    "cpJoinLse", "cpJoinLseOrig", "cpJoinLseOrigMean", "cpJoinTgOrig",
    "cpJoinTgOrigMean", "cpJoinLseOrigMeanTg"
  )) {
    cpTbl[[col]] <- gate
  }
  cpTbl
}
