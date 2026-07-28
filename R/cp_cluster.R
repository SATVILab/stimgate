# Get cutpoints using joint density clustering
#
# Samples are clustered from their paired unstimulated and stimulated densities
# on one common absolute-expression grid over the left-hand region. Only
# thresholds generated directly by the local-FDR procedure are donors. Within
# each cluster, direct thresholds are winsorised to the 15th and 85th
# percentiles when at least three direct thresholds are available. Every
# non-direct threshold is replaced by the 60th percentile whenever its cluster
# has at least one direct threshold. Clusters without a direct threshold retain
# their original high thresholds.
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

  direct <- gateTblStim$locGeneratedDirect %in% TRUE &
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
    gateTblStim = gateTblStim,
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
  leftUpperX <- .getCpClusterLocLeftUpperX(
    gateTblStim = gateTblStim,
    control = control
  )
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

  clusterInput <- featureTbl |>
    dplyr::left_join(
      gateTblStim |>
        dplyr::select(ind, locGeneratedDirect, gate),
      by = "ind"
    ) |>
    dplyr::mutate(
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE &
        is.finite(suppressWarnings(as.numeric(.data$gate)))
    )

  clusterObj <- .getCpClusterLocClusters(
    featureTbl = clusterInput,
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
      .data$locGeneratedDirect %in% TRUE,
      is.finite(.data$gateNumeric)
    ) |>
    dplyr::group_by(.data$grp) |>
    dplyr::summarise(
      locClusterNDirect = dplyr::n(),
      locClusterQ15 = dplyr::if_else(
        dplyr::n() >= control$minDirectForWinsor,
        .getCpClusterLocRqQuantile(
          .data$gateNumeric,
          tau = control$winsorLower
        ),
        NA_real_
      ),
      locClusterQ60 = .getCpClusterLocRqQuantile(
        .data$gateNumeric,
        tau = control$imputeQuantile
      ),
      locClusterQ85 = dplyr::if_else(
        dplyr::n() >= control$minDirectForWinsor,
        .getCpClusterLocRqQuantile(
          .data$gateNumeric,
          tau = control$winsorUpper
        ),
        NA_real_
      ),
      .groups = "drop"
    )

  locTbl |>
    dplyr::left_join(clusterSummary, by = "grp") |>
    dplyr::mutate(
      cpOrig = suppressWarnings(as.numeric(.data$gate)),
      isDirectDonor = .data$locGeneratedDirect %in% TRUE &
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
      locGeneratedDirect = .data$isDirectDonor,
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
  purrr::map_df(seq_len(nrow(gateTblStim)), function(i) {
    .getCpClusterLocRowUnchanged(
      row = gateTblStim[i, , drop = FALSE],
      reason = reason
    )
  })
}
