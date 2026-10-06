# Shared scalar local-FDR bandwidths
#
# `bwScope = "sample"` estimates the scalar local-FDR bandwidth separately for
# every stimulated sample. `"cytokine"` estimates the bandwidths of about
# `.bwSharedNTarget` tubes spread across the batches and uses their 10% trimmed
# mean for the whole channel. `"cluster"` clusters all tubes on their densities
# up to the right shoulder of the left modal complex, estimates bandwidths for
# about `.bwSharedNTarget` tubes spread across the clusters (at least one per
# cluster) and gives each tube its cluster's median bandwidth. Tubes with fewer
# than `minCell` cells are excluded throughout. Within the channel or cluster,
# tubes with at least `bwNcellMax` cells are preferred (see `.bwSharedSelect()`)
# and every selected tube's bandwidth is estimated on `bwNcellMax` cells, so the
# shared bandwidth matches one cell count. A sample then uses the smaller of its
# stimulated and unstimulated tube bandwidths, as in the per-sample calculation.
# Fixed (`bw`) and adaptive bandwidths are unaffected.

.bwSharedNTarget <- 100L
# Fraction of the left-complex peak height at which the clustering range ends
.bwSharedShoulderFrac <- 0.1
# Cells kept per tube for the clustering densities
.bwSharedNCellFeature <- 1e4L

#' @keywords internal
.completeChnlSettingsBwShared <- function(
  chnlSettings,
  indBatchList,
  .data,
  pathProject
) {
  chnlSettings$bwShared <- NULL
  chnlSettings$bwSharedTbl <- NULL
  if (
    identical(chnlSettings$bwScope, "sample") ||
      !is.null(chnlSettings$bw) ||
      .getCpUnsLocUseAdaptiveBw(chnlSettings)
  ) {
    return(chnlSettings)
  }

  batchCache <- new.env(parent = emptyenv())
  readBatch <- function(i) {
    key <- as.character(i)
    if (!exists(key, envir = batchCache, inherits = FALSE)) {
      batchCache[[key]] <- .bwSharedReadBatch(
        i = i,
        indBatchList = indBatchList,
        .data = .data,
        chnlSettings = chnlSettings,
        pathProject = pathProject
      )
    }
    batchCache[[key]]
  }
  # Estimate every selected tube's bandwidth on `bwNcellMax` cells: larger
  # tubes are downsampled and smaller ones upsampled (smoothed bootstrap).
  nCellPref <- suppressWarnings(as.numeric(chnlSettings$bwNcellMax))[1]
  if (length(nCellPref) == 0L || !is.finite(nCellPref)) nCellPref <- NULL
  chnlSettingsBw <- chnlSettings
  if (!is.null(nCellPref)) chnlSettingsBw$bwNcellMin <- nCellPref

  # Bandwidths of tubes `indSel`, reading only the batches that contain them
  estBw <- function(indSel) {
    batchSel <- which(vapply(
      indBatchList,
      function(indBatch) any(as.character(indBatch) %in% indSel),
      logical(1)
    ))
    bwEst <- purrr::map(batchSel, function(i) {
      xBatch <- readBatch(i)
      vapply(
        xBatch[intersect(names(xBatch), indSel)],
        .getCpUnsLocGetDensRawDensitiesBwInit,
        numeric(1),
        chnlSettings = chnlSettingsBw
      )
    }) |>
      unname() |>
      unlist()
    bwEst[!duplicated(names(bwEst))]
  }

  # Every tube with at least `minCell` cells. Tubes with fewer cells are not
  # gated by local-FDR and would add noisy bandwidths.
  tubeList <- purrr::map(seq_along(indBatchList), function(i) {
    readBatch(i) |>
      purrr::keep(function(x) length(x) >= chnlSettings$minCell)
  }) |>
    unname() |>
    purrr::flatten()
  tubeList <- tubeList[!duplicated(names(tubeList))]
  nCell <- lengths(tubeList)
  # Thinned expression for the clustering densities
  xList <- purrr::map(tubeList, function(x) {
    x <- sort(x)
    x[.spreadInd(length(x), .bwSharedNCellFeature)]
  })

  if (identical(chnlSettings$bwScope, "cytokine")) {
    bwVec <- estBw(.bwSharedSelect(nCell, .bwSharedNTarget, nCellPref))
    bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
    chnlSettings$bwShared <- if (length(bwVec) > 0L) {
      mean(bwVec, trim = 0.1)
    } else {
      chnlSettings$bwFallback
    }
    message(
      "shared bandwidth for ", chnlSettings$marker, ": ",
      signif(chnlSettings$bwShared, 3)
    )
    return(chnlSettings)
  }

  clusterTbl <- .bwSharedCluster(xList, bw = chnlSettings$bwFallback)
  if (nrow(clusterTbl) == 0L) {
    chnlSettings$bwShared <- chnlSettings$bwFallback
    return(chnlSettings)
  }

  # Tubes per cluster in proportion to cluster size, at least one each; the
  # tubes themselves are then chosen by cell count within each cluster.
  ord <- order(clusterTbl$grp)
  sel <- ord[.spreadInd(length(ord), .bwSharedNTarget)]
  sel <- union(sel, match(unique(clusterTbl$grp), clusterTbl$grp))
  nTargetGrp <- table(clusterTbl$grp[sel])
  indSel <- unlist(lapply(names(nTargetGrp), function(g) {
    indGrp <- clusterTbl$ind[clusterTbl$grp == g]
    .bwSharedSelect(nCell[indGrp], nTargetGrp[[g]], nCellPref)
  }), use.names = FALSE)
  bwEst <- estBw(indSel)

  clusterTbl$bwEst <- unname(bwEst[clusterTbl$ind])
  clusterTbl <- clusterTbl |>
    dplyr::group_by(.data$grp) |>
    dplyr::mutate(bw = stats::median(.data$bwEst, na.rm = TRUE)) |>
    dplyr::ungroup()

  # Fallback for tubes that could not be clustered, and for prejoined samples
  bwShared <- stats::median(bwEst, na.rm = TRUE)
  chnlSettings$bwShared <- if (is.finite(bwShared)) {
    bwShared
  } else {
    chnlSettings$bwFallback
  }
  chnlSettings$bwSharedTbl <- clusterTbl
  bwGrp <- unique(clusterTbl[c("grp", "bw")])
  message(
    "shared bandwidths for ", chnlSettings$marker, " by tube cluster: ",
    paste0(bwGrp$grp, " = ", signif(bwGrp$bw, 3), collapse = ", ")
  )
  chnlSettings
}

# Names of up to `nTarget` tubes for shared bandwidth estimation, from the
# named cell counts `nCell` (in batch order). Tubes with at least `nCellPref`
# cells are preferred, spread across the batches. If there are too few, the
# rest are drawn at random from bands one tenth of `nCellPref` wide below it,
# the highest band first, down to half of `nCellPref` (9-10k, ..., 5-6k for
# 10,000). If even that finds no tube, the bands continue down to the smallest
# tube. Without a finite `nCellPref`, tubes are spread across the batches.
#' @keywords internal
.bwSharedSelect <- function(nCell, nTarget, nCellPref) {
  ind <- names(nCell)
  if (length(ind) == 0L || nTarget < 1L) {
    return(character())
  }
  if (is.null(nCellPref)) {
    return(ind[.spreadInd(length(ind), nTarget)])
  }
  sel <- ind[nCell >= nCellPref]
  if (length(sel) >= nTarget) {
    return(sel[.spreadInd(length(sel), nTarget)])
  }
  width <- nCellPref / 10
  # Add random tubes from each band, highest first, until `nTarget` or `floor`
  fill <- function(sel, upper, floor) {
    while (length(sel) < nTarget && upper > floor) {
      band <- ind[nCell >= upper - width & nCell < upper]
      nTake <- min(length(band), nTarget - length(sel))
      if (nTake > 0L) {
        sel <- c(sel, band[sample.int(length(band), nTake)])
      }
      upper <- upper - width
    }
    sel
  }
  sel <- fill(sel, nCellPref, nCellPref / 2)
  if (length(sel) > 0L) {
    return(sel)
  }
  fill(sel, nCellPref / 2, min(nCell) - width)
}

# Named list of each tube's cut-channel expression in batch `i`
#' @keywords internal
.bwSharedReadBatch <- function(
  i,
  indBatchList,
  .data,
  chnlSettings,
  pathProject
) {
  exList <- .getExList(
    .data = .data,
    indBatch = indBatchList[[i]],
    batch = names(indBatchList)[i],
    pop = chnlSettings$popGate,
    chnlCut = chnlSettings$chnlCut,
    pathProject = pathProject
  )
  purrr::map(exList, function(ex) {
    x <- .getCut(ex)
    x <- x[is.finite(x)]
    if (isTRUE(chnlSettings$excMin) && length(x) > 0L) {
      x <- x[x > min(x)]
    }
    x
  })
}

# Cluster tubes on their densities from the low end of the pooled expression up
# to where the pooled density first falls to `.bwSharedShoulderFrac` of the
# right-most peak of its left modal complex. Returns a tibble of `ind` and `grp`
# for every tube with a usable density.
#' @keywords internal
.bwSharedCluster <- function(xList, bw) {
  empty <- tibble::tibble(ind = character(), grp = character())
  xList <- purrr::keep(
    xList,
    function(x) length(x) >= 3L && length(unique(x)) >= 3L
  )
  if (length(xList) == 0L) {
    return(empty)
  }

  pooled <- unlist(
    lapply(xList, function(x) x[.spreadInd(length(x), 1e3)]),
    use.names = FALSE
  )
  densPooled <- stats::density(pooled, bw = bw, n = 512L)
  iPeak <- .getPeakMainLeftIdx(densPooled$y)
  iEnd <- which(
    seq_along(densPooled$y) > iPeak &
      densPooled$y <= .bwSharedShoulderFrac * densPooled$y[iPeak]
  )[1]
  upper <- densPooled$x[if (is.na(iEnd)) length(densPooled$x) else iEnd]

  control <- .getCpClusterControlUpdate(list())
  densityGrid <- .getCpClusterLocDensityGrid(
    exprMin = stats::quantile(pooled, 0.0025, names = FALSE),
    leftUpperX = upper,
    nGrid = control$nGrid
  )
  featureMat <- t(vapply(
    xList,
    .getCpClusterLocDensityFeature,
    numeric(length(densityGrid)),
    densityGrid = densityGrid,
    bw = bw
  ))
  colnames(featureMat) <- sprintf("x%03d", seq_len(ncol(featureMat)))
  featureTbl <- tibble::as_tibble(featureMat) |>
    dplyr::mutate(ind = names(xList), .before = 1L)
  featureTbl <- featureTbl[stats::complete.cases(featureMat), , drop = FALSE]
  if (nrow(featureTbl) == 0L) {
    return(empty)
  }

  .getCpClusterLocClusters(featureTbl, control)$clusterTbl
}

# Shared bandwidth for one tube: its cluster's bandwidth, else the channel-level
# shared bandwidth (also used for prejoined samples, which have no single tube).
#' @keywords internal
.bwSharedGet <- function(chnlSettings, ind) {
  tbl <- chnlSettings$bwSharedTbl
  bw <- if (is.null(tbl)) {
    NA_real_
  } else {
    tbl$bw[match(as.character(ind), tbl$ind)]
  }
  if (length(bw) == 1L && is.finite(bw)) bw else chnlSettings$bwShared
}
