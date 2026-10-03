# Shared scalar local-FDR bandwidths
#
# `bwScope = "sample"` estimates the scalar local-FDR bandwidth separately for
# every stimulated sample. `"cytokine"` estimates the bandwidths of about
# `.bwSharedNTarget` tubes spread across the batches and uses their 10% trimmed
# mean for the whole channel. `"cluster"` clusters all tubes on their densities
# up to the right shoulder of the left modal complex, estimates bandwidths for
# about `.bwSharedNTarget` tubes spread across the clusters (at least one per
# cluster) and gives each tube its cluster's median bandwidth. Tubes with fewer
# than `minCell` cells are excluded throughout. A sample then uses the smaller of its stimulated and
# unstimulated tube bandwidths, as in the per-sample calculation. Fixed (`bw`)
# and adaptive bandwidths are unaffected.

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

  readBatch <- function(i) {
    .bwSharedReadBatch(i, indBatchList, .data, chnlSettings, pathProject)
  }
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
        chnlSettings = chnlSettings
      )
    }) |>
      unname() |>
      unlist()
    bwEst[!duplicated(names(bwEst))]
  }

  # Thinned expression of every tube with at least `minCell` cells. Tubes with
  # fewer cells are not gated by local-FDR and would add noisy bandwidths.
  xList <- purrr::map(seq_along(indBatchList), function(i) {
    readBatch(i) |>
      purrr::keep(function(x) length(x) >= chnlSettings$minCell) |>
      purrr::map(function(x) {
        x <- sort(x)
        x[.spreadInd(length(x), .bwSharedNCellFeature)]
      })
  }) |>
    unname() |>
    purrr::flatten()
  xList <- xList[!duplicated(names(xList))]

  if (identical(chnlSettings$bwScope, "cytokine")) {
    bwVec <- estBw(names(xList)[.spreadInd(length(xList), .bwSharedNTarget)])
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

  ord <- order(clusterTbl$grp)
  sel <- ord[.spreadInd(length(ord), .bwSharedNTarget)]
  sel <- union(sel, match(unique(clusterTbl$grp), clusterTbl$grp))
  bwEst <- estBw(clusterTbl$ind[sel])

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
