mbMod <- function(
  dataMod,
  chnlSettings,
  probCol,
  startX,
  lowerBoundX = NA_real_,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL,
  floorX = NA_real_
) {
  info <- list(
    applied = FALSE,
    reason = "marginal_filter_not_run",
    startX = startX,
    lowerBoundX = lowerBoundX,
    trimTbl = tibble::tibble(
      step = integer(), cutFrom = numeric(), cutTo = numeric(),
      nStim = integer(), nUns = integer(), contribution = numeric(),
      kept = logical()
    ),
    nStepsTrimmed = 0L,
    trimReason = "marginal_scan_not_run"
  )
  if (!is.data.frame(dataMod) || nrow(dataMod) < 4L || !is.finite(startX)) {
    info$reason <- "insufficient_data_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  x <- suppressWarnings(as.numeric(.getCut(dataMod)))

  # The main filtering routes already supply dataMod in expression order.
  # Only sort here when required, preserving this helper's existing behaviour
  # for direct or other callers.
  if (is.unsorted(x)) {
    ord <- order(x)
    dataMod <- .getCpUnsLocSubsetRows(dataMod, ord)
    x <- x[ord]
  }

  prob <- .getCpUnsLocProbability(dataMod, probCol)

  finite <- is.finite(x) & is.finite(prob)

  dm <- tibble::tibble(
    x = x[finite],
    prob = prob[finite]
  )

  if (nrow(dm) < 4L || diff(range(dm$x)) <= 0) {
    info$reason <- "insufficient_finite_data_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  grid <- .getCpUnsLocMarginalBreaks(dataMod, startX, nBin = 50L)
  if (is.null(grid)) {
    info$reason <- "insufficient_bins_for_marginal_filter"
    return(list(dataMod = dataMod, info = info))
  }

  breaks <- grid$breaks
  refIndex <- grid$refIndex
  leftBins <- rev(seq_len(refIndex - 1L))
  info$nLeftBinsBeforeDensityBound <- length(leftBins)
  info$nLeftBinsExcludedByDensityBound <- 0L

  if (is.finite(lowerBoundX) && length(leftBins) > 0L) {
    if (.partial) {
      crossing <- leftBins[breaks[leftBins] < lowerBoundX & breaks[leftBins + 1L] > lowerBoundX]
      if (length(crossing) == 1L) breaks[[crossing]] <- lowerBoundX
    }
    eligibleBin <- breaks[leftBins] >= lowerBoundX
    info$nLeftBinsExcludedByDensityBound <- sum(!eligibleBin)
    leftBins <- leftBins[eligibleBin]

    if (length(leftBins) == 0L) {
      keep <- is.finite(x) & x >= startX
      info$nDropped <- sum(!keep)
      info$applied <- info$nDropped > 0L
      info$reason <- if (info$applied) {
        "used_marginal_reference_because_density_bound_excluded_all_bins"
      } else {
        "marginal_reference_kept_all_values_after_density_bound"
      }
      info$stopReason <- "density_lower_bound_excluded_all_left_bins"
      info$finalStartX <- startX
      info$gridShift <- grid$shift
      info$gridSpacing <- grid$spacing
      info$refIndex <- refIndex
      info$scanTbl <- tibble::tibble()
      info$trimReason <- "no_acceptance_steps"
      return(list(
        dataMod = .getCpUnsLocSubsetRows(dataMod, keep),
        info = info
      ))
    }
  }

  refNBin <- length(breaks) - refIndex
  right <- dm$x >= startX
  refNCell <- sum(right)
  refExpectedResp <- sum(dm$prob[right])

  if (refNBin < 1L || refNCell < 1L || !is.finite(refExpectedResp)) {
    info$reason <- "invalid_right_reference_region"
    return(list(dataMod = dataMod, info = info))
  }

  refCellsPerBin <- refNCell / refNBin
  refPurity <- refExpectedResp / refNCell
  if (!is.finite(refPurity) || refPurity <= 0) {
    info$reason <- "invalid_right_reference_purity"
    return(list(dataMod = dataMod, info = info))
  }

  cellBinRatio <- suppressWarnings(as.numeric(
    .getCpUnsLocSetting(chnlSettings, "locMarginalCellBinRatio", 1.25)
  )[1])
  if (!is.finite(cellBinRatio) || cellBinRatio <= 0) {
    cellBinRatio <- 1.25
  }
  purityRel <- .getCpUnsLocUnitValue(
    .getCpUnsLocSetting(chnlSettings, "locMarginalPurityRel", 0.85),
    0.85,
    allowZero = TRUE
  )

  maxCells <- cellBinRatio * refCellsPerBin
  minPurity <- purityRel * refPurity
  if (length(leftBins) == 0L) {
    info$reason <- "no_bins_left_of_marginal_reference"
    return(list(dataMod = dataMod, info = info))
  }
  # Assign each cell to its marginal bin once, rather than rescanning all cells
  # separately for every bin.
  nBinTotal <- length(breaks) - 1L
  binIdx <- findInterval(dm$x, breaks)

  # Preserve the existing interval definition:
  # breaks[i] <= x < breaks[i + 1].
  validBin <- binIdx >= 1L & binIdx <= nBinTotal

  # Number of cells in every bin.
  binNCell <- tabulate(
    binIdx[validBin],
    nbins = nBinTotal
  )

  # Sum of fitted response probabilities in every bin.
  binExpectedResp <- numeric(nBinTotal)

  if (any(validBin)) {
    probSums <- rowsum(
      dm$prob[validBin],
      group = binIdx[validBin],
      reorder = FALSE
    )

    probSumIdx <- as.integer(rownames(probSums))
    binExpectedResp[probSumIdx] <- probSums[, 1L]
  }

  nzb <- if (.crowd == "ref") binNCell[seq.int(refIndex, nBinTotal)] else binNCell
  nzb <- sort(nzb[nzb > 0], decreasing = TRUE)
  if (.crowd != "none") maxCells <- .crowdFrac * mean(utils::head(nzb, 3L)) else maxCells <- Inf
  info$crowdTop3 <- mean(utils::head(nzb, 3L))

  currentCut <- startX
  scan <- vector("list", length(leftBins))
  steps <- list()
  stopReason <- "scanned_all_left_bins"
  consecutiveRejections <- 0L

  if (.multires) {
    fl <- if (is.finite(floorX)) floorX else min(dm$x)
    densLim <- if (.crowd == "none") Inf else .crowdFrac * mean(utils::head(nzb, 3L)) / grid$spacing
    info$floorX <- fl; info$densLim <- densLim
    xSr <- get(".getCut", asNamespace("stimgate"))(exTblStimOrig); xSr <- xSr[is.finite(xSr)]
    xUr <- get(".getCut", asNamespace("stimgate"))(exTblUnsOrig); xUr <- xUr[is.finite(xUr)]
    lefts <- if (startX > fl) c(startX - (startX - fl) / c(16, 8, 4, 2), fl) else numeric(0)
    stopReason <- "no_room_above_floor"
    repeat {
      cand <- sort(lefts[lefts < currentCut - 1e-12], decreasing = TRUE)
      if (length(cand) == 0L) { stopReason <- "reached_floor"; break }
      acc <- NA_real_
      for (l in cand) {
        w <- currentCut - l
        if (.purSrc == "count") {
          n <- sum(xSr >= l & xSr < currentCut); u <- sum(xUr >= l & xUr < currentCut)
          e <- n - u * length(xSr) / length(xUr); nRef <- sum(xSr >= startX)
        } else {
          inS <- dm$x >= l & dm$x < currentCut
          n <- sum(inS); e <- sum(dm$prob[inS]); nRef <- refNCell
        }
        barS <- minPurity * min(1, (n / nRef)^.purExp)
        okS <- n > 0L && e / n >= barS && n / w <= densLim
        scan[[length(scan) + 1L]] <- tibble::tibble(left = l, right = currentCut, nCell = n,
          expectedResp = e, purity = if (n > 0L) e / n else NA_real_, maxCells = densLim * w,
          minPurity = barS, accepted = okS, acceptedPair = NA, pending = n == 0L)
        if (okS) { acc <- l; break }
      }
      if (!is.finite(acc)) { stopReason <- "no_span_accepted"; break }
      steps[[length(steps) + 1L]] <- tibble::tibble(step = length(steps) + 1L, cutFrom = currentCut,
        cutTo = acc, nStim = NA_integer_, nUns = NA_integer_, contribution = NA_real_, kept = TRUE)
      currentCut <- acc
    }
  } else {
  ok <- function(n, e, w) n > 0L && n <= w * maxCells && e / n >= minPurity
  j <- 1L
  while (j <= length(leftBins)) {
    i <- leftBins[[j]]
    nCell <- binNCell[[i]]; expectedResp <- binExpectedResp[[i]]
    pending <- nCell == 0L
    acc1 <- !pending && ok(nCell, expectedResp, 1)
    acc2 <- FALSE
    if (!acc1 && .pair && j < length(leftBins)) {
      i2 <- leftBins[[j + 1L]]
      acc2 <- ok(nCell + binNCell[[i2]], expectedResp + binExpectedResp[[i2]], 2)
    }
    scan[[j]] <- tibble::tibble(left = breaks[[i]], right = breaks[[i + 1L]], nCell = nCell,
      expectedResp = expectedResp, purity = if (pending) NA_real_ else expectedResp / nCell,
      maxCells = maxCells, minPurity = minPurity, accepted = acc1, acceptedPair = acc2, pending = pending)
    if (acc1 || acc2) {
      newCut <- if (acc1) breaks[[i]] else breaks[[i2]]
      steps[[length(steps) + 1L]] <- tibble::tibble(step = length(steps) + 1L, cutFrom = currentCut,
        cutTo = newCut, nStim = NA_integer_, nUns = NA_integer_, contribution = NA_real_, kept = TRUE)
      consecutiveRejections <- 0L
      currentCut <- newCut
      j <- j + if (acc1) 1L else 2L
      next
    }
    if (!pending) {
      consecutiveRejections <- consecutiveRejections + 1L
      if (consecutiveRejections >= 3L) { stopReason <- "three_consecutive_rejections"; break }
    }
    j <- j + 1L
  }
  }

  if (length(steps) > 0L) {
    info$trimTbl <- dplyr::bind_rows(steps)
  }
  trimTbl <- info$trimTbl
  xStim <- if (!is.null(exTblStimOrig)) .getCut(exTblStimOrig) else NULL
  xUns <- if (!is.null(exTblUnsOrig)) .getCut(exTblUnsOrig) else NULL
  for (k in seq_len(nrow(trimTbl))) {
    if (!is.null(xStim)) {
      trimTbl$nStim[[k]] <- sum(
        xStim >= trimTbl$cutTo[[k]] & xStim < trimTbl$cutFrom[[k]],
        na.rm = TRUE
      )
    }
    if (!is.null(xUns)) {
      trimTbl$nUns[[k]] <- sum(
        xUns >= trimTbl$cutTo[[k]] & xUns < trimTbl$cutFrom[[k]],
        na.rm = TRUE
      )
    }
  }
  if (is.null(exTblUnsOrig)) {
    info$trimReason <- "unstimulated_expression_unavailable"
  } else if (is.null(exTblStimOrig)) {
    info$trimReason <- "stimulated_expression_unavailable"
  } else if (nrow(exTblStimOrig) == 0L || nrow(exTblUnsOrig) == 0L) {
    info$trimReason <- "empty_original_expression"
  } else if (nrow(trimTbl) == 0L) {
    info$trimReason <- "no_acceptance_steps"
  } else {
    trimTbl$contribution <- trimTbl$nStim / nrow(exTblStimOrig) -
      trimTbl$nUns / nrow(exTblUnsOrig)
    info$trimReason <- "trimmed_all_steps"
    for (k in rev(seq_len(nrow(trimTbl)))) {
      if (trimTbl$contribution[[k]] > 0) {
        info$trimReason <- "positive_span_contribution"
        break
      }
      trimTbl$kept[[k]] <- FALSE
      currentCut <- trimTbl$cutFrom[[k]]
      info$nStepsTrimmed <- info$nStepsTrimmed + 1L
    }
  }
  info$trimTbl <- trimTbl
  .debug("Marginal span trim reason", info$trimReason) # nolint
  .debug("Marginal steps trimmed", info$nStepsTrimmed) # nolint

  keep <- is.finite(x) & x >= currentCut
  info$applied <- sum(!keep) > 0L
  info$reason <- if (info$applied) {
    "dropped_values_left_of_marginal_boundary"
  } else {
    "marginal_filter_kept_all_values"
  }
  info$stopReason <- stopReason
  info$finalStartX <- currentCut
  info$gridShift <- grid$shift
  info$gridSpacing <- grid$spacing
  info$refIndex <- refIndex
  info$refNBin <- refNBin
  info$refNCell <- refNCell
  info$refExpectedResp <- refExpectedResp
  info$refCellsPerBin <- refCellsPerBin
  info$refPurity <- refPurity
  info$cellBinRatio <- cellBinRatio
  info$purityRel <- purityRel
  info$maxLeftBinCells <- maxCells
  info$minLeftBinPurity <- minPurity
  info$scanTbl <- dplyr::bind_rows(scan)
  if (nrow(info$scanTbl) > 0L) {
    info$scanTbl$retained <- info$scanTbl$left >= currentCut
  }

  list(dataMod = .getCpUnsLocSubsetRows(dataMod, keep), info = info)
}
