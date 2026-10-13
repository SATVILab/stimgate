# Plateau step (`locPlateauStep`)
#
# After a tube's local-FDR gate is chosen (after the cap), step it down in
# bandwidth-sized steps towards the left edge of the taut-string plateau that
# holds the gate, accepting each step only while the cells it adds are pure
# enough. The step is used only for tubes whose stimulated and unstimulated
# negative tails agree: a first gating pass over every stimulated tube of the
# marker records the two tails' cell counts below each tube's plateau edge, and
# the step is switched off for the marker, or for single tubes, when they do
# not agree. The marker is then gated again with the step allowed.

# Settings fixed by design, kept together here.
#' @keywords internal
.locPlateauConst <- function() {
  list(
    # alignment band, in bandwidths below the plateau edge
    bandFar = 8,
    bandNear = 4,
    # tubes with fewer unstimulated cells in the band are judged by the group
    bandMinCells = 20,
    # one-sided normal quantile for single-tube exclusion and the group bound
    z = stats::qnorm(0.95),
    # failed steps in a row that end the scan
    maxReject = 2L
  )
}

# Counts of stimulated and unstimulated cells in the alignment band below a
# tube's plateau edge, using unstimulated values without the bias.
#' @keywords internal
.getCpUnsLocPlateauBand <- function(xStim, xUns, edge, bw) {
  const <- .locPlateauConst()
  out <- list(
    edge = edge, nStim = NA_real_, nUns = NA_real_,
    nStimTotal = length(xStim), nUnsTotal = length(xUns)
  )
  if (!isTRUE(is.finite(edge) && is.finite(bw) && bw > 0)) {
    return(out)
  }
  lo <- edge - const$bandFar * bw
  hi <- edge - const$bandNear * bw
  out$nStim <- sum(xStim >= lo & xStim < hi)
  out$nUns <- sum(xUns >= lo & xUns < hi)
  out
}

# Decide, over all stimulated tubes of a marker, where the step may be used.
# `bandTbl` has one row per tube: ind, nStim, nUns, nStimTotal, nUnsTotal.
# The marker passes when the pooled stimulated/unstimulated ratio (random
# effects on the log scale, weighting tubes by counting precision) is at most
# `alignMax` and its one-sided upper bound for a new tube, allowing for excess
# between-tube spread, is at most `alignUpper`. A passing marker still skips
# tubes whose own ratio is clearly above `alignMax`; tubes with too few cells
# in the band follow the marker.
#' @keywords internal
.getCpUnsLocPlateauAlign <- function(bandTbl, alignMax, alignUpper) {
  const <- .locPlateauConst()
  tbl <- bandTbl
  tbl$ratio <- (tbl$nStim / tbl$nStimTotal) / (tbl$nUns / tbl$nUnsTotal)
  usable <- is.finite(tbl$nUns) & tbl$nUns >= const$bandMinCells &
    is.finite(tbl$nStim) & tbl$nStim > 0
  tbl$usable <- usable
  y <- log(tbl$ratio[usable])
  v <- 1 / tbl$nStim[usable] + 1 / tbl$nUns[usable]
  group <- list(
    nTube = nrow(tbl), nUsable = sum(usable), pooled = NA_real_,
    tau = NA_real_, upper = NA_real_, pass = FALSE
  )
  if (length(y) >= 2L) {
    w <- 1 / v
    yFix <- sum(w * y) / sum(w)
    q <- sum(w * (y - yFix)^2)
    tau2 <- max(0, (q - (length(y) - 1)) / (sum(w) - sum(w^2) / sum(w)))
    wR <- 1 / (v + tau2)
    yR <- sum(wR * y) / sum(wR)
    group$pooled <- exp(yR)
    group$tau <- sqrt(tau2)
    group$upper <- exp(yR + const$z * sqrt(tau2))
    group$pass <- group$pooled <= alignMax && group$upper <= alignUpper
  }
  tubeHigh <- rep(FALSE, nrow(tbl))
  tubeHigh[usable] <- (y - log(alignMax)) / sqrt(v) > const$z
  tbl$tubeHigh <- tubeHigh
  tbl$allowed <- group$pass & !tubeHigh
  list(tbl = tbl, group = group)
}

# Left edge of the taut-string plateau holding `gate`, fitted to the sorted
# stimulated values at or above `lowerX`. NA when there is no usable edge.
#' @keywords internal
.getCpUnsLocPlateauEdge <- function(xStim, gate, lowerX) {
  sub <- sort(xStim[is.finite(xStim) & xStim >= lowerX])
  den <- .getCpUnsLocAntimodeDensity(sub)
  if (is.null(den) || all(den$y == 0)) {
    return(NA_real_)
  }
  runs <- rle(den$y)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  j <- findInterval(gate, den$x)
  k <- which(starts <= j & ends >= j)[1L]
  if (is.na(k) || starts[k] <= 1L) {
    return(NA_real_)
  }
  edge <- sub[starts[k]]
  if (edge >= gate) NA_real_ else edge
}

# Step the gate down from `gate` towards `edge` in steps of `bw`. Each step is
# tested on the span from the last accepted cut down to the step's lower end;
# it is accepted when it holds at least `minCells` stimulated cells and its
# purity (stimulated excess over scaled unstimulated cells, as a share of
# stimulated cells) is at least `purityRel` times the purity above `gate`.
# The scan stops after two failures in a row,
# on reaching the edge, or after crossing `gapBw` bandwidths of steps without
# stimulated cells. `xUns` includes any unstimulated bias.
#' @keywords internal
.getCpUnsLocPlateauStepGate <- function(
  xStim, xUns, gate, edge, bw, purityRel, minCells = 2, gapBw = 1
) {
  const <- .locPlateauConst()
  out <- list(gate = gate, nAccepted = 0L, reason = "not_moved")
  bw <- suppressWarnings(as.numeric(bw)[1L])
  if (!isTRUE(is.finite(edge) && is.finite(bw) && bw > 0 && edge < gate)) {
    out$reason <- "no_plateau_edge"
    return(out)
  }
  ratio <- length(xStim) / length(xUns)
  nStimRef <- sum(xStim > gate)
  refPurity <- (nStimRef - sum(xUns > gate) * ratio) / nStimRef
  if (!is.finite(refPurity) || refPurity <= 0) {
    out$reason <- "no_reference_purity"
    return(out)
  }
  bar <- purityRel * refPurity
  cut <- gate
  hi <- gate
  nReject <- 0L
  emptyWidth <- 0
  while (nReject < const$maxReject) {
    lo <- max(hi - bw, edge)
    if (lo >= hi - 1e-12) break
    if (!any(xStim > lo & xStim <= hi)) {
      emptyWidth <- emptyWidth + (hi - lo)
      if (emptyWidth >= gapBw * bw - 1e-12) {
        out$reason <- "empty_gap"
        break
      }
    } else {
      emptyWidth <- 0
    }
    nS <- sum(xStim > lo & xStim <= cut)
    nU <- sum(xUns > lo & xUns <= cut) * ratio
    purity <- if (nS >= max(minCells, 1)) (nS - nU) / nS else -Inf
    if (purity >= bar) {
      cut <- lo
      nReject <- 0L
      out$nAccepted <- out$nAccepted + 1L
    } else {
      nReject <- nReject + 1L
    }
    hi <- lo
    if (lo <= edge) break
  }
  out$gate <- cut
  if (cut < gate && identical(out$reason, "not_moved")) {
    out$reason <- "moved"
  } else if (cut < gate) {
    out$reason <- paste0("moved_then_", out$reason)
  }
  out
}

# Apply the step to one tube's chosen gate. In the first gating pass
# (`chnlSettings$locPlateauRecord` is an environment) only record the tube's
# plateau edge and alignment-band counts; in the second, step the gate where
# `chnlSettings$locPlateauAllowed` allows it. Returns the (possibly lowered)
# gate with its provenance.
#' @keywords internal
.getCpUnsLocPlateauApply <- function(
  cpObj, ind, dataMod, exTblStimOrig, exTblUnsOrig, bias, chnlSettings
) {
  record <- chnlSettings$locPlateauRecord
  allowedVec <- chnlSettings$locPlateauAllowed
  allowed <- if (is.null(allowedVec)) NA else allowedVec[as.character(ind)]
  info <- list(
    applied = FALSE, allowed = as.logical(allowed)[1L], from = cpObj$cp,
    edge = NA_real_, bw = NA_real_, nAccepted = 0L, reason = "not_requested"
  )
  if (is.null(record) && !isTRUE(info$allowed)) {
    info$reason <- if (is.na(info$allowed)) "no_alignment_record" else "tails_not_aligned"
    return(list(cpObj = cpObj, info = info))
  }
  if (!isTRUE(cpObj$locGeneratedDirect) || !is.finite(cpObj$cp)) {
    info$reason <- "no_direct_gate"
    return(list(cpObj = cpObj, info = info))
  }
  bw <- suppressWarnings(as.numeric(unlist(attr(dataMod, "locDensityBw"))))
  bw <- stats::median(bw[is.finite(bw)])
  lowerX <- suppressWarnings(as.numeric(attr(dataMod, "minProbXPos"))[1L])
  info$bw <- bw
  if (!is.finite(bw) || !is.finite(lowerX)) {
    info$reason <- "no_bandwidth_or_lower_bound"
    return(list(cpObj = cpObj, info = info))
  }
  xStim <- .getCut(exTblStimOrig)
  xStim <- xStim[is.finite(xStim)]
  xUns <- .getCut(exTblUnsOrig)
  xUns <- xUns[is.finite(xUns)] + bias
  edge <- .getCpUnsLocPlateauEdge(xStim, cpObj$cp, lowerX)
  info$edge <- edge
  if (!is.null(record)) {
    # Raw unstimulated cells: the unstimulated bias (about one bandwidth by
    # default) would otherwise hide a heavier stimulated tail in the band.
    band <- .getCpUnsLocPlateauBand(xStim, xUns - bias, edge, bw)
    assign(as.character(ind), band, envir = record)
    info$reason <- "recorded"
    return(list(cpObj = cpObj, info = info))
  }
  step <- .getCpUnsLocPlateauStepGate(
    xStim = xStim,
    xUns = xUns,
    gate = cpObj$cp,
    edge = edge,
    bw = bw,
    purityRel = chnlSettings$locMarginalPurityRel %||% 0.5,
    minCells = chnlSettings$locPlateauMinCells %||% 2,
    gapBw = chnlSettings$locPlateauGapBw %||% 1
  )
  info$nAccepted <- step$nAccepted
  info$reason <- step$reason
  if (step$gate < cpObj$cp) {
    cpObj$cp <- step$gate
    info$applied <- TRUE
  }
  list(cpObj = cpObj, info = info)
}

# First gating pass for the plateau step: gate every tube with the step off,
# recording each stimulated tube's alignment-band counts below its plateau
# edge, then decide where the step may be used. The random-number state is
# restored afterwards, so the second pass draws exactly what a run without the
# step would. Returns chnlSettings with `locPlateauAllowed` (named by tube).
#' @keywords internal
.gateChnlPlateauAlign <- function(
  .data, indBatchList, chnlSettings, stage, pathProject
) {
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed <- if (hadSeed) get(".Random.seed", envir = globalenv())
  record <- new.env(parent = emptyenv())
  settingsRecord <- chnlSettings
  settingsRecord$locPlateauRecord <- record
  invisible(.gateChnlPreAdjGatesGate(
    indBatchList = indBatchList,
    .data = .data,
    chnlSettings = settingsRecord,
    stage = stage,
    pathProject = pathProject
  ))
  if (hadSeed) {
    assign(".Random.seed", seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  indStim <- unlist(lapply(indBatchList, function(x) as.character(x[-1L])))
  rows <- purrr::map_df(indStim, function(ind) {
    band <- record[[ind]]
    tibble::tibble(
      ind = ind,
      edge = band$edge %||% NA_real_,
      nStim = band$nStim %||% NA_real_,
      nUns = band$nUns %||% NA_real_,
      nStimTotal = band$nStimTotal %||% NA_real_,
      nUnsTotal = band$nUnsTotal %||% NA_real_
    )
  })
  align <- .getCpUnsLocPlateauAlign(
    rows,
    alignMax = chnlSettings$locPlateauAlignMax %||% 1.1,
    alignUpper = chnlSettings$locPlateauAlignUpper %||% 1.25
  )
  .intSaveNm(
    "locPlateauAlign",
    align,
    "all",
    file.path(stage, chnlSettings$chnlCut),
    pathProject
  )
  .debug("plateau step alignment", unlist(align$group)) # nolint
  chnlSettings$locPlateauAllowed <- stats::setNames(
    align$tbl$allowed, align$tbl$ind
  )
  chnlSettings
}
