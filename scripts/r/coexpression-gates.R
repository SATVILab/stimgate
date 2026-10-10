# Cytokine-positive gates lowered where two cytokines are co-expressed beyond
# independence (Analyses 11, 17 and 18). Each tube pair is a list with `stim`
# and `uns` data frames (one column per marker, values on the gating scale)
# and `gate`, the ordinary gates named by marker. For an ordered pair (a, b),
# b's gate is lowered for cells above a's gate (`.coexLowerGate()`); with
# several markers every ordered pair is evaluated (`.coexLowerGates()`) and a
# cell is b-positive if above b's ordinary gate or, for any a, above the
# raised a cut and the lowered b gate (`.coexPositive()`).

# Lowest value a lowered gate may take: the control tube's main negative peak
# plus `mult` robust SDs (left half-width at half maximum / 1.177, on a
# Gaussian KDE with the nrd0 bandwidth of the non-zero values), so a gate is
# never lowered into the bulk of the negatives.
.coexNegFloor <- function(x, mult = 1.5) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10L) {
    return(0)
  }
  d <- stats::density(x, bw = stats::bw.nrd0(x), n = 2048)
  i <- which.max(d$y)
  left <- which(d$x < d$x[[i]] & d$y <= d$y[[i]] / 2)
  sigma <- if (length(left)) (d$x[[i]] - d$x[[max(left)]]) / 1.177 else stats::mad(x)
  max(0, d$x[[i]] + mult * sigma)
}

# Purity of a set of stimulated cells: the estimated share that responded,
# 1 - (control share) / (stimulated share), with shares over each tube's cells
# (zeros included). Undefined (NA) with no cells in either tube and -Inf with
# control cells only.
.coexPurity <- function(s, u, nS, nU) {
  ifelse(s > 0, 1 - (u / nU) / (s / nS), ifelse(u > 0, -Inf, NA_real_))
}

# z of the net double-positive count (above both ordinary gates):
# (s - u nS / nU) / sqrt(s + u (nS / nU)^2).
.coexDpZ <- function(dat, a, b) {
  nS <- nrow(dat$stim)
  nU <- nrow(dat$uns)
  s <- sum(dat$stim[[a]] > dat$gate[[a]] & dat$stim[[b]] > dat$gate[[b]], na.rm = TRUE)
  u <- sum(dat$uns[[a]] > dat$gate[[a]] & dat$uns[[b]] > dat$gate[[b]], na.rm = TRUE)
  (s - u * nS / nU) / sqrt(max(s + u * (nS / nU)^2, 1e-12))
}

# Lowered gate for marker `b` among cells above `a`'s ordinary gate, and the
# raised `a` cut for the cells it adds.
#
# 1. Nothing moves unless the double-positive response is clear
#    (`.coexDpZ()` > `zMin`).
# 2. The distance from b's gate down to its floor (`.coexNegFloor()` of the
#    control tube) is cut into `nBin` bins. A set of bins passes when its a+
#    stimulated cells depart from independence (Pearson residual (O - E) /
#    sqrt(E), E from the stimulated tube's own a+ share and share of cells in
#    the bins) by at least `rMin` and its purity is at least `frac` times the
#    double positives'. Working back from the floor, the cut is the first bin
#    edge where that bin passes on its own and all bins from the gate down to
#    it pass together.
# 3. Band trim, purity only: the distance from a's gate down to its floor, in
#    `nBin` steps, sets slices above a's gate. Working up from a's gate, the a
#    cut is the first slice edge where that slice of the added band passes on
#    its own and the band above it passes together. If the band runs out, b's
#    gate is not lowered.
#
# A cell is then b-positive if above b's gate, or above `condCut` on a and
# `cut` on b.
.coexLowerGate <- function(dat, a, b, nBin = 20L, rMin = 3.5, frac = 0.75, zMin = 2) {
  nS <- nrow(dat$stim)
  nU <- nrow(dat$uns)
  gA <- dat$gate[[a]]
  gB <- dat$gate[[b]]
  floorA <- .coexNegFloor(dat$uns[[a]])
  floorB <- .coexNegFloor(dat$uns[[b]])
  z <- .coexDpZ(dat, a, b)
  out <- list(a = a, b = b, gateA = gA, gateB = gB, cut = gB, condCut = gA,
    floorA = floorA, floorB = floorB, z = z, purityDp = NA_real_)
  if (!is.finite(gA) || !is.finite(gB) || !(z > zMin) || floorB >= gB) {
    return(out)
  }
  sAv <- dat$stim[[a]]
  uAv <- dat$uns[[a]]
  sB <- dat$stim[[b]]
  uB <- dat$uns[[b]]
  sA <- (sAv > gA) %in% TRUE
  uA <- (uAv > gA) %in% TRUE
  pA <- mean(sA)
  pDp <- .coexPurity(sum(sA & sB > gB), sum(uA & uB > gB), nS, nU)
  out$purityDp <- pDp
  passBins <- function(lo, hi) {
    inS <- sB > lo & sB <= hi
    s <- sum(sA & inS)
    u <- sum(uA & uB > lo & uB <= hi)
    e <- nS * pA * mean(inS)
    s > 0 && (s - e) / sqrt(max(e, 1e-12)) >= rMin &&
      .coexPurity(s, u, nS, nU) >= frac * pDp
  }
  edges <- gB - (gB - floorB) * seq_len(nBin) / nBin
  upper <- c(gB, edges[-nBin])
  for (k in rev(seq_len(nBin))) {
    if (passBins(edges[[k]], upper[[k]]) && passBins(edges[[k]], gB)) {
      out$cut <- edges[[k]]
      break
    }
  }
  if (out$cut >= gB) {
    return(out)
  }
  bandS <- sB > out$cut & sB <= gB
  bandU <- uB > out$cut & uB <= gB
  passBand <- function(lo, hi) {
    s <- sum(bandS & sAv > lo & sAv <= hi)
    u <- sum(bandU & uAv > lo & uAv <= hi)
    s > 0 && .coexPurity(s, u, nS, nU) >= frac * pDp
  }
  step <- max(gA - floorA, 1e-6) / nBin
  maxA <- max(sAv[bandS & sA])
  lo <- gA
  while (lo < maxA && !(passBand(lo, lo + step) && passBand(lo, Inf))) {
    lo <- lo + step
  }
  if (lo >= maxA) {
    out$cut <- gB
    return(out)
  }
  out$condCut <- lo
  out
}

# Every ordered pair of `markers` with finite ordinary gates: one row per pair
# with the lowered gate (`cut`) for b among a+ cells and the raised a cut
# (`condCut`). Pairs with a missing gate keep b's ordinary gate.
.coexLowerGates <- function(dat, markers = names(dat$gate), ...) {
  pairs <- expand.grid(a = markers, b = markers, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$a != pairs$b, , drop = FALSE]
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    l <- .coexLowerGate(dat, a = pairs$a[[i]], b = pairs$b[[i]], ...)
    data.frame(a = l$a, b = l$b, gateA = l$gateA, gateB = l$gateB,
      cut = l$cut, condCut = l$condCut, floorA = l$floorA, floorB = l$floorB,
      z = l$z, purityDp = l$purityDp, lowered = isTRUE(l$cut < l$gateB),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Positivity of every cell in `ex` (a data frame with the markers as columns)
# for each marker, with the ordinary gates `gate` and the lowered gates
# `low` (`.coexLowerGates()`; NULL for the ordinary gates alone). Strict
# `x > gate`; a marker without a finite gate is negative throughout. Returns a
# logical matrix, one column per marker.
.coexPositive <- function(ex, gate, low = NULL, markers = names(gate)) {
  base <- vapply(markers, function(m) {
    (ex[[m]] > gate[[m]]) %in% TRUE
  }, logical(nrow(ex)))
  if (is.null(dim(base))) base <- matrix(base, nrow = nrow(ex), dimnames = list(NULL, markers))
  out <- base
  if (!is.null(low)) {
    low <- low[low$lowered & low$a %in% markers & low$b %in% markers, , drop = FALSE]
    for (i in seq_len(nrow(low))) {
      a <- low$a[[i]]
      b <- low$b[[i]]
      out[, b] <- out[, b] | ((ex[[a]] > low$condCut[[i]]) %in% TRUE & (ex[[b]] > low$cut[[i]]) %in% TRUE)
    }
  }
  out
}

# Counts of each combination of positive markers (all 2^k, in a fixed order),
# named like "IFNg+IL2-TNF+".
.coexCombnCounts <- function(pos) {
  markers <- colnames(pos)
  code <- as.vector(pos %*% (2L^(seq_along(markers) - 1L)))
  n <- tabulate(code + 1L, nbins = 2L^length(markers))
  lab <- vapply(seq_along(n) - 1L, function(k) {
    paste0(markers, ifelse(bitwAnd(k, 2L^(seq_along(markers) - 1L)) > 0, "+", "-"), collapse = "")
  }, character(1))
  stats::setNames(n, lab)
}
