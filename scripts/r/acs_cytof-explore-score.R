# Analysis 17 (17-explore-acs-cytof-coexpression.qmd): a co-expression score
# for each cell, read from Analysis 9's ACS caches. For markers A, the score is
# -log10 of the product of each marker's survival value under the control
# tube, S_j(x) = P_control(X_j >= x): large when a cell is jointly unusual for
# all markers in A. Survival values are conformal tail p-values, so the
# control tube's own scores give the null distribution, including any
# baseline dependence between markers.

# Tail p-values of `x` against the reference values `ref`: (#ref >= x + 1) /
# (n + 1) for new cells, #ref >= x / n for the reference cells themselves
# (`self = TRUE`, each counts itself). Exact zeros and other values at or
# below the reference minimum get 1.
.acsScoreSurv <- function(ref, x, self = FALSE) {
  ref <- sort(ref[is.finite(ref)])
  n <- length(ref)
  ge <- n - findInterval(x, ref, left.open = TRUE)
  p <- if (self) ge / n else (ge + 1) / (n + 1)
  p[x <= ref[[1]]] <- 1
  pmin(p, 1)
}

# Co-expression score of every cell in both tubes for the markers `markers`.
.acsScoreTbl <- function(dat, markers, donor) {
  one <- function(ex, tube, self) {
    s <- Reduce(`+`, lapply(markers, function(m) {
      -log10(.acsScoreSurv(dat$uns[[m]], ex[[m]], self = self))
    }))
    data.frame(donor = donor, tube = tube, score = s)
  }
  rbind(one(dat$stim, "stimulated", FALSE), one(dat$uns, "unstimulated", TRUE))
}

# Fraction of each tube's cells with a score at or above each threshold, and
# the net (stimulated minus control) fraction.
.acsScoreTail <- function(tbl, grid = seq(0, 8, by = 0.05)) {
  out <- lapply(split(tbl, list(tbl$donor, tbl$tube), drop = TRUE), function(d) {
    s <- sort(d$score)
    data.frame(
      donor = d$donor[[1]], tube = d$tube[[1]], threshold = grid,
      frac = (length(s) - findInterval(grid, s, left.open = TRUE)) / length(s)
    )
  })
  out <- do.call(rbind, out)
  wide <- stats::reshape(out, idvar = c("donor", "threshold"), timevar = "tube",
    direction = "wide")
  wide$net <- wide$frac.stimulated - wide$frac.unstimulated
  list(long = out, net = wide[, c("donor", "threshold", "net")])
}

# KL divergence (bits) of the stimulated from the control score distribution,
# on bins of `binwidth` with add-half smoothing.
.acsScoreKl <- function(tbl, binwidth = 0.25) {
  do.call(rbind, lapply(split(tbl, tbl$donor), function(d) {
    breaks <- seq(0, max(d$score) + binwidth, by = binwidth)
    p <- graphics::hist(d$score[d$tube == "stimulated"], breaks, plot = FALSE)$counts + 0.5
    q <- graphics::hist(d$score[d$tube == "unstimulated"], breaks, plot = FALSE)$counts + 0.5
    p <- p / sum(p)
    q <- q / sum(q)
    data.frame(donor = d$donor[[1]], klBits = sum(p * log2(p / q)))
  }))
}

.acsScoreTubeColours <- c(stimulated = "#B2182B", unstimulated = "#2166AC")

# Score histograms of both tubes as the fraction of each tube's cells per bin
# (log10 y). Cells with score 0 (not high on any marker) are left out of the
# plot but kept in each tube's denominator.
.acsScoreHistTbl <- function(tbl, binwidth = 0.25) {
  breaks <- seq(0, max(tbl$score) + binwidth, by = binwidth)
  do.call(rbind, lapply(split(tbl, list(tbl$donor, tbl$tube), drop = TRUE), function(d) {
    pos <- d$score[d$score > 0]
    counts <- if (length(pos)) graphics::hist(pos, breaks, plot = FALSE)$counts else 0 * breaks[-1]
    data.frame(donor = d$donor[[1]], tube = d$tube[[1]],
      mid = breaks[-1] - binwidth / 2, frac = counts / nrow(d))
  }))
}

.acsScoreHistPlot <- function(tbl, binwidth = 0.25) {
  h <- .acsScoreHistTbl(tbl, binwidth)
  h <- h[h$frac > 0, ]
  ggplot(h, aes(x = .data$mid, y = .data$frac, colour = .data$tube)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    scale_y_log10() +
    scale_colour_manual(values = .acsScoreTubeColours, name = "Tube") +
    facet_wrap(~donor, ncol = 4, axes = "all") +
    .analysis_theme() +
    theme(legend.position = "bottom") +
    labs(x = "Co-expression score (-log10 product of control-tube tail probabilities)",
      y = "Fraction of tube's cells per bin (log scale)")
}

# Fraction of cells at or above each score (log10 y), both tubes.
.acsScoreTailPlot <- function(tail) {
  d <- tail$long[tail$long$frac > 0, ]
  ggplot(d, aes(x = .data$threshold, y = .data$frac, colour = .data$tube)) +
    geom_line(linewidth = 0.6) +
    scale_y_log10() +
    scale_colour_manual(values = .acsScoreTubeColours, name = "Tube") +
    facet_wrap(~donor, ncol = 4, axes = "all") +
    .analysis_theme() +
    theme(legend.position = "bottom") +
    labs(x = "Score threshold", y = "Fraction of cells with score at or above it (log scale)")
}

# Net (stimulated minus control) fraction at or above each score, in percent,
# with the net frequency of cells above both ordinary gates (dashed).
.acsScoreNetPlot <- function(tail, rect) {
  d <- tail$net
  d$net <- 100 * d$net
  ggplot(d, aes(x = .data$threshold, y = .data$net)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_hline(data = rect, aes(yintercept = 100 * .data$net), linetype = "dashed",
      colour = "#D55E00") +
    geom_line(linewidth = 0.6) +
    facet_wrap(~donor, ncol = 4, scales = "free_y", axes = "all") +
    .analysis_theme() +
    labs(x = "Score threshold",
      y = "Net frequency at or above threshold (%)")
}

# Hexagon plots of two markers in both tubes with score contours (levels of
# -log10 S_x(x) - log10 S_y(y), using the control tube's survival values) and
# the ordinary gates.
.acsScoreHexPlot <- function(hexTbl, survGrid, gates, x, y, levels = c(2, 3, 4),
                             contourColour = "white") {
  hexTbl$panel <- paste0(hexTbl$donor, "\n", hexTbl$tube)
  survGrid <- merge(survGrid, unique(hexTbl[c("donor", "tube", "panel")]))
  lines <- do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  }))
  lines <- merge(lines, unique(hexTbl[c("donor", "tube", "panel")]))
  ggplot(hexTbl, aes(x = .data$x, y = .data$y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_contour(data = survGrid, aes(x = .data$gx, y = .data$gy, z = .data$z),
      breaks = levels, colour = contourColour, linewidth = 0.5, inherit.aes = FALSE) +
    geom_vline(data = lines[is.finite(lines$gx), ], aes(xintercept = .data$gx),
      colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines[is.finite(lines$gy), ], aes(yintercept = .data$gy),
      colour = "#D55E00", linewidth = 0.4) +
    facet_wrap(~panel, ncol = 6, axes = "all") +
    .analysis_theme() +
    labs(x = x, y = y)
}

# Score surface on a grid for the contours: the stimulated-cell formula under
# one donor's control tube.
.acsScoreGrid <- function(dat, x, y, donor, n = 120L) {
  gx <- seq(0, max(dat$stim[[x]], dat$uns[[x]]), length.out = n)
  gy <- seq(0, max(dat$stim[[y]], dat$uns[[y]]), length.out = n)
  sx <- -log10(.acsScoreSurv(dat$uns[[x]], gx))
  sy <- -log10(.acsScoreSurv(dat$uns[[y]], gy))
  g <- expand.grid(i = seq_len(n), j = seq_len(n))
  data.frame(donor = donor, gx = gx[g$i], gy = gy[g$j], z = sx[g$i] + sy[g$j])
}

# Net fraction of cells above both ordinary gates (strict x > gate).
.acsScoreRectNet <- function(dat, markers, donor) {
  bothPos <- function(ex) {
    Reduce(`&`, lapply(markers, function(m) (ex[[m]] > dat$gate[[m]]) %in% TRUE))
  }
  data.frame(donor = donor, net = mean(bothPos(dat$stim)) - mean(bothPos(dat$uns)))
}

# ---------------------------------------------------------------------------
# Two-marker co-expression beyond independence, and along its main direction
# ---------------------------------------------------------------------------

# Counts of each tube's cells in square bins of `width` over both markers.
.acsCoexBins <- function(dat, x, y, donor, width = 0.2) {
  maxX <- max(dat$stim[[x]], dat$uns[[x]])
  maxY <- max(dat$stim[[y]], dat$uns[[y]])
  bx <- seq(0, maxX + width, by = width)
  by <- seq(0, maxY + width, by = width)
  one <- function(ex, tube) {
    i <- findInterval(ex[[x]], bx, rightmost.closed = TRUE)
    j <- findInterval(ex[[y]], by, rightmost.closed = TRUE)
    tab <- table(factor(i, levels = seq_len(length(bx) - 1L)),
      factor(j, levels = seq_len(length(by) - 1L)))
    g <- expand.grid(i = seq_len(length(bx) - 1L), j = seq_len(length(by) - 1L))
    data.frame(donor = donor, tube = tube,
      x = bx[g$i] + width / 2, y = by[g$j] + width / 2,
      count = as.vector(tab[cbind(g$i, g$j)]))
  }
  rbind(one(dat$stim, "stimulated"), one(dat$uns, "unstimulated"))
}

# Pearson residuals (O - E) / sqrt(E) with E from each tube's own marginals.
.acsCoexResiduals <- function(bins) {
  do.call(rbind, lapply(split(bins, list(bins$donor, bins$tube), drop = TRUE), function(d) {
    n <- sum(d$count)
    px <- tapply(d$count, d$x, sum) / n
    py <- tapply(d$count, d$y, sum) / n
    d$expected <- as.numeric(n * px[as.character(d$x)] * py[as.character(d$y)])
    d$residual <- ifelse(d$expected > 0, (d$count - d$expected) / sqrt(d$expected), NA)
    d
  }))
}

# Co-expression beyond independence, from bin counts of two markers. In each
# tube, Poisson GAMs give an independence model (log rate = s(x) + s(y), the
# product of non-linear marginals) and an interaction model (adding
# ti(x, y)). Their difference, as a fraction of the tube's cells, is the
# co-expression the marginals do not explain. Only bins where cells are
# positive for at least one marker by its ordinary gate (strict, at the bin
# centre) are kept. Returns the bins of both tubes with both fits.
.acsCoexIndependence <- function(dat, x, y, donor, width = 0.2, k = 10L, kInt = 6L) {
  bins <- .acsCoexBins(dat, x, y, donor, width = width)
  # Few distinct bin centres (tiny ranges) cap the basis sizes.
  kx <- min(k, length(unique(bins$x)) - 1L)
  ky <- min(k, length(unique(bins$y)) - 1L)
  kIx <- min(kInt, kx)
  kIy <- min(kInt, ky)
  fits <- lapply(split(bins, bins$tube), function(d) {
    m0 <- mgcv::gam(count ~ s(x, k = kx) + s(y, k = ky),
      family = stats::poisson(), data = d, method = "REML")
    m1 <- mgcv::gam(count ~ s(x, k = kx) + s(y, k = ky) + ti(x, y, k = c(kIx, kIy)),
      family = stats::poisson(), data = d, method = "REML")
    d$fit0 <- as.numeric(stats::fitted(m0))
    d$fit1 <- as.numeric(stats::fitted(m1))
    d$excessFrac <- (d$fit1 - d$fit0) / sum(d$count)
    d$logRatio <- log2(d$fit1 / d$fit0)
    d
  })
  out <- do.call(rbind, fits)
  out$singlePos <- out$x > dat$gate[[x]] | out$y > dat$gate[[y]]
  out
}

# Per-donor summary: net (stimulated minus control) excess over independence
# in the single-positive region, and the cells called co-expressing: cells
# positive for at least one marker in bins where the stimulated tube has at
# least `minRatio` times the independence expectation and more excess than
# the control tube. Compared with the rectangle of both ordinary gates.
.acsCoexIndependenceSummary <- function(indep, dat, x, y, width = 0.2, minRatio = 2) {
  st <- indep[indep$tube == "stimulated", ]
  un <- indep[indep$tube == "unstimulated", ]
  un <- un[match(paste(st$x, st$y), paste(un$x, un$y)), ]
  net <- st$excessFrac - un$excessFrac
  region <- st$singlePos & st$logRatio >= log2(minRatio) & net > 0
  breaksX <- seq(0, max(st$x) + width, by = width)
  breaksY <- seq(0, max(st$y) + width, by = width)
  binKey <- function(a, b) paste(findInterval(a, breaksX), findInterval(b, breaksY))
  regionKeys <- binKey(st$x, st$y)[region]
  callTube <- function(ex) {
    single <- (ex[[x]] > dat$gate[[x]]) | (ex[[y]] > dat$gate[[y]])
    single & binKey(ex[[x]], ex[[y]]) %in% regionKeys
  }
  inS <- callTube(dat$stim)
  inU <- callTube(dat$uns)
  data.frame(
    donor = st$donor[[1]],
    netExcessPct = 100 * sum(pmax(net[st$singlePos], 0)),
    stimCalled = sum(inS), controlCalled = sum(inU),
    netCalledPct = 100 * (mean(inS) - mean(inU)),
    netRectanglePct = 100 * .acsScoreRectNet(dat, c(x, y), st$donor[[1]])$net
  )
}

# log2(interaction / independence) fit in the single-positive region of both
# tubes, with the ordinary gates.
.acsCoexIndependencePlot <- function(indep, gates, x, y, limit = 3) {
  d <- indep[indep$singlePos & (indep$count > 0 | indep$fit1 >= 0.5), ]
  d$panel <- paste0(d$donor, "\n", d$tube)
  d$logRatio <- pmax(pmin(d$logRatio, limit), -limit)
  lines <- do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  }))
  lines <- merge(lines, unique(d[c("donor", "tube", "panel")]))
  ggplot(d, aes(x = .data$x, y = .data$y, fill = .data$logRatio)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-limit, limit), name = "log2(interaction /\nindependence)") +
    geom_vline(data = lines, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4) +
    facet_wrap(~panel, ncol = 6, axes = "all") +
    .analysis_theme(grid = "none") +
    labs(x = x, y = y)
}

# Pearson residuals of each tube against independence (white = expected).
.acsCoexResidualPlot <- function(res, gates, x, y, limit = 5) {
  res <- res[is.finite(res$residual) & (res$count > 0 | res$expected >= 0.5), ]
  res$panel <- paste0(res$donor, "\n", res$tube)
  res$residual <- pmax(pmin(res$residual, limit), -limit)
  lines <- do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  }))
  lines <- merge(lines, unique(res[c("donor", "tube", "panel")]))
  ggplot(res, aes(x = .data$x, y = .data$y, fill = .data$residual)) +
    geom_tile() +
    geom_vline(data = lines, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-limit, limit), name = "(O - E) / sqrt(E)") +
    facet_wrap(~panel, ncol = 6, axes = "all") +
    .analysis_theme(grid = "none") +
    labs(x = x, y = y)
}

# Difference in departure from independence between the tubes: stimulated
# minus control excess (interaction fit minus independence fit, as a fraction
# of each tube's cells) per bin, in percentage points.
.acsCoexDifference <- function(indep) {
  st <- indep[indep$tube == "stimulated", ]
  un <- indep[indep$tube == "unstimulated", ]
  un <- un[match(paste(st$x, st$y), paste(un$x, un$y)), ]
  st$diffPct <- 100 * (st$excessFrac - un$excessFrac)
  st[c("donor", "x", "y", "singlePos", "diffPct")]
}

.acsCoexDifferencePlot <- function(diff, gates, x, y) {
  d <- diff[abs(diff$diffPct) > 0, ]
  limit <- max(stats::quantile(abs(d$diffPct), 0.99), 1e-3)
  d$diffPct <- pmax(pmin(d$diffPct, limit), -limit)
  lines <- do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  }))
  ggplot(d, aes(x = .data$x, y = .data$y, fill = .data$diffPct)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-limit, limit), name = "Stimulated minus control\nexcess over independence (pp)") +
    geom_vline(data = lines, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4) +
    facet_wrap(~donor, ncol = 4, axes = "all") +
    .analysis_theme(grid = "none") +
    labs(x = x, y = y)
}

# Lowest value a lowered gate may take: the control tube's main negative peak
# plus `mult` robust SDs (left half-width at half maximum / 1.177, on a
# Gaussian KDE with the nrd0 bandwidth of the non-zero values), so a gate is
# never lowered into the bulk of the negatives.
.acsCoexNegFloor <- function(x, mult = 1.5) {
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
.acsCoexPurity <- function(s, u, nS, nU) {
  ifelse(s > 0, 1 - (u / nU) / (s / nS), ifelse(u > 0, -Inf, NA_real_))
}

# z of the net double-positive count (above both ordinary gates):
# (s - u nS / nU) / sqrt(s + u (nS / nU)^2).
.acsCoexDpZ <- function(dat, a, b) {
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
#    (`.acsCoexDpZ()` > `zMin`).
# 2. The distance from b's gate down to its floor (`.acsCoexNegFloor()` of the
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
.acsCoexLowerGate <- function(dat, a, b, nBin = 20L, rMin = 3.5, frac = 0.75, zMin = 2) {
  nS <- nrow(dat$stim)
  nU <- nrow(dat$uns)
  gA <- dat$gate[[a]]
  gB <- dat$gate[[b]]
  floorA <- .acsCoexNegFloor(dat$uns[[a]])
  floorB <- .acsCoexNegFloor(dat$uns[[b]])
  z <- .acsCoexDpZ(dat, a, b)
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
  pDp <- .acsCoexPurity(sum(sA & sB > gB), sum(uA & uB > gB), nS, nU)
  out$purityDp <- pDp
  passBins <- function(lo, hi) {
    inS <- sB > lo & sB <= hi
    s <- sum(sA & inS)
    u <- sum(uA & uB > lo & uB <= hi)
    e <- nS * pA * mean(inS)
    s > 0 && (s - e) / sqrt(max(e, 1e-12)) >= rMin &&
      .acsCoexPurity(s, u, nS, nU) >= frac * pDp
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
    s > 0 && .acsCoexPurity(s, u, nS, nU) >= frac * pDp
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

# Cells positive for both markers: above both ordinary gates, or added by
# either lowered gate (`low[[x]]` lowers x among y+ cells).
.acsCoexLowerBothPos <- function(ex, dat, low, x, y) {
  base <- ex[[x]] > dat$gate[[x]] & ex[[y]] > dat$gate[[y]]
  addX <- ex[[y]] > low[[x]]$condCut & ex[[x]] > low[[x]]$cut
  addY <- ex[[x]] > low[[y]]$condCut & ex[[y]] > low[[y]]$cut
  (base | addX | addY) %in% TRUE
}

# Both lowered gates for one donor and the resulting double positives, against
# the rectangle of ordinary gates.
.acsCoexLowerSummary <- function(dat, x, y, donor, ...) {
  low <- stats::setNames(list(
    .acsCoexLowerGate(dat, a = y, b = x, ...),
    .acsCoexLowerGate(dat, a = x, b = y, ...)
  ), c(x, y))
  sPos <- .acsCoexLowerBothPos(dat$stim, dat, low, x, y)
  uPos <- .acsCoexLowerBothPos(dat$uns, dat, low, x, y)
  list(
    low = low,
    summary = data.frame(donor = donor, zDp = low[[x]]$z,
      gateX = dat$gate[[x]], floorX = low[[x]]$floorB, cutX = low[[x]]$cut,
      condCutY = low[[x]]$condCut,
      gateY = dat$gate[[y]], floorY = low[[y]]$floorB, cutY = low[[y]]$cut,
      condCutX = low[[y]]$condCut,
      stimCalled = sum(sPos), controlCalled = sum(uPos),
      netPct = 100 * (mean(sPos) - mean(uPos)),
      netRectanglePct = 100 * .acsScoreRectNet(dat, c(x, y), donor)$net,
      controlRectangle = sum(.acsCoexLowerBothPos(dat$uns, dat,
        lapply(low, function(l) list(cut = Inf, condCut = Inf)), x, y)))
  )
}

# Line segments for the lowered gates on x-against-y plots: x's lowered gate
# (vertical, from y's raised cut up) and y's (horizontal, from x's raised cut
# right), drawn only where a gate was lowered.
.acsCoexLowerSegments <- function(lower, x, y) {
  do.call(rbind, lapply(lower, function(l) {
    s <- l$summary
    rbind(
      if (s$cutX < s$gateX) data.frame(donor = s$donor, x = s$cutX, xend = s$cutX, y = s$condCutY, yend = Inf),
      if (s$cutY < s$gateY) data.frame(donor = s$donor, x = s$condCutX, xend = Inf, y = s$cutY, yend = s$cutY)
    )
  }))
}

# Hexagon plots with the ordinary gates (orange) and the lowered gates
# (blue), with their raised cuts on the other marker.
.acsCoexLowerHexPlot <- function(hexTbl, lower, gates, x, y) {
  hexTbl$panel <- paste0(hexTbl$donor, "\n", hexTbl$tube)
  panels <- unique(hexTbl[c("donor", "tube", "panel")])
  base <- merge(do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  })), panels)
  seg <- .acsCoexLowerSegments(lower, x, y)
  p <- ggplot(hexTbl, aes(x = .data$x, y = .data$y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_vline(data = base, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = base, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4)
  if (!is.null(seg) && nrow(seg)) {
    p <- p + geom_segment(data = merge(seg, panels),
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      colour = "#0072B2", linewidth = 0.6, inherit.aes = FALSE)
  }
  p +
    facet_wrap(~panel, ncol = 6, axes = "all") +
    .analysis_theme() +
    labs(x = x, y = y)
}

# Histograms of `b` among cells above `a`'s ordinary gate, in both tubes, as
# cells per 100,000 cells of the tube (so the bars compare directly: where
# the control bar is a quarter of the stimulated bar, purity is 0.75). Exact
# zeros are left out of the bars but counted in the tube sizes. Lines: b's
# ordinary gate (black dashed), the lowered gate (blue; drawn when lowered),
# and b's floor (grey dotted).
.acsCoexCondHistPlot <- function(dat, lower, a, b, binwidth = 0.1) {
  tbl <- do.call(rbind, lapply(names(dat), function(donor) {
    d <- dat[[donor]]
    one <- function(ex, tube) {
      v <- ex[[b]][(ex[[a]] > d$gate[[a]]) %in% TRUE]
      v <- v[v > 0]
      if (!length(v)) {
        return(NULL)
      }
      data.frame(donor = donor, tube = tube, x = v, w = 1e5 / nrow(ex))
    }
    rbind(one(d$stim, "stimulated"), one(d$uns, "unstimulated"))
  }))
  lines <- do.call(rbind, lapply(names(lower), function(donor) {
    l <- lower[[donor]]$low[[b]]
    stopifnot(identical(l$a, a))
    data.frame(donor = donor, gate = l$gateB, cut = l$cut, floor = l$floorB,
      lowered = l$cut < l$gateB)
  }))
  ggplot(tbl, aes(x = .data$x, weight = .data$w, fill = .data$tube, colour = .data$tube)) +
    geom_histogram(binwidth = binwidth, boundary = 0, position = "identity",
      alpha = 0.35, linewidth = 0.3) +
    geom_vline(data = lines, aes(xintercept = .data$gate), colour = "grey20",
      linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = lines[lines$lowered, ], aes(xintercept = .data$cut),
      colour = "#0072B2", linewidth = 0.6) +
    geom_vline(data = lines, aes(xintercept = .data$floor), colour = "grey50",
      linetype = "dotted", linewidth = 0.6) +
    scale_y_sqrt() +
    scale_colour_manual(values = .acsScoreTubeColours, name = "Tube") +
    scale_fill_manual(values = .acsScoreTubeColours, name = "Tube") +
    facet_wrap(~donor, ncol = 4, scales = "free_y", axes = "all") +
    .analysis_theme() +
    theme(legend.position = "bottom") +
    labs(x = paste0(b, " among ", a, "+ cells (non-zero values)"),
      y = "Cells per 100,000 in the tube (square-root scale)")
}
