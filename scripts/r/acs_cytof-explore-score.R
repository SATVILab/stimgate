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

# F-beta-style cytokine-positive gate: lower B's gate for cells above A's
# gate to the threshold that maximises an estimated F-beta of co-expression.
# For each candidate t from B's gate down to `floorB` (default
# `.acsCoexNegFloor()`, step `step`): S = stimulated A+ cells with B > t, and
# in each tube the excess over independence, A+ cells with B > t minus
# n * P(A+) * P(B > t) from that tube's marginals. TP = max(stimulated excess -
# control excess scaled to the stimulated tube's size, 0): co-expression
# beyond independence that the stimulation adds, so responders for A whose B
# is independent add nothing. Precision = TP / S, recall = TP / max_t TP. The
# highest t with the largest F-beta is chosen, and kept only if TP there is
# clear (z = TP / sqrt(S + U_ctrl (nS / nU)^2) > `zMin`, with U_ctrl the
# control A+ cells with B > t); otherwise B's gate is kept.
.acsCoexFbetaGate <- function(dat, a, b, beta = 1, step = 0.05, floorB = NULL, zMin = 2) {
  gateA <- dat$gate[[a]]
  gateB <- dat$gate[[b]]
  if (is.null(floorB)) floorB <- .acsCoexNegFloor(dat$uns[[b]])
  out <- data.frame(conditionOn = a, lowered = b, gateB = gateB, gateCytB = gateB,
    precision = NA_real_, fbeta = NA_real_, z = NA_real_)
  if (!is.finite(gateA) || !is.finite(gateB) || gateB < floorB) {
    return(out)
  }
  grid <- rev(seq(floorB, gateB, by = step))
  grid <- unique(c(gateB, grid))
  scale <- nrow(dat$stim) / nrow(dat$uns)
  excess <- function(ex) {
    aPos <- (ex[[a]] > gateA) %in% TRUE
    vapply(grid, function(t) {
      bPos <- ex[[b]] > t
      c(sum(aPos & bPos), sum(aPos & bPos) - nrow(ex) * mean(aPos) * mean(bPos))
    }, numeric(2))
  }
  eS <- excess(dat$stim)
  eU <- excess(dat$uns)
  S <- eS[1, ]
  tp <- pmax(eS[2, ] - eU[2, ] * scale, 0)
  if (!any(tp > 0)) {
    return(out)
  }
  precision <- ifelse(S > 0, tp / S, 0)
  recall <- tp / max(tp)
  f <- ifelse(precision + recall > 0,
    (1 + beta^2) * precision * recall / (beta^2 * precision + recall), 0)
  k <- which(f >= max(f) - 1e-12)[[1L]]
  z <- tp[[k]] / sqrt(max(S[[k]] + eU[1, k] * scale^2, 1e-12))
  out$precision <- precision[[k]]
  out$fbeta <- f[[k]]
  out$z <- z
  if (z > zMin) {
    out$gateCytB <- grid[[k]]
  }
  out
}

# Two-marker region growing. Bins (`width` wide) above both ordinary gates
# form the starting region (the rectangle). A bin above both markers' floors
# (`.acsCoexNegFloor()` of the control tube) and above at least one ordinary
# gate is acceptable when (a) the stimulated tube's residual against
# independence (from its own marginals) is at least `rMin` and (b) the
# background is low: its share of stimulated cells is at least `ratioMin`
# times the control tube's ((count + 0.5) / size) and it holds at most
# `maxCtrl` control cells, and holds at least `minStim` stimulated cells
# (sparser bins are neutral). An acceptable bin joins the region if it can be
# reached from the rectangle through neighbouring bins (8-connected, within
# the candidate area) crossing at most `maxRejected` rejected bins; neutral and
# acceptable bins cost nothing to cross. Returns the bins with their status,
# the floors and a per-donor summary of the cells called in both tubes.
.acsCoexGrow <- function(dat, x, y, donor, width = 0.2, rMin = 2, ratioMin = 2,
                         maxRejected = 2L, maxCtrl = 2L, minStim = 1L) {
  bins <- .acsCoexBins(dat, x, y, donor, width = width)
  res <- .acsCoexResiduals(bins)
  st <- res[res$tube == "stimulated", ]
  un <- bins[bins$tube == "unstimulated", ]
  st$ctrl <- un$count[match(paste(st$x, st$y), paste(un$x, un$y))]
  nS <- nrow(dat$stim)
  nU <- nrow(dat$uns)
  floorX <- .acsCoexNegFloor(dat$uns[[x]])
  floorY <- .acsCoexNegFloor(dat$uns[[y]])
  gx <- dat$gate[[x]]
  gy <- dat$gate[[y]]
  # Bins by their lower edges.
  st$x0 <- st$x - width / 2
  st$y0 <- st$y - width / 2
  st$rect <- st$x0 >= gx - 1e-9 & st$y0 >= gy - 1e-9
  # Candidates lie above both floors and above at least one ordinary gate, so
  # the region never grows into the negative bulk.
  st$domain <- st$x0 >= floorX - width & st$y0 >= floorY - width &
    (st$x0 + width > gx | st$y0 + width > gy)
  st$acceptable <- st$domain & st$count >= minStim & is.finite(st$residual) &
    st$residual >= rMin & st$count / nS >= ratioMin * (st$ctrl + 0.5) / nU &
    st$ctrl <= maxCtrl
  # Bins with fewer than `minStim` stimulated cells are neutral, like empty ones.
  st$rejected <- st$domain & st$count >= minStim & !st$acceptable & !st$rect
  # 0-1 shortest path (cost 1 to enter a rejected bin) from the rectangle.
  i <- round(st$x0 / width)
  j <- round(st$y0 / width)
  idx <- matrix(NA_integer_, max(i) + 1L, max(j) + 1L)
  idx[cbind(i + 1L, j + 1L)] <- seq_len(nrow(st))
  cost <- rep(Inf, nrow(st))
  start <- which(st$rect)
  cost[start] <- 0
  queue <- start
  while (length(queue)) {
    k <- queue[[1L]]
    queue <- queue[-1L]
    for (di in -1:1) for (dj in -1:1) {
      if (di == 0 && dj == 0) next
      ni <- i[[k]] + di + 1L
      nj <- j[[k]] + dj + 1L
      if (ni < 1L || nj < 1L || ni > nrow(idx) || nj > ncol(idx)) next
      m <- idx[ni, nj]
      if (is.na(m) || !st$domain[[m]] || st$rect[[m]]) next
      newCost <- cost[[k]] + as.integer(st$rejected[[m]])
      if (newCost < cost[[m]] && newCost <= maxRejected) {
        cost[[m]] <- newCost
        if (st$rejected[[m]]) queue <- c(queue, m) else queue <- c(m, queue)
      }
    }
  }
  st$region <- st$rect | (st$acceptable & cost <= maxRejected)
  breaksX <- seq(0, max(st$x) + width, by = width)
  breaksY <- seq(0, max(st$y) + width, by = width)
  binKey <- function(a, b) paste(findInterval(a, breaksX), findInterval(b, breaksY))
  rk <- binKey(st$x, st$y)[st$region]
  inS <- binKey(dat$stim[[x]], dat$stim[[y]]) %in% rk
  inU <- binKey(dat$uns[[x]], dat$uns[[y]]) %in% rk
  list(
    bins = st, floorX = floorX, floorY = floorY,
    summary = data.frame(donor = donor, floorX = floorX, floorY = floorY,
      binsAdded = sum(st$region & !st$rect), stimCalled = sum(inS),
      controlCalled = sum(inU), netCalledPct = 100 * (mean(inS) - mean(inU)),
      netRectanglePct = 100 * .acsScoreRectNet(dat, c(x, y), donor)$net)
  )
}

# Hexagon plots with the grown region: rectangle and added bins outlined
# (black tiles), the ordinary gates (orange) and the F-beta gates (dashed),
# each drawn only where it applies.
.acsCoexGrowHexPlot <- function(hexTbl, grow, gates, fbeta, x, y, width = 0.2) {
  hexTbl$panel <- paste0(hexTbl$donor, "\n", hexTbl$tube)
  panels <- unique(hexTbl[c("donor", "tube", "panel")])
  region <- do.call(rbind, lapply(grow, function(g) {
    b <- g$bins[g$bins$region & !g$bins$rect, c("donor", "x", "y")]
    if (nrow(b)) b else NULL
  }))
  base <- merge(do.call(rbind, lapply(names(gates), function(donor) {
    data.frame(donor = donor, gx = gates[[donor]][[x]], gy = gates[[donor]][[y]])
  })), panels)
  p <- ggplot(hexTbl, aes(x = .data$x, y = .data$y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10", name = "Cells")
  if (!is.null(region) && nrow(region)) {
    p <- p + geom_tile(data = merge(region, panels), aes(x = .data$x, y = .data$y),
      width = width, height = width, fill = NA, colour = "black", linewidth = 0.25,
      inherit.aes = FALSE)
  }
  xMax <- max(hexTbl$x)
  yMax <- max(hexTbl$y)
  lowY <- merge(fbeta[fbeta$gateCytY < fbeta$gateY, c("donor", "gateX", "gateCytY")], panels)
  lowX <- merge(fbeta[fbeta$gateCytX < fbeta$gateX, c("donor", "gateY", "gateCytX")], panels)
  if (nrow(lowY)) {
    p <- p + geom_segment(data = lowY,
      aes(x = .data$gateX, xend = xMax, y = .data$gateCytY, yend = .data$gateCytY),
      linetype = "dashed", colour = "grey30", linewidth = 0.5, inherit.aes = FALSE)
  }
  if (nrow(lowX)) {
    p <- p + geom_segment(data = lowX,
      aes(x = .data$gateCytX, xend = .data$gateCytX, y = .data$gateY, yend = yMax),
      linetype = "dashed", colour = "grey30", linewidth = 0.5, inherit.aes = FALSE)
  }
  p +
    geom_vline(data = base, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = base, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4) +
    facet_wrap(~panel, ncol = 6, axes = "all") +
    .analysis_theme() +
    labs(x = x, y = y)
}

# F-beta gates in both directions and the resulting net double-positive
# frequency (x > gate, or x > gateCyt when positive for the other marker).
.acsCoexFbetaSummary <- function(dat, x, y, donor, ...) {
  gx <- .acsCoexFbetaGate(dat, a = y, b = x, ...)
  gy <- .acsCoexFbetaGate(dat, a = x, b = y, ...)
  dp <- function(ex) {
    xBase <- (ex[[x]] > dat$gate[[x]]) %in% TRUE
    yBase <- (ex[[y]] > dat$gate[[y]]) %in% TRUE
    xPos <- xBase | (yBase & (ex[[x]] > gx$gateCytB) %in% TRUE)
    yPos <- yBase | (xBase & (ex[[y]] > gy$gateCytB) %in% TRUE)
    xPos & yPos
  }
  data.frame(donor = donor,
    gateX = dat$gate[[x]], gateCytX = gx$gateCytB, zX = gx$z,
    gateY = dat$gate[[y]], gateCytY = gy$gateCytB, zY = gy$z,
    stimCalled = sum(dp(dat$stim)), controlCalled = sum(dp(dat$uns)),
    netCalledPct = 100 * (mean(dp(dat$stim)) - mean(dp(dat$uns))))
}
