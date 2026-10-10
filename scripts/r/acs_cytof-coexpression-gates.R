# Analysis 18 (18-explore-acs-cytof-coexpression-gates.qmd): the co-expression
# gates (`scripts/r/coexpression-gates.R`) applied to every ACS population,
# stimulation and donor, from Analysis 9's ordinary gates. Read-only on
# Analysis 9's caches; results go to the projr cache `acs_cytof_coexpression`.

.acsCoexGatesSemantics <- "acs-coexpression-gates-v1"

# Rule settings shared by every tube (passed to `.coexLowerGates()`).
.acsCoexGatesSettings <- function() {
  list(nBin = 20L, rMin = 3.5, frac = 0.75, zMin = 2)
}

# One donor's tube pair: lowered gates for every ordered cytokine pair, the
# cells each lowered pair adds (stimulated and control), and the combination
# counts of both tubes under the ordinary and the co-expression rule.
.acsCoexGatesTube <- function(dat, pop, stim, donor, settings = .acsCoexGatesSettings()) {
  markers <- names(dat$gate)
  low <- do.call(.coexLowerGates, c(list(dat = dat, markers = markers), settings))
  added <- function(ex, i) {
    a <- low$a[[i]]
    b <- low$b[[i]]
    if (!low$lowered[[i]]) {
      return(0L)
    }
    sum(!((ex[[b]] > dat$gate[[b]]) %in% TRUE) &
      (ex[[a]] > low$condCut[[i]]) %in% TRUE & (ex[[b]] > low$cut[[i]]) %in% TRUE)
  }
  low$addedStim <- vapply(seq_len(nrow(low)), function(i) added(dat$stim, i), integer(1))
  low$addedUns <- vapply(seq_len(nrow(low)), function(i) added(dat$uns, i), integer(1))
  counts <- do.call(rbind, lapply(c("ordinary", "coexpr"), function(rule) {
    lowRule <- if (identical(rule, "coexpr")) low
    do.call(rbind, lapply(c(stim = "stim", uns = "uns"), function(tube) {
      n <- .coexCombnCounts(.coexPositive(dat[[tube]], dat$gate, lowRule, markers))
      data.frame(rule = rule, tube = tube, combn = names(n), n = as.integer(n),
        stringsAsFactors = FALSE)
    }))
  }))
  id <- data.frame(pop = pop, stim = stim, donor = donor, stringsAsFactors = FALSE)
  list(
    low = cbind(id[rep(1L, nrow(low)), ], low, row.names = NULL),
    counts = cbind(id[rep(1L, nrow(counts)), ], counts, row.names = NULL)
  )
}

# Every donor and stimulation of one population. `workers` > 1 forks
# (future multicore) over donors; the GatingSet is read in each fork.
.acsCoexGatesPopulation <- function(pop, paths, channelMap, settings = .acsCoexGatesSettings(),
                                    workers = 1L) {
  gs <- flowWorkspace::load_gs(paths$gs, backend_readonly = TRUE)
  sampleMap <- .acsCytofReadPreprocessing(paths$gs, gs)$sampleMap
  gates <- stimgate::getStimGates(paths$stimgate)
  stims <- setdiff(unique(as.character(sampleMap$stim)), "uns")
  jobs <- do.call(rbind, lapply(stims, function(s) {
    donors <- sort(unique(sampleMap$SampleID[sampleMap$stim == s]))
    donors <- donors[donors %in% sampleMap$SampleID[sampleMap$stim == "uns"]]
    data.frame(stim = rep(s, length(donors)), donor = donors, stringsAsFactors = FALSE)
  }))
  one <- function(i) {
    dat <- .acsCytposTubeData(gs, sampleMap, gates, jobs$donor[[i]], jobs$stim[[i]], channelMap)
    dat$gate <- dat$gate[is.finite(dat$gate)]
    keep <- names(dat$gate)
    dat$stim <- dat$stim[keep]
    dat$uns <- dat$uns[keep]
    .acsCoexGatesTube(dat, pop, jobs$stim[[i]], jobs$donor[[i]], settings)
  }
  out <- if (workers > 1L) {
    oplan <- future::plan(future::multicore, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    future.apply::future_lapply(seq_len(nrow(jobs)), one, future.seed = TRUE)
  } else {
    lapply(seq_len(nrow(jobs)), one)
  }
  list(
    low = do.call(rbind, lapply(out, `[[`, "low")),
    counts = do.call(rbind, lapply(out, `[[`, "counts"))
  )
}

# Background-subtracted frequencies (stimulated minus control, as fractions
# of each tube) of every combination, per tube and rule.
.acsCoexGatesFreq <- function(counts) {
  key <- c("pop", "stim", "donor", "rule", "combn")
  tot <- stats::aggregate(n ~ pop + stim + donor + rule + tube, counts, sum)
  names(tot)[names(tot) == "n"] <- "nTube"
  x <- merge(counts, tot, by = c("pop", "stim", "donor", "rule", "tube"))
  x$prop <- x$n / x$nTube
  s <- x[x$tube == "stim", c(key, "n", "prop")]
  u <- x[x$tube == "uns", c(key, "n", "prop")]
  names(s)[names(s) %in% c("n", "prop")] <- c("nStim", "propStim")
  names(u)[names(u) %in% c("n", "prop")] <- c("nUns", "propUns")
  out <- merge(s, u, by = key)
  out$propBs <- out$propStim - out$propUns
  out
}

# Number of positive markers in combination labels like "IFNg+IL2-TNF+".
.acsCoexGatesDegree <- function(combn) {
  lengths(regmatches(combn, gregexpr("\\+", combn)))
}

# Single- and multi-positive background-subtracted frequencies per tube and
# rule: overall (exactly one positive cytokine; two or more) and per cytokine
# (positive for it alone; positive for it and at least one other).
.acsCoexGatesSingleMulti <- function(freq) {
  freq <- freq[.acsCoexGatesDegree(freq$combn) > 0L, ]
  freq$degree <- .acsCoexGatesDegree(freq$combn)
  markers <- regmatches(freq$combn[[1]], gregexpr("[^+-]+(?=[+-])", freq$combn[[1]], perl = TRUE))[[1]]
  key <- c("pop", "stim", "donor", "rule")
  sumBy <- function(d, label, cyt) {
    if (!nrow(d)) {
      return(NULL)
    }
    a <- stats::aggregate(cbind(propBs, nStim, nUns) ~ pop + stim + donor + rule, d, sum)
    a$type <- label
    a$cytokine <- cyt
    a
  }
  rbind(
    sumBy(freq[freq$degree == 1L, ], "single", "any"),
    sumBy(freq[freq$degree >= 2L, ], "multi", "any"),
    do.call(rbind, lapply(markers, function(m) {
      has <- grepl(paste0(m, "+"), freq$combn, fixed = TRUE)
      rbind(
        sumBy(freq[has & freq$degree == 1L, ], "single", m),
        sumBy(freq[has & freq$degree >= 2L, ], "multi", m)
      )
    }))
  )
}

# Ordinary and co-expression values side by side, with their difference.
.acsCoexGatesPaired <- function(sm) {
  key <- c("pop", "stim", "donor", "type", "cytokine")
  o <- sm[sm$rule == "ordinary", c(key, "propBs", "nStim", "nUns")]
  cx <- sm[sm$rule == "coexpr", c(key, "propBs", "nStim", "nUns")]
  names(o)[-seq_along(key)] <- paste0(names(o)[-seq_along(key)], "Ordinary")
  names(cx)[-seq_along(key)] <- paste0(names(cx)[-seq_along(key)], "Coexpr")
  out <- merge(o, cx, by = key)
  out$diffBs <- out$propBsCoexpr - out$propBsOrdinary
  out
}

# Possible pathologies per tube: a gate lowered into the bottom bin (next to
# its floor); control cells called multi-positive at least doubled and up by
# ten or more; the net multi-positive frequency more than tripled (from at
# least 0.01%); and a net single-positive frequency driven below -0.05%.
.acsCoexGatesFlags <- function(low, paired, settings = .acsCoexGatesSettings()) {
  low$bottomBin <- low$lowered & (low$cut - low$floorB) < (low$gateB - low$floorB) / settings$nBin + 1e-9
  key <- c("pop", "stim", "donor")
  lowTube <- stats::aggregate(
    cbind(nLowered = lowered, nBottomBin = bottomBin, addedStim, addedUns) ~ pop + stim + donor,
    low, sum
  )
  multi <- paired[paired$type == "multi" & paired$cytokine == "any", ]
  single <- paired[paired$type == "single" & paired$cytokine == "any", ]
  out <- merge(lowTube, multi[c(key, "propBsOrdinary", "propBsCoexpr", "nUnsOrdinary", "nUnsCoexpr")], by = key)
  names(out)[names(out) %in% c("propBsOrdinary", "propBsCoexpr", "nUnsOrdinary", "nUnsCoexpr")] <-
    c("multiBsOrdinary", "multiBsCoexpr", "multiUnsOrdinary", "multiUnsCoexpr")
  out <- merge(out, single[c(key, "propBsOrdinary", "propBsCoexpr")], by = key)
  names(out)[names(out) %in% c("propBsOrdinary", "propBsCoexpr")] <- c("singleBsOrdinary", "singleBsCoexpr")
  out$flagBottomBin <- out$nBottomBin > 0L
  out$flagControl <- out$multiUnsCoexpr >= 2 * out$multiUnsOrdinary & out$multiUnsCoexpr - out$multiUnsOrdinary >= 10
  out$flagMultiJump <- out$multiBsCoexpr > 3 * pmax(out$multiBsOrdinary, 1e-4)
  out$flagSingleNegative <- out$singleBsCoexpr < -5e-4
  out
}

# COMPASS count matrices for one population and stimulation: rows are donors,
# columns every combination with at least one positive cytokine plus the
# all-negative one (last, as COMPASS requires), named in COMPASS form
# ("IFNg&!IL2&...").
.acsCoexGatesCompassInput <- function(counts, rule, keep = NULL) {
  d <- counts[counts$rule == rule, ]
  toCompass <- function(combn) {
    m <- regmatches(combn, gregexpr("[^+-]+[+-]", combn))
    vapply(m, function(v) {
      paste0(ifelse(endsWith(v, "+"), "", "!"), substr(v, 1L, nchar(v) - 1L), collapse = "&")
    }, character(1))
  }
  mat <- function(tube) {
    x <- d[d$tube == tube, ]
    w <- stats::xtabs(n ~ donor + combn, x)
    m <- matrix(as.numeric(w), nrow = nrow(w), dimnames = list(rownames(w), toCompass(colnames(w))))
    m
  }
  ns <- mat("stim")
  # COMPASS takes the all-negative category as the last column.
  deg <- lengths(regmatches(colnames(ns), gregexpr("(^|&)[^!]", colnames(ns))))
  ns <- ns[, order(deg == 0L, -deg), drop = FALSE]
  if (!is.null(keep)) ns <- ns[, colnames(ns) %in% keep, drop = FALSE]
  nu <- mat("uns")[rownames(ns), colnames(ns), drop = FALSE]
  list(n_s = ns, n_u = nu, meta = data.frame(donor = rownames(ns), stringsAsFactors = FALSE))
}

# Categories COMPASS is fitted on, shared by both rules: the all-negative one
# and those with at least `minCell` stimulated cells in at least `minDonor`
# donors under either rule (rare categories carry no information and slow the
# sampler).
.acsCoexGatesCompassKeep <- function(counts, minCell = 5L, minDonor = 3L) {
  keep <- unique(unlist(lapply(c("ordinary", "coexpr"), function(rule) {
    ns <- .acsCoexGatesCompassInput(counts, rule)$n_s
    colnames(ns)[colSums(ns >= minCell) >= minDonor]
  })))
  allNeg <- colnames(.acsCoexGatesCompassInput(counts, "ordinary")$n_s)
  union(keep, allNeg[length(allNeg)])
}

# SimpleCOMPASS fit for one population, stimulation and rule; returns the
# functionality and polyfunctionality scores per donor and the mean posterior
# response probability per subset.
.acsCoexGatesCompass <- function(counts, rule, keep = NULL, iterations = 10000L,
                                 replications = 8L, seed = 100L) {
  inp <- .acsCoexGatesCompassInput(counts, rule, keep)
  fit <- COMPASS::SimpleCOMPASS(
    n_s = inp$n_s, n_u = inp$n_u, meta = inp$meta, individual_id = "donor",
    iterations = iterations, replications = replications, verbose = FALSE, seed = seed
  )
  fs <- COMPASS::FunctionalityScore(fit)
  pfs <- COMPASS::PolyfunctionalityScore(fit)
  # The last column is the all-negative category, not a subset.
  gamma <- fit$fit$mean_gamma[, -ncol(fit$fit$mean_gamma), drop = FALSE]
  list(
    scores = data.frame(donor = names(fs), rule = rule, fs = as.numeric(fs),
      pfs = as.numeric(pfs[names(fs)]), stringsAsFactors = FALSE),
    subsets = data.frame(subset = colnames(gamma), rule = rule,
      meanGamma = colMeans(gamma), stringsAsFactors = FALSE, row.names = NULL)
  )
}

# Scatter of ordinary against co-expression background-subtracted
# frequencies (percent), one panel per stimulation, for single- and
# multi-positive cells.
.acsCoexGatesScatterPlot <- function(paired) {
  d <- paired[paired$cytokine == "any", ]
  d$type <- factor(d$type, levels = c("single", "multi"),
    labels = c("Single-positive (exactly one cytokine)", "Multi-positive (two or more)"))
  ggplot(d, aes(x = 100 * .data$propBsOrdinary, y = 100 * .data$propBsCoexpr)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(alpha = 0.5, size = 1) +
    facet_wrap(~ type + stim, scales = "free", nrow = 2, axes = "all",
      labeller = ggplot2::label_wrap_gen(width = 30, multi_line = FALSE)) +
    .analysis_theme() +
    labs(x = "Ordinary gates: background-subtracted frequency (%)",
      y = "Co-expression gates: background-subtracted frequency (%)")
}

# Change (co-expression minus ordinary, percentage points) in each cytokine's
# single- and multi-positive background-subtracted frequency, per donor.
.acsCoexGatesDiffPlot <- function(paired) {
  d <- paired[paired$cytokine != "any", ]
  d$type <- factor(d$type, levels = c("single", "multi"),
    labels = c("Positive for it alone", "Positive for it and another"))
  ggplot(d, aes(x = .data$cytokine, y = 100 * .data$diffBs, colour = .data$type)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2"), name = NULL) +
    facet_wrap(~ stim, scales = "free_y", axes = "all") +
    .analysis_theme() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Co-expression minus ordinary (percentage points)")
}

# COMPASS functionality and polyfunctionality scores, ordinary against
# co-expression gates, one panel per population and stimulation.
.acsCoexGatesCompassPlot <- function(scores, score = c("fs", "pfs")) {
  score <- match.arg(score)
  w <- stats::reshape(scores[c("pop", "stim", "donor", "rule", score)],
    idvar = c("pop", "stim", "donor"), timevar = "rule", direction = "wide")
  names(w) <- sub(paste0("^", score, "\\."), "", names(w))
  ggplot(w, aes(x = .data$ordinary, y = .data$coexpr)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(alpha = 0.5, size = 1) +
    facet_wrap(~ pop + stim, scales = "free", axes = "all",
      labeller = ggplot2::label_wrap_gen(multi_line = FALSE)) +
    .analysis_theme() +
    labs(
      x = paste0("Ordinary gates: ", if (score == "fs") "functionality" else "polyfunctionality", " score"),
      y = "Co-expression gates"
    )
}

# Hexagon plots of the lowered pair adding the most stimulated cells in each
# chosen tube: both tubes, the ordinary gates (orange) and the lowered gate
# (blue) from the raised cut on the other cytokine.
.acsCoexGatesTubeHexPlot <- function(tubes) {
  hex <- do.call(rbind, lapply(tubes, function(t) {
    rbind(
      data.frame(panel = paste0(t$label, "\nstimulated"), x = t$dat$stim[[t$a]], y = t$dat$stim[[t$b]]),
      data.frame(panel = paste0(t$label, "\ncontrol"), x = t$dat$uns[[t$a]], y = t$dat$uns[[t$b]])
    )
  }))
  lines <- do.call(rbind, lapply(tubes, function(t) {
    data.frame(panel = paste0(t$label, c("\nstimulated", "\ncontrol")),
      gx = t$dat$gate[[t$a]], gy = t$dat$gate[[t$b]], cut = t$cut, condCut = t$condCut)
  }))
  ggplot(hex, aes(x = .data$x, y = .data$y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_vline(data = lines, aes(xintercept = .data$gx), colour = "#D55E00", linewidth = 0.4) +
    geom_hline(data = lines, aes(yintercept = .data$gy), colour = "#D55E00", linewidth = 0.4) +
    geom_segment(data = lines, aes(x = .data$condCut, xend = Inf, y = .data$cut, yend = .data$cut),
      colour = "#0072B2", linewidth = 0.6, inherit.aes = FALSE) +
    facet_wrap(~ panel, ncol = 4, axes = "all") +
    .analysis_theme() +
    labs(x = "Cytokine the cells are positive for (a)", y = "Cytokine whose gate is lowered (b)")
}
