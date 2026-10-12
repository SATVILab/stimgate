# Analysis 18 (18-explore-acs-cytof-coexpression-gates.qmd): the co-expression
# gates (`scripts/r/coexpression-gates.R`) applied to every ACS population,
# stimulation and donor, from Analysis 9b's ordinary gates, over every ordered
# pair of all the cytokines Analysis 9b gates. Each tube is scored under three
# rules: the ordinary gates, the co-expression gates with the double-positive
# screen ("screened", `screen = TRUE`) and without it ("unscreened",
# `screen = FALSE`). Read-only on Analysis 9b's caches; results go to the projr
# cache `acs_cytof_coexpression`.

.acsCoexGatesSemantics <- "acs-coexpression-gates-v2"

# The three rules, ordinary gates first (the baseline).
.acsCoexGatesRules <- function() {
  c("ordinary", "screened", "unscreened")
}

# Display labels and colours of the two co-expression variants.
.acsCoexGatesRuleLabels <- c(
  ordinary = "Ordinary gates", screened = "Coexpression (screened)",
  unscreened = "Coexpression (unscreened)"
)
.acsCoexGatesRuleColours <- c(screened = "#0072B2", unscreened = "#CC79A7")

# Rule settings shared by every tube (passed to `.coexLowerGates()`); `screen`
# is set by the rule.
.acsCoexGatesSettings <- function() {
  list(nBin = 20L, rMin = 3.5, frac = 0.75, zMin = 2)
}

# One donor's tube pair: lowered gates for every ordered cytokine pair under
# each co-expression variant, the cells each lowered pair adds (stimulated and
# control), and the combination counts of both tubes under all three rules.
# Every marker in `dat$gate` is kept, so all tubes share the same
# combinations; a marker without a finite gate is negative throughout.
.acsCoexGatesTube <- function(dat, pop, stim, donor, settings = .acsCoexGatesSettings()) {
  markers <- names(dat$gate)
  screens <- c(screened = TRUE, unscreened = FALSE)
  low <- do.call(rbind, lapply(names(screens), function(rule) {
    l <- do.call(.coexLowerGates, c(list(dat = dat, markers = markers, screen = screens[[rule]]), settings))
    data.frame(rule = rule, l, stringsAsFactors = FALSE)
  }))
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
  counts <- do.call(rbind, lapply(.acsCoexGatesRules(), function(rule) {
    lowRule <- if (!identical(rule, "ordinary")) low[low$rule == rule, , drop = FALSE]
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

# Every donor and stimulation of one population, over all the cytokines in
# `channelMap` (those Analysis 9b gates). `workers` > 1 forks (future
# multicore) over donors; the GatingSet is read in each fork.
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
    .acsCoexGatesTube(dat[c("stim", "uns", "gate")], pop, jobs$stim[[i]], jobs$donor[[i]], settings)
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

# Positivity of each marker in combination labels like "IFNg+IL2-TNF+": a
# logical matrix, one row per label and one column per marker (read from the
# first label; every label must list the same markers in the same order).
.acsCoexGatesCombnMatrix <- function(combn) {
  tok <- regmatches(combn, gregexpr("[^+-]+[+-]", combn))
  markers <- sub("[+-]$", "", tok[[1]])
  out <- t(vapply(tok, function(v) {
    if (!identical(sub("[+-]$", "", v), markers)) {
      stop("Combination labels must list the same markers in the same order.")
    }
    endsWith(v, "+")
  }, logical(length(markers))))
  if (length(markers) == 1L) out <- t(out)
  dimnames(out) <- list(NULL, markers)
  out
}

# Single- and multi-positive background-subtracted frequencies per tube and
# rule: overall (exactly one positive cytokine; two or more) and per cytokine
# (positive for it alone; positive for it and at least one other).
.acsCoexGatesSingleMulti <- function(freq) {
  freq <- freq[.acsCoexGatesDegree(freq$combn) > 0L, ]
  freq$degree <- .acsCoexGatesDegree(freq$combn)
  pos <- .acsCoexGatesCombnMatrix(freq$combn)
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
    do.call(rbind, lapply(colnames(pos), function(m) {
      has <- pos[, m]
      rbind(
        sumBy(freq[has & freq$degree == 1L, ], "single", m),
        sumBy(freq[has & freq$degree >= 2L, ], "multi", m)
      )
    }))
  )
}

# Each co-expression variant beside the ordinary gates: one row per tube,
# quantity and variant (`rule`), the ordinary values (`*Ordinary`), the
# variant's (`*Coexpr`) and their difference.
.acsCoexGatesPaired <- function(sm) {
  key <- c("pop", "stim", "donor", "type", "cytokine")
  vals <- c("propBs", "nStim", "nUns")
  o <- sm[sm$rule == "ordinary", c(key, vals)]
  cx <- sm[sm$rule != "ordinary", c(key, "rule", vals)]
  names(o)[match(vals, names(o))] <- paste0(vals, "Ordinary")
  names(cx)[match(vals, names(cx))] <- paste0(vals, "Coexpr")
  out <- merge(o, cx, by = key)
  out$diffBs <- out$propBsCoexpr - out$propBsOrdinary
  out
}

# Possible pathologies per tube and co-expression variant: a gate lowered
# into the bottom bin (next to its floor); control cells called multi-positive
# at least doubled and up by ten or more; the net multi-positive frequency more
# than tripled (from at least 0.01%); and a net single-positive frequency
# driven below -0.05%.
.acsCoexGatesFlags <- function(low, paired, settings = .acsCoexGatesSettings()) {
  low$bottomBin <- low$lowered & (low$cut - low$floorB) < (low$gateB - low$floorB) / settings$nBin + 1e-9
  low$bottomBin <- low$bottomBin %in% TRUE
  key <- c("pop", "stim", "donor", "rule")
  lowTube <- stats::aggregate(
    cbind(nLowered = lowered, nBottomBin = bottomBin, addedStim, addedUns) ~ pop + stim + donor + rule,
    low, sum
  )
  multi <- paired[paired$type == "multi" & paired$cytokine == "any", ]
  single <- paired[paired$type == "single" & paired$cytokine == "any", ]
  out <- merge(lowTube, multi[c(key, "propBsOrdinary", "propBsCoexpr", "nUnsOrdinary", "nUnsCoexpr")], by = key)
  names(out)[match(c("propBsOrdinary", "propBsCoexpr", "nUnsOrdinary", "nUnsCoexpr"), names(out))] <-
    c("multiBsOrdinary", "multiBsCoexpr", "multiUnsOrdinary", "multiUnsCoexpr")
  out <- merge(out, single[c(key, "propBsOrdinary", "propBsCoexpr")], by = key)
  names(out)[match(c("propBsOrdinary", "propBsCoexpr"), names(out))] <- c("singleBsOrdinary", "singleBsCoexpr")
  out$flagBottomBin <- out$nBottomBin > 0L
  out$flagControl <- out$multiUnsCoexpr >= 2 * out$multiUnsOrdinary & out$multiUnsCoexpr - out$multiUnsOrdinary >= 10
  out$flagMultiJump <- out$multiBsCoexpr > 3 * pmax(out$multiBsOrdinary, 1e-4)
  out$flagSingleNegative <- out$singleBsCoexpr < -5e-4
  out
}

# Response tests ------------------------------------------------------------

# Counts per donor, rule and combination tested for a response, with both
# tubes side by side: every exclusive combination with at least one positive
# marker (`type = "exclusive"`, labelled like "IFNg+IL2-TNF+") and each
# marker's marginal (`type = "marginal"`, labelled like "IFNg+": positive for
# it, whatever the others). `NStim` and `NUns` are the tubes' cell counts.
.acsCoexGatesTestCounts <- function(counts) {
  key <- c("pop", "stim", "donor", "rule", "tube")
  tot <- stats::aggregate(n ~ pop + stim + donor + rule + tube, counts, sum)
  names(tot)[names(tot) == "n"] <- "nTube"
  excl <- counts[.acsCoexGatesDegree(counts$combn) > 0L, c(key, "combn", "n")]
  excl$type <- rep("exclusive", nrow(excl))
  pos <- .acsCoexGatesCombnMatrix(counts$combn)
  marg <- do.call(rbind, lapply(colnames(pos), function(m) {
    d <- counts[pos[, m], c(key, "n"), drop = FALSE]
    a <- stats::aggregate(n ~ pop + stim + donor + rule + tube, d, sum)
    a$combn <- rep(paste0(m, "+"), nrow(a))
    a$type <- rep("marginal", nrow(a))
    a[c(key, "combn", "n", "type")]
  }))
  x <- merge(rbind(excl, marg), tot, by = key)
  idKey <- c("pop", "stim", "donor", "rule", "type", "combn")
  s <- x[x$tube == "stim", c(idKey, "n", "nTube")]
  u <- x[x$tube == "uns", c(idKey, "n", "nTube")]
  names(s)[match(c("n", "nTube"), names(s))] <- c("nStim", "NStim")
  names(u)[match(c("n", "nTube"), names(u))] <- c("nUns", "NUns")
  out <- merge(s, u, by = idKey)
  out$propStim <- out$nStim / out$NStim
  out$propUns <- out$nUns / out$NUns
  out$propBs <- out$propStim - out$propUns
  out[order(out$pop, out$stim, out$rule, out$type, out$combn, out$donor), , drop = FALSE]
}

# One-sided Fisher's exact test that the stimulated tube has a larger share of
# cells in a combination than the control tube (`nS` of `NS` against `nU` of
# `NU`); vectorised, equal to `fisher.test(alternative = "greater")` on the
# 2x2 table (rows stimulated and control, columns in and out of the
# combination).
.acsCoexGatesFisherGreater <- function(nS, NS, nU, NU) {
  stats::phyper(nS - 1, nS + nU, NS + NU - nS - nU, NS, lower.tail = FALSE)
}

# One-sided paired Wilcoxon signed-rank test that the stimulated frequency
# exceeds the control frequency, from the per-donor differences `d`. Zero
# differences are dropped, as `wilcox.test()` does; with no non-zero
# difference the test is undefined (NA). Ties use the normal approximation.
.acsCoexGatesWilcoxGreater <- function(d) {
  d <- d[is.finite(d)]
  if (!any(d != 0)) {
    return(NA_real_)
  }
  suppressWarnings(stats::wilcox.test(d, mu = 0, alternative = "greater")$p.value)
}

# One-sided paired t-test of the same differences; NA with fewer than two
# donors or (essentially) constant differences, where `t.test()` is undefined.
.acsCoexGatesTGreater <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) < 2L) {
    return(NA_real_)
  }
  se <- stats::sd(d) / sqrt(length(d))
  if (!(se > 10 * .Machine$double.eps * abs(mean(d)))) {
    return(NA_real_)
  }
  stats::t.test(d, alternative = "greater")$p.value
}

# Response tests per population, stimulation, rule and combination, from
# `.acsCoexGatesTestCounts()`:
# - per donor, a one-sided Fisher's exact test of stimulated against control
#   counts at `alpha` (unadjusted; `nDonorFisher` donors significant);
# - across donors, one-sided paired Wilcoxon signed-rank and t-tests of
#   stimulated minus control frequency, with Benjamini-Hochberg adjustment
#   across the combinations (exclusive and marginal together) within each
#   population, stimulation and rule. Undefined tests (NA) are left out of the
#   adjustment and never called responding.
# A combination responds when its adjusted Wilcoxon p-value is below `alpha`.
.acsCoexGatesResponseTests <- function(tc, alpha = 0.05) {
  tc$pFisher <- .acsCoexGatesFisherGreater(tc$nStim, tc$NStim, tc$nUns, tc$NUns)
  idKey <- c("pop", "stim", "rule", "type", "combn")
  grp <- interaction(tc[idKey], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(tc, grp), function(d) {
    data.frame(
      d[1L, idKey], nDonor = nrow(d),
      nDonorFisher = sum(d$pFisher < alpha),
      nStim = sum(d$nStim), nUns = sum(d$nUns),
      meanPropBs = mean(d$propBs), medianPropBs = stats::median(d$propBs),
      nDonorNonZero = sum(d$propBs != 0),
      pWilcox = .acsCoexGatesWilcoxGreater(d$propBs),
      pT = .acsCoexGatesTGreater(d$propBs),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }))
  fam <- interaction(out[c("pop", "stim", "rule")], drop = TRUE)
  bh <- function(p) stats::p.adjust(p, method = "BH")
  out$pWilcoxAdj <- stats::ave(out$pWilcox, fam, FUN = bh)
  out$pTAdj <- stats::ave(out$pT, fam, FUN = bh)
  out$responding <- (out$pWilcoxAdj < alpha) %in% TRUE
  out$respondingT <- (out$pTAdj < alpha) %in% TRUE
  rownames(out) <- NULL
  out
}

# The three rules side by side for each population, stimulation and
# combination: whether each calls it responding (Wilcoxon, adjusted), its
# adjusted p-values, whether the call changes between rules and whether only
# the unscreened rule calls it responding.
.acsCoexGatesResponseCompare <- function(tests) {
  key <- c("pop", "stim", "type", "combn")
  rules <- .acsCoexGatesRules()
  out <- unique(tests[key])
  for (r in rules) {
    d <- tests[tests$rule == r, c(key, "responding", "pWilcoxAdj", "nDonorFisher")]
    names(d)[-seq_along(key)] <- paste0(names(d)[-seq_along(key)], "_", r)
    out <- merge(out, d, by = key, all.x = TRUE)
  }
  resp <- as.matrix(out[paste0("responding_", rules)])
  resp[is.na(resp)] <- FALSE
  out$nRulesResponding <- rowSums(resp)
  out$changed <- out$nRulesResponding > 0L & out$nRulesResponding < length(rules)
  out$unscreenedOnly <- resp[, "responding_unscreened"] &
    !resp[, "responding_screened"] & !resp[, "responding_ordinary"]
  out[order(out$pop, out$stim, out$type, out$combn), , drop = FALSE]
}

# Number of combinations each rule tests (a defined Wilcoxon p-value) and
# calls responding, per population, stimulation, rule and type of
# combination, and how many of those only this rule calls responding.
.acsCoexGatesResponseSummary <- function(tests, compare) {
  tests$one <- 1L
  tests$tested <- !is.na(tests$pWilcox)
  key <- c("pop", "stim", "rule", "type")
  out <- stats::aggregate(cbind(nCombn = one, nTested = tested, nResponding = responding,
    nRespondingT = respondingT) ~ pop + stim + rule + type, tests, sum)
  rules <- .acsCoexGatesRules()
  only <- do.call(rbind, lapply(rules, function(r) {
    hit <- compare[[paste0("responding_", r)]] %in% TRUE & compare$nRulesResponding == 1L
    if (!any(hit)) {
      return(NULL)
    }
    a <- stats::aggregate(list(nOnlyThisRule = hit[hit]), compare[hit, c("pop", "stim", "type")], sum)
    a$rule <- r
    a
  }))
  if (is.null(only)) {
    out$nOnlyThisRule <- 0L
  } else {
    out <- merge(out, only, by = key, all.x = TRUE)
    out$nOnlyThisRule[is.na(out$nOnlyThisRule)] <- 0L
  }
  out[order(out$pop, out$stim, out$type, match(out$rule, rules)), , drop = FALSE]
}

# Control-tube contamination: marker-positive calls each co-expression
# variant adds over the ordinary gates (summed over markers from the marginal
# counts, so a cell added for two markers counts twice), in the stimulated and
# the control tube, per population, stimulation and variant, per million cells
# of each tube. `controlShare` is the control rate over the stimulated rate:
# roughly the share of the stimulated tube's added calls the background
# explains (NA when nothing was added in the stimulated tube).
.acsCoexGatesContamination <- function(tc) {
  m <- tc[tc$type == "marginal", ]
  o <- m[m$rule == "ordinary", c("pop", "stim", "donor", "combn", "nStim", "nUns")]
  names(o)[5:6] <- c("nStimOrdinary", "nUnsOrdinary")
  x <- merge(m[m$rule != "ordinary", ], o, by = c("pop", "stim", "donor", "combn"))
  x$addedStim <- x$nStim - x$nStimOrdinary
  x$addedUns <- x$nUns - x$nUnsOrdinary
  added <- stats::aggregate(cbind(addedStim, addedUns) ~ pop + stim + rule, x, sum)
  cells <- unique(x[c("pop", "stim", "rule", "donor", "NStim", "NUns")])
  cells$one <- 1L
  cells <- stats::aggregate(cbind(nDonor = one, cellsStim = NStim, cellsUns = NUns) ~
    pop + stim + rule, cells, sum)
  out <- merge(cells, added, by = c("pop", "stim", "rule"))
  out$addedStimPerMillion <- 1e6 * out$addedStim / out$cellsStim
  out$addedUnsPerMillion <- 1e6 * out$addedUns / out$cellsUns
  out$controlShare <- ifelse(out$addedStim > 0, out$addedUnsPerMillion / out$addedStimPerMillion, NA_real_)
  out
}

# Short combination labels: the positive markers only ("IFNg+TNF+"), and
# marginals as "IFNg+ (any)".
.acsCoexGatesCombnLabel <- function(combn, type) {
  short <- vapply(regmatches(combn, gregexpr("[^+-]+\\+", combn)), paste0, character(1), collapse = "")
  ifelse(type == "marginal", paste0(combn, " (any)"), short)
}

# Heatmap of -log10 adjusted Wilcoxon p-values for one population: rows are
# the combinations any rule calls responding in any stimulation, plus every
# marginal; columns the rules; one panel per stimulation. Responding cells
# are marked with a dot; untestable ones are grey.
.acsCoexGatesResponseHeatmap <- function(tests) {
  hit <- unique(tests$combn[tests$responding & tests$type == "exclusive"])
  d <- tests[tests$type == "marginal" | tests$combn %in% hit, ]
  d$label <- .acsCoexGatesCombnLabel(d$combn, d$type)
  lev <- unique(d$label[order(d$type == "marginal", .acsCoexGatesDegree(d$combn), d$combn)])
  d$label <- factor(d$label, levels = rev(lev))
  d$rule <- factor(d$rule, levels = .acsCoexGatesRules(),
    labels = .acsCoexGatesRuleLabels[.acsCoexGatesRules()])
  d$score <- -log10(d$pWilcoxAdj)
  ggplot(d, aes(x = .data$rule, y = .data$label, fill = .data$score)) +
    geom_tile(colour = "white") +
    geom_point(data = d[d$responding, , drop = FALSE], size = 1, colour = "white") +
    scale_fill_viridis_c(name = "-log10 adjusted p", na.value = "grey85") +
    facet_wrap(~ stim, nrow = 1, axes = "all") +
    .analysis_theme(grid = "none") +
    labs(x = NULL, y = NULL)
}

# COMPASS -------------------------------------------------------------------

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

# Categories COMPASS is fitted on, shared by every rule: the all-negative one
# and those with at least `minCell` stimulated cells in at least `minDonor`
# donors under any rule (rare categories carry no information and slow the
# sampler).
.acsCoexGatesCompassKeep <- function(counts, minCell = 5L, minDonor = 3L) {
  keep <- unique(unlist(lapply(unique(counts$rule), function(rule) {
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

# Plots ---------------------------------------------------------------------

.acsCoexGatesRuleScale <- function() {
  scale_colour_manual(values = .acsCoexGatesRuleColours,
    labels = .acsCoexGatesRuleLabels[names(.acsCoexGatesRuleColours)], name = NULL)
}

# Scatter of ordinary against co-expression background-subtracted
# frequencies (percent), one panel per stimulation, for single- and
# multi-positive cells; colour marks the co-expression variant.
.acsCoexGatesScatterPlot <- function(paired) {
  d <- paired[paired$cytokine == "any", ]
  d$type <- factor(d$type, levels = c("single", "multi"),
    labels = c("Single-positive (exactly one cytokine)", "Multi-positive (two or more)"))
  ggplot(d, aes(x = 100 * .data$propBsOrdinary, y = 100 * .data$propBsCoexpr, colour = .data$rule)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(alpha = 0.5, size = 1) +
    .acsCoexGatesRuleScale() +
    facet_wrap(~ type + stim, scales = "free", nrow = 2, axes = "all",
      labeller = ggplot2::label_wrap_gen(width = 30, multi_line = FALSE)) +
    .analysis_theme() +
    labs(x = "Ordinary gates: background-subtracted frequency (%)",
      y = "Coexpression gates: background-subtracted frequency (%)")
}

# Change (co-expression minus ordinary, percentage points) in each cytokine's
# single- and multi-positive background-subtracted frequency, per donor, one
# panel per variant and stimulation.
.acsCoexGatesDiffPlot <- function(paired) {
  d <- paired[paired$cytokine != "any", ]
  d$type <- factor(d$type, levels = c("single", "multi"),
    labels = c("Positive for it alone", "Positive for it and another"))
  d$ruleLabel <- unname(.acsCoexGatesRuleLabels[d$rule])
  ggplot(d, aes(x = .data$cytokine, y = 100 * .data$diffBs, colour = .data$type)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2"), name = NULL) +
    facet_wrap(~ ruleLabel + stim, scales = "free_y", axes = "all",
      labeller = ggplot2::label_wrap_gen(width = 30, multi_line = FALSE)) +
    .analysis_theme() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Coexpression minus ordinary (percentage points)")
}

# COMPASS functionality and polyfunctionality scores, ordinary against each
# co-expression variant, one panel per population and stimulation.
.acsCoexGatesCompassPlot <- function(scores, score = c("fs", "pfs")) {
  score <- match.arg(score)
  key <- c("pop", "stim", "donor")
  o <- scores[scores$rule == "ordinary", c(key, score)]
  names(o)[names(o) == score] <- "ordinary"
  cx <- scores[scores$rule != "ordinary", c(key, "rule", score)]
  names(cx)[names(cx) == score] <- "coexpr"
  w <- merge(o, cx, by = key)
  ggplot(w, aes(x = .data$ordinary, y = .data$coexpr, colour = .data$rule)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(alpha = 0.5, size = 1) +
    .acsCoexGatesRuleScale() +
    facet_wrap(~ pop + stim, scales = "free", axes = "all",
      labeller = ggplot2::label_wrap_gen(multi_line = FALSE)) +
    .analysis_theme() +
    labs(
      x = paste0("Ordinary gates: ", if (score == "fs") "functionality" else "polyfunctionality", " score"),
      y = "Coexpression gates"
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
