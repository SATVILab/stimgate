suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
scen <- c(am_g002 = "Low sep 0.2% 100k", am_g0001 = "Low sep 0.01% 100k", am_g002_20k = "Low sep 0.2% 20k", am_g0005_20k = "Low sep 0.05% 20k",
  am_h002 = "High sep 0.2% 100k", am_h0001 = "High sep 0.01% 100k", am_h002_20k = "High sep 0.2% 20k", am_h0005_20k = "High sep 0.05% 20k",
  am_s002 = "Skew high sep 0.2% 100k")
out <- list()
for (t in names(scen)) { d <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  for (i in seq_along(d)) { e <- d[[i]]; s <- e$sim[[1]]; xs <- s$stim; xu <- s$uns; g <- e$threshold; lb <- e$filter[[1]]$info$minProbXPos
    if (!is.finite(g) || is.null(lb)) next
    sub <- sort(xs[xs >= lb]); den <- ns$.getCpUnsLocAntimodeDensity(sub); if (is.null(den)) next
    r <- rle(den$y); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    k <- which(starts <= findInterval(g, den$x) & ends >= findInterval(g, den$x))[1]; if (is.na(k)) k <- length(starts)
    cand <- if (starts[k] > 1) sub[starts[k]] else NA   # data point at left edge of the plateau holding the gate
    if (!is.finite(cand) || cand >= g) { cand <- NA }
    nS <- if (is.na(cand)) 0 else sum(xs > cand & xs <= g); nU <- if (is.na(cand)) 0 else sum(xu > cand & xu <= g) * length(xs) / length(xu)
    nSref <- sum(xs > g); nUref <- sum(xu > g) * length(xs) / length(xu)
    out[[length(out) + 1]] <- tibble(scenario = scen[[t]], i, gate = g, cand, nS, nU, pur = (nS - nU) / nS, purLo = (nS - nU - 1.28 * sqrt(nS + nU)) / nS,
      refPur = (nSref - nUref) / nSref, f1Gate = f1(xs, s$labelsStim, g), f1Cand = if (is.na(cand)) NA else f1(xs, s$labelsStim, cand)) } }
r <- bind_rows(out); saveRDS(r, file.path(P, "plateau.rds"))
rule <- function(r, bar, lo) { ok <- !is.na(r$cand) & r$nS > 0 & (if (lo) r$purLo else r$pur) >= bar * r$refPur; ifelse(ok, r$f1Cand, r$f1Gate) }
options(width = 220)
print(r |> group_by(scenario) |> summarise(n = n(), haveCand = sum(!is.na(cand)), gate = median(gate), cand = median(cand, na.rm = TRUE), nS = median(nS), nU = median(nU),
  pur = median(pur, na.rm = TRUE), f1Gate = median(f1Gate), f1CandAll = median(ifelse(is.na(cand), f1Gate, f1Cand)),
  candBetter = sum(f1Cand > f1Gate, na.rm = TRUE), candWorse = sum(f1Cand < f1Gate, na.rm = TRUE)) |> as.data.frame(), digits = 3)
cat("\nMedian F1 under 'move to plateau edge if span purity >= bar x purity above gate'\n")
tb <- r |> group_by(scenario) |> summarise(gate = median(f1Gate),
  pt0.5 = median(rule(pick(everything()), 0.5, FALSE)), pt0.7 = median(rule(pick(everything()), 0.7, FALSE)),
  lo0.35 = median(rule(pick(everything()), 0.35, TRUE)), lo0.5 = median(rule(pick(everything()), 0.5, TRUE)),
  moved_lo0.5 = sum(!is.na(cand) & purLo >= 0.5 * refPur), worse_lo0.5 = sum(!is.na(cand) & purLo >= 0.5 * refPur & f1Cand < f1Gate),
  moved_pt0.5 = sum(!is.na(cand) & pur >= 0.5 * refPur), worse_pt0.5 = sum(!is.na(cand) & pur >= 0.5 * refPur & f1Cand < f1Gate))
print(as.data.frame(tb), digits = 3)
