suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); options(width = 220)
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
bsf <- function(xs, xu, g) mean(xs > g) - mean(xu > g)
spanPure <- function(xs, xu, lo, hi, g) { r <- length(xs) / length(xu); nR <- sum(xs > g); ref <- (nR - sum(xu > g) * r) / nR
  nS <- sum(xs > lo & xs <= hi); nU <- sum(xu > lo & xu <= hi) * r; nS >= 2 && is.finite(ref) && ref > 0 && (nS - nU) / nS >= 0.5 * ref }
tags <- c(g002 = "Clean low 0.2% 100k", g002_20k = "Clean low 0.2% 20k", g0005_20k = "Clean low 0.05% 20k", g0001 = "Clean low 0.01% 100k",
  h002 = "Clean high 0.2% 100k", h002_20k = "Clean high 0.2% 20k", h0005_20k = "Clean high 0.05% 20k", h0001 = "Clean high 0.01% 100k",
  s002 = "Skew high 0.2% 100k", s0005 = "Skew high 0.05% 100k", s002_20k = "Skew high 0.2% 20k", s0005_20k = "Skew high 0.05% 20k")
out <- list()
for (t in names(tags)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { s <- A[[i]]$sim[[1]]; xs <- s$stim; l <- s$labelsStim; g <- A[[i]]$threshold
    bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA; lb <- A[[i]]$filter[[1]]$info$minProbXPos; if (!is.finite(g) || !is.finite(bw) || is.null(lb)) next
    xu <- s$uns + bw   # biasUns = one bandwidth, as in the analyses
    ed <- .getCpUnsLocPlateauEdge(xs, g, lb); truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    gStep <- .getCpUnsLocPlateauStepGate(xs, xu, g, ed, bw, 0.5, minCells = 2, gapBw = 1)$gate
    nBetween <- if (is.finite(ed)) sum(xs > ed & xs <= g) else 0
    half <- if (is.finite(ed)) (g + ed) / 2 else g
    gH <- if (!is.finite(ed)) g else if (nBetween < 20) half else gStep
    gH2 <- if (!is.finite(ed)) g else if (nBetween < 20) half else if (spanPure(xs, xu, ed, g, g)) ed else g
    row <- tibble(scenario = tags[[t]], i, nBetween)
    for (nm in c("g", "gStep", "gH", "gH2")) { gg <- get(nm); row[[paste0("f1_", nm)]] <- f1(xs, l, gg); row[[paste0("e_", nm)]] <- bsf(xs, s$uns, gg) / truth - 1 }
    out[[length(out) + 1]] <- row } }
r <- bind_rows(out)
cat("Median F1; cells between edge and gate (median); tubes worse/better than now\n")
print(r |> group_by(scenario) |> summarise(nBetween = median(nBetween), now = median(f1_g), step = median(f1_gStep), H = median(f1_gH), H2 = median(f1_gH2),
  stepWB = paste0(sum(f1_gStep < f1_g - 1e-9), "/", sum(f1_gStep > f1_g + 1e-9)), H_WB = paste0(sum(f1_gH < f1_g - 1e-9), "/", sum(f1_gH > f1_g + 1e-9)),
  H2_WB = paste0(sum(f1_gH2 < f1_g - 1e-9), "/", sum(f1_gH2 > f1_g + 1e-9)),
  q90now = round(100 * quantile(e_g, .9)), q90step = round(100 * quantile(e_gStep, .9)), q90H = round(100 * quantile(e_gH, .9)), q90H2 = round(100 * quantile(e_gH2, .9))) |> as.data.frame(), digits = 3)
