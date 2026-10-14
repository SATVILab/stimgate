suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); options(width = 230)
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:3]))   # f1, bsf
# Step with a choice of guard: z = 0 plain, z > 0 cautious bound; minS = minimum stim cells in an accepted span.
stepG <- function(xs, xu, g, edge, bw, z = 0, minS = 1, gapBw = 1) {
  if (!isTRUE(is.finite(edge) && edge < g)) return(g)
  ratio <- length(xs) / length(xu); nR <- sum(xs > g); ref <- (nR - sum(xu > g) * ratio) / nR; if (!is.finite(ref) || ref <= 0) return(g)
  cut <- g; hi <- g; rej <- 0; empty <- 0
  while (rej < 2) { lo <- max(hi - bw, edge); if (lo >= hi - 1e-12) break
    if (!any(xs > lo & xs <= hi)) { empty <- empty + (hi - lo); if (empty >= gapBw * bw - 1e-12) break } else empty <- 0
    nS <- sum(xs > lo & xs <= cut); nU <- sum(xu > lo & xu <= cut) * ratio
    pur <- if (nS >= minS) (nS - nU) / nS - z * sqrt(nS + nU * ratio) / nS else -Inf
    if (pur >= 0.5 * ref) { cut <- lo; rej <- 0 } else rej <- rej + 1
    hi <- lo; if (lo <= edge) break }
  cut }
# unstim shift that brings the band ratio (4-8 bw below the edge) down to target
shiftTo <- function(xs, xu, lo, hi, target) { fS <- mean(xs >= lo & xs < hi); if (fS == 0) return(0)
  for (d in seq(0, 1.5, by = 0.005)) if (fS / mean(xu + d >= lo & xu + d < hi) <= target) return(d); 1.5 }
scen <- c(g002 = "Clean low 0.2% 100k", g002_20k = "Clean low 0.2% 20k", g0005_20k = "Clean low 0.05% 20k", g0001 = "Clean low 0.01% 100k",
  h002 = "Clean high 0.2% 100k", h002_20k = "Clean high 0.2% 20k", h0005_20k = "Clean high 0.05% 20k", h0001 = "Clean high 0.01% 100k",
  s002 = "Skew high 0.2% 100k", s0005 = "Skew high 0.05% 100k", s002_20k = "Skew high 0.2% 20k", s0005_20k = "Skew high 0.05% 20k",
  g002_sd01 = "SD low 0.2% 100k", g002_20k_sd01 = "SD low 0.2% 20k", h002_20k_sd01 = "SD high 0.2% 20k",
  g002_ms02 = "Shift low 0.2% 100k", g002_20k_ms02 = "Shift low 0.2% 20k", h002_20k_ms02 = "Shift high 0.2% 20k")
out <- list()
for (t in names(scen)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim; g <- a$threshold
    bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA; lb <- a$filter[[1]]$info$minProbXPos
    if (!is.finite(g) || !is.finite(bw) || is.null(lb)) next
    ed <- .getCpUnsLocPlateauEdge(xs, g, lb); truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    # band ratio below the edge; failing tubes get a target shrinking with the mismatch
    r0 <- if (is.finite(ed)) mean(xs >= ed - 8 * bw & xs < ed - 4 * bw) / mean(xu >= ed - 8 * bw & xu < ed - 4 * bw) else NA
    fail <- grepl("SD|Shift", scen[[t]])
    tgtGrad <- if (fail && is.finite(r0)) max(0.7, 0.96 / r0) else 0.96
    dGrad <- if (is.finite(ed)) shiftTo(xs, xu, ed - 8 * bw, ed - 4 * bw, tgtGrad) else 0
    d70 <- if (is.finite(ed)) shiftTo(xs, xu, ed - 8 * bw, ed - 4 * bw, 0.70) else 0
    gates <- c(now = g,
      gap = if (fail) g else stepG(xs, xu, g, ed, bw),
      gapMin2 = if (fail) g else stepG(xs, xu, g, ed, bw, minS = 2),
      gapMin3 = if (fail) g else stepG(xs, xu, g, ed, bw, minS = 3),
      min2Graded = stepG(xs, xu + dGrad, g, ed, bw, minS = 2),
      min2Fail70 = if (fail) stepG(xs, xu + d70, g, ed, bw, minS = 2) else stepG(xs, xu, g, ed, bw, minS = 2))
    row <- tibble(scenario = scen[[t]], i)
    for (nm in names(gates)) { row[[paste0("f1_", nm)]] <- f1(xs, l, gates[[nm]]); row[[paste0("e_", nm)]] <- bsf(xs, xu, gates[[nm]]) / truth - 1 }
    out[[length(out) + 1]] <- row } }
r <- bind_rows(out); saveRDS(r, file.path(P, "guards.rds"))
v <- c("now", "gap", "gapMin2", "gapMin3", "min2Graded", "min2Fail70")
cat("Median F1 (batch scenarios: no step unless shifted)\n")
print(r |> group_by(scenario) |> summarise(across(all_of(paste0("f1_", v)), median)) |> rename_with(~ sub("f1_", "", .x)) |> as.data.frame(), digits = 3)
cat("\nTubes worse than now / better than now\n")
print(r |> group_by(scenario) |> summarise(across(all_of(paste0("f1_", v[-1])), ~ paste0(sum(.x < f1_now - 1e-9), "/", sum(.x > f1_now + 1e-9)))) |> rename_with(~ sub("f1_", "", .x)) |> as.data.frame())
cat("\n90th percentile of frequency error (%)\n")
print(r |> group_by(scenario) |> summarise(across(all_of(paste0("e_", v)), ~ round(100 * quantile(.x, 0.9)))) |> rename_with(~ sub("e_", "", .x)) |> as.data.frame())
