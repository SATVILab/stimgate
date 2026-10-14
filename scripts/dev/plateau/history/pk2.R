suppressMessages(library(dplyr)); P <- Sys.getenv("P"); options(width = 200)
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
bsf <- function(xs, xu, g) mean(xs > g) - mean(xu > g)
tags <- c(g002 = "Clean low 0.2% 100k", g002_20k = "Clean low 0.2% 20k", g0001 = "Clean low 0.01% 100k", h0001 = "Clean high 0.01% 100k",
  s0005 = "Skew high 0.05% 100k", s0005_20k = "Skew high 0.05% 20k", g002_20k_sd01 = "SD low 0.2% 20k", h002_20k_sd01 = "SD high 0.2% 20k", h002_20k_ms02 = "Shift high 0.2% 20k")
out <- list()
for (t in names(tags)) { fa <- file.path(P, paste0("dbg_am_", t, ".rds")); fk <- file.path(P, paste0("dbg_pk2_", t, ".rds")); if (!file.exists(fk)) next
  A <- readRDS(fa); K <- readRDS(fk); if (length(K) == 2 * length(A)) K <- K[seq_along(A) + length(A)]
  for (i in seq_along(A)) { s <- A[[i]]$sim[[1]]; stopifnot(identical(s$stim, K[[i]]$sim[[1]]$stim)); g0 <- A[[i]]$threshold; g1 <- K[[i]]$threshold
    if (!is.finite(g0) || !is.finite(g1)) next; truth <- mean(s$labelsStim == "gp") - mean(s$labelsUns == "gp")
    out[[length(out) + 1]] <- tibble(scenario = tags[[t]], i, lowered = g1 < g0 - 1e-9, raised = g1 > g0 + 1e-9, f0 = f1(s$stim, s$labelsStim, g0), f1 = f1(s$stim, s$labelsStim, g1),
      e0 = bsf(s$stim, s$uns, g0) / truth - 1, e1 = bsf(s$stim, s$uns, g1) / truth - 1) } }
r <- bind_rows(out)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
print(r |> group_by(scenario) |> summarise(n = n(), lowered = sum(lowered), raised = sum(raised), worse = sum(f1 < f0 - 1e-9), better = sum(f1 > f0 + 1e-9),
  F1now = median(f0), F1step = median(f1), errNow = q(e0), errStep = q(e1)) |> as.data.frame(), digits = 3)
