suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
bsf <- function(xs, xu, g) mean(xs > g) - mean(xu > g)
edge <- function(sub, g) { den <- ns$.getCpUnsLocAntimodeDensity(sub); if (is.null(den) || all(den$y == 0)) return(NA)
  r <- rle(den$y); ends <- cumsum(r$lengths); st <- ends - r$lengths + 1; j <- findInterval(g, den$x); k <- which(st <= j & ends >= j)[1]
  if (is.na(k) || st[k] <= 1) NA else sub[st[k]] }
stepDown <- function(xs, xu, g, fl, bw, k = 1, bar = 0.5, maxRej = 2) {
  ratio <- length(xs) / length(xu); nSr <- sum(xs > g); refPur <- (nSr - sum(xu > g) * ratio) / nSr
  if (!is.finite(fl) || fl >= g || !is.finite(refPur) || refPur <= 0) return(g)
  cur <- g; hi <- g; rej <- 0
  while (rej < maxRej) { lo <- max(hi - k * bw, fl); if (lo >= hi - 1e-9) break
    nS <- sum(xs > lo & xs <= cur); nU <- sum(xu > lo & xu <= cur) * ratio
    if (nS > 0 && (nS - nU) / nS >= bar * refPur) { cur <- lo; rej <- 0 } else rej <- rej + 1
    hi <- lo; if (lo <= fl) break }
  cur }
scen <- c(g002 = "Low sep 0.2% 100k", g002_20k = "Low sep 0.2% 20k", g0005_20k = "Low sep 0.05% 20k", g0001 = "Low sep 0.01% 100k",
  h002 = "High sep 0.2% 100k", h002_20k = "High sep 0.2% 20k", h0005_20k = "High sep 0.05% 20k", h0001 = "High sep 0.01% 100k", s002 = "Skew high sep 0.2% 100k", s0005 = "Skew high sep 0.05% 100k", s0001 = "Skew high sep 0.01% 100k", s002_20k = "Skew high sep 0.2% 20k", s0005_20k = "Skew high sep 0.05% 20k",
  g002_20k_ms02 = "Shift 0.2: low sep 0.2% 20k", h002_20k_ms02 = "Shift 0.2: high sep 0.2% 20k", g002_ms02 = "Shift 0.2: low sep 0.2% 100k",
  g002_20k_sd01 = "SD x1.1: low sep 0.2% 20k", h002_20k_sd01 = "SD x1.1: high sep 0.2% 20k", g002_sd01 = "SD x1.1: low sep 0.2% 100k")
out <- list()
for (t in names(scen)) { fa <- file.path(P, paste0("dbg_am_", t, ".rds")); fb <- file.path(P, paste0("dbg_pb_", t, ".rds")); if (!file.exists(fb) || !file.exists(fa)) next
  A <- readRDS(fa); B <- readRDS(fb)
  for (i in seq_along(A)) { a <- A[[i]]; b <- B[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    stopifnot(identical(xs, b$sim[[1]]$stim))
    g0 <- a$threshold; gb <- b$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- b$filter[[1]]$info$antimode$bw %||% NA
    ga <- if (is.finite(g0) && is.finite(bw) && !is.null(lb)) stepDown(xs, xu, g0, edge(sort(xs[xs >= lb]), g0), bw) else g0
    truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    out[[length(out) + 1]] <- tibble(scenario = scen[[t]], i, g0, ga, gb,
      f1_now = if (is.finite(g0)) f1(xs, l, g0) else NA, f1_after = if (is.finite(ga)) f1(xs, l, ga) else NA, f1_before = if (is.finite(gb)) f1(xs, l, gb) else NA,
      err_now = if (is.finite(g0)) bsf(xs, xu, g0) / truth - 1 else NA, err_after = if (is.finite(ga)) bsf(xs, xu, ga) / truth - 1 else NA,
      err_before = if (is.finite(gb)) bsf(xs, xu, gb) / truth - 1 else NA, capB = sub("local_fdr_", "", b$locReason)) } }
r <- bind_rows(out); saveRDS(r, file.path(P, "cmpab.rds")); options(width = 250)
cat("Per-tube gates (before = step in the antimode slot, before the cap; after = step applied to the per-tube gate)\n")
print(r |> group_by(scenario) |> summarise(n = n(), gateNow = median(g0, na.rm = TRUE), gateAfter = median(ga, na.rm = TRUE), gateBefore = median(gb, na.rm = TRUE),
  f1Now = median(f1_now, na.rm = TRUE), f1After = median(f1_after, na.rm = TRUE), f1Before = median(f1_before, na.rm = TRUE),
  f1NowMean = mean(f1_now, na.rm = TRUE), f1AfterMean = mean(f1_after, na.rm = TRUE), f1BeforeMean = mean(f1_before, na.rm = TRUE),
  capRaisedB = sum(capB == "cap_threshold_selected")) |> as.data.frame(), digits = 3)
cat("\nRelative error of raw background-subtracted frequency at the gate: median [10th, 90th]\n")
q <- function(v) sprintf("%+.0f%% [%+.0f, %+.0f]", 100 * median(v, na.rm = TRUE), 100 * quantile(v, 0.1, na.rm = TRUE), 100 * quantile(v, 0.9, na.rm = TRUE))
print(r |> group_by(scenario) |> summarise(now = q(err_now), after = q(err_after), before = q(err_before)) |> as.data.frame())
