suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15])); source(file.path(P, "scen_block.R"))
# Location-scale map of unstim, x' = m + b (x - m) + a with m the unstim median, chosen so the unstim share in two
# bands below the plateau edge matches the stim share (least squares on log shares, grid search).
lsMap <- function(xs, xu, bands) { m <- stats::median(xu); su <- sort(xu); nU <- length(xu)
  pS <- sapply(bands, function(bd) mean(xs >= bd[1] & xs < bd[2])); if (any(pS == 0)) return(c(a = 0, b = 1, m = m))
  best <- c(a = 0, b = 1, m = m); bestV <- Inf
  for (b in seq(0.8, 1.6, by = 0.01)) for (a in seq(-0.6, 0.6, by = 0.01)) {
    pU <- sapply(bands, function(bd) { lo <- (bd[1] - a - m) / b + m; hi <- (bd[2] - a - m) / b + m
      (findInterval(hi, su, left.open = TRUE) - findInterval(lo, su, left.open = TRUE)) / nU })
    if (any(pU == 0)) next
    v <- sum(log(pU / pS)^2); if (v < bestV) { bestV <- v; best <- c(a = a, b = b, m = m) } }
  best }
mapU <- function(xu, mp) mp[["m"]] + mp[["b"]] * (xu - mp[["m"]]) + mp[["a"]]
out <- list()
for (t in names(scenAll)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    g0 <- a$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    if (!is.finite(g0) || !is.finite(bw) || is.null(lb)) next
    fl <- edge(sort(xs[xs >= lb]), g0); if (!is.finite(fl)) fl <- g0
    truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    m1 <- lsMap(xs, xu, list(c(fl - 8 * bw, fl - 5 * bw), c(fl - 5 * bw, fl - 2 * bw)))
    m2 <- lsMap(xs, xu, list(c(fl - 12 * bw, fl - 8 * bw), c(fl - 8 * bw, fl - 4 * bw)))
    gS <- stepDown(xs, xu, g0, fl, bw)
    gL1 <- stepDown(xs, mapU(xu, m1), g0, fl, bw)
    gL2 <- stepDown(xs, mapU(xu, m2), g0, fl, bw)
    row <- tibble(scenario = scenAll[[t]], i, a1 = m1[["a"]], b1 = m1[["b"]], b2 = m2[["b"]])
    for (nm in c("g0", "gS", "gL1", "gL2")) { g <- get(nm); row[[paste0("f1_", nm)]] <- f1(xs, l, g); row[[paste0("e_", nm)]] <- bsf(xs, xu, g) / truth - 1 }
    out[[length(out) + 1]] <- row } }
r <- bind_rows(out); saveRDS(r, file.path(P, "ls.rds")); options(width = 250)
cat("Median F1. L1 = unstim widened/shifted about its median, bands 8-5 and 5-2 bw below the plateau edge; L2 = bands 12-8 and 8-4 bw below\n")
print(r |> group_by(scenario) |> summarise(a1med = median(a1), b1q10 = quantile(b1, .1), b1med = median(b1), b1q90 = quantile(b1, .9), b2med = median(b2),
  now = median(f1_g0), step = median(f1_gS), L1 = median(f1_gL1), L2 = median(f1_gL2)) |> as.data.frame(), digits = 3)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
cat("\nFrequency error %, median [10th, 90th]\n")
print(r |> group_by(scenario) |> summarise(now = q(e_g0), step = q(e_gS), L1 = q(e_gL1), L2 = q(e_gL2)) |> as.data.frame())
