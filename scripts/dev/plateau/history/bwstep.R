suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P")
r0 <- readRDS(file.path(P, "plateau.rds"))
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
map <- c("Low sep 0.2% 100k"="am_g002","Low sep 0.01% 100k"="am_g0001","Low sep 0.2% 20k"="am_g002_20k","Low sep 0.05% 20k"="am_g0005_20k",
  "High sep 0.2% 100k"="am_h002","High sep 0.01% 100k"="am_h0001","High sep 0.2% 20k"="am_h002_20k","High sep 0.05% 20k"="am_h0005_20k","Skew high sep 0.2% 100k"="am_s002")
step <- function(xs, xu, g, floor, bw, k, bar, maxRej) {
  ratio <- length(xs) / length(xu); nSr <- sum(xs > g); refPur <- (nSr - sum(xu > g) * ratio) / nSr
  cur <- g; rej <- 0; hi <- g
  while (rej < maxRej) { lo <- max(hi - k * bw, floor); if (lo >= hi - 1e-9) break
    nS <- sum(xs > lo & xs <= cur); nU <- sum(xu > lo & xu <= cur) * ratio   # span from accepted cut down to lo
    if (nS > 0 && (nS - nU) / nS >= bar * refPur) { cur <- lo; rej <- 0 } else rej <- rej + 1
    hi <- lo; if (lo <= floor) break }
  cur }
res <- list()
for (sc in unique(r0$scenario)) { d <- readRDS(file.path(P, paste0("dbg_", map[[sc]], ".rds"))); rr <- r0[r0$scenario == sc, ]
  for (j in seq_len(nrow(rr))) { e <- d[[rr$i[j]]]; s <- e$sim[[1]]; xs <- s$stim; xu <- s$uns; g <- rr$gate[j]; fl <- rr$cand[j]; if (is.na(fl)) fl <- g
    bw <- stats::bw.nrd0(xs[seq(1, length(xs), length.out = min(length(xs), 10000))])
    row <- tibble(scenario = sc, f1Gate = rr$f1Gate[j], f1Floor = f1(xs, s$labelsStim, fl))
    for (cfg in list(c(1, 0.5, 1), c(1, 0.5, 2), c(2, 0.5, 1), c(1, 0.35, 1), c(0.5, 0.5, 2))) {
      gg <- step(xs, xu, g, fl, bw, cfg[1], cfg[2], cfg[3]); nm <- sprintf("k%g_b%g_r%g", cfg[1], cfg[2], cfg[3])
      row[[nm]] <- f1(xs, s$labelsStim, gg); row[[paste0("mv_", nm)]] <- gg < g - 1e-9; row[[paste0("wr_", nm)]] <- row[[nm]] < row$f1Gate - 1e-9 }
    res[[length(res) + 1]] <- row } }
r <- bind_rows(res); options(width = 250)
print(r |> group_by(scenario) |> summarise(n = n(), gate = median(f1Gate), floor = median(f1Floor), across(starts_with("k"), median),
  across(starts_with("mv_"), sum), across(starts_with("wr_"), sum)) |> as.data.frame(), digits = 3)
