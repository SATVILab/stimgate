suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15]))   # f1, bsf, edge, stepDown
scen <- c(g002 = "Low sep 0.2% 100k", g002_20k = "Low sep 0.2% 20k", g0005_20k = "Low sep 0.05% 20k", g0001 = "Low sep 0.01% 100k",
  h002 = "High sep 0.2% 100k", h002_20k = "High sep 0.2% 20k", h0005_20k = "High sep 0.05% 20k", h0001 = "High sep 0.01% 100k", s002 = "Skew high 0.2% 100k",
  g002_ms02 = "Shift: low sep 0.2% 100k", g002_20k_ms02 = "Shift: low sep 0.2% 20k", h002_20k_ms02 = "Shift: high sep 0.2% 20k",
  g002_sd01 = "SD x1.1: low sep 0.2% 100k", g002_20k_sd01 = "SD x1.1: low sep 0.2% 20k", h002_20k_sd01 = "SD x1.1: high sep 0.2% 20k")
out <- list()
for (t in names(scen)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    g0 <- a$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    if (!is.finite(g0) || !is.finite(bw) || is.null(lb)) next
    eS <- edge(sort(xs[xs >= lb]), g0); eU <- edge(sort(xu[xu >= lb]), g0)
    d <- if (is.finite(eS) && is.finite(eU) && eS - eU > bw / 4) eS - eU else 0
    gBase <- stepDown(xs, xu, g0, eS, bw)
    gA <- stepDown(xs, xu, g0, eS + bw / 2, bw)
    gB <- stepDown(xs, xu + d, g0, eS, bw)
    gAB <- stepDown(xs, xu + d, g0, eS + bw / 2, bw)
    truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    row <- tibble(scenario = scen[[t]], i, mism = eS - eU, d)
    for (nm in c("g0", "gBase", "gA", "gB", "gAB")) { g <- get(nm); row[[paste0("f1_", nm)]] <- f1(xs, l, g); row[[paste0("err_", nm)]] <- bsf(xs, xu, g) / truth - 1 }
    out[[length(out) + 1]] <- row } }
r <- bind_rows(out); saveRDS(r, file.path(P, "bias.rds")); options(width = 250)
cat("Median F1 (current, step, A = edge + bw/2, B = unstim shifted by tail mismatch, A+B); mismatch = stim edge - unstim edge\n")
print(r |> group_by(scenario) |> summarise(mism = median(mism, na.rm = TRUE), nBiased = sum(d > 0), now = median(f1_g0), step = median(f1_gBase), A = median(f1_gA), B = median(f1_gB), AB = median(f1_gAB)) |> as.data.frame(), digits = 3)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
cat("\nFrequency error %: median [10th, 90th]\n")
print(r |> group_by(scenario) |> summarise(now = q(err_g0), step = q(err_gBase), A = q(err_gA), B = q(err_gB), AB = q(err_gAB)) |> as.data.frame())
