suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15]))
source(file.path(P, "scen_block.R"))
# smallest shift of unstim that brings the stim/unstim share in [lo, hi) to <= 1
findShift <- function(xs, xu, lo, hi) { fS <- mean(xs >= lo & xs < hi); if (fS == 0) return(0)
  for (d in seq(0, 1.5, by = 0.005)) if (mean(xu + d >= lo & xu + d < hi) >= fS) return(d); 1.5 }
out <- list()
for (t in names(scenAll)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    g0 <- a$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    if (!is.finite(g0) || !is.finite(bw) || is.null(lb)) next
    fl <- edge(sort(xs[xs >= lb]), g0); truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    dFix <- findShift(xs, xu, g0 - 2, g0 - 1.5)
    dEdge <- if (is.finite(fl)) findShift(xs, xu, fl - 8 * bw, fl - 4 * bw) else 0
    gS <- stepDown(xs, xu, g0, fl, bw)
    g1 <- stepDown(xs, xu + dEdge, g0, fl, bw)                              # all unstim shifted, edge-relative band
    g2 <- stepDown(xs, ifelse(xu < g0, xu + dEdge, xu), g0, fl, bw)         # only unstim below the gate shifted
    g3 <- stepDown(xs, xu + dFix, g0, fl, bw)                                # all shifted, fixed band
    row <- tibble(scenario = scenAll[[t]], i, dFix, dEdge)
    for (nm in c("g0", "gS", "g1", "g2", "g3")) { g <- get(nm); row[[paste0("f1_", nm)]] <- f1(xs, l, g); row[[paste0("e_", nm)]] <- bsf(xs, xu, g) / truth - 1 }
    out[[length(out) + 1]] <- row } }
r <- bind_rows(out); saveRDS(r, file.path(P, "shiftu.rds")); options(width = 250)
cat("Median F1. now = current gate; step; S1 = step with all unstim shifted (band 4-8 bw below plateau edge); S2 = only unstim below gate shifted; S3 = all shifted, band 1.5-2 below gate\n")
print(r |> group_by(scenario) |> summarise(dEdge = median(dEdge), dFix = median(dFix), now = median(f1_g0), step = median(f1_gS), S1 = median(f1_g1), S2 = median(f1_g2), S3 = median(f1_g3)) |> as.data.frame(), digits = 3)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
cat("\nFrequency error %, median [10th, 90th]\n")
print(r |> group_by(scenario) |> summarise(now = q(e_g0), step = q(e_gS), S1 = q(e_g1), S2 = q(e_g2), S3 = q(e_g3)) |> as.data.frame())
