suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15])); source(file.path(P, "scen_block.R"))
out <- list()
for (t in names(scenAll)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim
    g0 <- a$threshold; lb <- a$filter[[1]]$info$minProbXPos; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    if (!is.finite(g0) || !is.finite(bw) || is.null(lb)) next
    fl <- edge(sort(xs[xs >= lb]), g0); if (!is.finite(fl)) fl <- g0
    lo <- fl - 8 * bw; hi <- fl - 4 * bw; nSb <- sum(xs >= lo & xs < hi); nUb <- sum(xu >= lo & xu < hi)
    r0 <- (nSb / length(xs)) / (nUb / length(xu))
    gS <- stepDown(xs, xu, g0, fl, bw); truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    out[[length(out) + 1]] <- tibble(scenario = scenAll[[t]], i, r0, nUb, f1Now = f1(xs, l, g0), f1Step = f1(xs, l, gS),
      eNow = bsf(xs, xu, g0) / truth - 1, eStep = bsf(xs, xu, gS) / truth - 1) } }
r <- bind_rows(out); saveRDS(r, file.path(P, "gate_r0.rds")); options(width = 250)
q <- function(v) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(v), 100 * quantile(v, .1), 100 * quantile(v, .9))
res <- r |> group_by(scenario) |> summarise(cellsInBand = median(nUb), r0q10 = quantile(r0, .1), r0med = median(r0), r0q90 = quantile(r0, .9),
  on105 = sum(r0 <= 1.05), on110 = sum(r0 <= 1.10), on120 = sum(r0 <= 1.20),
  f1Now = median(f1Now), f1Step = median(f1Step), f1G105 = median(ifelse(r0 <= 1.05, f1Step, f1Now)), f1G110 = median(ifelse(r0 <= 1.10, f1Step, f1Now)), f1G120 = median(ifelse(r0 <= 1.2, f1Step, f1Now)),
  errNow = q(eNow), errStep = q(eStep), errG110 = q(ifelse(r0 <= 1.10, eStep, eNow)))
print(as.data.frame(res |> select(scenario:f1G120)), digits = 3); cat("\n"); print(as.data.frame(res |> select(scenario, errNow, errStep, errG110)))
