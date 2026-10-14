suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); options(width = 200)
for (t in c("g002_20k", "g0001", "g002", "h002")) {
  A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  band <- bind_rows(lapply(seq_along(A), function(i) { s <- A[[i]]$sim[[1]]; bw <- B[[i]]$filter[[1]]$info$antimode$bw
    b <- .getCpUnsLocPlateauBand(s$stim, s$uns, s$stim, bw); tibble(ind = as.character(i), flatX = b$flatX, nStim = b$nStim, nUns = b$nUns, nStimTotal = b$nStimTotal, nUnsTotal = b$nUnsTotal) }))
  al <- .getCpUnsLocPlateauAlign(band, 1.1, 1.25)
  cat("\n==", t, ": pooled", round(al$group$pooled, 3), "tau", round(al$group$tau, 3), "upper", round(al$group$upper, 3), "pass", al$group$pass, " allowed", sum(al$tbl$allowed), "of", nrow(al$tbl), " median band unstim", median(band$nUns), "\n")
  res <- bind_rows(lapply(seq_along(A), function(i) { a <- A[[i]]; s <- a$sim[[1]]; g <- a$threshold; bw <- B[[i]]$filter[[1]]$info$antimode$bw; lb <- a$filter[[1]]$info$minProbXPos
    if (!is.finite(g) || is.null(lb)) return(NULL); ed <- .getCpUnsLocPlateauEdge(s$stim, g, lb)
    f <- function(c, gap) .getCpUnsLocPlateauStepGate(s$stim, s$uns, g, ed, bw, 0.5, cautious = c, gapBw = gap)$gate
    tibble(i, plain = f(FALSE, Inf) < g, cautious = f(TRUE, Inf) < g, gap = f(FALSE, 1) < g, both = f(TRUE, 1) < g) }))
  cat("tubes moved: plain", sum(res$plain), " cautious only", sum(res$cautious), " gap only", sum(res$gap), " both", sum(res$both), "of", nrow(res), "\n")
}
