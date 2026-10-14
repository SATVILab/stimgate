suppressMessages(library(dplyr)); P <- Sys.getenv("P"); t <- commandArgs(TRUE)[1]
dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds"))); res <- readRDS(file.path(P, paste0("res_", t, ".rds"))) |> filter(method=="stimgate")
base <- readRDS(file.path(P, paste0("res_", commandArgs(TRUE)[2], ".rds"))) |> filter(method=="stimgate")
for (k in seq_along(dbg)) { s <- dbg[[k]]; smp <- as.integer(s$ind) %/% 2; r <- res[res$sample==smp,]; b <- base[base$sample==smp,]
 m <- s$filter[[1]]$info$marginal; st <- m$scanTbl
 acc <- if (!is.null(st) && nrow(st)) st[st$accepted %in% TRUE, ] else NULL
 cat(sprintf("s%2d est/true %.2f (base %.2f) gate %.2f (base %.2f) TP %d FP %d | start %.2f floor %.2f | accepted: %s\n", smp, r$propRespEst/r$propRespTruth, b$propRespEst/b$propRespTruth, r$threshold, b$threshold, r$nTruePos, r$nFalsePos, m$startX %||% NA, m$floorX %||% NA,
  if (is.null(acc) || !nrow(acc)) "none" else paste(sprintf("[%.2f-%.2f nS%d net%.1f pur%.2f bar%.2f]", acc$left, acc$right, acc$nCell, acc$expectedResp, acc$purity, acc$minPurity), collapse=" "))) }
