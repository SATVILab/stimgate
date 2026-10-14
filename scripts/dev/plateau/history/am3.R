suppressMessages(library(dplyr)); `%||%` <- function(a,b) if (is.null(a)) b else a
P <- Sys.getenv("P"); options(width=200)
for (t in c("am_g002","am_h002","am_h0001","am_s002")) {
  f <- file.path(P, paste0("dbg_", t, ".rds")); if (!file.exists(f)) next
  d <- readRDS(f)
  r <- bind_rows(lapply(d, function(e) { i <- e$filter[[1]]$info; a <- i$antimode
    tibble(reason=a$reason %||% NA, lb=i$minProbXPos %||% NA, xBase=i$final$xBase %||% NA, nX=a$nExpressionValuesForDensity %||% NA,
      nModes=NROW(a$allModes), nAnti=NROW(a$allAntimodes), antiX=paste(sprintf("%.2f", a$allAntimodes$x %||% numeric()), collapse=","),
      refX=if (NROW(a$referenceMode)) a$referenceMode$x else NA, ratio=if (NROW(a$candidateAntimodes) && NROW(a$referenceMode)) min(a$candidateAntimodes$height)/a$referenceModeHeight else NA,
      gate=e$threshold, tp=e$nTruePos, pos=e$nPosStim) }))
  cat("\n==", t, "\n"); print(table(r$reason)); print(as.data.frame(head(r, 8)), digits=3)
}
