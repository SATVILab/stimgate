suppressMessages(library(dplyr)); `%||%` <- function(a,b) if (is.null(a)) b else a
P <- Sys.getenv("P")
for (t in c("g01_base","g002_base","g0005_base","h002_base","h0.0005_base","h0.0001_base")) {
  d <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  r <- bind_rows(lapply(d, function(e) { f <- e$filter[[1]]$info; a <- f$antimode
    tibble(reason=a$reason %||% NA, nAnti=NROW(a$allAntimodes), nModes=NROW(a$allModes), nCand=NROW(a$candidateAntimodes), nRef=NROW(a$referenceModes),
      xBase=f$final$xBase %||% NA, lb=a$lowerBoundX %||% NA, nX=a$nExpressionValuesForDensity %||% NA,
      candMinRatio=if (NROW(a$candidateAntimodes) && NROW(a$referenceMode)) min(a$candidateAntimodes$height)/a$referenceModeHeight else NA) }))
  cat("\n==", t, "\n"); print(table(r$reason)); print(summarise(r, across(-reason, ~round(median(.x,na.rm=T),3))) |> as.data.frame())
}
