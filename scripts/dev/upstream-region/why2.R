suppressMessages(library(dplyr)); P <- Sys.getenv("P"); t <- commandArgs(TRUE)[1]
dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds"))); res <- readRDS(file.path(P, paste0("res_", t, ".rds"))) |> filter(method=="stimgate")
for (k in seq_along(dbg)) { m <- dbg[[k]]$filter[[1]]$info$marginal; s <- m$scanTbl; if (is.null(s) || !nrow(s)) { cat(k, "no scan:", m$reason %||% "-", "\n"); next }
 r <- s |> filter(right == min(right))
 cat(sprintf("%2d start %.2f cut %.2f floor %.2f | final-cut spans: %s\n", k, m$startX, m$finalStartX, m$floorX,
  paste(sprintf("[%.2f n%d p%.2f/bar%.2f lim%.0f]", r$left, r$nCell, r$purity, r$minPurity, r$maxCells), collapse=" "))) }
cat("stimgate reasons:\n"); print(table(res$locReason))
