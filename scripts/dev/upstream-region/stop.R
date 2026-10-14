P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
 r <- sapply(dbg, function(s) { m <- s$filter[[1]]$info$marginal; if (is.null(m)) NA else paste(m$stopReason %||% m$reason) })
 cat(t, "\n"); print(table(r, useNA="ifany"))
 s <- dbg[[which(!is.na(r))[1]]]$filter[[1]]$info$marginal; print(names(s)); if (!is.null(s$scanTbl)) print(head(as.data.frame(s$scanTbl), 4)) }
