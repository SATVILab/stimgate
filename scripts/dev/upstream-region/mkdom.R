P <- Sys.getenv("P"); a <- commandArgs(TRUE)
dbg <- readRDS(file.path(P, paste0("dbg_", a[1], ".rds")))
v <- sapply(dbg, function(s) { f <- s$filter[[1]]$info$final; if (is.null(f)) NA_real_ else as.numeric(f$xClear) })
names(v) <- sapply(dbg, `[[`, "ind"); print(v); saveRDS(v, file.path(P, paste0("xclear_", a[1], ".rds")))
