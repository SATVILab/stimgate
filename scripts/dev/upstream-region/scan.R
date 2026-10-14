P <- Sys.getenv("P"); dbg <- readRDS(file.path(P, "dbg_g002_refpair_t2.rds"))
for (k in 1:3) { m <- dbg[[k]]$filter[[1]]$info$marginal; cat("start", m$startX, "floor", m$floorX, "maxC", m$maxLeftBinCells, "minPur", m$minLeftBinPurity, "\n"); print(head(as.data.frame(m$scanTbl)[, c("left","nCell","purity","accepted","acceptedPair")], 6), digits=3) }
