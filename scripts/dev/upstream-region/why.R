suppressMessages(library(dplyr)); P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
 st <- bind_rows(lapply(seq_along(dbg), function(k){ m <- dbg[[k]]$filter[[1]]$info$marginal; s <- m$scanTbl; if (is.null(s) || !nrow(s)) return(NULL)
   s$k <- k; s$width <- s$right - s$left; s$failPur <- !(s$purity >= s$minPurity); s$failCrowd <- !(s$nCell <= s$maxCells); s }))
 last <- st |> group_by(k) |> filter(right == min(right)) |> ungroup()   # spans tried at the final cut (all rejected)
 cat(t, ": spans tried at final cut:", nrow(last), " fail purity", sum(last$failPur, na.rm=T), " fail crowding", sum(last$failCrowd, na.rm=T), " empty", sum(last$nCell==0), "\n")
 print(head(as.data.frame(last |> filter(k %in% unique(k)[1:2]) |> select(k,left,right,nCell,purity,minPurity,maxCells)), 10), digits=3) }
