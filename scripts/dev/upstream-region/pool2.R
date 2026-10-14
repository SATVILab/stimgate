suppressMessages(library(dplyr)); P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds"))); res <- readRDS(file.path(P, paste0("res_", t, ".rds"))) |> filter(method=="stimgate")
 init <- sapply(dbg, function(s) s$threshold); fin <- res$threshold[match(as.integer(sapply(dbg, `[[`, "ind")) %/% 2, res$sample)]
 ev <- function(g) median(sapply(seq_along(dbg), function(i){ sm <- dbg[[i]]$sim[[1]]; gp <- sm$labelsStim=="gp"; tp <- sum(sm$stim > g[i] & gp); 2*tp/(sum(sm$stim > g[i])+sum(gp)) }))
 cat(sprintf("%-14s tube gates: med %.2f sd %.2f F1 %.3f | final: sd %.2f F1 %.3f | pooled median of tube gates %.2f F1 %.3f\n", t, median(init,na.rm=T), sd(init,na.rm=T), ev(init), sd(fin), ev(fin), median(init,na.rm=T), ev(rep(median(init,na.rm=T), length(init))))) }
