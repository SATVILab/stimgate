suppressMessages(library(dplyr)); P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) {
  dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds"))); res <- readRDS(file.path(P, paste0("res_", t, ".rds"))) |> filter(method=="stimgate")
  init <- sapply(dbg, function(s) s$threshold); resp <- sapply(dbg, function(s) isTRUE(s$cpattr[[1]]$locGeneratedDirect))
  fin <- res$threshold[match(as.integer(sapply(dbg, `[[`, "ind")) %/% 2, res$sample)]
  ev <- function(g) sapply(seq_along(dbg), function(i){ sm <- dbg[[i]]$sim[[1]]; gp <- sm$labelsStim=="gp"; pos <- sm$stim > g[i]; tp<-sum(pos&gp); c(f1=2*tp/(sum(pos)+sum(gp)), est=(mean(pos)-mean(sm$uns>g[i]))/(mean(gp)-mean(sm$labelsUns=="gp"))-1) })
  med <- median(init[resp]); q40 <- quantile(init[resp], .4)
  for (nm in c("final","pooled_median","pooled_q40")) { g <- switch(nm, final=fin, pooled_median=rep(med,length(fin)), pooled_q40=rep(q40,length(fin))); e <- ev(g)
    cat(sprintf("%-22s %-14s gate %.2f  medF1 %.3f  meanF1 %.3f  relerr med %.2f q10 %.2f q90 %.2f\n", t, nm, median(g), median(e[1,]), mean(e[1,]), median(e[2,]), quantile(e[2,],.1), quantile(e[2,],.9))) }
}
