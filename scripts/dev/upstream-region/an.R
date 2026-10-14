suppressMessages(library(dplyr)); options(width=220)
P <- Sys.getenv("P"); tag <- commandArgs(TRUE)[1]; sid <- as.integer(commandArgs(TRUE)[2])
res <- readRDS(file.path(P, paste0("res_", tag, ".rds"))); dbg <- readRDS(file.path(P, paste0("dbg_", tag, ".rds")))
x <- readRDS("../stimgate-store/_tmp/sim/compare/freq_bs/current/collated/compare_raw.rds") |> filter(sim_id==sid, iter==1, method=="stimgate")
r <- res |> filter(method=="stimgate")
cat("parity max abs gate diff:", max(abs(r$threshold - x$threshold[match(r$sample, x$sample)])), "\n")
d <- bind_rows(lapply(dbg, function(s) { l <- s$lines[[1]]; tibble(ind=s$ind, initGate=s$threshold, xClearInit=l["xClearInit"], xDom=l["xDom"], xQual=l["xQual"], xAnti=l["xAntimode"], xSum=l["xSum"], bw=s$bw, bias=s$biasUns, reason=s$locReason) }))
d$finalGate <- r$threshold[match(as.integer(d$ind) %/% 2, r$sample)]
# posterior at gate from StimGate's own fitted prob
d$predAtGate <- sapply(dbg, function(s){ pf <- s$filter[[1]]$dataMod; NA })
print(as.data.frame(d), digits=3)
cat("which equals xSum: ", table(with(d, ifelse(abs(xSum-xQual)<1e-9 & (is.na(xAnti)|xAnti>xSum),"xQual", ifelse(!is.na(xAnti)&abs(xSum-xAnti)<1e-9,"xAnti","xClear")))), "\n")
s1 <- dbg[[1]]; str(s1$probFit[[1]], max.level=1); print(names(s1$densTblRaw[[1]]))
