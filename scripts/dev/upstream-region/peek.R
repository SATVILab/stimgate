S <- "../stimgate-store/_tmp/sim/compare"
for (a in c("freq_bs","null","freq_bs_batch")) { x <- readRDS(file.path(S,a,"current/collated/compare_raw.rds")); cat(a, nrow(x), "\n"); print(names(x)) }
x <- readRDS(file.path(S,"freq_bs","current/collated/compare_raw.rds")); print(head(as.data.frame(x),3))
