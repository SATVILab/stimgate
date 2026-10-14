x <- readRDS("../stimgate-store/_tmp/sim/compare/freq_bs/current/collated/compare_raw.rds")
print(dim(x)); print(names(x)); print(table(x$method))
m <- readRDS(list.files("../stimgate-store/_tmp/sim/compare/freq_bs/current", "manifest", recursive=TRUE, full.names=TRUE)[1])
str(m[c("analysis_key")]); print(m$params$comparison_semantics_version)
