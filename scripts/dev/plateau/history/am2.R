suppressMessages(devtools::load_all(quiet=TRUE)); P <- Sys.getenv("P")
for (t in c("g01_base","g002_base","h002_base","h0.0001_base")) for (i in 1:3) {
  e <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))[[i]]; x <- e$sim[[1]]$stim; xb <- e$filter[[1]]$info$final$xBase
  for (lb in c(-Inf, 3)) {
    den <- .getCpUnsLocAntimodeDensity(sort(x[x >= lb])); ex <- .getCpUnsLocTautStringExtremaCurrent(den)
    m <- ex$modes; a <- ex$antimodes; ref <- m[m$x > xb, ]
    cat(sprintf("%-13s s%d lb=%4s xBase=%.2f modes=%d antimodes=%d [%s] modesRightOfBase=%d [%s]\n", t, i, lb, xb, nrow(m), nrow(a),
      paste(sprintf("%.2f", a$x), collapse=","), nrow(ref), paste(sprintf("%.2f", ref$x), collapse=",")))
  }
}
