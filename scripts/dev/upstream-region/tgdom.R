suppressMessages(library(dplyr)); `%||%` <- function(a,b) if (is.null(a)) b else a
P <- Sys.getenv("P")
for (t in c("g01_nc_reg","g002_nc_reg","g0005_ce05_nc_reg","h002_nc_reg","h0.0005_nc_reg","h0.0001_nc_reg")) {
  d <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  r <- bind_rows(lapply(d, function(e) { f <- e$filter[[1]]$info; fin <- f$final; m <- f$marginal
    tibble(xDom = fin$xDom %||% NA, xTG = m$xTG %||% NA) }))
  cat(sprintf("%-18s n=%2d  TG<xDom in %2d  median xDom %.2f  median xTG %.2f  xTG-xDom range [%.2f, %.2f]\n", t, nrow(r),
    sum(r$xTG < r$xDom, na.rm=TRUE), median(r$xDom,na.rm=T), median(r$xTG,na.rm=T), min(r$xTG-r$xDom,na.rm=T), max(r$xTG-r$xDom,na.rm=T)))
}
