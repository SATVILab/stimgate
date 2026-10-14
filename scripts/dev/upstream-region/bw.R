P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
 o <- t(sapply(dbg[1:6], function(s) { m <- s$filter[[1]]$info$marginal; x <- s$sim[[1]]$stim; c(ours = m$bwTG %||% NA, hpi1 = ks::hpi(x, deriv.order = 1)) }))
 cat(t, "\n"); print(round(o, 3)) }
