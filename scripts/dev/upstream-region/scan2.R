P <- Sys.getenv("P"); t <- commandArgs(TRUE)[1]; dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
for (k in seq_along(dbg)) { f <- dbg[[k]]$filter[[1]]$info; m <- f$marginal; fin <- f$final; if (is.null(m$scanTbl) || nrow(m$scanTbl)==0) next
 st <- m$scanTbl; a <- ifelse(st$accepted, "A", ifelse(st$acceptedPair %in% TRUE, "P", ifelse(st$pending, ".", "r")))
 cat(sprintf("%2d start %.2f floor %.2f -> cut %.2f | xQual %.2f xSum %.2f thr %.2f | %s | trim %s %d\n", k, m$startX, m$floorX, m$finalStartX, fin$xQual, fin$xSum, dbg[[k]]$threshold,
   paste(sprintf("%.2f:%s%.2f", st$left, a, st$purity), collapse=" "), m$trimReason, m$nStepsTrimmed)) }
