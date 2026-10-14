suppressMessages(library(dplyr)); P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
d <- bind_rows(lapply(dbg, function(s){ f <- s$filter[[1]]$info; di <- f$clear$densityDominance; fin <- f$final
 if (is.null(fin)) return(NULL)
 tibble(clearInit=fin$xClearInit, onset=di$onsetX %||% NA, peak=di$peakX %||% NA, xDom=fin$xDom, xQual=fin$xQual, xSum=fin$xSum)}))
cat(t, "\n"); print(summarise(d, n=n(), medClearInit=median(clearInit), medOnset=median(onset,na.rm=T), medPeak=median(peak,na.rm=T), medDom=median(xDom,na.rm=T), medQual=median(xQual), nQualBelowDom=sum(xQual < xDom - 1e-9, na.rm=T), medQualDrop=median((xDom-xQual)[xQual < xDom - 1e-9], na.rm=T)) |> as.data.frame(), digits=3) }
