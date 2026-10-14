suppressMessages(library(dplyr)); options(width=220)
P <- Sys.getenv("P"); tags <- commandArgs(TRUE)
f1 <- function(tp,fp,fn) 2*tp/(2*tp+fp+fn)
o <- bind_rows(lapply(tags, function(t) readRDS(file.path(P, paste0("res_", t, ".rds"))) |>
  mutate(f1=f1(nTruePos,nFalsePos,nFalseNeg), variant=t)))
print(o |> filter(method %in% c("stimgate","fbeta","tailgate")) |> filter(method=="stimgate" | variant==tags[1]) |>
  group_by(variant, method) |> summarise(gate=median(threshold,na.rm=T), sdGate=sd(threshold,na.rm=T), F1=median(f1,na.rm=T),
  err=median(propRespEst/propRespTruth-1,na.rm=T), q10=quantile(propRespEst/propRespTruth-1,.1,na.rm=T), q90=quantile(propRespEst/propRespTruth-1,.9,na.rm=T), .groups="drop") |> as.data.frame(), digits=3)
for (t in tags) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  d <- bind_rows(lapply(dbg, function(s){ f <- s$filter[[1]]$info; fin <- f$final; m <- f$marginal; if (is.null(fin)) return(NULL)
    tibble(xClear=fin$xClear %||% NA, xDom=fin$xDom %||% NA, xQual=fin$xQual %||% NA, xSum=fin$xSum %||% NA, xTG=m$xTG %||% NA, floor=m$floorX %||% NA, stop=m$stopReason %||% NA, top3=m$crowdTop3 %||% NA, maxC=m$maxLeftBinCells %||% NA)}))
  cat(t, "\n"); print(summarise(d, n=n(), across(c(xClear,xDom,xTG,floor,xQual,xSum,top3,maxC), ~median(.x,na.rm=T)), nMoved=sum(xQual < xClear - 1e-9, na.rm=T)) |> as.data.frame(), digits=3)
  print(table(d$stop, useNA="ifany")) }
