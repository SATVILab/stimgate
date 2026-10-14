suppressMessages(library(dplyr)); P <- Sys.getenv("P")
for (t in commandArgs(TRUE)) { dbg <- readRDS(file.path(P, paste0("dbg_", t, ".rds"))); res <- readRDS(file.path(P, paste0("res_", t, ".rds")))
 tg <- res |> filter(method=="tailgate"); 
 d <- bind_rows(lapply(dbg, function(s){ m <- s$filter[[1]]$info$marginal; fin <- s$filter[[1]]$info$final; if (is.null(m)) return(NULL)
   tibble(sample=as.integer(s$ind) %/% 2, start=m$startX, xTG=m$xTG %||% NA, floor=m$floorX %||% NA, xDom=fin$xDom %||% NA, cut=m$finalStartX %||% NA, stop=m$stopReason %||% NA, nAcc=sum(m$scanTbl$accepted %in% TRUE))}))
 d$tailgate <- tg$threshold[match(d$sample, tg$sample)]
 cat(t, "\n"); print(summarise(d, across(c(xDom,start,xTG,tailgate,floor,cut,nAcc), ~median(.x,na.rm=T)), corTG=cor(xTG,tailgate,use="complete")) |> as.data.frame(), digits=3); print(table(d$stop)) }
