suppressMessages(library(dplyr)); options(width=200)
x <- readRDS(file.path(Sys.getenv("P"),"x.rds"))
s <- x |> group_by(prob_response,n_cell,method) |> summarise(f1=median(f1,na.rm=TRUE), sens=median(sens,na.rm=TRUE), prec=median(prec,na.rm=TRUE), gate=median(threshold,na.rm=TRUE), npos=median(nTruePos+nFalseNeg), .groups="drop") |> filter(method!="tailgate_default")
print(as.data.frame(s), digits=3)
w <- x |> filter(method %in% c("stimgate","fbeta","tailgate")) |> select(sim_id,iter,sample,prob_response,n_cell,method,threshold,f1) |>
  tidyr::pivot_wider(names_from=method, values_from=c(threshold,f1))
print(w |> group_by(prob_response,n_cell) |> summarise(dSG_FB=median(threshold_stimgate-threshold_fbeta,na.rm=T), dSG_TG=median(threshold_stimgate-threshold_tailgate,na.rm=T), pSGhighTG=mean(threshold_stimgate>threshold_tailgate,na.rm=T), sdGateSG=sd(threshold_stimgate,na.rm=T), sdGateTG=sd(threshold_tailgate,na.rm=T), .groups="drop") |> as.data.frame(), digits=3)
# within dataset spread of stimgate gates
print(x |> filter(method=="stimgate") |> group_by(prob_response,n_cell,sim_id,iter) |> summarise(sd=sd(threshold), nuniq=length(unique(round(threshold,4))), .groups="drop") |> group_by(prob_response,n_cell) |> summarise(medsd=median(sd), meduniq=median(nuniq)) |> as.data.frame(), digits=3)
