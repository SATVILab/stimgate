suppressMessages(library(dplyr))
x <- readRDS("../stimgate-store/_tmp/sim/compare/freq_bs/current/collated/compare_raw.rds")
print(names(x)[1:30])
x <- x |> filter(transformation=="gaussian", mean_pos_setting=="low", prob_response %in% c(0.01,0.002,0.0005)) |>
  mutate(prec=nTruePos/(nTruePos+nFalsePos), sens=nTruePos/(nTruePos+nFalseNeg), f1=2*nTruePos/(2*nTruePos+nFalsePos+nFalseNeg))
s <- x |> group_by(prob_response,n_cell,method) |> summarise(n=n(), f1=median(f1,na.rm=TRUE), sens=median(sens,na.rm=TRUE), prec=median(prec,na.rm=TRUE), gate=median(threshold,na.rm=TRUE), relerr=median(propRespEst/propRespTruth-1,na.rm=TRUE), .groups="drop") |> filter(method!="tailgate_default")
print(as.data.frame(s), digits=3)
print(x |> filter(method=="stimgate") |> count(prob_response, locSource, locReason) |> as.data.frame())
saveRDS(x, file.path(Sys.getenv("P"), "x.rds"))
