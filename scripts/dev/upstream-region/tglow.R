suppressMessages(library(dplyr)); options(width=250, dplyr.summarise.inform=FALSE)
S <- "../stimgate-store/_tmp/sim/compare"
sm <- function(d, ...) d |> filter(method %in% c("stimgate","tailgate","fbeta")) |>
  group_by(..., method) |>
  summarise(gate=median(threshold,na.rm=TRUE), prec=median(nTruePos/pmax(nTruePos+nFalsePos,1),na.rm=TRUE),
    fpGtTp=mean(nFalsePos>nTruePos,na.rm=TRUE), ratioQ50=median(propRespEst/propRespTruth,na.rm=TRUE),
    ratioQ90=quantile(propRespEst/propRespTruth,.9,na.rm=TRUE), over2=mean(propRespEst>2*propRespTruth,na.rm=TRUE)) |>
  tidyr::pivot_wider(names_from=method, values_from=c(gate,prec,fpGtTp,ratioQ90,over2)) |> select(-ratioQ50) |> as.data.frame()
x <- readRDS(file.path(S,"freq_bs/current/collated/compare_raw.rds"))
cat("Analysis 7, 100k cells, prob <= 0.002\n")
x7 <- x |> filter(n_cell==max(n_cell), prob_response <= 0.002)
r <- sm(x7, transformation, mean_pos_setting, prob_response); print(r |> mutate(across(where(is.numeric), ~signif(.x,3))))
cat("\nAnalysis 7 smaller cell counts: share of tubes where Tailgate over-calls >2x, by trans x sep (prob<=0.002)\n")
print(x |> filter(prob_response<=0.002, method %in% c("tailgate","stimgate")) |> group_by(n_cell, transformation, mean_pos_setting, method) |>
  summarise(over2=round(mean(propRespEst>2*propRespTruth,na.rm=TRUE),2)) |> tidyr::pivot_wider(names_from=c(method), values_from=over2) |> as.data.frame())
y <- readRDS(file.path(S,"null/current/collated/compare_raw.rds"))
cat("\nAnalysis 7b null: median and q90 of propRespEst (percentage points), share positive\n")
print(y |> filter(method %in% c("stimgate","tailgate","fbeta")) |> group_by(transformation, mean_pos_setting, prob_background, n_cell, method) |>
  summarise(q50=round(100*median(propRespEst,na.rm=TRUE),4), q90=round(100*quantile(propRespEst,.9,na.rm=TRUE),4), pos=round(mean(propRespEst>0,na.rm=TRUE),2)) |>
  tidyr::pivot_wider(names_from=method, values_from=c(q50,q90,pos)) |> as.data.frame())
z <- readRDS(file.path(S,"freq_bs_batch/current/collated/compare_raw.rds"))
cat("\nAnalysis 8: share over-calling >2x, prob<=0.002, by mismatch\n")
print(z |> filter(prob_response<=0.002, method %in% c("tailgate","stimgate")) |> group_by(mismatch_type, mismatch_val, transformation, mean_pos_setting, method) |>
  summarise(over2=round(mean(propRespEst>2*propRespTruth,na.rm=TRUE),2), q90=round(quantile(propRespEst/propRespTruth,.9,na.rm=TRUE),2)) |>
  tidyr::pivot_wider(names_from=method, values_from=c(over2,q90)) |> as.data.frame())
