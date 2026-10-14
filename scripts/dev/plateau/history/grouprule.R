suppressMessages(library(dplyr)); P <- Sys.getenv("P"); options(width = 220)
# Per tube: log ratio and its counting-noise variance (stim and unstim band counts; tube sizes treated as equal).
prep <- function(d, r, nU) d |> mutate(nS = pmax(.data[[r]] * .data[[nU]], 0.5), nUb = pmax(.data[[nU]], 0.5), y = log(nS / nUb), v = 1 / nS + 1 / nUb)
# Random-effects summary per group (DerSimonian-Laird) and the per-tube exclusion.
rule <- function(d, limit = log(1.10), tauMax = 0.15) {
  d |> mutate(w = 1 / v) |> group_by(grp) |> mutate(
    yFix = sum(w * y) / sum(w), Q = sum(w * (y - yFix)^2), k = n(),
    tau2 = pmax(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))), tau = sqrt(tau2),
    wR = 1 / (v + tau2), pooled = exp(sum(wR * y) / sum(wR)),
    groupOn = pooled <= exp(limit) & tau <= tauMax,
    tubeHigh = (y - limit) / sqrt(v) > 1.645,
    fires = groupOn & !tubeHigh) |> ungroup() }
summ <- function(d) d |> group_by(grp) |> summarise(tubes = n(), bandCells = median(nUb), pooled = first(pooled), tau = first(tau),
  groupOn = first(groupOn), tubesExcluded = sum(groupOn & tubeHigh), fires = sum(fires), .groups = "drop")
# Simulations: each scenario (20 tubes of one dataset) as a group
s <- readRDS(file.path(P, "gate_r0.rds")) |> filter(is.finite(r0)) |> mutate(grp = scenario) |> prep("r0", "nUb") |> rule()
cat("Simulations: one group per scenario\n"); print(as.data.frame(summ(s)), digits = 3)
f1 <- s |> group_by(grp) |> summarise(f1Now = median(f1Now), f1Step = median(f1Step), f1Rule = median(ifelse(fires, f1Step, f1Now)),
  q90Now = quantile(eNow, .9), q90Rule = quantile(ifelse(fires, eStep, eNow), .9))
cat("\nSimulations: median F1 and 90th-percentile frequency error under the rule\n"); print(as.data.frame(f1), digits = 3)
# Real data: unstim shifted by biasUns; groups are population x cytokine x stimulus (ACS) and strain/population x cytokine (OMIP-111)
r <- readRDS(file.path(P, "realratio.rds")) |> filter(is.finite(r0Bias), nUband >= 20) |> mutate(grp = paste(data, pop, cyt, stim)) |> prep("r0Bias", "nUband") |> rule()
g <- summ(r) |> mutate(data = sub(" .*", "", grp))
cat("\nReal data (unstim + biasUns, tubes with >= 20 unstim cells in the band)\n")
print(g |> group_by(data) |> summarise(groups = n(), groupsOn = sum(groupOn), tubes = sum(tubes), tubesFire = sum(fires), medTau = median(tau), medPooled = median(pooled)) |> as.data.frame(), digits = 3)
cat("\nGroups switched off (pooled ratio > 1.10 or tau > 0.15)\n"); print(as.data.frame(g |> filter(!groupOn) |> select(grp, tubes, bandCells, pooled, tau)), digits = 3)
cat("\nOMIP-111 groups\n"); print(as.data.frame(g |> filter(data == "OMIP-111") |> select(grp, tubes, bandCells, pooled, tau, groupOn, fires)), digits = 3)
