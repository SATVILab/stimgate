suppressMessages(library(dplyr)); options(width = 230, dplyr.summarise.inform = FALSE)
S <- "/scratch/01463748/projects/2024/stimgate/stimgate-store"
new <- readRDS(file.path(S, "_tmp/sim/compare/freq_bs/current/collated/compare_raw.rds")) |> filter(method == "stimgate")
r <- readRDS(file.path(P <- Sys.getenv("P"), "res7.rds"))$m
cat("plateau column present:", "locPlateauStep" %in% names(new), "\n")
sc <- new |> filter(mean_pos_setting == "low", (transformation == "skew" & prob_response <= 5e-4) | (transformation == "gaussian" & prob_response <= 5e-4))
print(as.data.frame(sc |> group_by(transformation, prob_response, n_cell) |> summarise(tubes = n(),
  noResponse = mean(grepl("No responding", locReason)), direct = mean(locGeneratedDirect %in% TRUE), shared = mean(locSource %in% "combined"),
  medPos = median(nTruePos + nFalseNeg), medTPabove = median(nTruePos), medFPabove = median(nFalsePos), medGate = median(threshold))), digits = 3)
cat("\nReasons:\n"); print(sc |> count(transformation, prob_response, n_cell, locReason) |> group_by(transformation, prob_response, n_cell) |> mutate(share = round(n / sum(n), 2)) |> as.data.frame())
cat("\nAmong tubes whose gate was lowered: F1 old -> new\n")
print(as.data.frame(r |> filter(mean_pos_setting == "low", prob_response <= 5e-4, transformation %in% c("gaussian", "skew")) |> group_by(transformation, prob_response, n_cell) |>
  summarise(lowered = sum(gate_new < gate_old - 1e-9), dGate = median((gate_new - gate_old)[gate_new < gate_old - 1e-9]), f1OldLow = median(F1_old[gate_new < gate_old - 1e-9]), f1NewLow = median(F1_new[gate_new < gate_old - 1e-9]),
  f1OldAll = median(F1_old), f1NewAll = median(F1_new), meanOld = mean(F1_old, na.rm = TRUE), meanNew = mean(F1_new, na.rm = TRUE))), digits = 3)
