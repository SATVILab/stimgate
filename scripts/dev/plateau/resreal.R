suppressMessages(library(dplyr)); options(width = 230, dplyr.summarise.inform = FALSE)
S <- "/scratch/01463748/projects/2024/stimgate/stimgate-store"; B <- file.path(S, "_backup/2026-10-13-pre-plateau-step")
old <- read.csv(file.path(B, "acs_cytof/manual-comparison/manual-comparison.csv")); new <- read.csv(file.path(S, "_tmp/acs_cytof/manual-comparison/manual-comparison.csv"))
k <- c("SampleID", "stim", "popCode", "cyt")
j <- inner_join(old |> filter(method == "stimgate") |> select(all_of(k), thrOld = threshold, bsOld = freq_bs_auto, unsOld = freq_uns_auto, man = freq_bs_man),
                new |> filter(method == "stimgate") |> select(all_of(k), thrNew = threshold, bsNew = freq_bs_auto, unsNew = freq_uns_auto), by = k)
ccc <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; 2 * cov(x, y) * (length(x) - 1) / length(x) / (var(x) * (length(x) - 1) / length(x) + var(y) * (length(y) - 1) / length(y) + (mean(x) - mean(y))^2) }
cat("ACS StimGate tubes:", nrow(j), " lowered:", sum(j$thrNew < j$thrOld - 1e-9, na.rm = TRUE), " raised:", sum(j$thrNew > j$thrOld + 1e-9, na.rm = TRUE), "\n")
cat("Overall vs manual  PCC old/new:", round(cor(j$bsOld, j$man, use = "complete"), 3), round(cor(j$bsNew, j$man, use = "complete"), 3),
    "  CCC old/new:", round(ccc(j$bsOld, j$man), 3), round(ccc(j$bsNew, j$man), 3),
    "  median abs diff old/new:", round(median(abs(j$bsOld - j$man), na.rm = TRUE), 4), round(median(abs(j$bsNew - j$man), na.rm = TRUE), 4), "\n")
cat("Unstim freq above gate (%), 90th pct old/new:", round(quantile(j$unsOld, .9, na.rm = TRUE), 3), round(quantile(j$unsNew, .9, na.rm = TRUE), 3), "  max old/new:", round(max(j$unsOld, na.rm = TRUE), 2), round(max(j$unsNew, na.rm = TRUE), 2), "\n")
print(as.data.frame(j |> group_by(popCode, cyt) |> summarise(n = n(), lowered = sum(thrNew < thrOld - 1e-9, na.rm = TRUE), pccOld = cor(bsOld, man, use = "complete"), pccNew = cor(bsNew, man, use = "complete"),
  medBsOld = median(bsOld, na.rm = TRUE), medBsNew = median(bsNew, na.rm = TRUE), medMan = median(man, na.rm = TRUE)) |> filter(lowered > 0) |> mutate(across(where(is.double), ~ round(.x, 3)))))
# OMIP-111
o1 <- read.csv(file.path(B, "omip111/current/comparison.csv")); n1 <- read.csv(file.path(S, "_tmp/omip111/current/comparison.csv"))
cat("\nOMIP-111 comparison columns:", paste(head(names(n1), 30), collapse = ", "), "\n")
