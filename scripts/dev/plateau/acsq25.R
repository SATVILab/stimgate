suppressMessages({library(dplyr)}); options(width = 200, dplyr.summarise.inform = FALSE)
S <- "/scratch/01463748/projects/2024/stimgate/stimgate-store/_tmp"; P <- Sys.getenv("P")
r <- readRDS(file.path(P, "acsedge.rds"))
acs <- read.csv(file.path(S, "acs_cytof/manual-comparison/manual-comparison.csv")) |> filter(method == "stimgate")
chMap <- c(IFNg = "Ho165Di", IL2 = "Gd158Di", TNF = "Nd146Di", TNFa = "Nd146Di", IL17 = "Dy164Di", IL6 = "Gd156Di", IL22 = "Nd150Di")
resp <- bind_rows(lapply(unique(acs$popCode), function(pc) bind_rows(lapply(unique(acs$cyt), function(cy) {
  g <- readRDS(file.path(S, "acs_cytof/scratch", pc, "stimgate/gates/poproot", paste0("chnl", chMap[[cy]]), "all/gateTbl.rds"))
  g <- g[g$gateName == "loc_minClust", ]; tibble(pop = pc, cyt = cy, ind = as.character(g$ind), responder = g$locResponder %in% TRUE, gateTbl = g$gate) }))))
a <- acs |> transmute(pop = popCode, cyt, stim, donor = SampleID, ind = as.character(ind))
r <- r |> left_join(a, by = c("pop", "cyt", "stim", "donor")) |> left_join(resp, by = c("pop", "cyt", "ind"))
cat("responder flag matched:", sum(!is.na(r$responder)), "of", nrow(r), "; responders:", sum(r$responder, na.rm = TRUE), "\n")
xsCache <- new.env()
readEx <- function(pp, ind, ch) readRDS(file.path(pp, "sampleData", "pop_root", paste0("ind_", ind), paste0("chnl_", ch, ".rds")))
r <- r |> group_by(pop, cyt) |> mutate(q25 = if (sum(responder %in% TRUE) >= 3) quantile(gate[responder %in% TRUE], 0.25) else NA_real_) |> ungroup() |>
  mutate(floor = pmax(q25, edge), gFloor = ifelse(is.finite(floor) & gate > floor, floor, gate))
# frequencies at the new gate
bl <- list(); out <- numeric(nrow(r)); unsOut <- numeric(nrow(r))
for (i in seq_len(nrow(r))) { pp <- file.path(S, "acs_cytof/scratch", r$pop[i], "stimgate"); ch <- chMap[[r$cyt[i]]]
  if (is.null(bl[[r$pop[i]]])) { b <- readRDS(file.path(pp, "metaData/batchList.rds")); bl[[r$pop[i]]] <- unlist(lapply(unname(b), function(z) stats::setNames(rep(z[1], length(z) - 1), z[-1]))) }
  xs <- readEx(pp, r$ind[i], ch); xu <- readEx(pp, bl[[r$pop[i]]][[r$ind[i]]], ch); g <- r$gFloor[i]
  out[i] <- 100 * (mean(xs > g) - mean(xu > g)); unsOut[i] <- 100 * mean(xu > g) }
r$bsFloor <- out; r$unsFloor <- unsOut
ccc <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; 2 * mean((x - mean(x)) * (y - mean(y))) / (mean((x - mean(x))^2) + mean((y - mean(y))^2) + (mean(x) - mean(y))^2) }
cat("moved:", sum(r$gFloor < r$gate - 1e-9), "of", nrow(r), " (responders moved:", sum(r$gFloor < r$gate - 1e-9 & r$responder %in% TRUE), ", non-responders moved:", sum(r$gFloor < r$gate - 1e-9 & !(r$responder %in% TRUE)), ")\n")
cat("all: PCC now/floor", round(cor(r$bsNow, r$man, use = "complete"), 3), round(cor(r$bsFloor, r$man, use = "complete"), 3), " CCC now/floor", round(ccc(r$bsNow, r$man), 3), round(ccc(r$bsFloor, r$man), 3),
  " unstim q90 now/floor", round(quantile(r$unsNow, .9), 3), round(quantile(r$unsFloor, .9), 3), "\n")
print(as.data.frame(r |> group_by(pop) |> summarise(moved = mean(gFloor < gate - 1e-9), pccNow = cor(bsNow, man, use = "complete"), pccFloor = cor(bsFloor, man, use = "complete"),
  maeNow = mean(abs(bsNow - man), na.rm = TRUE), maeFloor = mean(abs(bsFloor - man), na.rm = TRUE), unsQ90Now = quantile(unsNow, .9), unsQ90Floor = quantile(unsFloor, .9))), digits = 3)
s <- r |> filter(stim != "p4") |> group_by(pop, cyt, stim) |> summarise(pccNow = cor(bsNow, man, use = "complete"), pccFloor = cor(bsFloor, man, use = "complete"), .groups = "drop")
cat("\nMean PCC over non-p4 strata now/floor:", round(mean(s$pccNow, na.rm = TRUE), 3), round(mean(s$pccFloor, na.rm = TRUE), 3), "\n")
