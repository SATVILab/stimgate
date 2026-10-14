suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); options(width = 200, dplyr.summarise.inform = FALSE)
S <- "/scratch/01463748/projects/2024/stimgate/stimgate-store/_tmp"
acs <- read.csv(file.path(S, "acs_cytof/manual-comparison/manual-comparison.csv")) |> filter(method == "stimgate")
chMap <- c(IFNg = "Ho165Di", IL2 = "Gd158Di", TNF = "Nd146Di", TNFa = "Nd146Di", IL17 = "Dy164Di", IL6 = "Gd156Di", IL22 = "Nd150Di")
readEx <- function(pp, ind, ch) readRDS(file.path(pp, "sampleData", "pop_root", paste0("ind_", ind), paste0("chnl_", ch, ".rds")))
out <- list()
for (pc in unique(acs$popCode)) { pp <- file.path(S, "acs_cytof/scratch", pc, "stimgate"); bl <- readRDS(file.path(pp, "metaData/batchList.rds"))
  cs <- readRDS(file.path(pp, "metaData/chnlSettings.rds")); lab <- readRDS(file.path(pp, "metaData/chnlLab.rds"))
  unsOf <- unlist(lapply(unname(bl), function(b) stats::setNames(rep(b[1], length(b) - 1), b[-1])))
  a <- acs |> filter(popCode == pc)
  for (k in seq_len(nrow(a))) { ch <- chMap[[a$cyt[k]]]; bw <- cs[[lab[[ch]]]]$bwShared; g <- a$threshold[k]; if (!is.finite(g)) next
    xs <- readEx(pp, a$ind[k], ch); xu <- readEx(pp, unsOf[[as.character(a$ind[k])]], ch)
    nz <- xs[xs > 0]; if (length(nz) < 100) next
    d <- stats::density(nz, bw = bw, n = 1024); pk <- d$x[.getPeakMainLeftIdx(d$y)]
    e <- .getCpUnsLocPlateauEdge(xs, g, pk); ge <- if (is.finite(e)) min(g, e + bw) else g
    out[[length(out) + 1]] <- tibble(pop = pc, cyt = a$cyt[k], stim = a$stim[k], donor = a$SampleID[k], gate = g, edge = e, gEdge = ge, bw = bw,
      bsNow = 100 * (mean(xs > g) - mean(xu > g)), bsEdge = 100 * (mean(xs > ge) - mean(xu > ge)), unsNow = 100 * mean(xu > g), unsEdge = 100 * mean(xu > ge), man = a$freq_bs_man[k]) } }
r <- bind_rows(out); saveRDS(r, file.path(Sys.getenv("P"), "acsedge.rds"))
cat("tubes", nrow(r), " edge found", sum(is.finite(r$edge)), " gate moved", sum(r$gEdge < r$gate - 1e-9), " median move (bandwidths)", round(median(((r$gate - r$gEdge) / r$bw)[r$gEdge < r$gate]), 1), "\n")
ccc <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; 2 * mean((x - mean(x)) * (y - mean(y))) / (mean((x - mean(x))^2) + mean((y - mean(y))^2) + (mean(x) - mean(y))^2) }
cat("all tubes: PCC now/edge", round(cor(r$bsNow, r$man, use = "complete"), 3), round(cor(r$bsEdge, r$man, use = "complete"), 3), " CCC now/edge", round(ccc(r$bsNow, r$man), 3), round(ccc(r$bsEdge, r$man), 3),
  " 90th pct unstim above gate (%) now/edge", round(quantile(r$unsNow, .9), 3), round(quantile(r$unsEdge, .9), 3), "\n")
print(as.data.frame(r |> group_by(pop) |> summarise(tubes = n(), moved = mean(gEdge < gate - 1e-9), pccNow = cor(bsNow, man, use = "complete"), pccEdge = cor(bsEdge, man, use = "complete"),
  maeNow = mean(abs(bsNow - man), na.rm = TRUE), maeEdge = mean(abs(bsEdge - man), na.rm = TRUE), unsQ90Now = quantile(unsNow, .9), unsQ90Edge = quantile(unsEdge, .9))), digits = 3)
cat("\nTCRgd by cytokine:\n")
print(as.data.frame(r |> filter(pop == "tcrgd") |> group_by(cyt) |> summarise(moved = mean(gEdge < gate - 1e-9), pccNow = cor(bsNow, man, use = "complete"), pccEdge = cor(bsEdge, man, use = "complete"),
  medBsNow = median(bsNow), medBsEdge = median(bsEdge), medMan = median(man), unsQ90Now = quantile(unsNow, .9), unsQ90Edge = quantile(unsEdge, .9))), digits = 3)
