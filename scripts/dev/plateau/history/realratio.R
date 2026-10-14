# Baseline stim/unstim ratio (band 4-8 bandwidths below the taut-string plateau edge) on real data, read-only on saved caches.
suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
S <- "../stimgate-store/_tmp"
edgeAt <- function(sub, g) { den <- ns$.getCpUnsLocAntimodeDensity(sub); if (is.null(den) || all(den$y == 0)) return(NA_real_)
  r <- rle(den$y); ends <- cumsum(r$lengths); st <- ends - r$lengths + 1; j <- findInterval(g, den$x); k <- which(st <= j & ends >= j)[1]
  if (is.na(k) || st[k] <= 1) NA_real_ else sub[st[k]] }
tubeRatio <- function(xs, xu, g, bw, bias) {
  nz <- xs[xs > 0]; if (length(nz) < 100 || !is.finite(g)) return(NULL)
  d <- stats::density(nz, bw = bw, n = 1024); pk <- d$x[ns$.getPeakMainLeftIdx(d$y)]
  fl <- edgeAt(sort(xs[xs >= pk]), g); if (!is.finite(fl)) return(tibble(edge = NA_real_, r0 = NA_real_, r0Bias = NA_real_, nUband = NA_integer_))
  lo <- fl - 8 * bw; hi <- fl - 4 * bw
  sh <- function(v) mean(v >= lo & v < hi)
  tibble(edge = fl, gate = g, nUband = sum(xu >= lo & xu < hi), r0 = sh(xs) / sh(xu), r0Bias = sh(xs) / sh(xu + bias)) }
readEx <- function(pp, ind, chnl) readRDS(file.path(pp, "sampleData", "pop_root", paste0("ind_", ind), paste0("chnl_", chnl, ".rds")))
out <- list()
# ACS (Analysis 9b caches): tubes, stimuli and applied gates from the manual-comparison table
acs <- read.csv(file.path(S, "acs_cytof/manual-comparison/manual-comparison.csv")) |> filter(method == "stimgate")
chMap <- c(IFNg = "Ho165Di", IL2 = "Gd158Di", TNF = "Nd146Di", TNFa = "Nd146Di", IL17 = "Dy164Di", IL6 = "Gd156Di", IL22 = "Nd150Di")
for (pc in unique(acs$popCode)) { pp <- file.path(S, "acs_cytof/scratch", pc, "stimgate"); bl <- readRDS(file.path(pp, "metaData/batchList.rds"))
  cs <- readRDS(file.path(pp, "metaData/chnlSettings.rds")); lab <- readRDS(file.path(pp, "metaData/chnlLab.rds"))
  unsOf <- unlist(lapply(unname(bl), function(b) stats::setNames(rep(b[1], length(b) - 1), b[-1])))
  a <- acs |> filter(popCode == pc)
  for (k in seq_len(nrow(a))) { ch <- chMap[[a$cyt[k]]]; set <- cs[[lab[[ch]]]]; ind <- a$ind[k]; iu <- unsOf[[as.character(ind)]]
    rr <- tubeRatio(readEx(pp, ind, ch), readEx(pp, iu, ch), a$threshold[k], set$bwShared, set$biasUns %||% 0)
    if (!is.null(rr)) out[[length(out) + 1]] <- tibble(data = "ACS", pop = pc, cyt = a$cyt[k], stim = a$stim[k], donor = a$SampleID[k], responder = !grepl("No responding", a$locReason[k])) |> bind_cols(rr) } }
# OMIP-111 (Analysis 14 caches)
for (strain in c("C57", "BALB")) for (pop in c("CD4", "CD8")) { pp <- file.path(S, "omip111/current", strain, pop, "stimgate"); if (!dir.exists(pp)) next
  bl <- readRDS(file.path(pp, "metaData/batchList.rds")); cs <- readRDS(file.path(pp, "metaData/chnlSettings.rds"))
  for (cyt in names(cs)) { gt <- readRDS(file.path(pp, "gates/poproot", paste0("chnl", cyt), "all/gateTbl.rds")); gt <- gt[gt$gateName == "loc_minClust", ]
    for (b in names(bl)) for (ind in bl[[b]][-1]) { g <- gt$gate[gt$ind == ind][1]; rsn <- gt$locReason[gt$ind == ind][1]
      rr <- tubeRatio(readEx(pp, ind, cyt), readEx(pp, bl[[b]][1], cyt), g, cs[[cyt]]$bwShared, cs[[cyt]]$biasUns %||% 0)
      if (!is.null(rr)) out[[length(out) + 1]] <- tibble(data = "OMIP-111", pop = paste(strain, pop), cyt = cyt, stim = "stim", donor = b, responder = !grepl("No responding", rsn)) |> bind_cols(rr) } } }
r <- bind_rows(out); saveRDS(r, file.path(P, "realratio.rds")); options(width = 220)
cat("Tubes:", nrow(r), " with a ratio:", sum(is.finite(r$r0)), "\n\nPer data set: share of tubes aligned (ratio <= 1.10), raw unstim and unstim + biasUns\n")
print(r |> group_by(data) |> summarise(n = sum(is.finite(r0)), medBandCells = median(nUband, na.rm = TRUE), q10 = quantile(r0, .1, na.rm = TRUE), med = median(r0, na.rm = TRUE), q90 = quantile(r0, .9, na.rm = TRUE),
  aligned = mean(r0 <= 1.10, na.rm = TRUE), alignedBias = mean(r0Bias <= 1.10, na.rm = TRUE)) |> as.data.frame(), digits = 3)
cat("\nPer population x cytokine x stimulus (groups of donors)\n")
grp <- r |> filter(is.finite(r0)) |> group_by(data, pop, cyt, stim) |> summarise(n = n(), bandCells = median(nUband), med = median(r0), q10 = quantile(r0, .1), q90 = quantile(r0, .9),
  aligned = mean(r0 <= 1.10), alignedBias = mean(r0Bias <= 1.10), .groups = "drop")
print(as.data.frame(grp), digits = 3)
