suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
source(textConnection(paste(readLines(file.path(P, "tsplot.R"))[3:10], collapse = "\n")))
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l != "gn" & FALSE) + sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
out <- list()
for (t in c("am_g002","am_h002","am_s002")) { d <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  for (i in seq_along(d)) { e <- d[[i]]; s <- e$sim[[1]]; xs <- s$stim; n <- length(xs)
    bw <- stats::bw.nrd0(xs[seq(1, n, length.out = 10000)]); flat <- tgPoint(xs, bw) + bw/2
    ex <- ns$.getCpUnsLocTautStringExtremaCurrent(ns$.getCpUnsLocAntimodeDensity(sort(xs[xs >= flat])))
    am <- ex$antimodes$x; am <- if (length(am)) max(am[am < e$threshold], -Inf) else -Inf
    out[[length(out)+1]] <- tibble(t, i, gate = e$threshold, anti = am, f1Gate = f1(xs, s$labelsStim, e$threshold),
      f1Anti = if (is.finite(am)) f1(xs, s$labelsStim, am) else NA, f1Flat = f1(xs, s$labelsStim, flat)) } }
r <- bind_rows(out); print(r |> group_by(t) |> summarise(nAnti = sum(is.finite(anti)), gate = median(gate), anti = median(anti[is.finite(anti)]),
  f1Gate = median(f1Gate), f1AntiWhenFound = median(f1Anti, na.rm = TRUE), f1GateSameTubes = median(f1Gate[is.finite(anti)]), f1Flat = median(f1Flat)) |> as.data.frame(), digits = 3)
