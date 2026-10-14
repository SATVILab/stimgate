suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); P <- Sys.getenv("P"); options(width = 200)
# gate-free edge: taut string on stim cells above the main peak; left edge of the widest flat step
edgeFree <- function(xs, bw) { d <- stats::density(xs, bw = bw, n = 2048); pk <- d$x[.getPeakMainLeftIdx(d$y)]
  sub <- sort(xs[xs >= pk]); den <- .getCpUnsLocAntimodeDensity(sub); if (is.null(den) || all(den$y == 0)) return(NA)
  r <- rle(den$y); ends <- cumsum(r$lengths); st <- ends - r$lengths + 1; w <- sub[pmin(ends + 1, length(sub))] - sub[st]
  k <- which.max(w); sub[st[k]] }
tags <- c("g002", "g002_20k", "g0005_20k", "g0001", "h002", "h002_20k", "h0005_20k", "h0001", "s002", "s0005", "s002_20k", "s0005_20k", "g002_20k_sd01", "h002_20k_sd01", "h002_20k_ms02")
out <- list()
for (t in tags) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in seq_along(A)) { s <- A[[i]]$sim[[1]]; g <- A[[i]]$threshold; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA; lb <- A[[i]]$filter[[1]]$info$minProbXPos
    if (!is.finite(g) || !is.finite(bw) || is.null(lb)) next
    out[[length(out) + 1]] <- tibble(t, i, gate = g, eGate = .getCpUnsLocPlateauEdge(s$stim, g, lb), eFree = edgeFree(s$stim, bw), bw = bw) } }
r <- bind_rows(out)
print(r |> group_by(t) |> summarise(n = n(), gate = median(gate), eGate = median(eGate, na.rm = TRUE), eFree = median(eFree, na.rm = TRUE),
  same = sum(abs(eGate - eFree) < 1e-9, na.rm = TRUE), within1bw = sum(abs(eGate - eFree) <= bw, na.rm = TRUE), freeAboveGate = sum(eFree >= gate, na.rm = TRUE)) |> as.data.frame(), digits = 3)
