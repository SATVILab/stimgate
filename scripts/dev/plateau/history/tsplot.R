suppressMessages({devtools::load_all(quiet = TRUE); library(ggplot2); library(dplyr)})
P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
tgPoint <- function(x, bw, tol = 0.01) {
  d <- stats::density(x, bw = bw, n = 2048)
  dy <- diff(d$y) / diff(d$x); dx <- (d$x[-1] + d$x[-length(d$x)]) / 2
  pk <- d$x[.getPeakMainLeftIdx(d$y)]
  steep <- dx[dx > pk][which.min(dy[dx > pk])]
  cand <- dx > steep & abs(dy) < tol * max(abs(dy))
  if (!any(cand)) NA_real_ else dx[which(cand)[1]]
}
scen <- eval(parse(text = Sys.getenv("SCEN")))
f1 <- function(x, l, g) { if (!is.finite(g)) return(NA_real_); tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2*tp/(2*tp+fp+fn) }
outDir <- Sys.getenv("OUT", P)
tabs <- list()
for (t in names(scen)) {
  d <- readRDS(file.path(P, paste0("dbg_", t, ".rds")))
  curves <- list(); rugs <- list(); vl <- list()
  for (i in seq_along(d)) { if (i > 4) { e <- d[[i]]; s <- e$sim[[1]]; xs <- s$stim; n <- length(xs); bw <- stats::bw.nrd0(xs[seq(1, n, length.out = min(n, 10000))]); flat <- tgPoint(xs, bw) + bw / 2; sub <- sort(xs[xs >= flat]); ts <- ns$.getCpUnsLocAntimodeDensity(sub); am <- if (is.null(ts)) numeric() else ns$.getCpUnsLocTautStringExtremaCurrent(ts)$antimodes$x; amB <- am[am < e$threshold]
    tabs[[length(tabs) + 1]] <- tibble(scenario = scen[[t]], sample = i, flat = flat, nAbove = length(sub), nGpAbove = sum(s$labelsStim == "gp" & xs >= flat), modes = NA, antimodes = paste(sprintf("%.2f", am), collapse = ","), gate = e$threshold, f1Gate = f1(xs, s$labelsStim, e$threshold), f1Flat = f1(xs, s$labelsStim, flat), f1Anti = f1(xs, s$labelsStim, if (length(amB)) max(amB) else NA)); next }
    e <- d[[i]]; s <- e$sim[[1]]; xs <- s$stim; xu <- s$uns; n <- length(xs)
    bw <- stats::bw.nrd0(xs[seq(1, n, length.out = min(n, 10000))])
    flat <- tgPoint(xs, bw) + bw / 2
    lab <- paste0("Sample ", i)
    sub <- sort(xs[xs >= flat]); fr <- length(sub) / n
    ts <- ns$.getCpUnsLocAntimodeDensity(sub)
    ex <- if (is.null(ts)) list(modes = data.frame(x = numeric()), antimodes = data.frame(x = numeric())) else ns$.getCpUnsLocTautStringExtremaCurrent(ts)
    lbP <- e$filter[[1]]$info$minProbXPos
    tsP <- NULL; if (length(lbP) == 1L && is.finite(lbP)) { subP <- sort(xs[xs >= lbP]); tsP <- ns$.getCpUnsLocAntimodeDensity(subP) }
    rng <- c(-1, max(xs))
    ks <- stats::density(xs, bw = bw, n = 2048, from = rng[1], to = rng[2]); ku <- stats::density(xu, bw = bw, n = 2048, from = rng[1], to = rng[2])
    curves[[i]] <- bind_rows(
      tibble(x = ks$x, y = ks$y, curve = "Stim (kernel)"),
      tibble(x = ku$x, y = ku$y, curve = "Unstim (kernel)"),
      if (!is.null(ts)) tibble(x = ts$x, y = ts$y * fr, curve = "Taut string from flattening point"),
      if (!is.null(tsP)) tibble(x = tsP$x, y = tsP$y * length(subP) / n, curve = "Taut string as package fits it")) |> mutate(sample = lab)
    rugs[[i]] <- tibble(x = xs[s$labelsStim == "gp"], sample = lab)
    am <- ex$antimodes; md <- ex$modes
    vl[[i]] <- bind_rows(tibble(x = flat, line = "Flattening point"), tibble(x = e$threshold, line = "Current gate"),
                         tibble(x = am$x, line = rep("Antimode (from flattening point)", nrow(am)))) |> mutate(sample = lab)
    tabs[[length(tabs) + 1]] <- tibble(scenario = scen[[t]], sample = i, flat = flat, nAbove = length(sub), nGpAbove = sum(s$labelsStim == "gp" & xs >= flat),
      modes = paste(sprintf("%.2f", md$x), collapse = ","), antimodes = paste(sprintf("%.2f", am$x), collapse = ","), gate = e$threshold,
      f1Gate = f1(xs, s$labelsStim, e$threshold), f1Flat = f1(xs, s$labelsStim, flat),
      f1Anti = f1(xs, s$labelsStim, if (length(am$x[am$x < e$threshold])) max(am$x[am$x < e$threshold]) else NA))
  }
  cv <- bind_rows(curves) |> filter(y > 0); rg <- bind_rows(rugs); vv <- bind_rows(vl)
  ymin <- 0.3 / length(d[[1]]$sim[[1]]$stim)
  g <- ggplot(cv, aes(x, y, colour = curve, linetype = curve)) + geom_line(linewidth = 0.5) +
    geom_rug(data = rg, aes(x = x), inherit.aes = FALSE, sides = "b", colour = "firebrick", alpha = 0.5, length = unit(0.04, "npc")) +
    geom_vline(data = vv, aes(xintercept = x, linetype = line, colour = line), linewidth = 0.4) +
    scale_y_log10(limits = c(ymin, NA), oob = scales::squish) +
    scale_colour_manual(values = c("Stim (kernel)" = "black", "Unstim (kernel)" = "grey55", "Taut string from flattening point" = "firebrick",
      "Taut string as package fits it" = "steelblue", "Flattening point" = "darkorange", "Current gate" = "forestgreen", "Antimode (from flattening point)" = "firebrick")) +
    scale_linetype_manual(values = c("Stim (kernel)" = "solid", "Unstim (kernel)" = "solid", "Taut string from flattening point" = "solid",
      "Taut string as package fits it" = "dashed", "Flattening point" = "dashed", "Current gate" = "solid", "Antimode (from flattening point)" = "dotted")) +
    facet_wrap(~sample, ncol = 1) + labs(x = "Expression", y = "Density (share of all stim cells, log scale)", colour = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "bottom", strip.background = element_rect(fill = "white")) + guides(colour = guide_legend(ncol = 2))
  ggsave(file.path(outDir, paste0(t, ".png")), g, width = 9, height = 11, dpi = 110)
}
options(width = 200); tb <- bind_rows(tabs); saveRDS(tb, file.path(outDir, "per_tube.rds"))
print(tb |> group_by(scenario) |> summarise(n = n(), flat = median(flat), nAbove = median(nAbove), gpAbove = median(nGpAbove), gate = median(gate),
  nAnti = sum(is.finite(f1Anti)), f1Gate = median(f1Gate, na.rm = TRUE), f1Flat = median(f1Flat, na.rm = TRUE),
  f1GateAntiTubes = median(f1Gate[is.finite(f1Anti)]), f1Anti = median(f1Anti, na.rm = TRUE)) |> as.data.frame(), digits = 3)
