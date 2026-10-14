suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr); library(ggplot2)}); P <- Sys.getenv("P"); ns <- asNamespace("stimgate")
src <- readLines(file.path(P, "cmpab.R")); eval(parse(text = src[2:15]))
O <- "_tmp/taut-string-antimode/ratio"
scen <- c(g002 = "Clean, low sep, 100k", h002_20k = "Clean, high sep, 20k", h002 = "Clean, high sep, 100k",
  g002_sd01 = "SD x1.1, low sep, 100k", h002_20k_sd01 = "SD x1.1, high sep, 20k", g002_20k_sd01 = "SD x1.1, low sep, 20k",
  g002_ms02 = "Shift 0.2, low sep, 100k", h002_20k_ms02 = "Shift 0.2, high sep, 20k", g002_20k_ms02 = "Shift 0.2, low sep, 20k")
curves <- list(); marks <- list()
for (t in names(scen)) { A <- readRDS(file.path(P, paste0("dbg_am_", t, ".rds"))); B <- readRDS(file.path(P, paste0("dbg_pb_", t, ".rds")))
  for (i in 1:8) { a <- A[[i]]; s <- a$sim[[1]]; xs <- s$stim; xu <- s$uns; g <- a$threshold; bw <- B[[i]]$filter[[1]]$info$antimode$bw %||% NA
    lb <- a$filter[[1]]$info$minProbXPos; if (!is.finite(g) || !is.finite(bw)) next
    grid <- seq(g - 2, g + 1.5, length.out = 400)
    dS <- stats::density(xs, bw = bw, from = min(grid), to = max(grid), n = 400)$y; dU <- stats::density(xu, bw = bw, from = min(grid), to = max(grid), n = 400)$y
    eps <- 1 / (length(xs) * bw * sqrt(2 * pi)) / 2   # half of one cell's peak kernel height, to keep the ratio finite in empty tails
    gA <- stepDown(xs, xu, g, edge(sort(xs[xs >= lb]), g), bw)
    curves[[length(curves) + 1]] <- tibble(scenario = scen[[t]], tube = factor(i), x = grid - g, ratio = (dS + eps) / (dU + eps))
    marks[[length(marks) + 1]] <- tibble(scenario = scen[[t]], tube = factor(i), step = gA - g, f1Now = f1(xs, s$labelsStim, g), f1Step = f1(xs, s$labelsStim, gA)) } }
cv <- bind_rows(curves) |> mutate(scenario = factor(scenario, levels = scen)); mk <- bind_rows(marks) |> mutate(scenario = factor(scenario, levels = scen))
mk <- mk |> mutate(effect = case_when(f1Step > f1Now + 1e-9 ~ "Step helped", f1Step < f1Now - 1e-9 ~ "Step hurt", step < -1e-9 ~ "Moved, no change", TRUE ~ "Not moved"))
cv <- cv |> left_join(mk |> select(scenario, tube, effect), by = c("scenario", "tube"))
g <- ggplot(cv, aes(x, ratio, group = tube, colour = effect)) + geom_hline(yintercept = 1, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) + geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(data = mk |> filter(step < -1e-9) |> left_join(cv |> group_by(scenario, tube) |> summarise(.groups = "drop"), by = c("scenario", "tube")),
    aes(x = step, y = 0.5), shape = 17, size = 1.6, inherit.aes = FALSE, colour = "black") +
  scale_y_log10() + scale_colour_manual(values = c("Step helped" = "#1b7837", "Step hurt" = "#b2182b", "Moved, no change" = "#2166ac", "Not moved" = "grey60")) +
  facet_wrap(~scenario, ncol = 3) + coord_cartesian(ylim = c(0.3, 300)) +
  labs(x = "Distance from current gate (expression units)", y = "Smoothed stim / unstim density ratio (log scale)", colour = NULL) +
  theme_bw() + theme(legend.position = "bottom", strip.background = element_rect(fill = "white"))
ggsave(file.path(O, "stim_unstim_ratio.png"), g, width = 11, height = 10, dpi = 110)
saveRDS(list(curves = cv, marks = mk), file.path(O, "stim_unstim_ratio.rds"))
print(mk |> count(scenario, effect) |> tidyr::pivot_wider(names_from = effect, values_from = n, values_fill = 0) |> as.data.frame())
