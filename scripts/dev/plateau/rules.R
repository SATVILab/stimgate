suppressMessages({devtools::load_all(quiet = TRUE); library(dplyr)}); options(width = 250, dplyr.summarise.inform = FALSE)
P <- Sys.getenv("P"); S <- "/scratch/01463748/projects/2024/stimgate/stimgate-store"
on7 <- readRDS(file.path(S, "_tmp/sim/compare/freq_bs/current/collated/compare_raw.rds")) |> filter(method == "stimgate", iter == 1)
off7 <- readRDS(file.path(S, "_backup/2026-10-13-pre-plateau-step/sim/compare/freq_bs/current/collated/compare_raw.rds")) |> filter(method == "stimgate", iter == 1)
f1 <- function(x, l, g) { tp <- sum(x > g & l == "gp"); fp <- sum(x > g & l == "gn"); fn <- sum(x <= g & l == "gp"); 2 * tp / (2 * tp + fp + fn) }
# Richards et al. F-beta, as scripts/python/fbeta.py (common bins over both tubes, sqrt(n) bins, moving mean of 10 bins)
fbetaCurve <- function(neg, pos, beta = 0.8, theta = 2, width = 10) {
  nb <- floor(sqrt(max(length(neg), length(pos)))); br <- seq(min(c(neg, pos)), max(c(neg, pos)), length.out = nb + 1)
  dens <- function(v) { h <- graphics::hist(v, breaks = br, plot = FALSE); h$counts / (length(v) * diff(br)) }
  mm <- function(x) { cs <- cumsum(x); c(rep(NA, width - 1), (cs[width:length(cs)] - c(0, cs[seq_len(length(cs) - width)])) / width) }
  pn <- mm(dens(neg)); pp <- mm(dens(pos)); n <- length(pn)
  fpos <- ifelse(!is.na(pp) & !is.na(pn) & pp > theta * pn, pp - pn, 0)
  tp <- rev(cumsum(rev(fpos))); fn <- c(0, cumsum(fpos))[seq_len(n)]; fp <- rev(cumsum(rev(ifelse(is.na(pn), NA, pn))))
  prec <- ifelse(tp == 0, 0, tp / (tp + fp)); rec <- tp / (tp + fn); rec[is.na(rec)] <- 0
  fs <- (1 + beta^2) * prec * rec / (beta^2 * prec + rec); fs[!is.finite(fs)] <- 0
  list(x = (br[-1] + br[-length(br)]) / 2, f = fs) }
# Cumulative-count F-beta: true = stim minus scaled unstim above g; false = scaled unstim above g; missed = net above lower end minus true.
fbetaCum <- function(xs, xu, lo, hi, beta = 0.8, tol = 0.01) {
  r <- length(xs) / length(xu); cand <- sort(unique(c(hi, xs[xs >= lo & xs < hi])), decreasing = TRUE)
  netLo <- sum(xs > lo) - sum(xu > lo) * r
  if (netLo <= 0) return(hi)
  tp <- sapply(cand, function(g) sum(xs > g) - sum(xu > g) * r); fp <- sapply(cand, function(g) sum(xu > g) * r)
  prec <- ifelse(tp + fp > 0, pmax(tp, 0) / (pmax(tp, 0) + fp), 0); rec <- pmax(tp, 0) / netLo
  f <- (1 + beta^2) * prec * rec / (beta^2 * prec + rec); f[!is.finite(f)] <- 0
  cand[which(f >= max(f) - tol)[1]] }   # highest gate within tolerance of the best
out <- list()
for (fl in list.files(P, pattern = "^rr_off_.*\\.rds$", full.names = TRUE)) {
  x <- readRDS(fl); row <- x$row; res <- x$result |> filter(method == "stimgate")
  tb <- bind_rows(lapply(x$tubes, function(t) { s <- t$sim
    xs <- s$stim; xu <- s$uns; own <- t$gate
    tibble(ind = t$ind, own = own, direct = t$locGeneratedDirect, lowerX = t$lowerX, bw = t$bw,
      responder = isTRUE(t$locGeneratedDirect) && is.finite(own) && (mean(xs > own) - mean(xu > own)) > 0, sim = list(s)) }))
  tb <- tb |> left_join(res |> transmute(ind = as.character(ind), sample, fin = threshold), by = "ind")
  m <- function(d) d |> filter(transformation == row$transformation, mean_pos_setting == row$mean_pos_setting, abs(prob_response - row$prob_response) < 1e-12, n_cell == row$n_cell) |> transmute(sample, g = threshold)
  tb <- tb |> left_join(m(on7) |> rename(stepOn = g), by = "sample") |> left_join(m(off7) |> rename(check = g), by = "sample")
  nResp <- sum(tb$responder & is.finite(tb$fin))
  q25 <- unname(quantile(tb$fin[tb$responder], 0.25, na.rm = TRUE))
  q25min5 <- if (nResp >= 5) q25 else NA_real_
  for (i in seq_len(nrow(tb))) { s <- tb$sim[[i]]; xs <- s$stim; xu <- s$uns; l <- s$labelsStim; g <- tb$fin[i]
    truth <- mean(l == "gp") - mean(s$labelsUns == "gp")
    e <- if (is.finite(tb$lowerX[i]) && is.finite(g)) .getCpUnsLocPlateauEdge(xs, g, tb$lowerX[i]) else NA
    floorX <- max(q25, e, na.rm = TRUE); if (!is.finite(floorX)) floorX <- g
    gFloorAll <- if (is.finite(g) && g > floorX) floorX else g
    gFloorResp <- if (tb$responder[i]) gFloorAll else g
    # move to the 25th-percentile responder gate only when the edge lies below it
    gQ25All <- if (is.finite(g) && is.finite(q25) && is.finite(e) && e < q25 && g > q25) q25 else g
    gQ25Resp <- if (tb$responder[i]) gQ25All else g
    # the same, only when the cluster has at least 5 responders, for responders only
    gQ25R5 <- if (tb$responder[i] && is.finite(g) && is.finite(q25min5) && is.finite(e) && e < q25min5 && g > q25min5) q25min5 else g
    # combined with the current step: the lower of the two gates
    s1 <- tb$stepOn[i]
    gComboAll <- if (is.finite(s1)) min(s1, gQ25All) else gQ25All
    gComboR5 <- if (is.finite(s1)) min(s1, gQ25R5) else gQ25R5
    lo <- if (is.finite(e)) e else g
    fc <- fbetaCurve(xu, xs); inI <- fc$x >= lo & fc$x <= g
    gFbBin <- if (any(inI)) { xb <- fc$x[inI][which.max(fc$f[inI])]; min(g, xb) } else g
    gFbCum <- if (is.finite(e)) fbetaCum(xs, xu, e, g) else g
    gates <- c(off = g, stepOn = tb$stepOn[i], floorAll = gFloorAll, floorResp = gFloorResp, q25All = gQ25All, q25Resp = gQ25Resp, q25R5 = gQ25R5, stepPlusQ25 = gComboAll, stepPlusQ25R5 = gComboR5, fbBin = gFbBin, fbCum = gFbCum)
    rowOut <- tibble(nResp = nResp, scenario = sprintf("%s %s %g%% %gk", row$transformation, row$mean_pos_setting, 100 * row$prob_response, row$n_cell / 1000), sample = tb$sample[i],
      responder = tb$responder[i], checkOk = isTRUE(abs(tb$check[i] - g) < 1e-9))
    for (nm in names(gates)) { gg <- gates[[nm]]; rowOut[[paste0("f1_", nm)]] <- if (is.finite(gg)) f1(xs, l, gg) else NA; rowOut[[paste0("e_", nm)]] <- if (is.finite(gg)) (mean(xs > gg) - mean(xu > gg)) / truth - 1 else NA }
    out[[length(out) + 1]] <- rowOut } }
r <- bind_rows(out); saveRDS(r, file.path(P, "rules.rds"))
cat("rerun matches the pre-step final gates (iter 1):", sum(r$checkOk), "of", nrow(r), "\n\n")
v <- c("off", "stepOn", "q25Resp", "q25R5", "stepPlusQ25", "stepPlusQ25R5", "fbCum")
cat("Median F1\n"); print(as.data.frame(r |> group_by(scenario) |> summarise(n = n(), across(all_of(paste0("f1_", v)), ~ median(.x, na.rm = TRUE))) |> rename_with(~ sub("f1_", "", .x))), digits = 3)
cat("\nTubes worse / better than off (by > 0.02)\n"); print(as.data.frame(r |> group_by(scenario) |> summarise(across(all_of(paste0("f1_", v[-1])), ~ paste0(sum(.x < f1_off - 0.02, na.rm = TRUE), "/", sum(.x > f1_off + 0.02, na.rm = TRUE)))) |> rename_with(~ sub("f1_", "", .x))))
q <- function(e) sprintf("%+.0f [%+.0f,%+.0f]", 100 * median(e, na.rm = TRUE), 100 * quantile(e, .1, na.rm = TRUE), 100 * quantile(e, .9, na.rm = TRUE))
cat("\nFrequency error %, median [10th, 90th]\n"); print(as.data.frame(r |> group_by(scenario) |> summarise(across(all_of(paste0("e_", v)), q)) |> rename_with(~ sub("e_", "", .x))))
