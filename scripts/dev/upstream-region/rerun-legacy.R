Sys.setenv(RUN_SIMULATIONS="false", RUN_PLOTS="false")
P <- Sys.getenv("P")
suppressMessages(source(file.path(P, "q7setup.R")))
source(file.path(scripts_r_dir, "sim-debug-loc.R"))
devtools::load_all(quiet = TRUE)
ns <- asNamespace("stimgate")

# --- Proposed quality rule (QUAL_NEW=1) -----------------------------------
# Floor = min(xDom, own Tailgate-style point on raw unstim + bw/2).
# Crowding: bin cells <= CROWD_FRAC * mean of 3 largest non-zero bins
#   (CROWD=global: over all stim bins; CROWD=ref: over bins right of start).
if (identical(Sys.getenv("QUAL_NEW"), "1")) {
  tgTol <- as.numeric(Sys.getenv("TG_TOL", as.character(10^-2.5)))
  crowd <- Sys.getenv("CROWD", "ref")
  crowdFrac <- as.numeric(Sys.getenv("CROWD_FRAC", as.character(2/3)))
  .pair <- identical(Sys.getenv("PAIR"), "1")
  sys.source(file.path(P, "mb.R"), envir = (mbEnv <- new.env(parent = ns)))
  mbMod <- mbEnv$mbMod
  environment(mbMod) <- list2env(list(.crowd = crowd, .crowdFrac = crowdFrac, .pair = .pair, .partial = identical(Sys.getenv("PARTIAL"), "1"), .multires = identical(Sys.getenv("MULTIRES"), "1"), .purSrc = Sys.getenv("PUR_SRC", "fit"), .purExp = as.numeric(Sys.getenv("PUR_EXP", "0"))), parent = ns)

  tgPoint <- function(x, bw, tol) {
    x <- x[is.finite(x)]
    d <- stats::density(x, bw = bw, n = 2048)
    dy <- diff(d$y) / diff(d$x); dx <- (d$x[-1] + d$x[-length(d$x)]) / 2
    pk <- d$x[get(".getPeakMainLeftIdx", ns)(d$y)]
    right <- dx > pk
    steep <- dx[right][which.min(dy[right])]          # steepest descent after the peak
    cand <- dx > steep & abs(dy) < tol * max(abs(dy))
    if (!any(cand)) return(NA_real_)
    dx[which(cand)[1]]
  }

  qOrig <- get(".getCpUnsLocQualityBoundaryCurrent", ns)
  qNew <- function(dataMod, chnlSettings, probCol, xClear, lowerBoundX = NA_real_,
                   exTblStimOrig = NULL, exTblUnsOrig = NULL) {
    pf <- parent.frame()
    xDom <- tryCatch(get("xDom", pf), error = function(e) NA_real_)
    bw <- suppressWarnings(as.numeric(attr(get("dataModFull", pf), "locDensityBw")))
    bw <- stats::median(bw[is.finite(bw)])
    if (!is.finite(bw)) bw <- as.numeric(chnlSettings$bwShared)[1]
    xF <- get(".getCut", ns)(if (identical(Sys.getenv("FLOOR_TUBE"), "stim")) exTblStimOrig else exTblUnsOrig)
    xTG <- tgPoint(xF, bw, tgTol) + bw / 2
    floorX <- min(c(xDom, xTG), na.rm = TRUE)
    lb <- max(c(lowerBoundX, floorX), na.rm = TRUE)
    out <- mbMod(dataMod = dataMod, chnlSettings = chnlSettings, probCol = probCol,
                 startX = xClear, lowerBoundX = if (.mr) lowerBoundX else lb, floorX = lb,
                 exTblStimOrig = exTblStimOrig, exTblUnsOrig = exTblUnsOrig)
    xQual <- suppressWarnings(as.numeric(out$info$finalStartX)[1L])
    if (!is.finite(xQual)) xQual <- xClear
    out$info$thresholdX <- xClear; out$info$xQual <- xQual
    out$info$xTG <- xTG; out$info$floorX <- floorX; out$info$bwTG <- bw
    out$info$preliminaryLowerBoundX <- lowerBoundX
    list(thresholdX = xQual, dataMod = out$dataMod, info = out$info)
  }
  environment(qNew) <- list2env(list(mbMod = mbMod, tgPoint = tgPoint, tgTol = tgTol, ns = ns, .mr = identical(Sys.getenv("MULTIRES"), "1")), parent = ns)
  unlockBinding(".getCpUnsLocQualityBoundaryCurrent", ns)
  assign(".getCpUnsLocQualityBoundaryCurrent", qNew, ns)
  cat("new quality rule: tol", tgTol, "crowd", crowd, crowdFrac, "pair", .pair, "partial", Sys.getenv("PARTIAL"), "\n")
}

if (identical(Sys.getenv("DOM_ONSET_HALFBW"), "1")) {
  origD <- get(".getCpUnsLocDominanceBoundaryCurrent", ns)
  wrapD <- function(density, startX, densityBw = NULL, dominanceRatio = 2, onsetWeight = 2/3, lowerBoundX = NA_real_) {
    out <- origD(density, startX, densityBw, dominanceRatio, onsetWeight, lowerBoundX)
    on <- out$info$onsetX; bwD <- suppressWarnings(as.numeric(unlist(densityBw))); bwD <- stats::median(bwD[is.finite(bwD)])
    if (is.finite(out$startX) && length(on) == 1L && is.finite(on) && is.finite(bwD)) { out$info$origDominanceX <- out$startX; out$startX <- on + bwD / 2 }
    out
  }
  environment(wrapD) <- list2env(list(origD = origD), parent = ns)
  unlockBinding(".getCpUnsLocDominanceBoundaryCurrent", ns); assign(".getCpUnsLocDominanceBoundaryCurrent", wrapD, ns)
  cat("dominance: onset + bw/2\n")
}

if (identical(Sys.getenv("PLATEAU"), "1")) {
  # Replace the antimode step: step down from xBase in bandwidth-sized steps towards the taut-string plateau edge.
  wrapA <- function(dataMod, chnlSettings, xBase, lowerBoundX = NA_real_, heightFrac = 0.95) {
    pf <- parent.frame()
    info <- list(applied = FALSE, reason = "plateau_not_applied", xBase = xBase, lowerBoundX = lowerBoundX)
    xs <- .getCut(get("exTblStimOrig", pf)); xu <- .getCut(get("exTblUnsOrig", pf)); xs <- xs[is.finite(xs)]; xu <- xu[is.finite(xu)]
    bw <- suppressWarnings(as.numeric(unlist(attr(get("dataModFull", pf), "locDensityBw")))); bw <- stats::median(bw[is.finite(bw)])
    info$bw <- bw
    if (!is.finite(xBase) || !is.finite(lowerBoundX) || !is.finite(bw) || length(xu) == 0L) return(list(thresholdX = NA_real_, info = info))
    sub <- sort(xs[xs >= lowerBoundX]); den <- .getCpUnsLocAntimodeDensity(sub)
    if (is.null(den) || all(den$y == 0)) { info$reason <- "no_taut_string"; return(list(thresholdX = NA_real_, info = info)) }
    r <- rle(den$y); ends <- cumsum(r$lengths); st <- ends - r$lengths + 1; j <- findInterval(xBase, den$x)
    k <- which(st <= j & ends >= j)[1]
    fl <- if (is.na(k) || st[k] <= 1) NA_real_ else sub[st[k]]
    info$plateauEdge <- fl
    if (!is.finite(fl) || fl >= xBase) { info$reason <- "no_plateau_edge_below_xBase"; return(list(thresholdX = NA_real_, info = info)) }
    ratio <- length(xs) / length(xu); nSr <- sum(xs > xBase); refPur <- (nSr - sum(xu > xBase) * ratio) / nSr
    info$refPur <- refPur
    if (!is.finite(refPur) || refPur <= 0) { info$reason <- "no_reference_purity"; return(list(thresholdX = NA_real_, info = info)) }
    cur <- xBase; hi <- xBase; rej <- 0; steps <- list()
    while (rej < .maxRej) { lo <- max(hi - .k * bw, fl); if (lo >= hi - 1e-9) break
      nS <- sum(xs > lo & xs <= cur); nU <- sum(xu > lo & xu <= cur) * ratio
      ok <- nS > 0 && (nS - nU) / nS >= .bar * refPur
      steps[[length(steps) + 1]] <- c(lo = lo, top = cur, nS = nS, nU = nU, ok = ok)
      if (ok) { cur <- lo; rej <- 0 } else rej <- rej + 1
      hi <- lo; if (lo <= fl) break }
    info$steps <- do.call(rbind, steps); info$xPlateau <- cur
    if (cur >= xBase) { info$reason <- "no_step_accepted"; return(list(thresholdX = NA_real_, info = info)) }
    info$applied <- TRUE; info$reason <- "plateau_step_applied"
    list(thresholdX = cur, info = info)
  }
  environment(wrapA) <- list2env(list(.k = as.numeric(Sys.getenv("PL_K", "1")), .bar = as.numeric(Sys.getenv("PL_BAR", "0.5")), .maxRej = as.numeric(Sys.getenv("PL_REJ", "2"))), parent = ns)
  unlockBinding(".getCpUnsLocAntimodeBoundaryCurrent", ns); assign(".getCpUnsLocAntimodeBoundaryCurrent", wrapA, ns)
  cat("plateau step in antimode slot\n")
}

args <- commandArgs(TRUE); pr <- as.numeric(args[1]); nc <- as.numeric(args[2]); tag <- args[3]
ctrl_over <- if (length(args) > 3) eval(parse(text = args[4])) else list()
row <- sim_grid_full |> dplyr::filter(transformation==Sys.getenv("TRANS","gaussian"), mean_pos_setting==Sys.getenv("SEP","low"), abs(prob_response-pr)<1e-9, n_cell==nc)
stopifnot(nrow(row)==1)
if (length(args) > 4) { ro <- eval(parse(text = args[5])); for (k in names(ro)) row[[k]] <- ro[[k]] }
cat("sim_id", row$sim_id, "\n")
.simCompareFbetaEnvironment(pathFbeta = path_fbeta)
call_args <- list(row = row, nSample = n_sample_sim, nIter = 1L, nMarker = 1, nCondition = 2, nCluster = 2,
  probExact = TRUE, covEvMin = 1.5, covEvMax = 1.5, clusterGates = cluster_gates,
  locEnforceShapeThreshold = loc_enforce_shape_threshold, locThresholdMethod = stimgate_loc_threshold_method,
  calcCytPosGates = calc_cyt_pos_gates, pathFbeta = path_fbeta, fbetaBeta = fbeta_beta, fbetaTheta = fbeta_theta,
  fbetaWidth = fbeta_width, tailgateAdjust = tailgate_adjust, tailgateAutoTol = tailgate_auto_tol, tailgateBias = tailgate_bias,
  tailgateMethod = tailgate_method, tailgateTol = tailgate_tol, tailgateX = tailgate_x, resume = FALSE)
call_args[names(ctrl_over)] <- ctrl_over
rec <- function(d) {
  s <- .simDebugLocSummary(d)
  ln <- .simDebugLocLines(d)
  s$lines <- list(stats::setNames(ln$x, as.character(ln$line)))
  s$sim <- list(d$sim)
  s$filter <- list(list(info = list(final = d$filter$info$final, marginal = d$filter$info$marginal)))
  s
}
out <- .simDebugLoc(do.call(.simCompareRunScenario, call_args), sample = NULL, ind = NULL, stopAfter = FALSE, onRecord = rec)
saveRDS(out, file.path(P, paste0("dbg_", tag, ".rds")))
res <- attr(out, "result")
saveRDS(res, file.path(P, paste0("res_", tag, ".rds")))
cat("done\n")
