# Rerun one Analysis 7 scenario (all 20 tubes of one dataset) and record each tube's gating.
# Args: prob_response n_cell tag [call overrides as R list]; env TRANS, SEP.
P <- Sys.getenv("P")
suppressMessages(source(file.path(P, "q7setup_wt.R")))
source(file.path(scripts_r_dir, "sim-debug-loc.R"))
args <- commandArgs(TRUE); pr <- as.numeric(args[1]); nc <- as.numeric(args[2]); tag <- args[3]
over <- if (length(args) > 3) eval(parse(text = args[4])) else list()
row <- sim_grid_full |> dplyr::filter(transformation == Sys.getenv("TRANS", "gaussian"), mean_pos_setting == Sys.getenv("SEP", "low"), abs(prob_response - pr) < 1e-12, n_cell == nc)
stopifnot(nrow(row) == 1)
.simCompareFbetaEnvironment(pathFbeta = path_fbeta)
callArgs <- list(row = row, nSample = n_sample_sim, nIter = 1L, nMarker = 1, nCondition = 2, nCluster = 2, probExact = TRUE, covEvMin = 1.5, covEvMax = 1.5,
  clusterGates = cluster_gates, locEnforceShapeThreshold = loc_enforce_shape_threshold, locPlateauStep = loc_plateau_step, locThresholdMethod = stimgate_loc_threshold_method,
  calcCytPosGates = calc_cyt_pos_gates, pathFbeta = path_fbeta, fbetaBeta = fbeta_beta, fbetaTheta = fbeta_theta, fbetaWidth = fbeta_width,
  tailgateAdjust = tailgate_adjust, tailgateAutoTol = tailgate_auto_tol, tailgateBias = tailgate_bias, tailgateMethod = tailgate_method, tailgateTol = tailgate_tol, tailgateX = tailgate_x, resume = FALSE)
callArgs[names(over)] <- over
rec <- function(d) {
  bw <- suppressWarnings(as.numeric(unlist(attr(d$dataMod, "locDensityBw")))); bw <- stats::median(bw[is.finite(bw)])
  list(ind = as.character(d$ind), gate = d$cp$cp %||% NA_real_, locGeneratedDirect = isTRUE(d$cp$locGeneratedDirect), locReason = d$cp$locReason %||% NA_character_,
       lowerX = suppressWarnings(as.numeric(attr(d$dataMod, "minProbXPos"))[1]), bw = bw, sim = d$sim)
}
out <- .simDebugLoc(do.call(.simCompareRunScenario, callArgs), sample = NULL, ind = NULL, stopAfter = FALSE, onRecord = rec)
saveRDS(list(tubes = unclass(out)[names(out)], result = attr(out, "result"), row = row, over = over), file.path(P, paste0("rr_", tag, ".rds")))
cat("done\n")
