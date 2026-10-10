# StimGate, Tailgate and F-beta on the prepared OMIP-016 CD4 T cells, and
# scoring against the manual gates re-applied from the authors' workspace.
#
# Requires scripts/r/omip016-prepare.R, scripts/r/acs_cytof-helper.R,
# scripts/r/acs_cytof-gate.R (.acsCytofSetDebug()), scripts/r/acs_cytof-methods.R
# (.acsCytofThresholdOne()), scripts/r/sim-compare-freq_bs.R (the
# Tailgate/F-beta wrappers), scripts/r/sim-debug-loc.R (.simDebugLoc()) and
# scripts/r/acs_cytof-debug.R (.acsCytofDebugSlim()).

.omip016Methods <- function() c("stimgate", "tailgate", "fbeta")

# OMIP-016 is conventional flow data without the exact-zero spike that made
# the ACS CyTOF comparators need tuning, so both comparators keep their
# published defaults (Analysis 9's "_default" settings). Tuning them against
# these manual gates would make the comparison circular.
.omip016ComparatorSettings <- function(method = c("tailgate", "fbeta")) {
  method <- match.arg(method)
  params <- switch(
    method,
    tailgate = list(
      tailgateX = "stim", adjust = 1, bandwidth = NULL, numPeaks = 1L,
      refPeak = 1L, derivativeMethod = "firstDeriv", tol = 1e-2,
      side = "right", strict = FALSE, autoTol = TRUE, bias = 0,
      removeZero = FALSE
    ),
    fbeta = list(
      beta = 0.8, theta = 2, width = 10L, numBins = NULL,
      patchPy2Compat = TRUE, removeZero = FALSE
    )
  )
  list(cacheVersion = 2L, method = method, family = method, params = params)
}

.omip016StimGateSettings <- function(locThresholdMethod = "cap",
                                     clusterGates = TRUE,
                                     calcCytPosGates = TRUE) {
  list(
    biasUns = NULL, biasUnsFactor = 1, bwMtd = "nrd0", bwScope = "cytokine",
    bwNcellMax = 1e4, bwFallback = "auto", bwMin = "none", bwMax = "none",
    gateCombn = "min", clusterGates = clusterGates,
    calcCytPosGates = calcCytPosGates, minCell = 100,
    locThresholdMethod = locThresholdMethod,
    gateName = if (isTRUE(clusterGates)) "loc_minClust" else "loc"
  )
}

.omip016ReadPrepared <- function(pathOut) {
  paths <- .omip016Paths(pathOut)
  if (!dir.exists(paths$gs)) {
    stop("Prepared OMIP-016 GatingSet not found at ", paths$gs, "; run the preparation stage.")
  }
  pre <- readRDS(.omip016PreprocessingFile(paths$gs))
  gs <- flowWorkspace::load_gs(paths$gs, backend_readonly = TRUE)
  if (!identical(flowWorkspace::sampleNames(gs), pre$sampleMap$file)) {
    stop("OMIP-016 GatingSet samples do not match its preprocessing manifest.")
  }
  labels <- lapply(pre$sampleMap$file, function(file) {
    readRDS(file.path(paths$labels, paste0(tools::file_path_sans_ext(file), "-cd4-manual.rds")))
  })
  names(labels) <- pre$sampleMap$file
  list(gs = gs, pre = pre, labels = labels, paths = paths)
}

.omip016Expr <- function(gs, ind, channels) {
  ex <- flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[ind]], "root"))
  ex[, channels, drop = FALSE]
}

# The run is wrapped in .simDebugLoc(), which records every tube's initial
# local-FDR gating without changing the gates, so a gate can be explained
# afterwards from the same run (slimmed records in omip016-debug-records.rds).
.omip016RunStimGate <- function(prep, pathProject, settings) {
  restoreDebug <- .acsCytofSetDebug()
  on.exit(restoreDebug(), add = TRUE)
  sm <- prep$pre$sampleMap
  records <- .simDebugLoc(
    .omip016GateStim(prep, pathProject, settings),
    sample = NULL,
    tubeInfo = function(ind) {
      list(sample = 1L, label = sm$stim[[match(ind, sm$ind)]])
    },
    onRecord = .acsCytofDebugSlim
  )
  saveRDS(unclass(records), file.path(pathProject, "omip016-debug-records.rds"))
  saveRDS(
    list(settings = settings, preprocessing = prep$pre[c("version", "inputContentHash", "settings")]),
    file.path(pathProject, "omip016-manifest.rds")
  )
  invisible(pathProject)
}

.omip016GateStim <- function(prep, pathProject, settings) {
  stimgate::gateStim(
    pathProject = pathProject,
    .data = prep$gs,
    popGate = "root",
    batchList = prep$pre$batchList,
    chnl = names(prep$pre$responseChannels),
    biasUns = settings$biasUns,
    control = stimgate::stimControl(
      biasUnsFactor = settings$biasUnsFactor,
      bwMtd = settings$bwMtd,
      bwScope = settings$bwScope,
      bwNcellMax = settings$bwNcellMax,
      bwFallback = settings$bwFallback,
      bwMin = settings$bwMin,
      bwMax = settings$bwMax,
      gateCombn = settings$gateCombn,
      clusterGates = settings$clusterGates,
      calcCytPosGates = settings$calcCytPosGates,
      minCell = settings$minCell,
      locThresholdMethod = settings$locThresholdMethod
    )
  )
}

# The recorded initial local-FDR gating of one stimulated tube and channel,
# drawn as .simDebugLocPlots() panels with the final StimGate gates and the
# manual threshold added. Returns the plot grid and a one-row summary.
.omip016DebugGate <- function(pathProject, ind, chnl, scored) {
  records <- readRDS(file.path(pathProject, "omip016-debug-records.rds"))
  nm <- paste0("dataset1_ind", ind, "_", chnl)
  rec <- records[[nm]]
  if (is.null(rec) || !isTRUE(rec$found)) {
    stop("No recorded local-FDR gating for tube ", ind, ", channel ", chnl, ".")
  }
  rows <- scored[as.character(scored$ind) == as.character(ind) & scored$chnl == chnl, ]
  sg <- rows[rows$method == "stimgate", ]
  man <- rows[rows$method == "manual", ]
  extra <- data.frame(
    # The line key appends each value.
    line = c("final gate", "final cyt+ gate", "manual"),
    x = c(sg$gate, sg$gateCyt, man$gate),
    colour = c("#D55E00", "#D55E00", "#CC0000"),
    linetype = c("solid", "dotted", "solid"),
    linewidth = c(1.1, 0.9, 0.9)
  )
  plots <- .simDebugLocPlots(rec, extraLines = extra)
  keep <- plots[c("density", "prob", "respCells")]
  attr(keep, "lines") <- attr(plots, "lines")
  summary <- .simDebugLocSummary(rec)
  list(
    plot = .simDebugLocPlotGrid(keep),
    summary = data.frame(
      ownGate = summary$threshold,
      ownReason = summary$locReason,
      ownFilterReason = summary$filterReason,
      ownPropRespEst = summary$propRespEst,
      smoothing = attr(rec$dataMod, "locProbSmoothMethod"),
      finalGate = sg$gate, finalGateCyt = sg$gateCyt,
      finalReason = sg$thresholdReason, manualGate = man$gate,
      stringsAsFactors = FALSE
    )
  )
}

# StimGate's final gates as one row per stimulated tube and channel.
.omip016StimGateGates <- function(prep, pathProject) {
  manifest <- readRDS(file.path(pathProject, "omip016-manifest.rds"))
  gateName <- manifest$settings$gateName
  gates <- stimgate::getStimGates(pathProject)
  gates <- gates[gates$gateName == gateName, , drop = FALSE]
  stimInd <- prep$pre$batchList[[1]][-1]
  channels <- names(prep$pre$responseChannels)
  out <- expand.grid(ind = stimInd, chnl = channels, stringsAsFactors = FALSE)
  k <- match(paste(out$ind, out$chnl), paste(gates$ind, gates$chnl))
  out$method <- "stimgate"
  out$gate <- unname(gates$gate[k])
  out$gateCyt <- unname(gates$gateCyt[k])
  out$thresholdOrigin <- as.character(gates$locSource[k])
  out$thresholdReason <- as.character(gates$locReason[k])
  out$thresholdFallbackUsed <- !(gates$locGenerated[k] %in% TRUE)
  out
}

.omip016RunComparator <- function(prep, method, pathFbeta = NULL) {
  settings <- .omip016ComparatorSettings(method)
  fbetaEnv <- if (method == "fbeta") {
    .simCompareFbetaEnvironment(pathFbeta = pathFbeta, patchPy2Compat = TRUE)
  }
  batch <- prep$pre$batchList[[1]]
  channels <- names(prep$pre$responseChannels)
  xUns <- .omip016Expr(prep$gs, batch[[1]], channels)
  rows <- lapply(batch[-1], function(ind) {
    xStim <- .omip016Expr(prep$gs, ind, channels)
    do.call(rbind, lapply(channels, function(chnl) {
      res <- .acsCytofThresholdOne(
        method = method, xUns = xUns[, chnl], xStim = xStim[, chnl],
        settings = settings, pathFbeta = pathFbeta, fbetaEnv = fbetaEnv
      )
      data.frame(
        ind = ind, chnl = chnl, method = method, gate = res$threshold,
        gateCyt = NA_real_,
        thresholdOrigin = res$thresholdOrigin, thresholdReason = NA_character_,
        thresholdFallbackUsed = res$thresholdFallbackUsed,
        stringsAsFactors = FALSE
      )
    }))
  })
  list(settings = settings, gates = do.call(rbind, rows))
}

# Manual gates as rows of the same shape, on the GatingSet scale.
.omip016ManualGates <- function(prep, thresholds) {
  stimInd <- prep$pre$batchList[[1]][-1]
  out <- expand.grid(ind = stimInd, chnl = thresholds$channel, stringsAsFactors = FALSE)
  out$method <- "manual"
  out$gate <- thresholds$threshold_trans[match(out$chnl, thresholds$channel)]
  out$gateCyt <- NA_real_
  out$thresholdOrigin <- "flowjo_workspace"
  out$thresholdReason <- NA_character_
  out$thresholdFallbackUsed <- FALSE
  out
}

# Positivity of every response marker for one tube and one method's gates,
# as the package classifies cells: x > gate (strict), or x > gateCyt for a
# cell that is ordinarily positive for another marker (StimGate's
# cytokine-positive gates; NA gateCyt means no such gate). A non-finite gate
# gives NA calls for that marker.
.omip016Classify <- function(x, gate, gateCyt) {
  base <- sweep(x, 2L, gate, FUN = ">")
  base[, !is.finite(gate)] <- NA
  out <- base
  for (j in which(is.finite(gateCyt) & is.finite(gate))) {
    other <- rowSums(base[, -j, drop = FALSE], na.rm = TRUE) > 0
    out[, j] <- base[, j] | (x[, j] > gateCyt[[j]] & other)
  }
  out
}

# Frequencies and per-cell agreement with the manual labels for each gate.
# The unstimulated tube is classified with the stimulated tube's gates, and
# the manual labels are the re-applied workspace gates (not the 1D manual
# threshold).
.omip016Score <- function(prep, gates) {
  channels <- names(prep$pre$responseChannels)
  markers <- unname(prep$pre$responseChannels)
  sm <- prep$pre$sampleMap
  indUns <- prep$pre$batchList[[1]][[1]]
  xUns <- .omip016Expr(prep$gs, indUns, channels)
  labUns <- prep$labels[[sm$file[[indUns]]]]
  if (!"gateCyt" %in% names(gates)) gates$gateCyt <- NA_real_
  groups <- split(seq_len(nrow(gates)), paste(gates$method, gates$ind))
  rows <- lapply(groups, function(k) {
    g <- gates[k, , drop = FALSE]
    if (anyDuplicated(g$chnl) || !setequal(g$chnl, channels)) {
      stop("Expected one gate per response channel for ", g$method[[1]], " tube ", g$ind[[1]], ".")
    }
    g <- g[match(channels, g$chnl), , drop = FALSE]
    ind <- g$ind[[1]]
    x <- .omip016Expr(prep$gs, ind, channels)
    lab <- prep$labels[[sm$file[[ind]]]]
    if (nrow(lab) != nrow(x)) stop("Manual labels do not match the GatingSet cells.")
    pos <- .omip016Classify(x, g$gate, g$gateCyt)
    posUns <- .omip016Classify(xUns, g$gate, g$gateCyt)
    labStim <- as.matrix(lab[, markers])
    labU <- as.matrix(labUns[, markers])
    data.frame(
      g,
      stim = sm$stim[[ind]], cyt = markers,
      nCellStim = nrow(x), countStim = colSums(pos),
      nCellUns = nrow(xUns), countUns = colSums(posUns),
      tp = colSums(pos & labStim), fp = colSums(pos & !labStim),
      fn = colSums(!pos & labStim),
      nManualStim = colSums(labStim), nManualUns = colSums(labU),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  out <- do.call(rbind, rows)
  out$freq_stim <- 100 * out$countStim / out$nCellStim
  out$freq_uns <- 100 * out$countUns / out$nCellUns
  out$freq_bs <- pmax(out$freq_stim - out$freq_uns, 0)
  out$freq_stim_man <- 100 * out$nManualStim / out$nCellStim
  out$freq_uns_man <- 100 * out$nManualUns / out$nCellUns
  out$freq_bs_man <- pmax(out$freq_stim_man - out$freq_uns_man, 0)
  # Undefined proportions stay NA.
  out$fdp <- ifelse(out$tp + out$fp > 0, out$fp / (out$tp + out$fp), NA_real_)
  out$sensitivity <- ifelse(out$tp + out$fn > 0, out$tp / (out$tp + out$fn), NA_real_)
  out$f1 <- ifelse(2 * out$tp + out$fp + out$fn > 0,
                   2 * out$tp / (2 * out$tp + out$fp + out$fn), NA_real_)
  out
}

# StimGate's own single-marker counts must equal the counts from
# .omip016Classify() used in the score table.
.omip016CheckStimGateCounts <- function(prep, pathProject, scored) {
  stats <- stimgate::getStimStats(pathProject)
  manifest <- readRDS(file.path(pathProject, "omip016-manifest.rds"))
  stats <- stats[stats$gateName == manifest$settings$gateName, , drop = FALSE]
  sg <- scored[scored$method == "stimgate", , drop = FALSE]
  ok <- vapply(seq_len(nrow(sg)), function(i) {
    # Combination labels concatenate "<channel>~+~" / "<channel>~-~" tokens.
    token <- paste0(sg$chnl[[i]], "~+~")
    tokens <- regmatches(stats$cytCombn, gregexpr("[^~]+~[+-]~", stats$cytCombn))
    hasPos <- vapply(tokens, function(x) token %in% x, logical(1))
    st <- stats[as.character(stats$ind) == as.character(sg$ind[[i]]) & hasPos, ]
    isTRUE(sum(st$countStim) == sg$countStim[[i]] &&
      sum(st$countUns) == sg$countUns[[i]])
  }, logical(1))
  if (!all(ok)) {
    stop("StimGate statistics do not reproduce the strict x > gate counts.")
  }
  invisible(TRUE)
}
