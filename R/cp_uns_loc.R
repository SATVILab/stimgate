# Local-FDR gating orchestration
#
# Entry points, sample preparation, gate-combination handling, and condition-
# level control flow. Density estimation, smoothing, filtering, threshold
# selection, and output assembly are kept in separate files.

# Calculate local FDR-based cutpoint for each level of bias
#' @keywords internal
.getCpUnsLoc <- function(
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  purrr::map(chnlSettings$biasUns, function(bias) {
    .debug("biasUns", bias) # nolint

    exListPrep <- .prepareDataWithBiasAndNoise(
      exList = exList,
      bias = bias,
      excMin = chnlSettings$excMin
    )

    # get gates for given level of bias across gate combination methods
    list("loc" = .getCpUnsLocGateCombn(
      exListOrig = exListPrep[["exListOrig"]],
      exListNoMin = exListPrep[["exListNoMin"]],
      exTblUnsBias = exListPrep[["exTblUnsBias"]],
      chnlSettings = chnlSettings,
      bias = bias,
      pathProject = pathProject,
      stage = stage
    ))
  }) |>
    purrr::flatten()
}

#' @keywords internal
.prepareDataWithBiasAndNoise <- function(
  exList,
  bias,
  excMin
) {
  # Keep the complete original data.
  exListOrig <- .prepareExListWithBiasAndNoise(
    exList = exList,
    ind = names(exList),
    excMin = FALSE,
    bias = 0
  ) |>
    .arrangeSamplesByExpr()

  # Remove minimum-expression values where requested. Row subsetting keeps the
  # ascending order of exListOrig.
  if (isTRUE(excMin)) {
    exListNoMin <- .prepareExListWithBiasAndNoise(
      exList = exListOrig,
      ind = names(exListOrig),
      excMin = TRUE,
      bias = 0
    )
  } else {
    exListNoMin <- exListOrig
  }

  # The biased unstimulated data differ from the already-prepared unstimulated
  # data only by a constant expression shift. Adding a constant preserves its
  # existing ordering.
  exTblUnsBias <- exListNoMin[[1L]]

  chnl <- attr(exTblUnsBias, "chnlCut")
  exTblUnsBias[[chnl]] <- .getCut(exTblUnsBias) + bias

  list(
    "exListOrig" = exListOrig,
    "exListNoMin" = exListNoMin,
    "exTblUnsBias" = exTblUnsBias
  )
}

# Get the unstim-based local fdr-method
# cutpoint for a given bias across gate combination methods
#' @keywords internal
.getCpUnsLocGateCombn <- function(
  exListOrig,
  exListNoMin,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage
) {
  .debug("getting gateCombn") # nolint
  cp <- list()

  # gate using prejoined stim data
  if ("prejoin" %in% chnlSettings$gateCombn) {
    .debug("prejoin") # nolint
    # join marker expression for stim samples and sort into ascending order
    cp[["prejoin"]] <- .getCpUnsLocSample(
      exListOrig = .prepareDataForPrejoinInd(exListOrig),
      exListNoMinStim = .prepareDataForPrejoinInd(exListNoMin)[-1],
      exTblUnsBias = exTblUnsBias,
      chnlSettings = chnlSettings,
      bias = bias,
      stage = stage,
      pathProject = pathProject,
      indStim = names(exListNoMin)[-1],
      exListOrigOutput = exListOrig,
      prejoin = TRUE
    )[["loc"]]
  }

  # gate each sample individually, then combine
  nonPrejoinCombn <- setdiff(chnlSettings$gateCombn, "prejoin")
  if (length(nonPrejoinCombn) > 0L) {
    .debug("non-prejoin") # nolint
    cpUnsListNonjoin <- .getCpUnsLocSample(
      exListOrig = exListOrig,
      exListNoMinStim = exListNoMin[-1],
      exTblUnsBias = exTblUnsBias,
      chnlSettings = chnlSettings,
      bias = bias,
      stage = stage,
      pathProject = pathProject
    )
    .debug("Combining cutpoints") # nolint
    cpCombn <- .getCpUnsLocCombineCpWithMeta(
      cp = cpUnsListNonjoin[["loc"]],
      gateCombn = nonPrejoinCombn,
      exListOrig = exListOrig,
      shareCap = chnlSettings$locShareCap %||% 1.5,
      cellCap = chnlSettings$locShareCellCap %||% 0.5
    )
    .intSaveNm(
      "locDetailBatchShare",
      .getCpUnsLocShareDetailTbl(cpCombn),
      .createCombinedIdentifier(names(exListNoMin)[-1]),
      file.path(stage, .getCpUnsLocGetChnl(exListOrig[[1L]])),
      pathProject
    )
    cp <- c(cp, cpCombn)
  }

  .debug("done getting gateCombn") # nolint
  list("cp" = cp)
}

#' @keywords internal
.prepareDataForPrejoinInd <- function(exList) {
  indUns <- names(exList)[1]
  indStim <- names(exList)[-1]
  exTblStim <- exList[seq.int(2, length(exList))] |>
    dplyr::bind_rows()
  exTblStim <- exTblStim[order(.getCut(exTblStim)), ]
  list(
    exList[[1]],
    exTblStim
  ) |>
    stats::setNames(c(indUns, .createCombinedIdentifier(indStim)))
}

# Sort each sample's rows into ascending order of the cut channel
#' @keywords internal
.arrangeSamplesByExpr <- function(exList) {
  purrr::map(exList, function(x) {
    cut <- .getCut(x)

    if (is.unsorted(cut)) {
      x[order(cut), , drop = FALSE]
    } else {
      x
    }
  })
}

# Share gates within a batch. Only responders (see .getCpShareResponder())
# donate: their gates are combined, and each tube accepts the combined gate
# within the limits of .getCpShareApply().
#' @keywords internal
.getCpUnsLocCombineCpWithMeta <- function(
  cp,
  gateCombn,
  exListOrig,
  shareCap = 1.5,
  cellCap = 0.5
) {
  meta <- .getCpUnsLocMetaFromCp(cp)
  cpNum <- suppressWarnings(as.numeric(cp))
  stimRow <- !(meta$locSource %in% "unstim_summary")
  xUns <- .getCut(exListOrig[[1L]])
  xStimList <- lapply(meta$ind, function(i) {
    if (i %in% names(exListOrig)[-1L]) .getCut(exListOrig[[i]])
  })
  freqOwn <- vapply(seq_along(cpNum), function(i) {
    if (is.null(xStimList[[i]])) {
      return(NA_real_)
    }
    .getCpShareFreq(cpNum[[i]], xStimList[[i]], xUns)
  }, numeric(1))
  # Frequency at each tube's own gate before sharing; the cluster step uses
  # it for its donors' median frequency.
  meta$locOwnFreq <- ifelse(stimRow, freqOwn, NA_real_)
  meta$locResponder <- stimRow & .getCpShareResponder(
    locGeneratedDirect = meta$locGeneratedDirect,
    cp = cpNum,
    freq = freqOwn
  )
  stimGenerated <- meta$locGenerated & stimRow

  purrr::map(gateCombn, function(gateCombnCurr) {
    if (is.null(gateCombnCurr) || gateCombnCurr %in% c("no", "prejoin")) {
      return(.getCpUnsLocCpAttachMeta(cp, meta))
    }

    # No generated gate: each tube keeps its own fallback. Combining them
    # could give a tube another tube's lower fallback, below its own cells.
    if (!any(stimGenerated)) {
      metaOut <- meta
      metaOut$locGenerated[] <- FALSE
      metaOut$locGeneratedDirect[] <- FALSE
      metaOut$locSource[] <- "not_calculated"
      metaOut$locReason[] <- "no_generated_local_fdr_threshold_to_combine"
      return(.getCpUnsLocCpAttachMeta(cp, metaOut))
    }

    # Generated gates, but none from a responder: nothing to share.
    if (!any(meta$locResponder)) {
      return(.getCpUnsLocCpAttachMeta(cp, meta))
    }

    cpShared <- .combineCp(
      cp = cpNum[meta$locResponder],
      gateCombn = gateCombnCurr
    )[[1]][[1]]
    freqDonor <- stats::median(freqOwn[meta$locResponder])
    cpOut <- stats::setNames(rep(cpShared, length(cpNum)), meta$ind)
    metaOut <- meta
    metaOut$locShareProposed[stimRow] <- cpShared
    for (i in which(stimRow & !vapply(xStimList, is.null, logical(1)))) {
      res <- .getCpShareApply(
        gs = cpShared,
        gc = cpNum[[i]],
        responder = meta$locResponder[[i]],
        propBsEst = meta$propBsEst[[i]],
        freqDonor = freqDonor,
        xStim = xStimList[[i]],
        xUns = xUns,
        shareCap = shareCap,
        cellCap = cellCap
      )
      cpOut[[i]] <- res$gate
      metaOut$locShareLimit[[i]] <- res$limit
    }

    metaOut$locGenerated[stimRow] <- TRUE
    metaOut$locGeneratedDirect[stimRow] <- (meta$locGeneratedDirect &
      abs(cpNum - cpOut) < 1e-7)[stimRow] %in% TRUE
    shared <- stimRow & !metaOut$locGeneratedDirect
    limited <- shared & !(metaOut$locShareLimit %in% "none")
    metaOut$locSource[shared] <- "combined"
    metaOut$locReason[shared] <- "combined_from_generated_local_fdr_thresholds"
    metaOut$locReason[limited] <- paste0(
      "combined_limited_by_", metaOut$locShareLimit[limited]
    )
    metaOut$locGenerated[!stimRow] <- TRUE
    metaOut$locGeneratedDirect[!stimRow] <- FALSE
    metaOut$locSource[!stimRow] <- "unstim_summary"
    metaOut$locReason[!stimRow] <-
      "summary_of_combined_generated_local_fdr_thresholds"
    .getCpUnsLocCpAttachMeta(cpOut, metaOut)
  }) |>
    stats::setNames(gateCombn)
}

# Background-subtracted frequency: proportion of stimulated cells strictly
# above each gate minus that of unstimulated cells.
#' @keywords internal
.getCpShareFreq <- function(gate, xStim, xUns) {
  vapply(gate, function(g) mean(xStim > g) - mean(xUns > g), numeric(1))
}

# A responder's own gate was generated directly by local FDR, is finite and
# has a positive background-subtracted frequency.
#' @keywords internal
.getCpShareResponder <- function(locGeneratedDirect, cp, freq) {
  (locGeneratedDirect %in% TRUE) & is.finite(cp) & ((freq > 0) %in% TRUE)
}

# Lowest gate at or above `gs` whose background-subtracted frequency is at
# most `limit`. Candidates are stimulated cell values above `gs`; the gate is
# placed below the selected cell (.getCpUnsLocGateBelowCell(), without a
# bandwidth) so that strict `x > gate` counts it. Above every stimulated cell
# the frequency is at most zero.
#' @keywords internal
.getCpShareLimitGate <- function(gs, limit, xStim, xUns) {
  if (!is.finite(gs) || isTRUE(.getCpShareFreq(gs, xStim, xUns) <= limit)) {
    return(gs)
  }
  cand <- sort(unique(xStim[xStim > gs]))
  freq <- .getCpUnsLocTailPropAtThresholds(xStim, cand, length(xStim)) -
    .getCpUnsLocTailPropAtThresholds(xUns, cand, length(xUns))
  ok <- which(freq <= limit)
  if (length(ok) == 0L) {
    return(max(xStim))
  }
  max(gs, .getCpUnsLocGateBelowCell(cand[[ok[[1L]]]], c(xStim, xUns)))
}

# Gate a tube accepts from a shared gate `gs`, given its current gate `gc`.
# Responders accept a higher gate, and a lower one only while their frequency
# is at most `shareCap * propBsEst` (never above `gc`; `gc` is kept when
# `propBsEst` is unavailable). Other tubes accept it only while their
# frequency is at most `cellCap / nStim` and the donors' median frequency
# `freqDonor`. Inf limits switch the rule off.
#' @keywords internal
.getCpShareApply <- function(
  gs,
  gc,
  responder,
  propBsEst,
  freqDonor,
  xStim,
  xUns,
  shareCap,
  cellCap
) {
  out <- function(gate, limit) {
    list(gate = gate, limit = if (isTRUE(gate == gs)) "none" else limit)
  }
  if (isTRUE(responder)) {
    if (!isTRUE(gs < gc) || is.infinite(shareCap)) {
      return(out(gs, "none"))
    }
    if (!is.finite(propBsEst) || propBsEst <= 0) {
      return(out(gc, "responder_no_estimate"))
    }
    gate <- .getCpShareLimitGate(gs, shareCap * propBsEst, xStim, xUns)
    return(out(min(gc, gate), "responder_cap"))
  }
  if (is.infinite(cellCap)) {
    return(out(gs, "none"))
  }
  limitCell <- cellCap / length(xStim)
  limitDonor <- if (is.finite(freqDonor)) max(0, freqDonor) else Inf
  gate <- .getCpShareLimitGate(gs, min(limitCell, limitDonor), xStim, xUns)
  out(
    gate,
    if (limitCell <= limitDonor) {
      "nonresponder_cell_cap"
    } else {
      "nonresponder_donor_cap"
    }
  )
}

#' @keywords internal
.getCpUnsLocMetaFromCp <- function(cp) {
  n <- length(cp)
  nm <- names(cp)
  if (is.null(nm)) {
    nm <- rep(NA_character_, n)
  }
  tibble::tibble(
    ind = as.character(nm),
    locGenerated = attr(cp, "locGenerated") %||% rep(FALSE, n),
    locGeneratedDirect = attr(cp, "locGeneratedDirect") %||% rep(FALSE, n),
    locSource = attr(cp, "locSource") %||% rep("not_calculated", n),
    locReason = attr(cp, "locReason") %||% rep(NA_character_, n),
    locResponder = attr(cp, "locResponder") %||% rep(FALSE, n),
    propBsEst = attr(cp, "propBsEst") %||% rep(NA_real_, n),
    locOwnFreq = attr(cp, "locOwnFreq") %||% rep(NA_real_, n),
    locShareLimit = attr(cp, "locShareLimit") %||% rep("none", n),
    locShareProposed = attr(cp, "locShareProposed") %||% rep(NA_real_, n)
  ) |>
    dplyr::mutate(
      locGenerated = .data$locGenerated %in% TRUE,
      locGeneratedDirect = .data$locGeneratedDirect %in% TRUE,
      locResponder = .data$locResponder %in% TRUE
    )
}

#' @keywords internal
.getCpUnsLocCpAttachMeta <- function(cp, meta) {
  n <- length(cp)
  if (nrow(meta) != n) {
    stop("Local-FDR metadata length does not match cutpoint vector length")
  }
  attr(cp, "locGenerated") <- meta$locGenerated %in% TRUE
  attr(cp, "locGeneratedDirect") <- meta$locGeneratedDirect %in% TRUE
  attr(cp, "locSource") <- as.character(meta$locSource)
  attr(cp, "locReason") <- as.character(meta$locReason)
  attr(cp, "locResponder") <-
    (meta[["locResponder"]] %||% rep(FALSE, n)) %in% TRUE
  attr(cp, "propBsEst") <- as.numeric(meta[["propBsEst"]] %||% rep(NA_real_, n))
  attr(cp, "locOwnFreq") <- as.numeric(
    meta[["locOwnFreq"]] %||% rep(NA_real_, n)
  )
  attr(cp, "locShareLimit") <- as.character(
    meta[["locShareLimit"]] %||% rep("none", n)
  )
  attr(cp, "locShareProposed") <- as.numeric(
    meta[["locShareProposed"]] %||% rep(NA_real_, n)
  )
  cp
}


# ------------------------------------------
# get cutpoints for a range of samples, and then individual samples
# ------------------------------------------

# Get cutpoint for a range of samples given the q-value and fdr
#' @keywords internal
.getCpUnsLocSample <- function(
  exListOrig,
  exListNoMinStim,
  exTblUnsBias,
  chnlSettings,
  bias,
  pathProject,
  stage,
  indStim = NULL,
  exListOrigOutput = NULL,
  prejoin = FALSE
) {
  .debug("getting loc gate at sample level") # nolint
  indStim <- indStim %||% names(exListNoMinStim)
  exListOrigOutput <- exListOrigOutput %||% exListOrig

  # get cutpoints for each sample
  cpUnsLocObjList <- purrr::map(
    seq_along(exListNoMinStim),
    function(i) {
      .debug("sample", i) # nolint

      exTblNoMinStim <- exListNoMinStim[[i]]
      ind <- .getInd(exTblNoMinStim)
      .debug("ind", ind) # nolint
      exTblUnsOrig <- exListOrig[[1]]
      exTblStimOrig <- exListOrig[[i + 1]]
      chnl <- chnlSettings$chnlCut %||%
        .getCpUnsLocGetChnl(exTblNoMinStim)
      stageChnl <- file.path(stage, chnl)
      .intSave(
        ind,
        stageChnl,
        pathProject,
        exTblNoMinStim,
        exTblUnsOrig,
        exTblStimOrig
      )

      # return early if there are too few cells
      tooFewCellsLgl <- .getCpUnsLocSampleCheckCellNumber(
        exTblStimNoMin = exTblNoMinStim,
        minCell = chnlSettings$minCell,
        exTblUnsBias = exTblUnsBias
      )
      if (tooFewCellsLgl) {
        objOut <- .getCpUnsLocTooFew(
          label = "too_few_cells_sample_fn",
          stage = stage,
          pathProject = pathProject,
          exTblNoMinStim = exTblNoMinStim,
          exTblUnsBias = exTblUnsBias,
          cpMin = chnlSettings$cpMin
        )
        return(objOut)
      }

      # Small samples widen a shared bandwidth; an automatic biasUns scales
      # with it, using the smaller of the stimulated and unstimulated tubes.
      sampleScale <- .getCpUnsLocSampleScale(
        nCell = min(nrow(exTblNoMinStim), nrow(exTblUnsBias)),
        chnlSettings = chnlSettings
      )
      chnlSettings$sampleScale <- sampleScale
      if (isTRUE(chnlSettings$biasUnsAuto) && sampleScale != 1) {
        biasSample <- bias * sampleScale
        chnlUns <- attr(exTblUnsBias, "chnlCut")
        exTblUnsBias[[chnlUns]] <- .getCut(exTblUnsBias) + (biasSample - bias)
        bias <- biasSample
      }

      # remove any cytokine-positive cells from unstim using gates from
      # sample for which gates are required
      exTblUnsBiasRm <- .getCpUnsLocSampleUnsRmCytPos(
        exTblUnsOrig = exTblUnsOrig,
        chnlSettings = chnlSettings,
        exTblStimNoMin = exTblNoMinStim,
        bias = bias,
        exTblUnsBias = exTblUnsBias,
        stage = stage
      )
      # removal can leave too few unstim cells, so check again
      if (nrow(exTblUnsBiasRm) < chnlSettings$minCell) {
        objOut <- .getCpUnsLocTooFew(
          label = "too_few_cells_sample_fn",
          stage = stage,
          pathProject = pathProject,
          exTblNoMinStim = exTblNoMinStim,
          exTblUnsBias = exTblUnsBias,
          cpMin = chnlSettings$cpMin
        )
        return(objOut)
      }
      exTblUnsBias <- exTblUnsBiasRm
      .intSave(ind, stageChnl, pathProject, exTblUnsBias)

      .getCpUnsLocCondition(
        exTblUnsBias = exTblUnsBias,
        exTblStimNoMin = exTblNoMinStim,
        chnlSettings = chnlSettings,
        exTblStimOrig = exTblStimOrig,
        exTblUnsOrig = exTblUnsOrig,
        bias = bias,
        pathProject = pathProject,
        stage = stage
      )
    }
  ) |>
    stats::setNames(names(exListNoMinStim))

  chnl <- .getCpUnsLocGetChnl(exListNoMinStim[[1]])
  .getCpUnsLocOutput(
    cpUnsLocObjList = cpUnsLocObjList,
    indUns = names(exListOrig)[1],
    indStim = indStim,
    stage = stage,
    pathProject = pathProject,
    chnl = chnl,
    exListOrig = exListOrigOutput,
    prejoin = prejoin
  )
}

#' @keywords internal
.getCpUnsLocSampleCheckCellNumber <- function(
  exTblStimNoMin,
  minCell,
  exTblUnsBias
) {
  nrow(exTblStimNoMin) < minCell ||
    nrow(exTblUnsBias) < minCell
}

# Fallback result when there are too few cells to gate; label names the
# intermediate save marking where the check failed
#' @keywords internal
.getCpUnsLocTooFew <- function(
  label,
  stage,
  pathProject,
  exTblNoMinStim,
  exTblUnsBias,
  cpMin
) {
  chnl <- .getCpUnsLocGetChnl(exTblNoMinStim)
  stageChnl <- file.path(stage, chnl)
  .intSaveNm(
    label,
    NULL,
    .getInd(exTblNoMinStim),
    stageChnl,
    pathProject
  ) # nolint
  objOut <- .getCpUnsLocConditionCheckOut(
    cpMin = cpMin,
    exTblStimNoMin = exTblNoMinStim,
    exTblUnsBias = exTblUnsBias,
    stage = stage,
    msg = "Too few cells"
  )
  .intSaveNm(
    "cpInd",
    objOut$cp,
    .getInd(exTblNoMinStim),
    stageChnl,
    pathProject
  )
  objOut
}

#' @keywords internal
.getCpUnsLocSampleUnsRmCytPos <- function(
  exTblUnsOrig,
  chnlSettings,
  exTblStimNoMin,
  bias,
  exTblUnsBias,
  stage
) {
  if (stage == "init") {
    return(exTblUnsBias)
  }
  .debug("Removing cytokine-positive cells from unstim") # nolint

  # first filter
  gateTblGnInd <- chnlSettings$gateTbl |>
    dplyr::filter(
      ind == exTblStimNoMin$ind[1], # nolint
      .data$gateName == .env$chnlSettings$gateNameCurr # nolint
    )

  posIndVecButSinglePosCurr <-
    .getPosIndButSinglePosForOneCyt(
      ex = exTblUnsOrig,
      gateTbl = gateTblGnInd,
      chnlSingleExc = chnlSettings$chnlCut,
      chnl = NULL,
      gateTypeCytPos = ifelse(
        chnlSettings$calcCytPosGates,
        "cyt",
        "base"
      )
    )

  exTblUnsOrig <- exTblUnsOrig[
    !posIndVecButSinglePosCurr, ,
    drop = FALSE
  ]

  # re-apply bias, noise and exclude minimum after removing cytokine-positive cells
  .prepareExListWithBiasAndNoise(
    exList = stats::setNames(
      list(exTblUnsOrig),
      attr(exTblUnsOrig, "indUns")
    ),
    ind = attr(exTblUnsOrig, "indUns"),
    excMin = chnlSettings$excMin,
    bias = bias
  )[[1]]
}

#' @keywords internal
.getCpUnsLocGetChnl <- function(exTbl) {
  if (is.null(exTbl) || !is.data.frame(exTbl)) {
    return("unknown_chnl")
  }
  chnl <- attr(exTbl, "chnlCut")
  if (!is.null(chnl) && nzchar(chnl)) {
    return(chnl)
  }
  cutVec <- .getCut(exTbl)
  cnVec <- colnames(exTbl)
  if (is.null(cnVec) || length(cnVec) == 0L) {
    return("unknown_chnl")
  }
  matchInd <- vapply(
    cnVec,
    function(nm) identical(exTbl[[nm]], cutVec),
    logical(1)
  )
  if (any(matchInd)) {
    return(cnVec[matchInd][1])
  }
  cnVec[1]
}

#' @keywords internal
.getCpUnsLocCondition <- function(
  exTblUnsBias,
  exTblStimNoMin,
  chnlSettings,
  exTblStimOrig,
  exTblUnsOrig,
  bias,
  pathProject,
  stage
) {
  .debug("getting loc gate for single sample") # nolint
  ind <- .getInd(exTblStimNoMin)
  chnl <- chnlSettings$chnlCut %||%
    .getCpUnsLocGetChnl(exTblStimNoMin)
  stageChnl <- file.path(stage, chnl)
  .debug("ind", ind) # nolint

  # return early if almost all stim expression is below cpMin (cell numbers
  # were already checked by .getCpUnsLocSample)
  cutStim <- .getCut(exTblStimNoMin)
  if (
    (stats::quantile(cutStim, 0.9) + 3 * stats::sd(cutStim)) <=
      chnlSettings$cpMin
  ) {
    objOut <- .getCpUnsLocTooFew(
      label = "too_few_cells_ind_fn",
      stage = stage,
      pathProject = pathProject,
      exTblNoMinStim = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias,
      cpMin = chnlSettings$cpMin
    )
    return(objOut)
  }

  exTblStimThreshold <- exTblStimNoMin
  exTblUnsThreshold <- exTblUnsBias
  .intSave(
    .getInd(exTblStimNoMin),
    stageChnl,
    pathProject,
    exTblStimThreshold,
    exTblUnsThreshold
  )

  # get smoothed probabilities
  dataMod <- .getCpUnsLocGetProb(
    exTblStimNoMin = exTblStimNoMin,
    exTblStimThreshold = exTblStimThreshold,
    exTblUnsThreshold = exTblUnsThreshold,
    exTblUnsBias = exTblUnsBias,
    stage = stage,
    bias = bias,
    exTblUnsOrig = exTblUnsOrig,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
  .intSave(
    .getInd(exTblStimNoMin),
    stageChnl,
    pathProject,
    dataMod
  )

  # get threshold
  .getCpUnsLocGetCp(
    dataMod = dataMod,
    exTblStimNoMin = exTblStimNoMin,
    exTblStimOrig = exTblStimOrig,
    exTblUnsOrig = exTblUnsOrig,
    exTblUnsBias = exTblUnsBias,
    bias = bias,
    cpMin = chnlSettings$cpMin,
    stage = stage,
    pathProject = pathProject,
    chnlSettings = chnlSettings
  )
}

#' @keywords internal
.getCpUnsLocConditionCheckOut <- function(
  cpMin,
  exTblStimNoMin,
  exTblUnsBias,
  stage,
  msg
) {
  .debug(msg) # nolint
  .getCpUnsLocConditionOut(
    cp = .getCpUnsLocConditionCpNonLoc(
      cpMin = cpMin,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsBias = exTblUnsBias
    ),
    locGenerated = FALSE,
    locGeneratedDirect = FALSE,
    locSource = "not_calculated",
    locReason = msg
  )
}


#' @keywords internal
.getCpUnsLocConditionOut <- function(
  cp,
  locGenerated,
  locGeneratedDirect,
  locSource,
  locReason
) {
  list(
    cp = cp,
    locGenerated = isTRUE(locGenerated),
    locGeneratedDirect = isTRUE(locGeneratedDirect),
    locSource = locSource,
    locReason = locReason
  )
}


#' @keywords internal
.getCpUnsLocConditionCpNonLoc <- function(
  cpMin,
  exTblStimNoMin,
  exTblUnsBias
) {
  rangeVecStim <- range(.getCut(exTblStimNoMin))
  rangeVecUns <- range(.getCut(exTblUnsBias))
  rangeLen <- max(diff(rangeVecStim), diff(rangeVecUns))
  max(
    cpMin,
    rangeVecStim[[2]] + rangeLen / 5,
    rangeVecUns[[2]] + rangeLen / 3
  )
}
