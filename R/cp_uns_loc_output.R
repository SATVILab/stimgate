# Local-FDR diagnostics, metadata, and output assembly
#
# Builds condition- and sample-level diagnostic tables, propagates threshold
# provenance, combines sample thresholds, and assembles the final output.

.getCpUnsLocThresholdOrigin <- function(
  locGenerated,
  locGeneratedDirect,
  locSource,
  locReason = NA_character_
) {
  locGenerated <- locGenerated %in% TRUE
  locGeneratedDirect <- locGeneratedDirect %in% TRUE
  locSource <- as.character(locSource %||% NA_character_)
  locReason <- as.character(locReason %||% NA_character_)

  if (locGenerated && locGeneratedDirect && locSource %in% "direct") {
    return("condition_detected_response")
  }
  if (locGenerated && locSource %in% "combined") {
    return("sample_imputed_from_other_stim_conditions")
  }
  if (locGenerated && locSource %in% "prejoin") {
    return("prejoin_generated_from_joined_stim_conditions")
  }
  if (locGenerated && locSource %in% "unstim_summary") {
    return("unstim_summary_from_generated_thresholds")
  }
  if (!locGenerated) {
    return("not_generated_fallback")
  }
  if (!is.na(locReason) && nzchar(locReason)) {
    return(paste0("generated_", locReason))
  }
  "generated_unknown_source"
}

#' @keywords internal
.getCpUnsLocPropBsAtCp <- function(cp, exTblStim, exTblUns) {
  cp <- suppressWarnings(as.numeric(cp))[1]
  if (
    !is.finite(cp) ||
      is.null(exTblStim) ||
      !is.data.frame(exTblStim) ||
      is.null(exTblUns) ||
      !is.data.frame(exTblUns) ||
      nrow(exTblStim) == 0L ||
      nrow(exTblUns) == 0L
  ) {
    return(tibble::tibble(
      nCellStim = if (is.data.frame(exTblStim)) {
        nrow(exTblStim)
      } else {
        NA_integer_
      },
      nCellUns = if (is.data.frame(exTblUns)) nrow(exTblUns) else NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propBs = NA_real_
    ))
  }
  propStim <- sum(.getCut(exTblStim) >= cp, na.rm = TRUE) / nrow(exTblStim)
  propUns <- sum(.getCut(exTblUns) >= cp, na.rm = TRUE) / nrow(exTblUns)
  tibble::tibble(
    nCellStim = nrow(exTblStim),
    nCellUns = nrow(exTblUns),
    propStim = propStim,
    propUns = propUns,
    propBs = propStim - propUns
  )
}

#' @keywords internal
.getCpUnsLocSelectedThresholdRow <- function(dataThreshold, cp) {
  if (!is.data.frame(dataThreshold) || nrow(dataThreshold) == 0L) {
    return(NULL)
  }
  cp <- suppressWarnings(as.numeric(cp))[1]
  if (!is.finite(cp)) {
    return(NULL)
  }
  cutVec <- suppressWarnings(as.numeric(.getCut(dataThreshold)))
  idx <- which(is.finite(cutVec))
  if (length(idx) == 0L) {
    return(NULL)
  }
  idx <- idx[which.min(abs(cutVec[idx] - cp))]
  dataThreshold[idx, , drop = FALSE]
}

#' @keywords internal
.getCpUnsLocConditionDetailRow <- function(
  cpObj,
  dataThreshold,
  exTblStimOrig,
  exTblUnsOrig,
  exTblStimNoMin,
  bias,
  stage,
  chnl
) {
  cp <- suppressWarnings(as.numeric(cpObj$cp))[1]
  selectedRow <- .getCpUnsLocSelectedThresholdRow(dataThreshold, cp)
  freqTbl <- .getCpUnsLocPropBsAtCp(
    cp = cp,
    exTblStim = exTblStimOrig,
    exTblUns = exTblUnsOrig
  )

  propBsEst <- NA_real_
  propBsDiff <- NA_real_
  if (!is.null(selectedRow) && "propBsDiff" %in% names(selectedRow)) {
    propBsDiff <- suppressWarnings(as.numeric(selectedRow$propBsDiff[1]))
    if (is.finite(propBsDiff) && "propBs" %in% names(selectedRow)) {
      propBsEst <- suppressWarnings(as.numeric(selectedRow$propBs[1])) -
        propBsDiff
    }
  }

  tibble::tibble(
    detailLevel = "condition",
    stage = stage,
    chnl = chnl,
    ind = as.character(.getInd(exTblStimNoMin)),
    threshold = cp,
    thresholdOrigin = .getCpUnsLocThresholdOrigin(
      locGenerated = cpObj$locGenerated,
      locGeneratedDirect = cpObj$locGeneratedDirect,
      locSource = cpObj$locSource,
      locReason = cpObj$locReason
    ),
    locGenerated = cpObj$locGenerated %in% TRUE,
    locGeneratedDirect = cpObj$locGeneratedDirect %in% TRUE,
    locSource = as.character(cpObj$locSource %||% NA_character_),
    locReason = as.character(cpObj$locReason %||% NA_character_),
    bias = suppressWarnings(as.numeric(bias))[1],
    propBsEst = propBsEst,
    propBsDiff = propBsDiff
  ) |>
    dplyr::bind_cols(freqTbl)
}

#' @keywords internal
.getCpUnsLocExByInd <- function(exList, ind, pos = NULL) {
  ind <- as.character(ind)
  if (!is.null(names(exList)) && ind %in% names(exList)) {
    return(exList[[ind]])
  }
  if (!is.null(pos) && is.finite(pos) && pos >= 1L && pos <= length(exList)) {
    return(exList[[pos]])
  }
  NULL
}

#' @keywords internal
.getCpUnsLocSampleDetailTbl <- function(
  cpVec,
  exListOrig,
  indUns,
  indStim,
  stage,
  chnl
) {
  meta <- .getCpUnsLocMetaFromCp(cpVec)
  exTblUns <- .getCpUnsLocExByInd(exListOrig, indUns, pos = 1L)

  purrr::map_df(seq_along(indStim), function(i) {
    indCurr <- as.character(indStim[[i]])
    cp <- suppressWarnings(as.numeric(cpVec[indCurr]))[1]
    metaRow <- meta[meta$ind == indCurr, , drop = FALSE]
    if (nrow(metaRow) == 0L) {
      metaRow <- tibble::tibble(
        ind = indCurr,
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "not_calculated",
        locReason = NA_character_
      )
    }
    exTblStim <- .getCpUnsLocExByInd(exListOrig, indCurr, pos = i + 1L)
    freqTbl <- .getCpUnsLocPropBsAtCp(
      cp = cp,
      exTblStim = exTblStim,
      exTblUns = exTblUns
    )
    tibble::tibble(
      detailLevel = "sample",
      stage = stage,
      chnl = chnl,
      ind = indCurr,
      threshold = cp,
      thresholdOrigin = .getCpUnsLocThresholdOrigin(
        locGenerated = metaRow$locGenerated[1],
        locGeneratedDirect = metaRow$locGeneratedDirect[1],
        locSource = metaRow$locSource[1],
        locReason = metaRow$locReason[1]
      ),
      locGenerated = metaRow$locGenerated[1] %in% TRUE,
      locGeneratedDirect = metaRow$locGeneratedDirect[1] %in% TRUE,
      locSource = as.character(metaRow$locSource[1] %||% NA_character_),
      locReason = as.character(metaRow$locReason[1] %||% NA_character_),
      propBsEst = NA_real_,
      propBsDiff = NA_real_
    ) |>
      dplyr::bind_cols(freqTbl)
  })
}

# get cp
#' @keywords internal


.createCombinedIdentifier <- function(indStim) {
  if (is.null(indStim) || length(indStim) == 0) {
    "empty_batch"
  } else {
    paste(indStim, collapse = "_")
  }
}

#' @keywords internal
.getCpUnsLocOutput <- function(
  cpUnsLocObjList,
  indUns,
  indStim,
  stage,
  pathProject,
  chnl,
  exListOrig
) {
  stageChnl <- file.path(stage, chnl)
  cpVec <- .getCpUnsLocSampleCpRep(
    stage = stage,
    cpUnsLocObjList = cpUnsLocObjList,
    indUns = indUns,
    indStim = indStim,
    pathProject = pathProject,
    chnl = chnl
  )
  indCombined <- .createCombinedIdentifier(indStim)
  .intSave(indCombined, stageChnl, pathProject, cpVec)
  .intSaveNm(
    "cpUnsLocMeta",
    .getCpUnsLocMetaFromCp(cpVec),
    indCombined,
    stageChnl,
    pathProject
  )
  locDetailSample <- .getCpUnsLocSampleDetailTbl(
    cpVec = cpVec,
    exListOrig = exListOrig,
    indUns = indUns,
    indStim = indStim,
    stage = stage,
    chnl = chnl
  )
  .intSaveNm(
    "locDetailSample",
    locDetailSample,
    indCombined,
    stageChnl,
    pathProject
  )
  .debug("done getting loc gate at sample level") # nolint
  list(
    "loc" = cpVec,
    "pList" = list()
  )
}


#' @keywords internal
.getCpUnsLocSampleCpRep <- function(
  stage,
  cpUnsLocObjList,
  indUns,
  indStim,
  pathProject,
  chnl
) {
  .debug("Possibly re-using calculated cutpoints") # nolint
  indCombined <- .createCombinedIdentifier(indStim)
  stageChnl <- file.path(stage, chnl)

  cpVec <- purrr::map_dbl(cpUnsLocObjList, ~ .x[["cp"]])
  meta <- .getCpUnsLocSampleMeta(
    cpUnsLocObjList = cpUnsLocObjList,
    ind = names(cpUnsLocObjList)
  )
  .intSaveNm(
    "cpVecBeforeRep",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaBeforeRep",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  if (length(cpVec) != length(indStim)) {
    .intSaveNm(
      "prejoinedCpUsed",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- stats::setNames(
      rep(cpVec, length(indStim)),
      indStim
    )
    meta <- .getCpUnsLocSampleMetaRep(
      meta = meta,
      indStim = indStim,
      source = "prejoin"
    )
  } else {
    .intSaveNm(
      "individualCpUsed",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- stats::setNames(cpVec, indStim)
    meta$ind <- as.character(indStim)
  }
  .intSaveNm(
    "cpVecAfterRep",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaAfterRep",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  locGenerated <- meta$locGenerated %in% TRUE
  if (any(locGenerated, na.rm = TRUE)) {
    .intSaveNm(
      "addingUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpUns <- mean(cpVec[locGenerated], na.rm = TRUE)
    metaUns <- tibble::tibble(
      ind = as.character(indUns),
      locGenerated = TRUE,
      locGeneratedDirect = FALSE,
      locSource = "unstim_summary",
      locReason = "mean_of_generated_local_fdr_thresholds"
    )
    cpVec <- c(cpVec, stats::setNames(cpUns, indUns))
    meta <- dplyr::bind_rows(meta, metaUns)
  } else {
    .intSaveNm(
      "addNaForUnsCp",
      NULL,
      indCombined,
      stageChnl,
      pathProject
    )
    cpVec <- c(cpVec, stats::setNames(NA_real_, indUns))
    meta <- dplyr::bind_rows(
      meta,
      tibble::tibble(
        ind = as.character(indUns),
        locGenerated = FALSE,
        locGeneratedDirect = FALSE,
        locSource = "unstim_summary",
        locReason = "no_generated_local_fdr_thresholds"
      )
    )
  }
  .intSaveNm(
    "cpVecAfterUnsAdded",
    cpVec,
    indCombined,
    stageChnl,
    pathProject
  )
  .intSaveNm(
    "cpMetaAfterUnsAdded",
    meta,
    indCombined,
    stageChnl,
    pathProject
  )

  .getCpUnsLocCpAttachMeta(cpVec, meta)
}

#' @keywords internal
.getCpUnsLocSampleMeta <- function(cpUnsLocObjList, ind) {
  tibble::tibble(
    ind = as.character(ind),
    locGenerated = purrr::map_lgl(
      cpUnsLocObjList,
      ~ .x[["locGenerated"]] %||% FALSE
    ),
    locGeneratedDirect = purrr::map_lgl(
      cpUnsLocObjList,
      ~ .x[["locGeneratedDirect"]] %||% FALSE
    ),
    locSource = purrr::map_chr(
      cpUnsLocObjList,
      ~ .x[["locSource"]] %||% "not_calculated"
    ),
    locReason = purrr::map_chr(
      cpUnsLocObjList,
      ~ .x[["locReason"]] %||% NA_character_
    )
  )
}

#' @keywords internal
.getCpUnsLocSampleMetaRep <- function(meta, indStim, source) {
  if (nrow(meta) == 0L) {
    return(tibble::tibble(
      ind = as.character(indStim),
      locGenerated = FALSE,
      locGeneratedDirect = FALSE,
      locSource = "not_calculated",
      locReason = "no_local_fdr_objects"
    ))
  }
  metaOne <- meta[1, , drop = FALSE]
  metaOut <- metaOne[rep(1, length(indStim)), , drop = FALSE]
  metaOut$ind <- as.character(indStim)
  metaOut$locSource <- ifelse(
    metaOut$locGenerated %in% TRUE,
    source,
    metaOut$locSource
  )
  metaOut
}
