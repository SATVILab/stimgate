# Complete marker parameter list with all required settings
# Ensures all parameters for each marker requiring a gate are properly defined

#' @keywords internal
.completeChnlSettings <- function(
  chnl,
  markerControl,
  control,
  biasUns,
  bw,
  .data,
  popGate,
  indBatchList,
  pathProject
) {
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  chnlSettings <- .resolveMarkerControl(
    markerControl = markerControl,
    chnl = chnl,
    chnlLab = chnlLab
  )
  chnlSettingsCommon <- c(
    list(popGate = popGate, biasUns = biasUns, bw = bw),
    unclass(control)
  )
  chnlList <- purrr::map(chnl, function(chnlCurr) {
    chnlSettingsSpecCurr <- list(
      marker = chnlLab[[chnlCurr]],
      chnlCut = chnlCurr
    ) |>
      append(chnlSettings[[chnlCurr]])

    chnlOut <- .completeChnlSettingsInd(
      chnl = chnlCurr,
      chnlSettingsCommon = chnlSettingsCommon,
      chnlSettingsSpec = chnlSettingsSpecCurr,
      .data = .data,
      indBatchList = indBatchList,
      pathProject = pathProject
    )
    .verifyChnlSettingsChnl(chnlCurr, chnlOut)
    chnlOut
  }) |>
    stats::setNames(chnlLab[chnl])

  .completeChnlSettingsSave(
    chnlList = chnlList,
    pathProject = pathProject
  )

  chnlList
}

#' @keywords internal
.completeChnlSettingsInd <- function(
  chnlSettingsCommon,
  chnlSettingsSpec,
  chnl,
  .data,
  indBatchList,
  pathProject
) {
  chnlSettings <- .completeChnlSettingsAddCommon(
    chnlSettingsCommon = chnlSettingsCommon,
    chnlSettings = chnlSettingsSpec
  )
  if (
    !is.logical(chnlSettings$locEnforceShapeThreshold) ||
      length(chnlSettings$locEnforceShapeThreshold) != 1L ||
      is.na(chnlSettings$locEnforceShapeThreshold)
  ) {
    stop(
      "`locEnforceShapeThreshold` must be TRUE or FALSE for channel `",
      chnl,
      "`"
    )
  }

  # Bandwidth-method settings shared by every automatic bandwidth below
  bwArgs <- lapply(
    stats::setNames(nm = c(
      "bwMtd", "bwAdj", "normPeakMinRel", "normExtraFrac",
      "normExtraMax", "normLambda", "normDensityN",
      "normExcessBwMtd", "normExcessNcell", "normAdaptiveNcell", "normMtd"
    )),
    function(nm) chnlSettings[[nm]]
  )

  needBwMin <- .completeChnlSettingsBwLimitIsAuto(chnlSettings$bwMin) &&
    !.completeChnlSettingsBwLimitIsNone(chnlSettings$bwMin)
  needBwMax <- .completeChnlSettingsBwLimitIsAuto(chnlSettings$bwMax) &&
    !.completeChnlSettingsBwLimitIsNone(chnlSettings$bwMax)
  needBwFallback <- .completeChnlSettingsBwLimitIsAuto(chnlSettings$bwFallback)
  needCpMin <- is.null(chnlSettings$cpMin)

  exListByBatch <- if (needBwMin || needBwMax || needBwFallback || needCpMin) {
    purrr::map(
      .completeChnlSettingsBatchInd(indBatchList),
      function(i) {
        .getExList(
          .data = .data,
          indBatch = indBatchList[[i]],
          pop = chnlSettings$popGate,
          chnlCut = chnl,
          batch = names(indBatchList)[i],
          pathProject = pathProject
        )
      }
    )
  } else {
    NULL
  }

  xList <- if (needBwMin || needBwMax || needBwFallback) {
    .completeChnlSettingsGetBwExprList(
      indBatchList = indBatchList,
      .data = .data,
      popGate = chnlSettings$popGate,
      chnlCut = chnl,
      pathProject = pathProject,
      exListByBatch = exListByBatch
    )
  } else {
    NULL
  }

  chnlSettings$bwMin <- .completeChnlSettingsBwLimit(
    bwLimit = chnlSettings$bwMin,
    noneValue = -Inf,
    nSampleBw = 1e5,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwArgs = bwArgs,
    xList = xList
  )

  chnlSettings$bwMax <- .completeChnlSettingsBwLimit(
    bwLimit = chnlSettings$bwMax,
    noneValue = Inf,
    nSampleBw = 1e2,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwArgs = bwArgs,
    xList = xList
  )

  chnlSettings$bwFallback <- .completeChnlSettingsBwFallback(
    bwFallback = chnlSettings$bwFallback,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwArgs = bwArgs,
    xList = xList
  )

  chnlSettings$biasUns <- .completeChnlSettingsBiasUns(
    biasUns = chnlSettings$biasUns,
    biasUnsFactor = chnlSettings$biasUnsFactor,
    bwMin = chnlSettings$bwMin,
    bwMax = chnlSettings$bwMax,
    bwFallback = chnlSettings$bwFallback
  )

  chnlSettings <- .completeChnlSettingsBwShared(
    chnlSettings = chnlSettings,
    indBatchList = indBatchList,
    .data = .data,
    pathProject = pathProject
  )

  chnlSettings$cpMin <- .completeChnlSettingsCpMin(
    cpMin = chnlSettings$cpMin,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    indBatchList = indBatchList,
    pathProject = pathProject,
    exListByBatch = exListByBatch
  )

  chnlSettings
}

#' @keywords internal
.completeChnlSettingsAddCommon <- function(
  chnlSettingsCommon,
  chnlSettings
) {
  chnlSettings |>
    append(chnlSettingsCommon[
      setdiff(names(chnlSettingsCommon), names(chnlSettings))
    ])
}

#' @keywords internal
.completeChnlSettingsBiasUns <- function(
  biasUns,
  biasUnsFactor,
  bwMin,
  bwMax,
  bwFallback
) {
  if (!is.null(biasUns)) {
    return(biasUns)
  }
  if (!is.null(bwFallback)) {
    return(0.25 * bwFallback * biasUnsFactor)
  }

  bwRef <- c(bwMin, bwMax)
  bwRef <- bwRef[is.finite(bwRef) & bwRef > 0]

  if (length(bwRef) == 0L) {
    return(0)
  }

  0.25 * mean(bwRef) * biasUnsFactor
}


#' @keywords internal
.completeChnlSettingsBwLimitIsAuto <- function(x) {
  is.null(x) ||
    (is.character(x) &&
      length(x) == 1L &&
      tolower(x) == "auto")
}

#' @keywords internal
.completeChnlSettingsBwLimitIsNone <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    tolower(x) == "none"
}

# Evenly spaced indices: deterministic, and spread across the batch order
#' @keywords internal
.spreadInd <- function(n, size) {
  unique(round(seq(1, n, length.out = min(size, n))))
}

# Batches used to estimate automatic channel settings
#' @keywords internal
.completeChnlSettingsBatchInd <- function(indBatchList) {
  .spreadInd(length(indBatchList), 5)
}

#' @keywords internal
.completeChnlSettingsGetBwExprList <- function(
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  exListByBatch = NULL
) {
  if (is.null(exListByBatch)) {
    exListByBatch <- purrr::map(
      .completeChnlSettingsBatchInd(indBatchList),
      function(i) {
        .getExList(
          .data = .data,
          indBatch = indBatchList[[i]],
          pop = popGate,
          chnlCut = chnlCut,
          batch = names(indBatchList)[i],
          pathProject = pathProject
        )
      }
    )
  }
  purrr::map(
    exListByBatch,
    function(exList) {
      purrr::map(exList, function(ex) {
        xVec <- .getCut(ex)
        xVec <- xVec[is.finite(xVec)]
        xVec <- xVec[xVec > min(xVec, na.rm = TRUE)]
        xVec
      })
    }
  ) |>
    purrr::flatten() |>
    purrr::keep(function(x) length(x) >= 2L && length(unique(x)) >= 2L)
}

#' @keywords internal
.completeChnlSettingsBwLimit <- function(
  bwLimit,
  noneValue,
  nSampleBw,
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwArgs,
  xList = NULL
) {
  if (.completeChnlSettingsBwLimitIsNone(bwLimit)) {
    return(noneValue)
  }

  if (!.completeChnlSettingsBwLimitIsAuto(bwLimit)) {
    return(bwLimit)
  }

  if (is.null(xList)) {
    xList <- .completeChnlSettingsGetBwExprList(
      indBatchList = indBatchList,
      .data = .data,
      popGate = popGate,
      chnlCut = chnlCut,
      pathProject = pathProject
    )
  }

  bwVec <- purrr::map_dbl(xList, function(xVec) {
    as.numeric(do.call(.bwCalcOne, c(
      list(x = xVec, bwNcellMin = nSampleBw, bwNcellMax = nSampleBw),
      bwArgs
    )))[1]
  })

  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
  if (length(bwVec) == 0L) {
    return(.Machine$double.eps)
  }

  mean(bwVec, trim = 0.1, na.rm = TRUE)
}

#' @keywords internal
.completeChnlSettingsBwFallback <- function(
  bwFallback,
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwArgs,
  xList = NULL
) {
  if (!.completeChnlSettingsBwLimitIsAuto(bwFallback)) {
    return(bwFallback)
  }

  if (is.null(xList)) {
    xList <- .completeChnlSettingsGetBwExprList(
      indBatchList = indBatchList,
      .data = .data,
      popGate = popGate,
      chnlCut = chnlCut,
      pathProject = pathProject
    )
  }

  if (length(xList) == 0L) {
    return(.Machine$double.eps)
  }

  nCellFallback <- stats::median(purrr::map_int(xList, length), na.rm = TRUE)
  nCellFallback <- max(2L, as.integer(round(nCellFallback)))

  xListFallback <- xList[.spreadInd(
    length(xList),
    max(1L, ceiling(sqrt(length(xList))))
  )]

  bwVec <- purrr::map_dbl(xListFallback, function(xVec) {
    as.numeric(do.call(.bwCalcOne, c(
      list(x = xVec, bwNcellMin = nCellFallback, bwNcellMax = nCellFallback),
      bwArgs
    )))[1]
  })
  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]

  if (length(bwVec) == 0L) {
    stop(
      "Failed to calculate fallback bandwidth for channel ",
      chnlCut,
      ". Specify bwFallback manually."
    )
  }

  stats::median(bwVec, na.rm = TRUE)
}

#' @keywords internal
.completeChnlSettingsCpMin <- function(
  cpMin,
  .data,
  popGate,
  chnlCut,
  indBatchList,
  pathProject,
  exListByBatch = NULL
) {
  if (!is.null(cpMin)) {
    return(cpMin)
  }
  .debug("calculating cpMin automatically") # nolint
  if (is.null(exListByBatch)) {
    exListByBatch <- purrr::map(
      .completeChnlSettingsBatchInd(indBatchList),
      function(i) {
        .getExList(
          # nolint
          .data = .data,
          indBatch = indBatchList[[i]],
          pop = popGate,
          chnlCut = chnlCut,
          batch = names(indBatchList)[i],
          pathProject = pathProject
        )
      }
    )
  }
  purrr::map(
    exListByBatch,
    function(exList) {
      purrr::map_dbl(exList, function(ex) {
        cutVals <- .getCut(ex)
        stats::median(cutVals[cutVals > min(cutVals)], na.rm = TRUE)[[1]]
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}

# Get all cutpoint type names
# Returns character vector of all available cutpoint names
#' @keywords internal
.completeChnlSettingsSave <- function(chnlList, pathProject) {
  pathSave <- file.path(pathProject, "metaData", "chnlSettings.rds")
  if (file.exists(pathSave)) {
    invisible(file.remove(pathSave))
  }
  if (!dir.exists(pathProject)) {
    dir.create(pathProject, recursive = TRUE)
  }
  saveRDS(
    chnlList,
    file = pathSave
  )
}

#' @title Read saved gating settings
#' @description Read settings saved by [gateStim()].
#' @param pathProject character Project directory.
#' @return A list of settings per marker, named by the saved keys (marker labels
#'   for projects created by [gateStim()]).
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(list(IFNg = list(bw = 0.1)),
#'   file.path(pathProject, "metaData", "chnlSettings.rds")
#' )
#' stimgateMetaReadSettingsChnls(pathProject)
#' @export
stimgateMetaReadSettingsChnls <- function(pathProject) {
  pathChnlList <- file.path(pathProject, "metaData", "chnlSettings.rds")
  if (!file.exists(pathChnlList)) {
    stop("Channel list file not found in project metaData folder")
  }
  readRDS(pathChnlList)
}

#' @title Relabel saved settings
#' @description Read [stimgateMetaReadSettingsChnls()] and map its names through
#'   [stimgateMetaReadChnlLab()]. Use only when the saved keys are channel names.
#' @param pathProject character Project directory.
#' @return A list with names replaced by marker labels; unmatched keys become NA.
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(list(BC1 = list(bw = 0.1)),
#'   file.path(pathProject, "metaData", "chnlSettings.rds")
#' )
#' saveRDS(c(BC1 = "IFNg"), file.path(pathProject, "metaData", "chnlLab.rds"))
#' stimgateMetaReadSettingsMarkers(pathProject)
#' @export
stimgateMetaReadSettingsMarkers <- function(pathProject) {
  markerList <- stimgateMetaReadSettingsChnls(pathProject)
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  names(markerList) <- chnlLab[names(markerList)]
  markerList
}

#' @title Read settings by saved key
#' @description Extract one entry from [stimgateMetaReadSettingsChnls()].
#'   The key must match exactly; channel names are not converted to marker labels.
#' @param pathProject character Project directory.
#' @param chnl character Exact saved key, usually a marker label.
#' @return A list of settings for the key; an unknown key raises an error.
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(list(IFNg = list(bw = 0.1)),
#'   file.path(pathProject, "metaData", "chnlSettings.rds")
#' )
#' stimgateMetaReadSettingsChnl(pathProject, "IFNg")
#' @export
stimgateMetaReadSettingsChnl <- function(pathProject, chnl) {
  chnlList <- stimgateMetaReadSettingsChnls(pathProject)
  if (!chnl %in% names(chnlList)) {
    stop(sprintf("Channel %s not found in marker list", chnl))
  }
  chnlList[[chnl]]
}

#' @title Read settings for a marker
#' @description Extract one entry from [stimgateMetaReadSettingsChnls()].
#' @param pathProject character Project directory.
#' @param marker character Exact saved marker key.
#' @return A list of marker settings; an unknown key raises an error.
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(list(IFNg = list(bw = 0.1)),
#'   file.path(pathProject, "metaData", "chnlSettings.rds")
#' )
#' stimgateMetaReadSettingsMarker(pathProject, "IFNg")
#' @export
stimgateMetaReadSettingsMarker <- function(pathProject, marker) {
  markerList <- stimgateMetaReadSettingsChnls(pathProject)
  if (!marker %in% names(markerList)) {
    stop(sprintf("Marker %s not found in marker list", marker))
  }
  markerList[[marker]]
}

#' @rdname stimgateMetaReadLab
#' @title Read channel and marker mappings
#' @description Read channel-to-marker labels with `stimgateMetaReadChnlLab()`;
#'   read the reverse mapping with `stimgateMetaReadMarkerLab()`.
#' @param pathProject character Project directory from [gateStim()].
#' @return A named character vector: channel names to marker labels for
#'   `stimgateMetaReadChnlLab()`, marker labels to channels for
#'   `stimgateMetaReadMarkerLab()`.
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(c(BC1 = "IFNg"), file.path(pathProject, "metaData", "chnlLab.rds"))
#' stimgateMetaReadChnlLab(pathProject)
#' stimgateMetaReadMarkerLab(pathProject)
#' @export
stimgateMetaReadChnlLab <- function(pathProject) {
  pathChnlLab <- file.path(pathProject, "metaData", "chnlLab.rds")
  if (!file.exists(pathChnlLab)) {
    stop("Channel label file not found in project metaData folder")
  }
  readRDS(pathChnlLab)
}

#' @rdname stimgateMetaReadLab
#' @export
stimgateMetaReadMarkerLab <- function(pathProject) {
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  stats::setNames(names(chnlLab), chnlLab)
}

#' @keywords internal
.saveMetaData <- function(.data, batchList, pathProject) {
  pathDirMetaData <- file.path(pathProject, "metaData")
  if (!dir.exists(pathDirMetaData)) {
    dir.create(pathDirMetaData, recursive = TRUE)
  }
  .saveMetaDataChnlLab(.data, pathDirMetaData)
  .saveMetaDataBatchList(batchList, pathDirMetaData)
}

#' @keywords internal
.saveMetaDataChnlLab <- function(.data, pathDir) {
  chnlLab <- chnlLab(.data)
  saveRDS(
    chnlLab,
    file = file.path(pathDir, "chnlLab.rds")
  )
}

#' @keywords internal
.saveMetaDataBatchList <- function(batchList, pathDir) {
  saveRDS(
    batchList,
    file = file.path(pathDir, "batchList.rds")
  )
}

#' @title Read saved batches
#' @description Read the sample grouping saved by [gateStim()].
#' @param pathProject character Project directory from [gateStim()].
#' @return A list of sample indices by batch, with unstimulated controls first.
#' @examples
#' pathProject <- tempfile("stimgate_meta_")
#' dir.create(file.path(pathProject, "metaData"), recursive = TRUE)
#' saveRDS(list(batch1 = c(1, 2)),
#'   file.path(pathProject, "metaData", "batchList.rds")
#' )
#' stimgateMetaReadBatchList(pathProject)
#' @export
stimgateMetaReadBatchList <- function(pathProject) {
  pathBatchList <- file.path(pathProject, "metaData", "batchList.rds")
  if (!file.exists(pathBatchList)) {
    stop("Batch list file not found in project metaData folder")
  }
  readRDS(pathBatchList)
}

.extractChnl <- function(chnl, marker, pathProject) {
  if (!is.null(chnl)) {
    if (!length(chnl) == length(unique(chnl))) {
      stop(
        "Duplicate channel names found in `chnl`. Please ensure that each channel is unique."
      )
    }
    return(chnl)
  }
  markerLab <- stimgateMetaReadMarkerLab(pathProject)
  chnlVec <- markerLab[marker] |> stats::setNames(NULL)
  if (!length(chnlVec) == length(unique(chnlVec))) {
    stop(
      "Duplicate channel labels found for the specified markers. ",
      "Please ensure that each marker has a unique channel label. ",
      "Otherwise, simply specify `chnl` instead of `marker`. "
    )
  }
  chnlVec
}

#' @keywords internal
.resolveMarkerControl <- function(markerControl, chnl, chnlLab) {
  if (is.null(markerControl)) {
    return(stats::setNames(lapply(chnl, function(x) list()), chnl))
  }
  if (!is.list(markerControl)) {
    stop("`markerControl` must be NULL or a named list.")
  }
  nms <- names(markerControl)
  if (
    length(markerControl) > 0L &&
      (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))
  ) {
    stop("`markerControl` elements must have non-empty names.")
  }
  allowed <- c(
    setdiff(
      names(formals(stimControl)),
      c("locEnforceShapeThreshold", "calcCytPosGates")
    ),
    "biasUns", "bw", "popGate"
  )
  resolved <- vapply(seq_along(markerControl), function(i) {
    nm <- nms[[i]]
    matches <- chnl[chnl == nm | chnlLab[chnl] == nm]
    if (length(matches) != 1L) {
      stop("Unknown or ambiguous marker/channel in `markerControl`: ", nm)
    }
    settings <- markerControl[[i]]
    if (!is.list(settings)) {
      stop("`markerControl` setting for '", nm, "' must be a list.")
    }
    settingNames <- names(settings)
    if (
      length(settings) > 0L &&
        (is.null(settingNames) || anyNA(settingNames) ||
          any(!nzchar(settingNames)) || anyDuplicated(settingNames) > 0L)
    ) {
      stop("`markerControl` settings for '", nm, "' must have unique, non-empty names.")
    }
    invalid <- setdiff(settingNames, allowed)
    if (length(invalid) > 0L) {
      stop(
        "Invalid settings for marker/channel '", nm, "': ",
        paste(invalid, collapse = ", ")
      )
    }
    .verifyChnlSettingsChnl(nm, settings)
    matches[[1]]
  }, character(1))
  if (anyDuplicated(resolved) > 0L) {
    stop("`markerControl` entries must resolve to distinct channels.")
  }
  stats::setNames(markerControl, resolved)
}
