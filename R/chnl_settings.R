# Complete marker parameter list with all required settings
# Ensures all parameters for each marker requiring a gate are properly defined

#' @keywords internal
.completeChnlSettings <- function(
  chnl,
  marker,
  chnlSettings,
  markerSettings,
  biasUns,
  biasUnsFactor,
  excMin,
  .data,
  popGate,
  indBatchList,
  bw,
  bwMin,
  bwMax,
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwCluster,
  cpMin,
  minCell,
  tolClust,
  locProbCol,
  locMinPeakProb,
  locDipAlpha,
  locAntimodeHeightFrac,
  locAntimodeLowRel,
  locAntimodeLowAbs,
  locFlatDerivFrac,
  locFlatHardDerivFrac,
  locLeftLowRel,
  locLeftLowAbs,
  locLeftCellFrac,
  locLeftLengthFrac,
  locMarginalPurityRel,
  locMarginalCellBinRatio,
  locMarginalRefQuantile,
  locTolRefPeak,
  maxPosProbX,
  gateCombn,
  gateQuant,
  pathProject
) {
  chnlSettings <- .extractChnlSettings(
    chnl = chnl,
    marker = marker,
    chnlSettings = chnlSettings,
    markerSettings = markerSettings,
    pathProject = pathProject
  )
  chnlSettingsCommon <- list(
    popGate = popGate,
    biasUns = biasUns,
    biasUnsFactor = biasUnsFactor,
    excMin = excMin,
    cpMin = cpMin,
    bwMin = bwMin,
    bw = bw,
    bwMax = bwMax,
    bwFallback = bwFallback,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax,
    bwCluster = bwCluster,
    minCell = minCell,
    tolClust = tolClust,
    locProbCol = locProbCol,
    locMinPeakProb = locMinPeakProb,
    locDipAlpha = locDipAlpha,
    locAntimodeHeightFrac = locAntimodeHeightFrac,
    locAntimodeLowRel = locAntimodeLowRel,
    locAntimodeLowAbs = locAntimodeLowAbs,
    locFlatDerivFrac = locFlatDerivFrac,
    locFlatHardDerivFrac = locFlatHardDerivFrac,
    locLeftLowRel = locLeftLowRel,
    locLeftLowAbs = locLeftLowAbs,
    locLeftCellFrac = locLeftCellFrac,
    locLeftLengthFrac = locLeftLengthFrac,
    locMarginalPurityRel = locMarginalPurityRel,
    locMarginalCellBinRatio = locMarginalCellBinRatio,
    locMarginalRefQuantile = locMarginalRefQuantile,
    locTolRefPeak = locTolRefPeak,
    gateCombn = gateCombn,
    maxPosProbX = maxPosProbX,
    gateQuant = gateQuant
  )
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
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

  chnlSettings$bwMin <- .completeChnlSettingsBwMin(
    bwMin = chnlSettings$bwMin,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax
  )

  chnlSettings$bwMax <- .completeChnlSettingsBwMax(
    bwMax = chnlSettings$bwMax,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax
  )

  chnlSettings$bwFallback <- .completeChnlSettingsBwFallback(
    bwFallback = chnlSettings$bwFallback,
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax
  )

  chnlSettings$biasUns <- .completeChnlSettingsBiasUns(
    biasUns = chnlSettings$biasUns,
    biasUnsFactor = chnlSettings$biasUnsFactor,
    bwMin = chnlSettings$bwMin,
    bwMax = chnlSettings$bwMax,
    bwFallback = chnlSettings$bwFallback
  )

  chnlSettings$bwCluster <- .completeChnlSettingsBwCluster(
    indBatchList = indBatchList,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    pathProject = pathProject,
    bwCluster = chnlSettings$bwCluster,
    bwMtd = chnlSettings$bwMtd,
    bwAdj = chnlSettings$bwAdj,
    bwNcellMin = chnlSettings$bwNcellMin,
    bwNcellMax = chnlSettings$bwNcellMax,
    bwFallback = chnlSettings$bwFallback
  )

  chnlSettings$cpMin <- .completeChnlSettingsCpMin(
    cpMin = chnlSettings$cpMin,
    .data = .data,
    popGate = chnlSettings$popGate,
    chnlCut = chnl,
    indBatchList = indBatchList,
    pathProject = pathProject
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
    return(bwFallback * biasUnsFactor)
  }

  bwRef <- c(bwMin, bwMax)
  bwRef <- bwRef[is.finite(bwRef) & bwRef > 0]

  if (length(bwRef) == 0L) {
    return(0)
  }

  mean(bwRef) * biasUnsFactor
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

#' @keywords internal
.completeChnlSettingsBwFallbackIsAuto <- function(x) {
  is.null(x) ||
    (is.character(x) &&
      length(x) == 1L &&
      tolower(x) == "auto")
}

#' @keywords internal
.completeChnlSettingsGetBwExprList <- function(
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject
) {
  batchInd <- seq_along(indBatchList)
  batchInd <- sample(batchInd, size = min(5, length(batchInd)))

  purrr::map(
    batchInd,
    function(i) {
      exList <- .getExList(
        .data = .data,
        indBatch = indBatchList[[i]],
        pop = popGate,
        chnlCut,
        batch = names(indBatchList)[i],
        pathProject = pathProject
      )

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
.completeChnlSettingsBwCalcOne <- function(
  x,
  bwMtd,
  bwAdj,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  .bwCalcOne(
    x = x,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    bwNcellMin = bwNcellMin,
    bwNcellMax = bwNcellMax
  )
}

#' @keywords internal
.completeChnlSettingsBwLimitAuto <- function(
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwMtd,
  bwAdj,
  nSampleBw
) {
  xList <- .completeChnlSettingsGetBwExprList(
    indBatchList = indBatchList,
    .data = .data,
    popGate = popGate,
    chnlCut = chnlCut,
    pathProject = pathProject
  )

  bwVec <- purrr::map_dbl(xList, function(xVec) {
    .completeChnlSettingsBwCalcOne(
      x = xVec,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = nSampleBw,
      bwNcellMax = nSampleBw
    )
  })

  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]
  if (length(bwVec) == 0L) {
    return(.Machine$double.eps)
  }

  mean(bwVec, trim = 0.1, na.rm = TRUE)
}

#' @keywords internal
.completeChnlSettingsBwMax <- function(
  bwMax,
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwMtd,
  bwAdj,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  if (.completeChnlSettingsBwLimitIsNone(bwMax)) {
    return(Inf)
  }

  if (!.completeChnlSettingsBwLimitIsAuto(bwMax)) {
    return(bwMax)
  }

  .completeChnlSettingsBwLimitAuto(
    indBatchList = indBatchList,
    .data = .data,
    popGate = popGate,
    chnlCut = chnlCut,
    pathProject = pathProject,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    nSampleBw = 1e2
  )
}

#' @keywords internal
.completeChnlSettingsBwMin <- function(
  bwMin,
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwMtd,
  bwAdj,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  if (.completeChnlSettingsBwLimitIsNone(bwMin)) {
    return(-1)
  }

  if (!.completeChnlSettingsBwLimitIsAuto(bwMin)) {
    return(bwMin)
  }

  .completeChnlSettingsBwLimitAuto(
    indBatchList = indBatchList,
    .data = .data,
    popGate = popGate,
    chnlCut = chnlCut,
    pathProject = pathProject,
    bwMtd = bwMtd,
    bwAdj = bwAdj,
    nSampleBw = 1e5
  )
}

#' @keywords internal
.completeChnlSettingsBwFallback <- function(
  bwFallback,
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwMtd,
  bwAdj,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  if (!.completeChnlSettingsBwFallbackIsAuto(bwFallback)) {
    return(bwFallback)
  }

  xList <- .completeChnlSettingsGetBwExprList(
    indBatchList = indBatchList,
    .data = .data,
    popGate = popGate,
    chnlCut = chnlCut,
    pathProject = pathProject
  )

  if (length(xList) == 0L) {
    return(.Machine$double.eps)
  }

  nCellFallback <- stats::median(purrr::map_int(xList, length), na.rm = TRUE)
  nCellFallback <- max(2L, as.integer(round(nCellFallback)))

  xListFallback <- sample(
    xList,
    size = min(length(xList), max(1L, ceiling(sqrt(length(xList))))),
    replace = FALSE
  )

  bwVec <- purrr::map_dbl(xListFallback, function(xVec) {
    .completeChnlSettingsBwCalcOne(
      x = xVec,
      bwMtd = bwMtd,
      bwAdj = bwAdj,
      bwNcellMin = nCellFallback,
      bwNcellMax = nCellFallback
    )
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
.completeChnlSettingsBwCluster <- function(
  indBatchList,
  .data,
  popGate,
  chnlCut,
  pathProject,
  bwCluster,
  bwMtd,
  bwAdj,
  bwFallback,
  bwNcellMin = NULL,
  bwNcellMax = NULL
) {
  if (!is.null(bwCluster)) {
    if (
      is.numeric(bwCluster) &&
        length(bwCluster) == 1L &&
        is.finite(bwCluster) &&
        bwCluster > 0
    ) {
      return(bwCluster)
    }
    return(bwFallback)
  }

  bwVec <- purrr::map(
    seq_len(min(5, length(indBatchList))),
    function(i) {
      exList <- try(
        .getExList(
          .data = .data,
          indBatch = indBatchList[[i]],
          pop = popGate,
          chnlCut,
          batch = names(indBatchList)[i],
          pathProject = pathProject
        ),
        silent = TRUE
      )

      if (inherits(exList, "try-error")) {
        return(bwFallback)
      }

      purrr::map_dbl(exList, function(ex) {
        xVec <- .getCut(ex)
        xVec <- xVec[is.finite(xVec)]

        if (length(xVec) < 2L) {
          return(bwFallback)
        }

        xVec <- xVec[xVec > min(xVec, na.rm = TRUE)]

        if (length(xVec) < 2L || length(unique(xVec)) < 2L) {
          return(bwFallback)
        }

        bwOut <- .completeChnlSettingsBwCalcOne(
          x = xVec,
          bwMtd = bwMtd,
          bwAdj = bwAdj,
          bwNcellMin = 1e4,
          bwNcellMax = 1e4
        )

        if (!is.finite(bwOut) || bwOut <= 0) {
          return(bwFallback)
        }

        bwOut
      })
    }
  ) |>
    unlist()

  bwVec <- bwVec[is.finite(bwVec) & bwVec > 0]

  if (length(bwVec) == 0L) {
    return(bwFallback)
  }

  mean(bwVec, trim = 0.1, na.rm = TRUE)
}

#' @keywords internal
.completeChnlSettingsCpMin <- function(
  cpMin,
  .data,
  popGate,
  chnlCut,
  indBatchList,
  pathProject
) {
  if (!is.null(cpMin)) {
    return(cpMin)
  }
  .debug("calculating cpMin automatically") # nolint
  purrr::map(
    seq_len(min(5, length(indBatchList))),
    function(i) {
      exList <- .getExList(
        # nolint
        .data = .data,
        indBatch = indBatchList[[i]],
        pop = popGate,
        chnlCut,
        batch = names(indBatchList)[i],
        pathProject = pathProject
      )
      purrr::map_dbl(exList, function(ex) {
        median(.getCut(ex)[.getCut(ex) > min(.getCut(ex))], na.rm = TRUE)[[
          1
        ]] # nolint
      })
    }
  ) |>
    unlist() |>
    mean(trim = 0.1)
}

# Get all cutpoint type names
# Returns character vector of all available cutpoint names
#' @keywords internal
.getFullCpTypeVec <- function(fdr) {
  c(
    "man",
    "tg",
    "dcp",
    "midp",
    "scp",
    "uns",
    "unsr",
    "loc"
  )
}

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

#' @keywords internal
.chnlLab <- function(.data) {
  adf <- switch(
    class(.data)[1],
    "GatingSet" = {
      gh <- .data[[1]]
      fr <- flowWorkspace::gh_pop_get_data(gh)
      flowCore::parameters(fr)@data
    },
    "GatingHierarchy" = {
      fr <- flowWorkspace::gh_pop_get_data(.data)
      flowCore::parameters(fr)@data
    },
    "flowFrame" = flowCore::parameters(.data)@data,
    "flowSet" = flowCore::parameters(.data[[1]])@data,
    "cytoframe" = flowCore::parameters(.data)@data,
    "cytoset" = flowCore::parameters(.data[[1]])@data,
    stop(paste0("Unsupported data class: ", class(.data)[1]))
  )

  labVec <- setNames(adf$desc, adf$name)
  for (i in seq_along(labVec)) {
    if (is.na(labVec[i])) {
      labVec[i] <- names(labVec)[i]
    }
  }

  labVec
}

#' @title Read marker settings from project
#' @description Read the saved marker settings list from the project's metaData folder.
#' @param pathProject character Path to project.
#' @return A named list of marker settings (as saved by .completeChnlSettingsSave()).
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
#' stimgateMetaReadSettingsChnls(tmp)
#' }
#' @export
stimgateMetaReadSettingsChnls <- function(pathProject) {
  pathChnlList <- file.path(pathProject, "metaData", "chnlSettings.rds")
  if (!file.exists(pathChnlList)) {
    stop("Channel list file not found in project metaData folder")
  }
  readRDS(pathChnlList)
}

#' @title Read marker list with channel labels
#' @description Read the project's marker list and return it with element names
#'   replaced by channel labels (from chnlLab).
#' @param pathProject character Path to project.
#' @return A named list of marker settings where names are channel labels.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
#' stimgateMetaReadSettingsChnls(tmp)
#' }
#' @export
stimgateMetaReadSettingsMarkers <- function(pathProject) {
  markerList <- stimgateMetaReadSettingsChnls(pathProject)
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  names(markerList) <- chnlLab[names(markerList)]
  markerList
}

#' @title Get marker settings for a single channel
#' @description Retrieve the marker settings for a single channel. The function
#'   accepts either a channel label (as returned by stimgateMetaReadChnlLab)
#'   or the original channel name/key used in markerList.
#' @param pathProject character Path to project.
#' @param chnl character Channel label or channel name.
#' @return A list of settings for the requested channel.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
#' stimgateMetaReadSettingsChnl(tmp, "BC1 label")
#' stimgateMetaReadSettingsChnl(tmp, "BC1")
#' }
#' @export
stimgateMetaReadSettingsChnl <- function(pathProject, chnl) {
  chnlList <- stimgateMetaReadSettingsChnls(pathProject)
  if (!chnl %in% names(chnlList)) {
    stop(sprintf("Channel %s not found in marker list", chnl))
  }
  chnlList[[chnl]]
}

#' @title Get settings for a named marker
#' @description Retrieve the settings for a marker by its original name/key.
#' @param pathProject character Path to project.
#' @param marker character Marker name/key as stored in markerList.
#' @return A list of settings for the requested marker.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(list(BC1 = list(a = 1)), file.path(tmp, "metaData", "markerList.rds"))
#' stimgateMetaReadSettingsMarker(tmp, "BC1")
#' }
#' @export
stimgateMetaReadSettingsMarker <- function(pathProject, marker) {
  markerList <- stimgateMetaReadSettingsChnls(pathProject)
  if (!marker %in% names(markerList)) {
    stop(sprintf("Marker %s not found in marker list", marker))
  }
  markerList[[marker]]
}

#' @rdname stimgateMetaReadLab
#' @title Read channel or marker label mapping
#' @description Read the saved channel label mapping (chnlLab.rds) from the project's metaData folder.
#' @param pathProject character Path to project.
#' @return Named character vector mapping channel names to labels.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(c(BC1 = "BC1 label"), file.path(tmp, "metaData", "chnlLab.rds"))
#' stimgateMetaReadChnlLab(tmp)
#' }
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
  chnlLab <- .chnlLab(.data)
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

#' @title Read batch list from project
#' @description Read the saved batchList object from the project's metaData folder.
#' @param pathProject character Path to project.
#' @return A list describing sample grouping into batches (as saved by .saveMetaDataBatchList()).
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "metaData"), showWarnings = FALSE)
#' saveRDS(list(batch1 = c(1, 2)), file.path(tmp, "metaData", "batchList.rds"))
#' stimgateMetaReadBatchList(tmp)
#' }
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
      "Otherwise, simply specify `chnl` instead of `marker` ",
      "(and then also `chnlSettings` instead of `markerSettings` if applicable). "
    )
  }
  chnlVec
}

.extractChnlSettings <- function(
  chnlSettings,
  markerSettings,
  chnl,
  marker,
  pathProject
) {
  .verifyChnlSettings(
    chnlSettings = chnlSettings,
    markerSettings = markerSettings,
    chnl = chnl,
    marker = marker
  )
  if (!is.null(chnlSettings)) {
    return(chnlSettings)
  }
  if (is.null(markerSettings)) {
    return(lapply(chnl, function(x) list()) |> stats::setNames(chnl))
  }
  markerLab <- stimgateMetaReadMarkerLab(pathProject)
  chnlSettings <- markerSettings
  names(chnlSettings) <- markerLab[names(markerSettings)]
  for (i in seq_along(chnlSettings)) {
    chnlSettings[[i]]$chnlCut <- markerLab[[chnlSettings[[i]]$markerCut]]
    chnlSettings[[i]]$markerCut <- NULL
  }
  chnlSettings
}
