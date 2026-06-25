#' @keywords internal
.verifyGateInputs <- function(
  pathProject,
  .data,
  batchList,
  popGate,
  chnl,
  marker,
  chnlSettings,
  markerSettings,
  calcCytPosGates,
  calcSinglePosGates,
  bwMtd,
  gateCombn,
  gateQuant,
  bwNcellMin,
  bwNcellMax,
  bwCluster,
  minCell,
  tolClust
) {
  # 1. Channel / Marker Mutual Exclusivity Checks
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Specify only one of 'chnl' or 'marker', not both.")
  }
  if (is.null(chnl) && is.null(marker)) {
    stop("Must specify one of 'chnl' or 'marker'.")
  }
  if (!is.null(chnl) && !is.null(markerSettings)) {
    stop("When 'chnl' is specified, 'markerSettings' must be NULL.")
  }
  if (!is.null(marker) && !is.null(chnlSettings)) {
    stop("When 'marker' is specified, 'chnlSettings' must be NULL.")
  }

  # 2. Structural & Type Checks
  if (
    !is.character(pathProject) ||
      length(pathProject) != 1 ||
      pathProject == ""
  ) {
    stop(
      "`pathProject` must be a single, non-empty character string specifying a directory."
    )
  }
  if (!is.character(popGate) || length(popGate) != 1) {
    stop("`popGate` must be a single character string (e.g., 'root').")
  }
  if (
    !inherits(
      .data,
      c(
        "GatingSet",
        "GatingHierarchy",
        "flowFrame",
        "flowSet",
        "cytoframe",
        "cytoset"
      )
    )
  ) {
    stop(
      "`.data` must be a valid flow core/workspace object (e.g., GatingSet, flowFrame)."
    )
  }
  if (!is.list(batchList) || length(batchList) == 0) {
    stop("`batchList` must be a non-empty list of sample indices.")
  }

  # 3. Logical & Hyperparameter Checks
  if (!is.logical(calcCytPosGates) || length(calcCytPosGates) != 1) {
    stop("`calcCytPosGates` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(calcSinglePosGates) || length(calcSinglePosGates) != 1) {
    stop("`calcSinglePosGates` must be a single logical value (TRUE/FALSE).")
  }

  validBwMtds <- c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")
  if (!bwMtd %in% validBwMtds) {
    stop(sprintf(
      "`bwMtd` must be one of: %s",
      paste(validBwMtds, collapse = ", ")
    ))
  }

  validGateCombns <- c("min", "median", "max")
  if (!gateCombn %in% validGateCombns) {
    stop(sprintf(
      "`gateCombn` must be one of: %s",
      paste(validGateCombns, collapse = ", ")
    ))
  }

  if (
    !is.numeric(gateQuant) ||
      length(gateQuant) != 2 ||
      any(gateQuant < 0 | gateQuant > 1)
  ) {
    stop(
      "`gateQuant` must be a numeric vector of two probabilities between 0 and 1."
    )
  }
  if (!is.numeric(bwNcellMin) || bwNcellMin <= 0) {
    stop("`bwNcellMin` must be a positive number.")
  }
  if (!is.numeric(bwNcellMax) || bwNcellMax < bwNcellMin) {
    stop(
      "`bwNcellMax` must be a positive number greater than or equal to `bwNcellMin`."
    )
  }
  if (!is.numeric(minCell) || minCell <= 0) {
    stop("`minCell` must be a positive number.")
  }
  if (!is.numeric(tolClust) || tolClust <= 0) {
    stop("`tolClust` must be a positive numeric threshold.")
  }
  if (!is.null(bwCluster) && (!is.numeric(bwCluster) || bwCluster <= 0)) {
    stop("`bwCluster` must be a positive numeric value if specified.")
  }
  chnlLab <- .chnlLab(.data)
  if (!is.null(chnl)) {
    allChnlValid <- all(chnl %in% names(chnlLab))
    if (!allChnlValid) {
      stop(
        "The following channels are not found in the GatingSet: ",
        paste(setdiff(chnl, names(chnlLab)), collapse = ", ")
      )
    }
  }
  if (!is.null(marker)) {
    allMarkerValid <- all(marker %in% chnlLab)
    if (!allMarkerValid) {
      stop(
        "The following markers are not found in the GatingSet: ",
        paste(setdiff(marker, chnlLab), collapse = ", ")
      )
    }
  }

  invisible(TRUE)
}

#' @keywords internal
.verifyChnlSettings <- function(chnlSettings, chnl, markerSettings, marker) {
  permissibleSettings <- c(
    "biasUns",
    "biasUnsFactor",
    "excMin",
    "cpMin",
    "bw",
    "bwMax",
    "bwMtd",
    "bwAdj",
    "bwNcellMin",
    "bwNcellMax",
    "minCell",
    "tolClust",
    "maxPosProbX",
    "bwCluster"
  )
  if (
    !is.null(chnlSettings) && !all(names(chnlSettings) %in% permissibleSettings)
  ) {
    stop(
      "Invalid channel settings detected. The invalid settings are: ",
      paste(setdiff(names(chnlSettings), permissibleSettings), collapse = ", ")
    )
  }
  if (
    !is.null(markerSettings) &&
      !all(names(markerSettings) %in% permissibleSettings)
  ) {
    stop(
      "Invalid marker settings detected. The invalid settings are: ",
      paste(
        setdiff(names(markerSettings), permissibleSettings),
        collapse = ", "
      )
    )
  }
  if (!is.null(chnlSettings) && !is.null(markerSettings)) {
    stop("Specify only one of `chnlSettings` or `markerSettings`, not both.")
  }
  if (!is.null(chnlSettings) && !is.list(chnlSettings)) {
    stop("`chnlSettings` must be a list of channel-specific settings.")
  }
  if (!is.null(markerSettings) && !is.list(markerSettings)) {
    stop("`markerSettings` must be a list of marker-specific settings.")
  }
  if (!is.null(chnlSettings) && length(chnlSettings) > 0L) {
    chnlVecFromSettings <- names(chnlSettings)
    if (length(chnlVecFromSettings) == 0L) {
      stop("`chnlSettings` must have named elements corresponding to channels.")
    }
    if (!length(chnlVecFromSettings) == length(unique(chnlVecFromSettings))) {
      stop("`chnlSettings` must have unique channel names.")
    }
    if (!length(chnlVecFromSettings) == length(chnlSettings)) {
      stop("`chnlSettings` must have a name for each element.")
    }
    if (!all(chnlVecFromSettings %in% chnl)) {
      stop(
        "All channels specified in `chnlSettings` must be included in `chnl`"
      )
    }
  } else if (!is.null(markerSettings) && length(markerSettings) > 0L) {
    markerVecFromSettings <- names(markerSettings)
    if (length(markerVecFromSettings) == 0L) {
      stop(
        "`markerSettings` must have named elements corresponding to markers."
      )
    }
    if (
      !length(markerVecFromSettings) == length(unique(markerVecFromSettings))
    ) {
      stop("`markerSettings` must have unique marker names.")
    }
    if (!length(markerVecFromSettings) == length(markerSettings)) {
      stop("`markerSettings` must have a name for each element.")
    }
    if (!all(markerVecFromSettings %in% marker)) {
      stop(
        "All markers specified in `markerSettings` must be included in `marker`"
      )
    }
  }
  if (is.null(chnlSettings) && is.null(markerSettings)) {
    return(invisible(TRUE))
  }
  if (is.list(chnlSettings) && length(chnlSettings) == 0L) {
    return(invisible(TRUE))
  }
  if (is.list(markerSettings) && length(markerSettings) == 0L) {
    return(invisible(TRUE))
  }
  if (!is.null(chnlSettings)) {
    lapply(
      names(chnlSettings),
      function(chnlCurr) {
        .verifyChnlSettingsChnl(chnlCurr, chnlSettings[[chnlCurr]])
      }
    )
  }
  if (!is.null(markerSettings)) {
    lapply(
      names(markerSettings),
      function(markerCurr) {
        .verifyChnlSettingsChnl(markerCurr, markerSettings[[markerCurr]])
      }
    )
  }
  invisible(TRUE)
}

#' @keywords internal
.verifyChnlSettingsChnl <- function(chnlCurr, settings) {
  prefix <- sprintf("Channel '%s' setting error: ", chnlCurr)

  # Validate logic flags if present
  if (
    !is.null(settings$excMin) &&
      (!is.logical(settings$excMin) || length(settings$excMin) != 1)
  ) {
    stop(paste0(prefix, "`excMin` must be a single logical value."))
  }

  # Validate numeric limits and hyper-parameters
  if (
    !is.null(settings$biasUns) &&
      (!is.numeric(settings$biasUns) || length(settings$biasUns) != 1)
  ) {
    stop(paste0(prefix, "`biasUns` must be a single numeric value."))
  }
  if (
    !is.null(settings$biasUnsFactor) &&
      (!is.numeric(settings$biasUnsFactor) || settings$biasUnsFactor <= 0)
  ) {
    stop(paste0(prefix, "`biasUnsFactor` must be a positive number."))
  }
  if (!is.null(settings$bw) && (!is.numeric(settings$bw) || settings$bw <= 0)) {
    stop(paste0(prefix, "`bw` must be a positive numeric value."))
  }
  if (
    !is.null(settings$bwAdj) &&
      (!is.numeric(settings$bwAdj) || settings$bwAdj <= 0)
  ) {
    stop(paste0(prefix, "`bwAdj` must be a positive numeric multiplier."))
  }
  if (
    !is.null(settings$cpMin) &&
      (!is.numeric(settings$cpMin) || length(settings$cpMin) != 1)
  ) {
    stop(paste0(prefix, "`cpMin` must be a single numeric value."))
  }
  if (
    !is.null(settings$maxPosProbX) &&
      (!is.numeric(settings$maxPosProbX) ||
        length(settings$maxPosProbX) != 1)
  ) {
    stop(paste0(prefix, "`maxPosProbX` must be a single numeric value."))
  }
  if (!"popGate" %in% names(settings)) {
    stop(paste0(prefix, "`popGate` must be specified in the settings."))
  }
  if (!is.character(settings$popGate) || length(settings$popGate) != 1) {
    stop(paste0(prefix, "`popGate` must be a single character string."))
  }

  invisible(TRUE)
}
