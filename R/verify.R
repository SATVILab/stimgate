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
  biasUns,
  biasUnsFactor,
  excMin,
  cpMin,
  bw,
  bwMin,
  bwMax,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwCluster,
  minCell,
  tolClust,
  maxPosProbX,
  gateCombn,
  gateQuant
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
    !is.character(pathProject) || length(pathProject) != 1 || pathProject == ""
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

  # 3. Global Hyperparameter & Logic Checks
  if (!is.logical(calcCytPosGates) || length(calcCytPosGates) != 1) {
    stop("`calcCytPosGates` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(calcSinglePosGates) || length(calcSinglePosGates) != 1) {
    stop("`calcSinglePosGates` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.logical(excMin) || length(excMin) != 1) {
    stop("`excMin` must be a single logical value (TRUE/FALSE).")
  }
  if (!is.null(biasUns) && (!is.numeric(biasUns) || length(biasUns) != 1)) {
    stop("`biasUns` must be a single numeric value, or NULL.")
  }
  if (
    !is.numeric(biasUnsFactor) ||
      length(biasUnsFactor) != 1 ||
      biasUnsFactor <= 0
  ) {
    stop("`biasUnsFactor` must be a positive single numeric value.")
  }
  if (!is.null(cpMin) && (!is.numeric(cpMin) || length(cpMin) != 1)) {
    stop("`cpMin` must be a single numeric value, or NULL.")
  }
  if (!is.numeric(maxPosProbX) || length(maxPosProbX) != 1) {
    stop("`maxPosProbX` must be a single numeric value.")
  }

  # Bandwidth checks
  if (!is.null(bw) && (!is.numeric(bw) || length(bw) != 1 || bw <= 0)) {
    stop("`bw` must be a single positive numeric value, or NULL.")
  }
  if (
    !is.null(bwMin) && (!is.numeric(bwMin) || length(bwMin) != 1 || bwMin <= 0)
  ) {
    stop("`bwMin` must be a single positive numeric value, or NULL.")
  }
  if (
    !is.null(bwMax) && (!is.numeric(bwMax) || length(bwMax) != 1 || bwMax <= 0)
  ) {
    stop("`bwMax` must be a single positive numeric value, or NULL.")
  }
  if (!is.null(bwMin) && !is.null(bwMax) && bwMax < bwMin) {
    stop("`bwMax` must be greater than or equal to `bwMin`.")
  }
  if (!is.numeric(bwAdj) || length(bwAdj) != 1 || bwAdj <= 0) {
    stop("`bwAdj` must be a single positive numeric multiplier.")
  }

  validBwMtds <- c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")
  if (is.null(bw) && !bwMtd %in% validBwMtds) {
    stop(sprintf(
      "`bwMtd` must be one of: %s",
      paste(validBwMtds, collapse = ", ")
    ))
  }

  # Cell limits
  if (is.null(bw) && !is.null(bwNcellMin) && !is.numeric(bwNcellMin)) {
    stop("`bwNcellMin` must be numeric.")
  }
  if (is.null(bw) && !is.null(bwNcellMax)) {
    if (!is.numeric(bwNcellMax)) {
      stop("`bwNcellMax` must be numeric.")
    }
    if (!is.null(bwNcellMin) && bwNcellMax < bwNcellMin) {
      stop("`bwNcellMax` must be >= `bwNcellMin`.")
    }
  }
  if (!is.numeric(minCell) || length(minCell) != 1 || minCell <= 0) {
    stop("`minCell` must be a positive number.")
  }

  if (
    !is.null(bwCluster) &&
      (!is.numeric(bwCluster) || length(bwCluster) != 1 || bwCluster <= 0)
  ) {
    stop("`bwCluster` must be a single positive numeric value, or NULL.")
  }
  if (
    !is.null(tolClust) &&
      (!is.numeric(tolClust) || length(tolClust) != 1 || tolClust <= 0)
  ) {
    stop("`tolClust` must be a positive numeric threshold.")
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

  # Channel presence
  chnlLab <- .chnlLab(.data)
  if (!is.null(chnl)) {
    if (!all(chnl %in% names(chnlLab))) {
      stop(
        "Channels not found in GatingSet: ",
        paste(setdiff(chnl, names(chnlLab)), collapse = ", ")
      )
    }
  }
  if (!is.null(marker)) {
    if (!all(marker %in% chnlLab)) {
      stop(
        "Markers not found in GatingSet: ",
        paste(setdiff(marker, chnlLab), collapse = ", ")
      )
    }
  }

  invisible(TRUE)
}

#' @keywords internal
.verifyChnlSettings <- function(chnlSettings, chnl, markerSettings, marker) {
  # Added 'bwMin' which was previously missing
  permissibleSettings <- c(
    "biasUns",
    "biasUnsFactor",
    "excMin",
    "cpMin",
    "bw",
    "bwMin",
    "bwMax",
    "bwMtd",
    "bwAdj",
    "bwNcellMin",
    "bwNcellMax",
    "minCell",
    "tolClust",
    "maxPosProbX",
    "bwCluster",
    "popGate",
    "gateCombn",
    "gateQuant"
  )

  if (
    !is.null(chnlSettings) && !all(names(chnlSettings) %in% permissibleSettings)
  ) {
    stop(
      "Invalid channel settings detected: ",
      paste(setdiff(names(chnlSettings), permissibleSettings), collapse = ", ")
    )
  }
  if (
    !is.null(markerSettings) &&
      !all(names(markerSettings) %in% permissibleSettings)
  ) {
    stop(
      "Invalid marker settings detected: ",
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
    chnlVec <- names(chnlSettings)
    if (length(chnlVec) == 0L || length(chnlVec) != length(chnlSettings)) {
      stop("`chnlSettings` elements must be named.")
    }
    if (length(chnlVec) != length(unique(chnlVec))) {
      stop("`chnlSettings` must have unique channel names.")
    }
    if (!all(chnlVec %in% chnl)) {
      stop("All channels in `chnlSettings` must be included in `chnl`")
    }
  } else if (!is.null(markerSettings) && length(markerSettings) > 0L) {
    markerVec <- names(markerSettings)
    if (
      length(markerVec) == 0L || length(markerVec) != length(markerSettings)
    ) {
      stop("`markerSettings` elements must be named.")
    }
    if (length(markerVec) != length(unique(markerVec))) {
      stop("`markerSettings` must have unique marker names.")
    }
    if (!all(markerVec %in% marker)) {
      stop("All markers in `markerSettings` must be included in `marker`")
    }
  }

  if (!is.null(chnlSettings)) {
    purrr::walk(names(chnlSettings), function(chnlCurr) {
      if (!is.list(chnlSettings[[chnlCurr]])) {
        stop(sprintf("Channel '%s' setting must be a list.", chnlCurr))
      }
      .verifyChnlSettingsChnl(chnlCurr, chnlSettings[[chnlCurr]])
    })
  }
  if (!is.null(markerSettings)) {
    purrr::walk(names(markerSettings), function(markerCurr) {
      if (!is.list(markerSettings[[markerCurr]])) {
        stop(sprintf("Marker '%s' setting must be a list.", markerCurr))
      }
      .verifyChnlSettingsChnl(markerCurr, markerSettings[[markerCurr]])
    })
  }
  invisible(TRUE)
}

#' @keywords internal
.verifyChnlSettingsChnl <- function(chnlCurr, settings) {
  prefix <- sprintf("Channel '%s' setting error: ", chnlCurr)

  # Check logical flags
  if (
    !is.null(settings$excMin) &&
      (!is.logical(settings$excMin) || length(settings$excMin) != 1)
  ) {
    stop(paste0(prefix, "`excMin` must be a single logical value."))
  }
  # Check numeric scalars
  if (
    !is.null(settings$biasUns) &&
      (!is.numeric(settings$biasUns) || length(settings$biasUns) != 1)
  ) {
    stop(paste0(prefix, "`biasUns` must be a single numeric value."))
  }
  if (
    !is.null(settings$biasUnsFactor) &&
      (!is.numeric(settings$biasUnsFactor) ||
        length(settings$biasUnsFactor) != 1 ||
        settings$biasUnsFactor <= 0)
  ) {
    stop(paste0(prefix, "`biasUnsFactor` must be a single positive number."))
  }
  if (
    !is.null(settings$bw) &&
      (!is.numeric(settings$bw) || length(settings$bw) != 1 || settings$bw <= 0)
  ) {
    stop(paste0(prefix, "`bw` must be a single positive numeric value."))
  }
  if (
    !is.null(settings$bwMin) &&
      (!is.numeric(settings$bwMin) ||
        length(settings$bwMin) != 1 ||
        settings$bwMin <= 0)
  ) {
    stop(paste0(prefix, "`bwMin` must be a single positive numeric value."))
  }
  if (
    !is.null(settings$bwMax) &&
      (!is.numeric(settings$bwMax) ||
        length(settings$bwMax) != 1 ||
        settings$bwMax <= 0)
  ) {
    stop(paste0(prefix, "`bwMax` must be a single positive numeric value."))
  }
  if (
    !is.null(settings$bwMin) &&
      !is.null(settings$bwMax) &&
      settings$bwMax < settings$bwMin
  ) {
    stop(paste0(prefix, "`bwMax` must be >= `bwMin`."))
  }
  if (
    !is.null(settings$bwAdj) &&
      (!is.numeric(settings$bwAdj) ||
        length(settings$bwAdj) != 1 ||
        settings$bwAdj <= 0)
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
      (!is.numeric(settings$maxPosProbX) || length(settings$maxPosProbX) != 1)
  ) {
    stop(paste0(prefix, "`maxPosProbX` must be a single numeric value."))
  }

  # Check Character strings
  if (
    "popGate" %in%
      names(settings) &&
      (!is.character(settings$popGate) || length(settings$popGate) != 1)
  ) {
    stop(paste0(prefix, "`popGate` must be a single character string."))
  }

  invisible(TRUE)
}
