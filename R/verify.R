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
  bwFallback,
  bwMtd,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwCluster,
  bwAdaptive,
  bwAdaptiveDensityN,
  bwAdaptivePadFrac,
  normPeakFrac,
  normPeakMinRel,
  normExtraFrac,
  normExtraMax,
  normExtraJitterFrac,
  normLambda,
  normDensityN,
  normExcessBwMtd,
  normExcessNcell,
  normAdaptiveNcell,
  normMtd,
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

  .verifyBwLimit <- function(x, nm, allow_none = TRUE) {
    if (is.null(x)) {
      return(invisible(TRUE))
    }
    if (is.character(x) && length(x) == 1L) {
      validChar <- if (allow_none) c("auto", "none") else "auto"
      if (tolower(x) %in% validChar) {
        return(invisible(TRUE))
      }
    }
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(sprintf(
        "`%s` must be a single positive numeric value%s.",
        nm,
        if (allow_none) ", `auto`, `none`, or NULL" else " or `auto`"
      ))
    }
    invisible(TRUE)
  }

  .verifyBwLimit(bwMin, "bwMin", allow_none = TRUE)
  .verifyBwLimit(bwMax, "bwMax", allow_none = TRUE)
  .verifyBwLimit(bwFallback, "bwFallback", allow_none = FALSE)

  if (is.numeric(bwMin) && is.numeric(bwMax) && bwMax < bwMin) {
    stop("`bwMax` must be greater than or equal to `bwMin`.")
  }
  if (!is.numeric(bwAdj) || length(bwAdj) != 1 || bwAdj <= 0) {
    stop("`bwAdj` must be a single positive numeric multiplier.")
  }

  validBwMtds <- c(
    "nrd0",
    "sj",
    "hpi0",
    "hpi1",
    "hpi2",
    "hpi3",
    "nrd0Norm",
    "sjNorm",
    "hpi0Norm",
    "hpi1Norm",
    "hpi2Norm",
    "hpi3Norm"
  )
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
    stop("`tolClust` must be a positive numeric value, or NULL.")
  }

  .verifyNormBwSettings(
    settings = list(
      bwAdaptive = bwAdaptive,
      bwAdaptiveDensityN = bwAdaptiveDensityN,
      bwAdaptivePadFrac = bwAdaptivePadFrac,
      normPeakFrac = normPeakFrac,
      normPeakMinRel = normPeakMinRel,
      normExtraFrac = normExtraFrac,
      normExtraMax = normExtraMax,
      normExtraJitterFrac = normExtraJitterFrac,
      normLambda = normLambda,
      normDensityN = normDensityN,
      normExcessBwMtd = normExcessBwMtd,
      normExcessNcell = normExcessNcell,
      normAdaptiveNcell = normAdaptiveNcell,
      normMtd = normMtd
    ),
    prefix = ""
  )

  .verifyLocSettings(
    settings = list(
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
      locTolRefPeak = locTolRefPeak
    ),
    prefix = ""
  )

  validGateCombns <- c("no", "min", "median", "max", "prejoin")
  if (
    !is.character(gateCombn) ||
      length(gateCombn) == 0L ||
      !all(gateCombn %in% validGateCombns)
  ) {
    stop(sprintf(
      "`gateCombn` must contain only: %s",
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
    "bwFallback",
    "bwMtd",
    "bwAdj",
    "bwNcellMin",
    "bwNcellMax",
    "minCell",
    "tolClust",
    "locProbCol",
    "locMinPeakProb",
    "locDipAlpha",
    "locAntimodeHeightFrac",
    "locAntimodeLowRel",
    "locAntimodeLowAbs",
    "locFlatDerivFrac",
    "locFlatHardDerivFrac",
    "locLeftLowRel",
    "locLeftLowAbs",
    "locLeftCellFrac",
    "locLeftLengthFrac",
    "locMarginalPurityRel",
    "locMarginalCellBinRatio",
    "locMarginalRefQuantile",
    "locTolRefPeak",
    "maxPosProbX",
    "bwCluster",
    "bwAdaptive",
    "bwAdaptiveDensityN",
    "bwAdaptivePadFrac",
    "normPeakFrac",
    "normPeakMinRel",
    "normExtraFrac",
    "normExtraMax",
    "normExtraJitterFrac",
    "normLambda",
    "normDensityN",
    "normExcessBwMtd",
    "normExcessNcell",
    "normAdaptiveNcell",
    "normMtd",
    "popGate",
    "gateCombn",
    "gateQuant"
  )

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
      settingsCurr <- chnlSettings[[chnlCurr]]
      if (!is.list(settingsCurr)) {
        stop(sprintf("Channel '%s' setting must be a list.", chnlCurr))
      }
      invalidSettings <- setdiff(names(settingsCurr), permissibleSettings)
      if (length(invalidSettings) > 0L) {
        stop(
          sprintf("Invalid settings for channel '%s': ", chnlCurr),
          paste(invalidSettings, collapse = ", ")
        )
      }
      .verifyChnlSettingsChnl(chnlCurr, settingsCurr)
    })
  }
  if (!is.null(markerSettings)) {
    purrr::walk(names(markerSettings), function(markerCurr) {
      settingsCurr <- markerSettings[[markerCurr]]
      if (!is.list(settingsCurr)) {
        stop(sprintf("Marker '%s' setting must be a list.", markerCurr))
      }
      invalidSettings <- setdiff(names(settingsCurr), permissibleSettings)
      if (length(invalidSettings) > 0L) {
        stop(
          sprintf("Invalid settings for marker '%s': ", markerCurr),
          paste(invalidSettings, collapse = ", ")
        )
      }
      .verifyChnlSettingsChnl(markerCurr, settingsCurr)
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

  .verifyBwLimitSetting <- function(x, nm, allow_none = TRUE) {
    if (is.null(x)) {
      return(invisible(TRUE))
    }
    if (is.character(x) && length(x) == 1L) {
      validChar <- if (allow_none) c("auto", "none") else "auto"
      if (tolower(x) %in% validChar) {
        return(invisible(TRUE))
      }
    }
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(paste0(
        prefix,
        "`",
        nm,
        "` must be a positive numeric value",
        if (allow_none) ", `auto`, `none`, or NULL." else " or `auto`."
      ))
    }
    invisible(TRUE)
  }

  .verifyBwLimitSetting(settings$bwMin, "bwMin", allow_none = TRUE)
  .verifyBwLimitSetting(settings$bwMax, "bwMax", allow_none = TRUE)
  .verifyBwLimitSetting(settings$bwFallback, "bwFallback", allow_none = FALSE)

  # Check Character strings
  if (
    "popGate" %in%
      names(settings) &&
      (!is.character(settings$popGate) || length(settings$popGate) != 1)
  ) {
    stop(paste0(prefix, "`popGate` must be a single character string."))
  }

  validBwMtds <- c(
    "nrd0",
    "sj",
    "hpi0",
    "hpi1",
    "hpi2",
    "hpi3",
    "nrd0Norm",
    "sjNorm",
    "hpi0Norm",
    "hpi1Norm",
    "hpi2Norm",
    "hpi3Norm"
  )
  if (
    !is.null(settings$bwMtd) &&
      (!is.character(settings$bwMtd) ||
        length(settings$bwMtd) != 1 ||
        !settings$bwMtd %in% validBwMtds)
  ) {
    stop(paste0(
      prefix,
      "`bwMtd` must be one of: ",
      paste(validBwMtds, collapse = ", "),
      "."
    ))
  }

  if (
    !is.null(settings$bwCluster) &&
      (!is.numeric(settings$bwCluster) ||
        length(settings$bwCluster) != 1 ||
        settings$bwCluster <= 0)
  ) {
    stop(paste0(prefix, "`bwCluster` must be a single positive numeric value."))
  }

  if (
    !is.null(settings$tolClust) &&
      (!is.numeric(settings$tolClust) ||
        length(settings$tolClust) != 1 ||
        settings$tolClust <= 0)
  ) {
    stop(paste0(
      prefix,
      "`tolClust` must be a positive numeric value, or NULL."
    ))
  }

  validGateCombns <- c("no", "min", "median", "max", "prejoin")
  if (
    !is.null(settings$gateCombn) &&
      (!is.character(settings$gateCombn) ||
        length(settings$gateCombn) == 0L ||
        !all(settings$gateCombn %in% validGateCombns))
  ) {
    stop(paste0(
      prefix,
      "`gateCombn` must contain only: ",
      paste(validGateCombns, collapse = ", "),
      "."
    ))
  }

  if (
    !is.null(settings$gateQuant) &&
      (!is.numeric(settings$gateQuant) ||
        length(settings$gateQuant) != 2 ||
        any(settings$gateQuant < 0 | settings$gateQuant > 1))
  ) {
    stop(paste0(
      prefix,
      "`gateQuant` must be two probabilities between 0 and 1."
    ))
  }

  .verifyNormBwSettings(settings = settings, prefix = prefix)

  .verifyLocSettings(settings = settings, prefix = prefix)

  invisible(TRUE)
}


#' @keywords internal
.verifyNormBwSettings <- function(settings, prefix = "") {
  if (
    !is.null(settings$bwAdaptive) &&
      (!is.logical(settings$bwAdaptive) || length(settings$bwAdaptive) != 1L)
  ) {
    stop(paste0(prefix, "`bwAdaptive` must be a single logical value."))
  }

  .check_positive_n <- function(nm, allow_null = TRUE, allow_inf = FALSE) {
    val <- settings[[nm]]
    if (is.null(val) && isTRUE(allow_null)) {
      return(invisible(TRUE))
    }
    if (
      !is.numeric(val) ||
        length(val) != 1L ||
        (!allow_inf && !is.finite(val)) ||
        (allow_inf && !is.finite(val) && !is.infinite(val)) ||
        val <= 0
    ) {
      stop(paste0(
        prefix,
        "`",
        nm,
        "` must be a single positive numeric value."
      ))
    }
    invisible(TRUE)
  }

  .check_prob <- function(nm) {
    val <- settings[[nm]]
    if (is.null(val)) {
      return(invisible(TRUE))
    }
    if (
      !is.numeric(val) ||
        length(val) != 1L ||
        !is.finite(val) ||
        val < 0 ||
        val > 1
    ) {
      stop(paste0(
        prefix,
        "`",
        nm,
        "` must be a single numeric value between 0 and 1."
      ))
    }
    invisible(TRUE)
  }

  .check_prob("normPeakFrac")
  .check_prob("normPeakMinRel")
  .check_prob("normExtraFrac")
  .check_prob("normExtraJitterFrac")

  .check_positive_n("normExtraMax", allow_inf = TRUE)
  .check_positive_n("normDensityN")
  .check_positive_n("normExcessNcell")
  .check_positive_n("normAdaptiveNcell")
  .check_positive_n("bwAdaptiveDensityN")

  if (!is.null(settings$bwAdaptivePadFrac)) {
    val <- settings$bwAdaptivePadFrac
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val < 0) {
      stop(paste0(
        prefix,
        "`bwAdaptivePadFrac` must be a single non-negative numeric value."
      ))
    }
  }

  if (!is.null(settings$normLambda)) {
    val <- settings$normLambda
    if (!is.numeric(val) || length(val) == 0L || any(!is.finite(val))) {
      stop(paste0(prefix, "`normLambda` must be a finite numeric vector."))
    }
  }

  if (!is.null(settings$normExcessBwMtd)) {
    validOrdinaryBwMtds <- c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")
    if (
      !is.character(settings$normExcessBwMtd) ||
        length(settings$normExcessBwMtd) != 1L ||
        !settings$normExcessBwMtd %in% validOrdinaryBwMtds
    ) {
      stop(paste0(
        prefix,
        "`normExcessBwMtd` must be one of: ",
        paste(validOrdinaryBwMtds, collapse = ", "),
        "."
      ))
    }
  }

  if (!is.null(settings$normMtd)) {
    if (
      !is.character(settings$normMtd) ||
        length(settings$normMtd) != 1L ||
        !settings$normMtd %in% c("moments", "boxcox")
    ) {
      stop(paste0(prefix, "`normMtd` must be either 'moments' or 'boxcox'."))
    }
  }

  if (isTRUE(settings$bwAdaptive) && identical(settings$normMtd, "boxcox")) {
    stop(paste0(
      prefix,
      "`bwAdaptive = TRUE` currently requires `normMtd = 'moments'`."
    ))
  }

  invisible(TRUE)
}

#' @keywords internal
.verifyLocSettings <- function(settings, prefix = "") {
  if (
    !is.null(settings$locProbCol) &&
      (!is.character(settings$locProbCol) ||
        length(settings$locProbCol) != 1 ||
        !settings$locProbCol %in% c("pred", "probSmooth"))
  ) {
    stop(paste0(prefix, "`locProbCol` must be either 'pred' or 'probSmooth'."))
  }

  probSettings <- c(
    "locMinPeakProb",
    "locDipAlpha",
    "locAntimodeHeightFrac",
    "locAntimodeLowRel",
    "locAntimodeLowAbs",
    "locFlatDerivFrac",
    "locFlatHardDerivFrac",
    "locLeftLowRel",
    "locLeftLowAbs",
    "locLeftCellFrac",
    "locLeftLengthFrac",
    "locMarginalPurityRel",
    "locMarginalRefQuantile"
  )
  purrr::walk(probSettings, function(nm) {
    val <- settings[[nm]]
    if (is.null(val)) {
      return(invisible(TRUE))
    }
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop(paste0(prefix, "`", nm, "` must be a single finite numeric value."))
    }
    if (val < 0 || val > 1) {
      stop(paste0(prefix, "`", nm, "` must be between 0 and 1."))
    }
    invisible(TRUE)
  })

  if (!is.null(settings$locMarginalCellBinRatio)) {
    val <- settings$locMarginalCellBinRatio
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0) {
      stop(paste0(
        prefix,
        "`locMarginalCellBinRatio` must be a single positive finite numeric value."
      ))
    }
  }

  if (
    !is.null(settings$locTolRefPeak) &&
      (!is.character(settings$locTolRefPeak) ||
        length(settings$locTolRefPeak) != 1 ||
        !settings$locTolRefPeak %in% c("highest", "first"))
  ) {
    stop(paste0(prefix, "`locTolRefPeak` must be either 'highest' or 'first'."))
  }

  invisible(TRUE)
}
