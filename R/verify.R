#' @keywords internal
.verifyIsNullOrNa <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  }
  if (length(x) != 1L || !is.atomic(x)) {
    return(FALSE)
  }
  isTRUE(is.na(x))
}

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
    bwAdaptiveCore,
    bwAdaptiveExtra,
    bwAdaptiveCrossover,
    bwAdaptiveTransitionWidth,
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
    gateQuant) {
  # Snapshot of all arguments, taken before any local variable exists.
  settings <- as.list(environment())

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

  # 3. Global-only checks
  if (!is.logical(calcCytPosGates) || length(calcCytPosGates) != 1) {
    stop("`calcCytPosGates` must be a single logical value (TRUE/FALSE).")
  }
  if (
    .verifyIsNullOrNa(bw) &&
      !.verifyIsNullOrNa(bwNcellMin) &&
      !is.numeric(bwNcellMin)
  ) {
    stop("`bwNcellMin` must be numeric.")
  }
  if (.verifyIsNullOrNa(bw) && !.verifyIsNullOrNa(bwNcellMax)) {
    if (!is.numeric(bwNcellMax)) {
      stop("`bwNcellMax` must be numeric.")
    }
    if (
      !.verifyIsNullOrNa(bwNcellMin) &&
        is.numeric(bwNcellMin) &&
        length(bwNcellMin) == 1L &&
        is.finite(bwNcellMin) &&
        is.finite(bwNcellMax) &&
        bwNcellMax < bwNcellMin
    ) {
      stop("`bwNcellMax` must be >= `bwNcellMin`.")
    }
  }
  if (!is.numeric(minCell) || length(minCell) != 1 || minCell <= 0) {
    stop("`minCell` must be a positive number.")
  }

  # 4. Settings shared with per-channel validation. Per channel, NULL means
  # "inherit the global value", so the global value must itself be supplied.
  # `bwMtd` is only used (and so only checked) when `bw` is not fixed.
  if (!.verifyIsNullOrNa(bw)) {
    settings[["bwMtd"]] <- NULL
  }
  required <- c(
    "excMin", "biasUnsFactor", "maxPosProbX", "bwAdj", "gateCombn",
    "gateQuant", if (.verifyIsNullOrNa(bw)) "bwMtd"
  )
  isMissing <- vapply(settings[required], .verifyIsNullOrNa, logical(1))
  if (any(isMissing)) {
    stop(
      "Must be supplied (not NULL or NA): ",
      paste0("`", required[isMissing], "`", collapse = ", ")
    )
  }
  .verifyChnlSettingsChnl(settings = settings, prefix = "")

  # Channel presence
  chnlLab <- chnlLab(.data)
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
  # `.verifyGateInputs()` (the only caller's precondition) guarantees that at
  # most one of `chnlSettings` and `markerSettings` is supplied.
  isChnl <- !is.null(chnlSettings)
  settingsList <- if (isChnl) chnlSettings else markerSettings
  if (is.null(settingsList)) {
    return(invisible(TRUE))
  }
  arg <- if (isChnl) "chnlSettings" else "markerSettings"
  type <- if (isChnl) "channel" else "marker"
  allowed <- if (isChnl) chnl else marker

  if (!is.list(settingsList)) {
    stop(sprintf("`%s` must be a list of %s-specific settings.", arg, type))
  }
  nms <- names(settingsList)
  if (length(settingsList) > 0L) {
    if (is.null(nms)) {
      stop(sprintf("`%s` elements must be named.", arg))
    }
    if (anyDuplicated(nms) > 0L) {
      stop(sprintf("`%s` must have unique %s names.", arg, type))
    }
    if (!all(nms %in% allowed)) {
      stop(sprintf(
        "All %ss in `%s` must be included in `%s`",
        type, arg, if (isChnl) "chnl" else "marker"
      ))
    }
  }

  # Every per-channel-capable `gateStim()` argument. `locEnforceShapeThreshold`
  # is legacy and deliberately global only.
  permissibleSettings <- setdiff(
    names(formals(gateStim)),
    c(
      "pathProject", ".data", "batchList", "chnl", "marker", "chnlSettings",
      "markerSettings", "calcCytPosGates", "locEnforceShapeThreshold"
    )
  )
  purrr::walk(nms, function(nm) {
    settingsCurr <- settingsList[[nm]]
    if (!is.list(settingsCurr)) {
      stop(sprintf("%s '%s' setting must be a list.", type, nm))
    }
    invalidSettings <- setdiff(names(settingsCurr), permissibleSettings)
    if (length(invalidSettings) > 0L) {
      stop(
        sprintf("Invalid settings for %s '%s': ", type, nm),
        paste(invalidSettings, collapse = ", ")
      )
    }
    .verifyChnlSettingsChnl(nm, settingsCurr)
  })
  invisible(TRUE)
}

.verifyBwMtdsOrdinary <- c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")
.verifyBwMtds <- c(.verifyBwMtdsOrdinary, paste0(.verifyBwMtdsOrdinary, "Norm"))
.verifyGateCombns <- c("no", "min", "median", "max", "prejoin")

#' @keywords internal
.verifyChnlSettingsChnl <- function(
    chnlCurr,
    settings,
    prefix = sprintf("Channel '%s' setting error: ", chnlCurr)) {
  if (
    !.verifyIsNullOrNa(settings[["excMin"]]) &&
      (!is.logical(settings[["excMin"]]) || length(settings[["excMin"]]) != 1)
  ) {
    stop(paste0(prefix, "`excMin` must be a single logical value."))
  }
  for (nm in c("biasUns", "cpMin", "maxPosProbX")) {
    val <- settings[[nm]]
    if (!.verifyIsNullOrNa(val) && (!is.numeric(val) || length(val) != 1)) {
      stop(paste0(prefix, "`", nm, "` must be a single numeric value."))
    }
  }
  for (nm in c("biasUnsFactor", "bwAdj")) {
    .check_positive_n(
      nm,
      allow_inf = TRUE, settings = settings, prefix = prefix
    )
  }
  for (nm in c("bw", "bwCluster", "tolClust")) {
    .check_positive_n(nm, settings = settings, prefix = prefix)
  }

  bwMin <- settings[["bwMin"]]
  bwMax <- settings[["bwMax"]]
  .verifyBwLimitSetting(
    bwMin, "bwMin",
    allow_neg = TRUE, allow_inf = TRUE, prefix = prefix
  )
  .verifyBwLimitSetting(bwMax, "bwMax", allow_inf = TRUE, prefix = prefix)
  .verifyBwLimitSetting(
    settings[["bwFallback"]], "bwFallback",
    allow_none = FALSE, prefix = prefix
  )
  if (
    is.numeric(bwMin) && is.numeric(bwMax) && all(is.finite(c(bwMin, bwMax))) &&
      bwMax < bwMin
  ) {
    stop(paste0(prefix, "`bwMax` must be >= `bwMin`."))
  }

  if (
    "popGate" %in% names(settings) &&
      (!is.character(settings[["popGate"]]) ||
        length(settings[["popGate"]]) != 1)
  ) {
    stop(paste0(prefix, "`popGate` must be a single character string."))
  }

  bwMtd <- settings[["bwMtd"]]
  if (
    !.verifyIsNullOrNa(bwMtd) &&
      (!is.character(bwMtd) || length(bwMtd) != 1 || !bwMtd %in% .verifyBwMtds)
  ) {
    stop(paste0(
      prefix, "`bwMtd` must be one of: ",
      paste(.verifyBwMtds, collapse = ", "), "."
    ))
  }

  gateCombn <- settings[["gateCombn"]]
  if (
    !.verifyIsNullOrNa(gateCombn) &&
      (!is.character(gateCombn) ||
        length(gateCombn) == 0L ||
        !all(gateCombn %in% .verifyGateCombns))
  ) {
    stop(paste0(
      prefix, "`gateCombn` must contain only: ",
      paste(.verifyGateCombns, collapse = ", "), "."
    ))
  }

  gateQuant <- settings[["gateQuant"]]
  if (
    !.verifyIsNullOrNa(gateQuant) &&
      (!is.numeric(gateQuant) ||
        length(gateQuant) != 2 ||
        any(gateQuant < 0 | gateQuant > 1))
  ) {
    stop(paste0(
      prefix, "`gateQuant` must be two probabilities between 0 and 1."
    ))
  }

  .verifyNormBwSettings(settings = settings, prefix = prefix)
  .verifyLocSettings(settings = settings, prefix = prefix)

  invisible(TRUE)
}

#' @keywords internal
.verifyNormBwSettings <- function(settings, prefix = "") {
  if (
    !.verifyIsNullOrNa(settings[["bwAdaptive"]]) &&
      (!is.logical(settings[["bwAdaptive"]]) ||
        length(settings[["bwAdaptive"]]) != 1L)
  ) {
    stop(paste0(prefix, "`bwAdaptive` must be a single logical value."))
  }

  .check_prob("normPeakFrac", settings = settings, prefix = prefix)
  .check_prob("normPeakMinRel", settings = settings, prefix = prefix)
  .check_prob("normExtraFrac", settings = settings, prefix = prefix)
  .check_prob("normExtraJitterFrac", settings = settings, prefix = prefix)

  .check_positive_n(
    "normExtraMax",
    allow_inf = TRUE,
    settings = settings,
    prefix = prefix
  )
  .check_positive_n("normDensityN", settings = settings, prefix = prefix)
  .check_positive_n("normExcessNcell", settings = settings, prefix = prefix)
  .check_positive_n("normAdaptiveNcell", settings = settings, prefix = prefix)
  .check_positive_n("bwAdaptiveDensityN", settings = settings, prefix = prefix)
  .check_positive_n(
    "bwAdaptivePadFrac",
    allow_zero = TRUE,
    settings = settings,
    prefix = prefix
  )

  .check_positive_n("bwAdaptiveCore", settings = settings, prefix = prefix)
  .check_positive_n("bwAdaptiveExtra", settings = settings, prefix = prefix)

  if (!.verifyIsNullOrNa(settings[["bwAdaptiveCrossover"]])) {
    val <- settings[["bwAdaptiveCrossover"]]
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop(paste0(
        prefix,
        "`bwAdaptiveCrossover` must be a single finite numeric value or NULL."
      ))
    }
  }

  .check_positive_n(
    "bwAdaptiveTransitionWidth",
    allow_zero = TRUE,
    settings = settings,
    prefix = prefix
  )

  if (!.verifyIsNullOrNa(settings[["normLambda"]])) {
    val <- settings[["normLambda"]]
    if (!is.numeric(val) || length(val) == 0L || any(!is.finite(val))) {
      stop(paste0(prefix, "`normLambda` must be a finite numeric vector."))
    }
  }

  if (!.verifyIsNullOrNa(settings[["normExcessBwMtd"]])) {
    if (
      !is.character(settings[["normExcessBwMtd"]]) ||
        length(settings[["normExcessBwMtd"]]) != 1L ||
        !settings[["normExcessBwMtd"]] %in% .verifyBwMtdsOrdinary
    ) {
      stop(paste0(
        prefix,
        "`normExcessBwMtd` must be one of: ",
        paste(.verifyBwMtdsOrdinary, collapse = ", "),
        "."
      ))
    }
  }

  if (!.verifyIsNullOrNa(settings[["normMtd"]])) {
    if (
      !is.character(settings[["normMtd"]]) ||
        length(settings[["normMtd"]]) != 1L ||
        !settings[["normMtd"]] %in% c("moments", "boxcox")
    ) {
      stop(paste0(prefix, "`normMtd` must be either 'moments' or 'boxcox'."))
    }
  }

  if (
    isTRUE(settings[["bwAdaptive"]]) &&
      identical(settings[["normMtd"]], "boxcox")
  ) {
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
    !.verifyIsNullOrNa(settings[["locProbCol"]]) &&
      (!is.character(settings[["locProbCol"]]) ||
        length(settings[["locProbCol"]]) != 1 ||
        !settings[["locProbCol"]] %in% c("pred", "probSmooth"))
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
    if (.verifyIsNullOrNa(val)) {
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

  .check_positive_n(
    "locMarginalCellBinRatio",
    settings = settings,
    prefix = prefix
  )

  if (
    !.verifyIsNullOrNa(settings[["locTolRefPeak"]]) &&
      (!is.character(settings[["locTolRefPeak"]]) ||
        length(settings[["locTolRefPeak"]]) != 1 ||
        !settings[["locTolRefPeak"]] %in% c("highest", "first"))
  ) {
    stop(paste0(prefix, "`locTolRefPeak` must be either 'highest' or 'first'."))
  }

  invisible(TRUE)
}


.verifyBwLimitSetting <- function(
    x,
    nm,
    allow_none = TRUE,
    allow_neg = FALSE,
    allow_inf = FALSE,
    prefix = "") {
  if (.verifyIsNullOrNa(x)) {
    return(invisible(TRUE))
  }
  if (
    is.character(x) &&
      length(x) == 1L &&
      tolower(x) %in% c("auto", if (allow_none) "none")
  ) {
    return(invisible(TRUE))
  }
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      (!allow_inf && is.infinite(x)) ||
      (!allow_neg && x <= 0)
  ) {
    stop(paste0(
      prefix,
      "`",
      nm,
      "` must be a single ",
      if (allow_neg) "" else "positive ",
      "numeric value",
      if (allow_none) ", `auto`, `none`, or NULL." else " or `auto`."
    ))
  }
  invisible(TRUE)
}

.check_positive_n <- function(
    nm,
    allow_inf = FALSE,
    allow_zero = FALSE,
    settings,
    prefix = "") {
  val <- settings[[nm]]
  if (.verifyIsNullOrNa(val)) {
    return(invisible(TRUE))
  }
  if (
    !is.numeric(val) ||
      length(val) != 1L ||
      (!allow_inf && !is.finite(val)) ||
      val < 0 ||
      (!allow_zero && val == 0)
  ) {
    stop(paste0(
      prefix,
      "`",
      nm,
      "` must be a single ",
      if (allow_zero) "non-negative" else "positive",
      " numeric value."
    ))
  }
  invisible(TRUE)
}

.check_prob <- function(nm, settings, prefix = "") {
  val <- settings[[nm]]
  if (.verifyIsNullOrNa(val)) {
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
