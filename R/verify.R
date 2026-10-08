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
  biasUns,
  bw,
  control,
  markerControl
) {
  # 1. Control class check: no other argument is read before this succeeds.
  if (!inherits(control, "stimControl")) {
    stop("`control` must be a stimControl object; use `stimControl()`.")
  }

  # 2. Channel / Marker Mutual Exclusivity Checks
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Specify only one of 'chnl' or 'marker', not both.")
  }
  if (is.null(chnl) && is.null(marker)) {
    stop("Must specify one of 'chnl' or 'marker'.")
  }

  # 3. Structural & Type Checks
  if (
    !is.character(pathProject) || length(pathProject) != 1 || pathProject == ""
  ) {
    stop(
      "`pathProject` must be a single, non-empty character string specifying a directory."
    )
  }
  if (!inherits(.data, "GatingSet")) {
    stop("`.data` must be a GatingSet after input conversion.")
  }
  if (!is.list(batchList) || length(batchList) == 0) {
    stop("`batchList` must be a non-empty list of sample indices or names.")
  }

  # 4. Global-only checks
  if (!is.null(markerControl) && !is.list(markerControl)) {
    stop("`markerControl` must be NULL or a named list.")
  }
  .verifyChnlSettingsChnl(
    settings = list(popGate = popGate, biasUns = biasUns, bw = bw),
    prefix = ""
  )

  # Channel presence
  chnlLab <- chnlLab(.data)
  if (!is.null(chnl)) {
    if (!all(chnl %in% names(chnlLab))) {
      stop(
        "Channels not found in GatingSet: ",
        paste(setdiff(chnl, names(chnlLab)), collapse = ", ")
      )
    }
    chnlSel <- chnl
  }
  if (!is.null(marker)) {
    if (!all(marker %in% chnlLab)) {
      stop(
        "Markers not found in GatingSet: ",
        paste(setdiff(marker, chnlLab), collapse = ", ")
      )
    }
    chnlSel <- names(chnlLab)[chnlLab %in% marker]
  }

  # Validate `markerControl` here, before any project directory is created. The
  # resolved settings are recomputed in `.completeChnlSettings()`.
  .resolveMarkerControl(
    markerControl = markerControl,
    chnl = chnlSel,
    chnlLab = chnlLab
  )

  invisible(TRUE)
}

.verifyBwMtdsOrdinary <- c("nrd0", "sj", "hpi0", "hpi1", "hpi2", "hpi3")
.verifyBwMtds <- c(.verifyBwMtdsOrdinary, paste0(.verifyBwMtdsOrdinary, "Norm"))
.verifyGateCombns <- c("no", "min", "median", "max", "prejoin")

#' @keywords internal
.verifyChnlSettingsChnl <- function(
  chnlCurr,
  settings,
  prefix = sprintf("Channel '%s' setting error: ", chnlCurr)
) {
  if (
    !.verifyIsNullOrNa(settings[["excMin"]]) &&
      (!is.logical(settings[["excMin"]]) || length(settings[["excMin"]]) != 1)
  ) {
    stop(paste0(prefix, "`excMin` must be a single logical value."))
  }
  for (nm in c("biasUns", "cpMin")) {
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
  for (nm in c("bw", "bwCluster")) {
    .check_positive_n(nm, settings = settings, prefix = prefix)
  }

  if (
    !.verifyIsNullOrNa(settings[["minCell"]]) &&
      (!is.numeric(settings[["minCell"]]) ||
        length(settings[["minCell"]]) != 1L ||
        settings[["minCell"]] <= 0)
  ) {
    stop(paste0(prefix, "`minCell` must be a positive number."))
  }

  # The bandwidth-selector settings are only used when no fixed bandwidth is
  # supplied for this channel.
  if (.verifyIsNullOrNa(settings[["bw"]])) {
    bwNcellMin <- settings[["bwNcellMin"]]
    bwNcellMax <- settings[["bwNcellMax"]]
    if (!.verifyIsNullOrNa(bwNcellMin) && !is.numeric(bwNcellMin)) {
      stop(paste0(prefix, "`bwNcellMin` must be numeric."))
    }
    if (!.verifyIsNullOrNa(bwNcellMax)) {
      if (!is.numeric(bwNcellMax)) {
        stop(paste0(prefix, "`bwNcellMax` must be numeric."))
      }
      if (
        !.verifyIsNullOrNa(bwNcellMin) &&
          is.numeric(bwNcellMin) &&
          length(bwNcellMin) == 1L &&
          is.finite(bwNcellMin) &&
          is.finite(bwNcellMax) &&
          bwNcellMax < bwNcellMin
      ) {
        stop(paste0(prefix, "`bwNcellMax` must be >= `bwNcellMin`."))
      }
    }
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
    .verifyIsNullOrNa(settings[["bw"]]) &&
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

  if (
    "clusterGates" %in% names(settings) &&
      (!is.logical(settings[["clusterGates"]]) ||
        length(settings[["clusterGates"]]) != 1L ||
        is.na(settings[["clusterGates"]]))
  ) {
    stop(paste0(prefix, "`clusterGates` must be TRUE or FALSE."))
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

  .check_prob("normPeakMinRel", settings = settings, prefix = prefix)
  .check_prob("normExtraFrac", settings = settings, prefix = prefix)

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
    !.verifyIsNullOrNa(settings[["bwScope"]]) &&
      (!is.character(settings[["bwScope"]]) ||
        length(settings[["bwScope"]]) != 1L ||
        !settings[["bwScope"]] %in% c("cytokine", "cluster", "sample"))
  ) {
    stop(paste0(
      prefix,
      "`bwScope` must be one of 'cytokine', 'cluster' or 'sample'."
    ))
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

  if (
    !.verifyIsNullOrNa(settings[["locThresholdMethod"]]) &&
      (!is.character(settings[["locThresholdMethod"]]) ||
        length(settings[["locThresholdMethod"]]) != 1L ||
        !settings[["locThresholdMethod"]] %in% c("region", "match"))
  ) {
    stop(paste0(
      prefix, "`locThresholdMethod` must be either 'region' or 'match'."
    ))
  }

  probSettings <- c(
    "locMinPeakProb",
    "locDipAlpha",
    "locAntimodeHeightFrac",
    "locAntimodeLowRel",
    "locAntimodeLowAbs",
    "locFlatDerivFrac",
    "locFlatHardDerivFrac",
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

  invisible(TRUE)
}


.verifyBwLimitSetting <- function(
  x,
  nm,
  allow_none = TRUE,
  allow_neg = FALSE,
  allow_inf = FALSE,
  prefix = ""
) {
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
  prefix = ""
) {
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

#' @keywords internal
.verifyBatchList <- function(batchList) {
  isValidBatch <- vapply(batchList, function(x) {
    (is.numeric(x) || is.character(x)) && length(x) >= 2L && !anyNA(x)
  }, logical(1))
  if (!all(isValidBatch)) {
    stop(
      "Each `batchList` element must hold at least two non-missing sample ",
      "indices: the unstimulated sample first, then the stimulated samples."
    )
  }
  # an unstim may be shared across batches, but must be first in each;
  # a stim sample must appear exactly once
  indUns <- unique(as.character(lapply(batchList, `[[`, 1L)))
  indStim <- unlist(lapply(batchList, function(x) as.character(x[-1])))
  indStimUns <- intersect(indStim, indUns)
  if (length(indStimUns) > 0L) {
    stop(
      "Sample(s) ", paste0(indStimUns, collapse = ", "),
      " are both unstimulated (first in a `batchList` element) and ",
      "stimulated."
    )
  }
  indStimDup <- unique(indStim[duplicated(indStim)])
  if (length(indStimDup) > 0L) {
    stop(
      "Stimulated sample(s) ", paste0(indStimDup, collapse = ", "),
      " appear more than once in `batchList`."
    )
  }
  invisible(TRUE)
}
