# Global variable bindings to avoid R CMD check notes
# These are primarily used in dplyr and ggplot2 contexts
globalVariables(c(
  # Variables used in dplyr operations
  "marker",
  "batch",
  "ind",
  "gate",
  "gateCyt",
  "gateName",
  "chnl",
  "gateUse",
  "gateType",
  "gateCombn",
  "popGate",
  "gateTbl",
  "chnlCut",
  "tol",
  "indInBatchGate",
  "tolClustSingle",
  "indBatchGate",
  "cp",
  "grp",
  "cpJoinLseOrigMeanTg",
  "cpOrigQuantMin",
  "cpJoin",
  "cpJoinLse",
  "cpJoinLseOrig",
  "cpJoinLseOrigMean",
  "cpJoinTgOrig",
  "cpJoinTgOrigMean",
  "propBsOrig",
  "propBsCpDiff",
  "propBsCpDiffSd",
  "propBsCp",
  "propL1se",
  "pred",
  "der",
  "cpOrig",
  "maxExpr",
  "gate05",
  "propBsCpDiffSdMax",
  "grpLevel",
  "indVec",
  "x1",
  "x",
  "y",
  "xInd",
  "countStim",
  "nCellStim",
  "countUns",
  "nCellUns",
  "propStimPos",
  "propUnsPos",
  "propStimSd",
  "propUnsSd",
  "x512",
  "xVec",
  "freqBs",
  "freqStim",
  "popGateCurr",
  "cpJoinTg",
  "lseOrig",
  "cpTgCtrl",
  "chnlPos",
  "dirSave",
  "isNullGateTbl",
  "pathProject",
  "excMin",
  "propBsDiff",
  "propStim",
  "propUns",
  "propBs",
  "probSmooth",
  "nRow",
  "yStim",
  "yUns",
  "stim",
  "xStim",
  "prob",
  "xUns",
  "propLab",
  "type",
  "lineId",
  "dens",
  "no",
  "yes",
  "probStim",
  "probStimNorm",
  "propPos",
  # Variables from .get_cp_uns_loc_prob_tbl_filter
  "minorResponseInd",
  "moderateResponseInd",
  "nRemaining",
  "probLargerCount",
  "probLargerProp",
  # Variables from .get_prop_bs_by_cp_tbl_ind_calc
  "countStimCp",
  "countUnsCp",
  "propStimCp",
  "propUnsCp",
  "propBsSd",
  "propStimPosCp",
  "propUnsPosCp",
  "propStimSdCp",
  "propUnsSdCp",
  "propBsSdCp",
  # Variables from other functions
  "cytCombn",
  "freqUns",
  "V1",
  "V2",
  "i",
  "tolGateSingle",
  # Variables used in plots and ggplot2 context
  ".debug",
  "rbeta",
  # Variables from fcs_write.R
  "concat",
  "gateConcat",
  # Additional variables for R CMD check
  ".env",
  "everything",
  "approx",
  "as.formula",
  "binomial",
  "density",
  "glm",
  "median",
  "optim",
  "predict",
  "quantile",
  "read.csv",
  "rnorm",
  "sd",
  "setNames",
  "locGeneratedDirect",
  "calc_skew",
  "calc_gamma",
  "_stimgate_stimgate_cpPmden"
))

.debugState <- new.env(parent = emptyenv())
.debugState$file <- NULL
.debugState$initialized <- FALSE

#' Reset internal debug state
#'
#' @return invisible(NULL)
#' @keywords internal
.debugStateReset <- function() {
  .debugState$file <- NULL
  .debugState$initialized <- FALSE
  invisible(NULL)
}

#' Initialise textual debug state for a StimGate run
#'
#' @param pathProject character Path to project directory.
#' @param reset logical Whether to reset existing debug directory.
#'   Default: TRUE.
#' @return logical TRUE if debug is active and initialized, FALSE otherwise.
#' @keywords internal
.debugInit <- function(pathProject, reset = TRUE) {
  mustDebug <- tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")
  if (!mustDebug) {
    .debugStateReset()
    return(FALSE)
  }
  if (
    !is.character(pathProject) ||
      length(pathProject) != 1L ||
      is.na(pathProject) ||
      !nzchar(pathProject)
  ) {
    .debugStateReset()
    return(FALSE)
  }
  tryCatch(
    {
      dirDebug <- file.path(pathProject, "debug")
      if (isTRUE(reset) && dir.exists(dirDebug)) {
        unlink(dirDebug, recursive = TRUE, force = TRUE)
      }
      if (!dir.exists(dirDebug)) {
        dir.create(dirDebug, recursive = TRUE, showWarnings = FALSE)
      }
      pathDebugFile <- file.path(dirDebug, "debug.txt")
      file.create(pathDebugFile, showWarnings = FALSE)
      .debugState$file <- pathDebugFile
      .debugState$initialized <- TRUE
      TRUE
    },
    error = function(e) {
      .debugStateReset()
      FALSE
    }
  )
}

#' Print debug message conditionally
#'
#' Writes debug output directly to pathProject/debug/debug.txt when
#' STIMGATE_DEBUG is enabled.
#'
#' @param msg character Message to print.
#' @param val object Optional value to append to message. Default: NULL.
#' @return logical invisibly TRUE if message was written, FALSE otherwise.
#' @keywords internal
.debug <- function(msg, val = NULL) {
  mustDebug <- tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")
  if (!mustDebug) {
    return(invisible(FALSE))
  }
  tryCatch(
    {
      if (!is.null(val)) {
        msg <- paste0(msg, ": ", val)
      }
      pathDebug <- .debugState$file
      if (is.null(pathDebug) || !is.character(pathDebug) || !nzchar(pathDebug)) {
        return(invisible(FALSE))
      }
      if (!dir.exists(dirname(pathDebug))) {
        dir.create(dirname(pathDebug), recursive = TRUE, showWarnings = FALSE)
      }
      cat(msg, file = pathDebug, sep = "\n", append = TRUE)
      invisible(TRUE)
    },
    error = function(e) {
      invisible(FALSE)
    }
  )
}

#' @keywords internal
.intSaveNm <- function(name, obj, ind, stage, pathProject) {
  if (!.intSaveCheck(ind)) {
    return(invisible(FALSE))
  }
  pathSave <- .intSavePathSave(
    pathProject = pathProject,
    stage = stage,
    ind = ind,
    name = name
  )
  saveRDS(obj, pathSave)
  invisible(TRUE)
}

#' @keywords internal
.intSave <- function(ind, stage, pathProject, ...) {
  if (!.intSaveCheck(ind)) {
    return(invisible(FALSE))
  }

  dots <- list(...)
  dotNames <- names(dots)

  callNames <- as.list(substitute(list(...)))[-1]
  callNames <- vapply(
    callNames,
    function(x) paste(deparse(x), collapse = ""),
    character(1)
  )

  if (is.null(dotNames)) {
    dotNames <- callNames
  } else {
    dotNames[dotNames == ""] <- callNames[dotNames == ""]
  }

  for (i in seq_along(dots)) {
    .intSaveNm(
      name = dotNames[[i]],
      obj = dots[[i]],
      ind = ind,
      stage = stage,
      pathProject = pathProject
    )
  }

  invisible(TRUE)
}

#' @keywords internal
.isInvalidInd <- function(ind) {
  is.null(ind) || length(ind) == 0 || all(is.na(ind))
}

#' @keywords internal
.intSaveCheck <- function(ind) {
  if (.isInvalidInd(ind)) {
    return(FALSE)
  }

  envVar <- Sys.getenv("STIMGATE_INTERMEDIATE") |>
    trimws() |>
    tolower()
  if (is.null(envVar) || length(envVar) == 0 || envVar == "") {
    return(FALSE)
  }
  if (envVar %in% c("y", "true", "yes", "all")) {
    return(TRUE)
  }
  envVarSplit <- strsplit(envVar, ",|;") |>
    unlist() |>
    trimws()
  any(as.character(ind) %in% envVarSplit)
}

#' @keywords internal
.intSavePathSave <- function(pathProject, stage, ind, name) {
  name <- paste0(name, ".rds")
  name <- gsub("\\.rds(\\.rds)*$", ".rds", name, ignore.case = TRUE)
  pathSave <- file.path(
    pathProject,
    "intermediateData",
    stage,
    "ind",
    paste0(as.character(ind), collapse = "_"),
    name
  )
  if (!dir.exists(dirname(pathSave))) {
    dir.create(dirname(pathSave), recursive = TRUE, showWarnings = FALSE)
  }
  pathSave
}

#' @keywords internal
.browse <- function(ind) {
  if (!.browseCheck(ind)) {
    return(invisible(FALSE))
  }
  eval(quote(browser()), envir = parent.frame())
  invisible(TRUE)
}

#' @keywords internal
.browseCheck <- function(ind) {
  if (.isInvalidInd(ind)) {
    return(FALSE)
  }

  envVar <- Sys.getenv("STIMGATE_BROWSE") |>
    trimws() |>
    tolower()
  if (is.null(envVar) || length(envVar) == 0 || envVar == "") {
    return(FALSE)
  }
  if (envVar %in% c("y", "true", "yes", "all")) {
    return(TRUE)
  }
  envVarSplit <- strsplit(envVar, ",|;") |>
    unlist() |>
    trimws()
  any(as.character(ind) %in% envVarSplit)
}
