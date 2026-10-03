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
  "pred",
  "cpOrig",
  "indVec",
  "x1",
  "x",
  "y",
  "countStim",
  "nCellStim",
  "countUns",
  "nCellUns",
  "xVec",
  "freqBs",
  "freqStim",
  "chnlPos",
  "dirSave",
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
  "type",
  "dens",
  "no",
  "yes",
  "probStim",
  "probStimNorm",
  # Variables from other functions
  "cytCombn",
  "freqUns",
  "V1",
  "V2",
  "i",
  "tolGateSingle",
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
#' @return logical TRUE if debug is active and initialized, FALSE otherwise.
#' @keywords internal
.debugInit <- function(pathProject) {
  if (!.profileEnabled()) {
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
      if (dir.exists(dirDebug)) {
        unlink(dirDebug, recursive = TRUE, force = TRUE)
      }
      if (!dir.exists(dirDebug)) {
        dir.create(dirDebug, recursive = TRUE, showWarnings = FALSE)
      }
      pathDebugFile <- file.path(dirDebug, "debug.txt")
      if (!file.create(pathDebugFile, showWarnings = FALSE)) {
        .debugStateReset()
        return(FALSE)
      }
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
  if (!.profileEnabled()) {
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
    saveRDS(dots[[i]], .intSavePathSave(
      pathProject = pathProject,
      stage = stage,
      ind = ind,
      name = dotNames[[i]]
    ))
  }

  invisible(TRUE)
}

#' @keywords internal
.intSaveCheck <- function(ind) {
  if (is.null(ind) || length(ind) == 0L || all(is.na(ind))) {
    return(FALSE)
  }

  envVar <- Sys.getenv("STIMGATE_INTERMEDIATE") |>
    trimws() |>
    tolower()
  if (envVar == "") {
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
