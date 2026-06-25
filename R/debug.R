#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom  utils head

# Global variable bindings to avoid R CMD check notes
# These are primarily used in dplyr and ggplot2 contexts
globalVariables(c(
  # Variables used in dplyr operations
  "marker",
  "batch",
  "ind",
  "gate",
  "gateCyt",
  "gateSingle",
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
  "gateConcat"
))

#' Create a debug file in tempdir()
#'
#' @return character Path to the created debug file (invisibly).
#' @keywords internal
.debugFileCreate <- function() {
  dirPath <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dirPath)) {
    dir.create(dirPath, recursive = TRUE, showWarnings = FALSE)
  }
  timestamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
  filePath <- file.path(dirPath, paste0("stimgate_", timestamp, ".txt"))
  file.create(filePath, showWarnings = FALSE)
  invisible(filePath)
}

#' Get the most recent debug file path
#'
#' @return character Path to the most recent debug file, or NULL if none.
#' @keywords internal
.debugFileGetPath <- function() {
  dirPath <- file.path(tempdir(), "stimgate")
  if (!dir.exists(dirPath)) {
    return(NULL)
  }
  files <- list.files(
    dirPath,
    pattern = "^stimgate_\\d{4}-\\d{2}-\\d{2}-\\d{6}\\.txt$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(NULL)
  }
  files[which.max(file.info(files)$mtime)]
}

#' Print debug message conditionally
#'
#' Writes debug output to the most recent stimgate debug file when
#' STIMGATE_DEBUG is enabled.
#'
#' @param msg character Message to print
#' @param val object Optional value to append to message. Default is NULL.
#' @return logical invisibly TRUE if message was written, FALSE otherwise
#' @keywords internal
.debug <- function(msg, val = NULL) {
  mustDebug <- tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")
  if (!mustDebug) {
    return(invisible(FALSE))
  }
  if (!is.null(val)) {
    msg <- paste0(msg, ": ", val)
  }
  pathDebug <- .debugFileGetPath()
  if (is.null(pathDebug)) {
    pathDebug <- .debugFileCreate()
  }
  cat(msg, file = pathDebug, sep = "\n", append = TRUE)
  invisible(TRUE)
}

#' Copy the latest stimgate debug file to the working directory
#'
#' @description Copies the most recent debug file created by stimgate
#' to the current working directory. The copied file uses the same
#' filename as the source.
#' @param pathDir character. Directory to copy the debug file to.
#' Default is `getwd()` (i.e. the working directory).
#'
#' @return character Path to the copied file (invisibly), or NULL if no
#'   debug file exists.
#' @export
stimgate_debug_copy <- function(pathDir = getwd()) {
  src <- .debugFileGetPath()
  if (is.null(src)) {
    message("No stimgate debug file found.")
    return(invisible(NULL))
  }
  if (!dir.exists(pathDir)) {
    dir.create(pathDir, recursive = TRUE)
  }
  dest <- file.path(pathDir, basename(src))
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) {
    message("Failed to copy debug file to working directory.")
    return(invisible(NULL))
  }
  invisible(dest)
}

#' Print the latest stimgate debug file to the console
#'
#' @description Prints the most recent debug file created by stimgate
#' to the console.
#'
#' @return character The debug text (invisibly), or NULL if no debug file exists.
#' @export
stimgate_debug_print <- function() {
  src <- .debugFileGetPath()
  if (is.null(src)) {
    message("No stimgate debug file found.")
    return(invisible(NULL))
  }
  txt <- readLines(src, warn = FALSE)
  cat(txt, sep = "\n")
  invisible(txt)
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
    "intermediate_data",
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
