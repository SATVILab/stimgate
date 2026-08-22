#' @keywords internal
#' @aliases stimgate-package
#' @useDynLib stimgate, .registration = TRUE
"_PACKAGE"

.simcyto_compat_on_load <- function() {
  if (!requireNamespace("simcyto", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  ns <- asNamespace("simcyto")
  current <- tryCatch(get("simCytExperiment", mode = "function", envir = ns), error = function(e) NULL)
  if (is.null(current)) {
    return(invisible(FALSE))
  }

  formals_current <- names(formals(current))
  if (
    "stimMeanShiftClusters" %in% formals_current &&
      "stimSdMultiplierClusters" %in% formals_current
  ) {
    return(invisible(TRUE))
  }

  orig_fn <- current
  compat_fn <- function(
      ...,
      stimMeanShift = 0,
      stimSdMultiplier = 1,
      stimMeanShiftClusters = NULL,
      stimSdMultiplierClusters = NULL,
      scenario = NULL) {
    call_args <- list(...)
    call_args$stimMeanShift <- stimMeanShift
    call_args$stimSdMultiplier <- stimSdMultiplier
    if (!is.null(stimMeanShiftClusters)) {
      call_args$stimMeanShiftClusters <- stimMeanShiftClusters
    }
    if (!is.null(stimSdMultiplierClusters)) {
      call_args$stimSdMultiplierClusters <- stimSdMultiplierClusters
    }
    if (!is.null(scenario)) {
      call_args$scenario <- scenario
    }

    out <- tryCatch(
      do.call(orig_fn, call_args),
      error = function(e) {
        msg <- conditionMessage(e)
        if (
          grepl("unused argument|unused arguments|unexpected argument|did not match", msg, ignore.case = TRUE)
        ) {
          call_args$stimMeanShiftClusters <- NULL
          call_args$stimSdMultiplierClusters <- NULL
          return(do.call(orig_fn, call_args))
        }
        stop(e)
      }
    )

    if (!is.null(out[["flowFrameList"]])) {
      out[["flowFrameList"]] <- unname(out[["flowFrameList"]])
    }
    if (!is.null(out[["labelsList"]])) {
      out[["labelsList"]] <- unname(out[["labelsList"]])
    }

    if (
      !is.null(stimMeanShiftClusters) ||
        !is.null(stimSdMultiplierClusters)
    ) {
      if (exists(".simCompareApplyClusterMismatch", mode = "function")) {
        out <- .simCompareApplyClusterMismatch(
          outListExperiment = out,
          stimMeanShift = stimMeanShift,
          stimSdMultiplier = stimSdMultiplier,
          stimMeanShiftClusters = stimMeanShiftClusters,
          stimSdMultiplierClusters = stimSdMultiplierClusters
        )
      }
    }

    out
  }

  assignInNamespace("simCytExperiment", compat_fn, ns = "simcyto")
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  .simcyto_compat_on_load()
  invisible(TRUE)
}
