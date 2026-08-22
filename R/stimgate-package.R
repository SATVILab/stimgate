#' @keywords internal
#' @aliases stimgate-package
#' @useDynLib stimgate, .registration = TRUE
"_PACKAGE"

.simCompareResolveClusterMismatchNames <- function(clusterSpec) {
  if (is.null(clusterSpec) || length(clusterSpec) == 0L) {
    return(character(0))
  }

  spec <- if (length(clusterSpec) == 1L && is.character(clusterSpec)) {
    strsplit(as.character(clusterSpec), ",", fixed = TRUE)[[1]]
  } else {
    as.character(clusterSpec)
  }

  spec <- trimws(spec)
  spec <- spec[nzchar(spec)]
  unique(spec)
}

.simCompareApplyClusterMismatch <- function(
    outListExperiment,
    stimMeanShift = 0,
    stimSdMultiplier = 1,
    stimMeanShiftClusters = NULL,
    stimSdMultiplierClusters = NULL) {
  if (is.null(outListExperiment) || is.null(outListExperiment[["flowFrameList"]])) {
    return(outListExperiment)
  }

  shiftLabels <- .simCompareResolveClusterMismatchNames(stimMeanShiftClusters)
  sdLabels <- .simCompareResolveClusterMismatchNames(stimSdMultiplierClusters)
  if (length(shiftLabels) == 0L && length(sdLabels) == 0L) {
    return(outListExperiment)
  }

  flowFrameList <- outListExperiment[["flowFrameList"]]
  labelsList <- outListExperiment[["labelsList"]]
  nCondition <- if (!is.null(outListExperiment[["nCondition"]])) {
    outListExperiment[["nCondition"]]
  } else {
    2L
  }

  for (idx in seq_along(flowFrameList)) {
    condIndex <- ((idx - 1L) %% nCondition) + 1L
    if (condIndex == 1L) {
      next
    }

    labelVec <- labelsList[[idx]]
    expr <- flowCore::exprs(flowFrameList[[idx]])
    if (ncol(expr) == 0L || length(labelVec) != nrow(expr)) {
      next
    }

    if (length(shiftLabels) > 0L) {
      shiftMask <- labelVec %in% shiftLabels
      if (any(shiftMask)) {
        expr[shiftMask, 1L] <- expr[shiftMask, 1L] + stimMeanShift
        flowCore::exprs(flowFrameList[[idx]]) <- expr
      }
    }

    if (length(sdLabels) > 0L) {
      sdMask <- labelVec %in% sdLabels
      if (any(sdMask)) {
        negVals <- expr[sdMask, 1L]
        center <- mean(negVals)
        expr[sdMask, 1L] <- center + (negVals - center) * stimSdMultiplier
        flowCore::exprs(flowFrameList[[idx]]) <- expr
      }
    }
  }

  outListExperiment[["flowFrameList"]] <- flowFrameList
  outListExperiment
}

.simcyto_compat_on_load <- function() {
  if (!requireNamespace("simcyto", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  ns <- asNamespace("simcyto")
  current <- tryCatch(get("simCytExperiment", mode = "function", envir = ns), error = function(e) NULL)
  if (is.null(current)) {
    return(invisible(FALSE))
  }

  orig_fn <- current
  compat_fn <- function(
      ...,
      stimMeanShift = 0,
      stimSdMultiplier = 1,
      stimMeanShiftClusters = NULL,
      stimSdMultiplierClusters = NULL,
      scenario = NULL) {
    cluster_shift <- stimMeanShiftClusters
    cluster_sd <- stimSdMultiplierClusters
    global_shift <- stimMeanShift
    global_sd <- stimSdMultiplier

    if (!is.null(cluster_shift) || !is.null(cluster_sd)) {
      stimMeanShift <- 0
      stimSdMultiplier <- 1
    }

    call_args <- list(...)
    call_args$stimMeanShift <- stimMeanShift
    call_args$stimSdMultiplier <- stimSdMultiplier
    if (!is.null(cluster_shift)) {
      call_args$stimMeanShiftClusters <- cluster_shift
    }
    if (!is.null(cluster_sd)) {
      call_args$stimSdMultiplierClusters <- cluster_sd
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
      !is.null(cluster_shift) ||
        !is.null(cluster_sd)
    ) {
      if (exists(".simCompareApplyClusterMismatch", mode = "function")) {
        out <- .simCompareApplyClusterMismatch(
          outListExperiment = out,
          stimMeanShift = global_shift,
          stimSdMultiplier = global_sd,
          stimMeanShiftClusters = cluster_shift,
          stimSdMultiplierClusters = cluster_sd
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
