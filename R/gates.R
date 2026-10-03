#' @title Get gates
#'
#' @description Get all the gates for each of the markers gated.
#'
#' @param pathProject character. Path to the project directory.
#' @param pop character. Optional population name(s) to filter gates by. Default is NULL (all populations).
#' @param marker character. Optional marker name(s) to filter gates by. Default is NULL (all markers).
#' @param chnl character. Optional channel name(s) to filter gates by. Default is NULL (all channels).
#'
#' @return Gate table with gates for each sample for each marker.
#' @examples
#' # Get example dataset
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#'
#' # Run the stimgate pipeline
#' pathProject <- gateStim(
#'   pathProject = file.path(tempdir(), "getGateExample"),
#'   .data = gs,
#'   batchList = exampleData$batchList,
#'   marker = exampleData$marker,
#'   popGate = "root"
#' )
#'
#' # Get identified gates
#' gates <- getStimGates(pathProject)
#' @export
getStimGates <- function(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL
) {
  pop <- pop %|c|% .gateGetPop(pathProject)

  markerChnl <- NULL
  if (!is.null(marker) && length(pop) > 0L) {
    markerLab <- stimgateMetaReadMarkerLab(pathProject)
    unknown <- marker[!marker %in% names(markerLab)]
    if (length(unknown) > 0L) {
      stop("Unknown marker: ", paste(unknown, collapse = ", "))
    }
    markerChnl <- unname(markerLab[marker])
  }

  purrr::map_df(pop, function(popCurr) {
    chnlVec <- if (!is.null(marker)) {
      markerChnl
    } else {
      chnl %|c|% .gateGetChnl(pathProject, popCurr)
    }
    chnlLab <- if (length(chnlVec) > 0L) {
      stimgateMetaReadChnlLab(pathProject)
    }

    purrr::map_df(chnlVec, function(chnlCurr) {
      markerCurr <- chnlLab[chnlCurr] |>
        stats::setNames(NULL)

      .gatesGetPathAll(
        pathProject = pathProject,
        pop = popCurr,
        chnlCut = chnlCurr,
        init = FALSE
      ) |>
        readRDS() |>
        dplyr::mutate(marker = markerCurr) |>
        dplyr::mutate(pop = popCurr) |>
        dplyr::select(pop, dplyr::everything())
    })
  })
}

#' @keywords internal
.gateGetPop <- function(pathProject) {
  .gateGetDirs(file.path(pathProject, "gates"), "pop")
}

#' @keywords internal
.gateGetChnl <- function(pathProject, pop) {
  .gateGetDirs(file.path(pathProject, "gates", paste0("pop", pop)), "chnl")
}

#' @keywords internal
.gateGetDirs <- function(pathDir, prefix) {
  dirVec <- list.dirs(pathDir, full.names = FALSE, recursive = FALSE)
  dirVec <- dirVec[nzchar(dirVec) & startsWith(dirVec, prefix)]
  unique(sub(paste0("^", prefix), "", dirVec))
}

#' @keywords internal
.gatesGetPathAll <- function(pathProject, pop, chnlCut, init) {
  file.path(
    pathProject,
    "gates",
    paste0("pop", pop),
    paste0("chnl", chnlCut),
    "all",
    if (init) "gateTblInit.rds" else "gateTbl.rds"
  )
}


#' @title Get detailed gate diagnostics
#'
#' @description
#' Read detailed local-FDR threshold diagnostics saved during gating, including
#' condition-level, sample-level and final cluster-level thresholds with the
#' corresponding background-subtracted frequencies.
#'
#' @param pathProject character. Path to the project directory.
#' @param pop character. Optional population name(s) to retain. Population is
#'   currently recorded as `NA` for intermediate diagnostics.
#' @param marker character. Optional marker name(s) to retain.
#' @param chnl character. Optional channel name(s) to retain.
#' @param save logical. If TRUE, save the detailed table as an RDS file.
#' @param pathSave character. Optional path for the saved RDS file. Defaults to
#'   `file.path(pathProject, "gatesDetailed.rds")`.
#'
#' @return A tibble with one row per saved threshold diagnostic.
#' @export
getStimGatesDetailed <- function(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL,
  save = FALSE,
  pathSave = NULL
) {
  detailTbl <- .gateGetDetailedIntermediate(pathProject)

  if (nrow(detailTbl) > 0L) {
    if (!"chnl" %in% names(detailTbl)) {
      detailTbl$chnl <- detailTbl$detailPathChnl
    } else if ("detailPathChnl" %in% names(detailTbl)) {
      detailTbl$chnl <- dplyr::coalesce(
        detailTbl$chnl, detailTbl$detailPathChnl
      )
    }

    chnlLab <- try(stimgateMetaReadChnlLab(pathProject), silent = TRUE)
    if (!inherits(chnlLab, "try-error") && length(chnlLab) > 0L) {
      detailTbl$marker <- unname(chnlLab[detailTbl$chnl])
    } else {
      detailTbl$marker <- NA_character_
    }
    detailTbl$pop <- NA_character_

    if (!is.null(chnl)) {
      detailTbl <- detailTbl |>
        dplyr::filter(.data$chnl %in% .env$chnl)
    }
    if (!is.null(marker)) {
      detailTbl <- detailTbl |>
        dplyr::filter(.data$marker %in% .env$marker)
    }
    detailTbl <- detailTbl |>
      dplyr::select(
        pop,
        marker,
        chnl,
        dplyr::everything()
      )
  }

  if (isTRUE(save)) {
    pathSave <- pathSave %||% file.path(pathProject, "gatesDetailed.rds")
    saveRDS(detailTbl, pathSave)
  }

  detailTbl
}

#' @keywords internal
.gateGetDetailedIntermediate <- function(pathProject) {
  pathInt <- file.path(pathProject, "intermediateData")
  if (!dir.exists(pathInt)) {
    return(tibble::tibble())
  }

  pathVec <- list.files(
    pathInt,
    pattern = "^(locDetail.*|locClusterQuantileTbl)\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(pathVec) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_df(pathVec, function(pathCurr) {
    obj <- try(readRDS(pathCurr), silent = TRUE)
    if (inherits(obj, "try-error") || !is.data.frame(obj)) {
      return(NULL)
    }
    meta <- .gateGetDetailedPathMeta(pathCurr, pathInt)
    obj <- .gateGetDetailedNormaliseObject(obj, meta$detailObject)
    obj |>
      dplyr::mutate(
        detailObject = meta$detailObject,
        detailPathStage = meta$stage,
        detailPathChnl = meta$chnl,
        detailPathInd = meta$ind,
        detailSourceFile = pathCurr
      )
  })
}


#' @keywords internal
.gateGetDetailedNormaliseObject <- function(obj, detailObject) {
  if (detailObject == "locClusterQuantileTbl") {
    if (!"detailLevel" %in% names(obj)) {
      obj$detailLevel <- "cluster_final"
    }
    if (!"threshold" %in% names(obj) && "cpJoinTgOrig" %in% names(obj)) {
      obj$threshold <- obj$cpJoinTgOrig
    }
  }
  if (detailObject == "locDetailClusterFinal") {
    if (!"detailLevel" %in% names(obj)) {
      obj$detailLevel <- "cluster_final"
    }
    columns <- c(
      threshold = "locFinalThreshold",
      thresholdOrigin = "locFinalThresholdOrigin",
      nCellStim = "locFinalNCellStim",
      nCellUns = "locFinalNCellUns",
      propStim = "locFinalPropStim",
      propUns = "locFinalPropUns",
      propBs = "locFinalPropBs"
    )
    for (column in names(columns)) {
      source <- columns[[column]]
      if (!column %in% names(obj) && source %in% names(obj)) {
        obj[[column]] <- obj[[source]]
      }
    }
  }
  obj
}

#' @keywords internal
.gateGetDetailedPathMeta <- function(pathCurr, pathInt) {
  pathCurrNorm <- normalizePath(pathCurr, winslash = "/", mustWork = FALSE)
  pathIntNorm <- normalizePath(pathInt, winslash = "/", mustWork = FALSE)
  rel <- substring(pathCurrNorm, nchar(pathIntNorm) + 2L)
  parts <- strsplit(dirname(rel), "/", fixed = TRUE)[[1]]
  if (identical(parts, ".")) {
    parts <- character(0)
  }
  detailObject <- sub("\\.rds$", "", basename(pathCurr))
  list(
    stage = if (length(parts) >= 1L) parts[[1]] else NA_character_,
    chnl = if (length(parts) >= 2L) parts[[2]] else NA_character_,
    ind = if (length(parts) >= 4L && parts[[3]] == "ind") {
      parts[[4]]
    } else {
      NA_character_
    },
    detailObject = detailObject
  )
}
