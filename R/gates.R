#' @title Read stimulation gates
#' @description Read final gates saved by [gateStim()], optionally selecting
#'   populations and markers. For threshold diagnostics, use [getStimGatesDetailed()].
#' @param pathProject character Project directory from [gateStim()].
#' @param pop character or NULL Populations to retain; NULL selects all.
#'   Default: NULL.
#' @param marker character or NULL Marker labels to retain; takes precedence
#'   over `chnl`. Default: NULL (all markers).
#' @param chnl character or NULL Channels to retain. Default: NULL (all channels).
#' @return A tibble of stimulated-sample gates with identifiers `pop`, `marker`,
#'   `chnl`, `batch`, `ind`, `gateName`, threshold `gate`, and refinement and
#'   threshold-provenance columns when available.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' getStimGates(pathProject)
#' @export
getStimGates <- function(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL
) {
  pop <- as.character(pop %||% .gateGetPop(pathProject))

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
      as.character(chnl %||% .gateGetChnl(pathProject, popCurr))
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
.gateGetChnlPopMap <- function(pathProject) {
  popMap <- character(0)
  pops <- try(.gateGetPop(pathProject), silent = TRUE)
  if (!inherits(pops, "try-error") && length(pops) > 0L) {
    for (p in pops) {
      chnls <- try(.gateGetChnl(pathProject, p), silent = TRUE)
      if (!inherits(chnls, "try-error") && length(chnls) > 0L) {
        popMap[chnls] <- p
      }
    }
  }
  if (length(popMap) == 0L) {
    settings <- try(stimgateMetaReadSettingsChnls(pathProject), silent = TRUE)
    if (!inherits(settings, "try-error") && is.list(settings)) {
      for (s in settings) {
        if (!is.null(s$chnlCut) && !is.null(s$popGate)) {
          popMap[s$chnlCut] <- s$popGate
        }
      }
    }
  }
  popMap
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


#' @title Read threshold diagnostics
#' @description Read saved condition, sample and cluster threshold diagnostics
#'   with background-subtracted frequencies. Use [getStimGates()] for final gates.
#' @param pathProject character Project directory from [gateStim()].
#' @param pop character or NULL Populations to retain. Default: NULL (all).
#' @param marker character or NULL Marker labels to retain. Default: NULL (all).
#' @param chnl character or NULL Channels to retain. Default: NULL (all).
#' @param save logical Save the returned table as RDS. Default: FALSE.
#' @param pathSave character or NULL Output path; NULL uses
#'   `file.path(pathProject, "gatesDetailed.rds")`. Default: NULL.
#' @details
#' Enable diagnostic saving with `Sys.setenv(STIMGATE_INTERMEDIATE = "all")`
#' before running [gateStim()]. Without saved diagnostics, returns an empty tibble.
#' @return A tibble with one row per saved diagnostic, including `pop`, `marker`,
#'   `chnl`, threshold and frequency columns, and source-file metadata.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' Sys.setenv(STIMGATE_INTERMEDIATE = "all")
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' Sys.unsetenv("STIMGATE_INTERMEDIATE")
#' getStimGatesDetailed(pathProject)
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

    if (!"pop" %in% names(detailTbl)) {
      detailTbl$pop <- NA_character_
    }
    if ("popGate" %in% names(detailTbl)) {
      detailTbl$pop <- dplyr::coalesce(detailTbl$pop, detailTbl$popGate)
    }
    if ("detailPathPop" %in% names(detailTbl)) {
      detailTbl$pop <- dplyr::coalesce(detailTbl$pop, detailTbl$detailPathPop)
    }

    popMap <- .gateGetChnlPopMap(pathProject)
    if (length(popMap) > 0L) {
      detailTbl$pop <- dplyr::coalesce(
        detailTbl$pop, unname(popMap[detailTbl$chnl])
      )
    }
    pops <- try(.gateGetPop(pathProject), silent = TRUE)
    if (!inherits(pops, "try-error") && length(pops) == 1L) {
      detailTbl$pop <- dplyr::coalesce(detailTbl$pop, pops)
    }

    if (!is.null(pop)) {
      detailTbl <- detailTbl |>
        dplyr::filter(.data$pop %in% .env$pop)
    }
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
        detailPathPop = meta$pop,
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
  popParts <- parts[startsWith(parts, "pop")]
  popFromPath <- if (length(popParts) > 0L) {
    sub("^pop", "", popParts[[1]])
  } else {
    NA_character_
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
    pop = popFromPath,
    detailObject = detailObject
  )
}
