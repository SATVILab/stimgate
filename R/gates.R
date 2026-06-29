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
#' @examples{
#' # Get example dataset
#' exampleData <- get_example_data()
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
#' # Get statistics for the identified gates
#' gates <- getStimGates(pathProject)
#' }
#' @export
getStimGates <- function(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL
) {
  pop <- pop %|c|% .gateGetPop(pathProject)

  purrr::map_df(pop, function(popCurr) {
    chnlVec <- if (!is.null(marker)) {
      markerLab <- stimgateMetReadMarkerLab(pathProject)
      markerLab[marker] |> stats::setNames(NULL)
    } else {
      chnl %|c|% .gateGetChnl(pathProject, popCurr)
    }

    purrr::map_df(chnlVec, function(chnlCurr) {
      markerCurr <- stimgateMetReadChnlLab(pathProject)[chnlCurr] |>
        stats::setNames(NULL)

      .gatesGetPathAll(
        pathProject = pathProject,
        pop = popCurr,
        chnlCut = chnlCurr,
        init = FALSE
      ) |>
        readRDS() |>
        dplyr::filter(chnl == chnlCurr) |>
        dplyr::mutate(marker = markerCurr) |>
        dplyr::mutate(pop = popCurr) |>
        dplyr::select(pop, dplyr::everything())
    })
  })
}

#' @keywords internal
.gateGetPop <- function(pathProject) {
  pathDir <- file.path(pathProject, "gates")
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  dirVec <- list.dirs(pathDir, full.names = FALSE, recursive = FALSE)
  popVec <- unique(sub("^pop(.*)$", "\\1", dirVec))
  popVec
}

#' @keywords internal
.gateGetChnl <- function(pathProject, pop) {
  pathDir <- file.path(pathProject, "gates", paste0("pop", pop))
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  dirVec <- list.dirs(pathDir, full.names = FALSE, recursive = FALSE)
  chnlVec <- unique(sub("^chnl(.*)$", "\\1", dirVec))
  chnlVec
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

  if (nrow(detailTbl) == 0L) {
    if (isTRUE(save)) {
      pathSave <- pathSave %||% file.path(pathProject, "gatesDetailed.rds")
      saveRDS(detailTbl, pathSave)
    }
    return(detailTbl)
  }

  if (!"chnl" %in% names(detailTbl)) {
    detailTbl$chnl <- detailTbl$detailPathChnl
  } else if ("detailPathChnl" %in% names(detailTbl)) {
    detailTbl$chnl <- dplyr::coalesce(detailTbl$chnl, detailTbl$detailPathChnl)
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
  if (!is.null(pop)) {
    detailTbl <- detailTbl |>
      dplyr::filter(is.na(.data$pop) | .data$pop %in% .env$pop)
  }

  detailTbl <- detailTbl |>
    dplyr::select(
      pop,
      marker,
      chnl,
      dplyr::everything()
    )

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
    pattern = "^locDetail.*\\.rds$",
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
  if (detailObject %in% "locDetailClusterFinal") {
    if (!"detailLevel" %in% names(obj)) {
      obj$detailLevel <- "cluster_final"
    }
    if (!"threshold" %in% names(obj) && "locFinalThreshold" %in% names(obj)) {
      obj$threshold <- obj$locFinalThreshold
    }
    if (
      !"thresholdOrigin" %in% names(obj) &&
        "locFinalThresholdOrigin" %in% names(obj)
    ) {
      obj$thresholdOrigin <- obj$locFinalThresholdOrigin
    }
    if (!"nCellStim" %in% names(obj) && "locFinalNCellStim" %in% names(obj)) {
      obj$nCellStim <- obj$locFinalNCellStim
    }
    if (!"nCellUns" %in% names(obj) && "locFinalNCellUns" %in% names(obj)) {
      obj$nCellUns <- obj$locFinalNCellUns
    }
    if (!"propStim" %in% names(obj) && "locFinalPropStim" %in% names(obj)) {
      obj$propStim <- obj$locFinalPropStim
    }
    if (!"propUns" %in% names(obj) && "locFinalPropUns" %in% names(obj)) {
      obj$propUns <- obj$locFinalPropUns
    }
    if (!"propBs" %in% names(obj) && "locFinalPropBs" %in% names(obj)) {
      obj$propBs <- obj$locFinalPropBs
    }
  }
  obj
}

#' @keywords internal
.gateGetDetailedPathMeta <- function(pathCurr, pathInt) {
  pathCurrNorm <- normalizePath(pathCurr, winslash = "/", mustWork = FALSE)
  pathIntNorm <- normalizePath(pathInt, winslash = "/", mustWork = FALSE)
  rel <- sub(
    paste0("^", gsub("([\\\\.\"])", "\\\\\\1", pathIntNorm), "/?"),
    "",
    pathCurrNorm
  )
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  detailObject <- sub("\\.rds$", "", basename(pathCurr))
  list(
    stage = parts[[1]] %||% NA_character_,
    chnl = parts[[2]] %||% NA_character_,
    ind = if (length(parts) >= 4L && parts[[3]] == "ind") {
      parts[[4]]
    } else {
      NA_character_
    },
    detailObject = detailObject
  )
}
