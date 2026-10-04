#' @keywords internal
.getExList <- function(
  .data,
  indBatch,
  batch,
  pop,
  chnlCut,
  extraChnl = NULL,
  pathProject
) {
  isPathGiven <- is.character(pathProject) && nzchar(pathProject)
  if (!isPathGiven) {
    stop("pathProject must be a non-empty character string.")
  }
  # get expression .data for each batch
  lapply(indBatch, function(i) {
    .getEx(
      .data = if (is.null(.data)) NULL else .data[[i]],
      pop = pop,
      chnlCut = chnlCut,
      extraChnl = extraChnl,
      ind = i,
      # specify corresponding unstim
      indUns = indBatch[[1]],
      batch = batch,
      pathProject = pathProject
    )
  }) |>
    stats::setNames(as.character(indBatch))
}


#' @keywords internal
.getEx <- function(
  .data,
  pop,
  chnlCut,
  ind,
  indUns,
  batch,
  extraChnl = NULL,
  pathProject,
  addAttributes = TRUE
) {
  chnl <- c(chnlCut, extraChnl)
  ex <- if (.getExCheckChnlSaved(chnl, ind, pop, pathProject)) {
    lapply(chnl, \(x) readRDS(.getExChnlPath(x, ind, pop, pathProject))) |>
      stats::setNames(chnl) |>
      tibble::as_tibble()
  } else {
    if (is.null(.data)) {
      stop(
        "Incomplete expression cache for sample ", ind,
        ", population ", pop, ", channel(s) ", paste(chnl, collapse = ", "),
        ". Parallel channel workers require a complete, accessible disk cache."
      )
    }
    fr <- flowWorkspace::gh_pop_get_data(.data, y = pop)
    exNew <- flowCore::exprs(fr)[, chnl, drop = FALSE] |>
      tibble::as_tibble()
    # cache each channel so later reads skip the GatingSet
    for (chnlCurr in chnl) {
      pathChnl <- .getExChnlPath(chnlCurr, ind, pop, pathProject)
      if (!file.exists(pathChnl)) {
        dir.create(dirname(pathChnl), recursive = TRUE, showWarnings = FALSE)
        saveRDS(exNew[[chnlCurr]], pathChnl)
      }
    }
    exNew
  }
  if (!addAttributes) {
    return(ex)
  }
  attr(ex, "ind") <- ind |> as.character()
  attr(ex, "indUns") <- indUns |> as.character()
  attr(ex, "isUns") <- ind == indUns
  attr(ex, "chnlCut") <- chnlCut
  attr(ex, "batch") <- batch
  attr(ex, "popGate") <- pop
  ex
}

#' @keywords internal
.getExCheckChnlSaved <- function(chnl, ind, pop, pathProject) {
  pathChnlDir <- .getExChnlPathDir(ind, pop, pathProject)
  if (!dir.exists(pathChnlDir)) {
    return(FALSE)
  }
  fnVec <- list.files(pathChnlDir)
  reqVec <- paste0("chnl_", chnl, ".rds")
  all(reqVec %in% fnVec)
}

#' @keywords internal
.getExChnlPath <- function(chnl, ind, pop, pathProject) {
  file.path(
    .getExChnlPathDir(ind, pop, pathProject),
    paste0("chnl_", chnl, ".rds")
  )
}
#' @keywords internal
.getExChnlPathDir <- function(ind, pop, pathProject) {
  file.path(
    pathProject,
    "sampleData",
    paste0("pop_", pop),
    paste0("ind_", ind)
  )
}

.getInd <- function(ex) {
  attr(ex, "ind")
}

#' @keywords internal
.getCut <- function(ex) {
  ex[[attr(ex, "chnlCut")]]
}

#' @keywords internal
.getIndUns <- function(ind, indBatchList) {
  # the unstim is the first sample of the batch containing `ind`;
  # a shared unstim is first in every batch it appears in
  indBatch <- Filter(function(x) ind %in% x, indBatchList)[[1]]
  indBatch[[1]]
}

#' @keywords internal
.getBatch <- function(ind, indBatchList) {
  hasInd <- vapply(
    indBatchList,
    function(x) ind %in% x,
    logical(1)
  )
  names(indBatchList)[hasInd]
}

#' @title Read cell expression values
#' @description Read expression saved by [gateStim()], optionally selecting
#'   stimulation-positive cells. Supply marker labels or channel names.
#' @param pathProject character Project directory from [gateStim()].
#' @param .data GatingSet, other input accepted by [gateStim()], or NULL Data
#'   passed to [gateStim()], in the same sample order; used only when
#'   the expression values were not saved. Default: NULL.
#' @param pop character or NULL Population names; NULL selects all saved
#'   populations. Default: NULL.
#' @param ind character or numeric vector or NULL Sample indices; NULL selects
#'   all saved samples, including controls. Default: NULL.
#' @param chnl character or NULL Channels to return; NULL selects all saved
#'   channels unless `marker` is supplied. Default: NULL.
#' @param marker character or NULL Marker labels to return; cannot be combined
#'   with `chnl`. Default: NULL.
#' @param bias logical Add the saved `biasUns` shift to controls. Default: FALSE.
#' @param excMin logical Exclude cells at the minimum of any requested channel.
#'   Default: FALSE.
#' @param chnlGate character or NULL Channels used to select positive cells;
#'   include these in the requested expression columns. Default: NULL.
#' @param markerGate character or NULL Marker labels used to select positive
#'   cells; cannot be combined with `chnlGate`. Default: NULL.
#' @param gateTypeCytPos character Positivity rule: "base" uses the main gate;
#'   "cyt" also admits cells above a refined gate when another marker clears
#'   its main gate. Default: "cyt".
#' @param mult logical Require positivity for at least two gating markers.
#'   Applies only when `chnlGate` or `markerGate` is supplied. Default: FALSE.
#' @param combnExc list or NULL Channel combinations to exclude: each vector
#'   specifies positive channels, with other gating channels negative.
#'   Applies only when gating channels are supplied. Default: NULL.
#' @param transFn function or NULL Transformation applied to the expression
#'   tibble before adding metadata columns. Default: NULL.
#' @param transChnl character or NULL Columns to transform when using channels;
#'   NULL transforms all expression columns. Default: NULL.
#' @param transMarker character or NULL Columns to transform when using markers;
#'   NULL transforms all expression columns. Default: NULL.
#' @return A tibble with one row per selected cell, `pop`, `ind`, and expression
#'   columns named by channel (or marker when `marker` is supplied). Empty
#'   selections have zero rows. The `nCellPos` attribute is a tibble with `pop`,
#'   `ind`, `nCellPos` for every requested population/sample pair; `probGMin`
#'   records the fraction of cells kept after removing minimum-expression cells.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' getStimExpr(pathProject, marker = exampleData$marker)
#' @export
getStimExpr <- function(
  pathProject,
  .data = NULL,
  pop = NULL,
  ind = NULL,
  chnl = NULL,
  marker = NULL,
  bias = FALSE,
  excMin = FALSE,
  combnExc = NULL,
  chnlGate = NULL,
  markerGate = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE,
  transFn = NULL,
  transChnl = NULL,
  transMarker = NULL
) {
  .assertString(pathProject)
  pop <- as.character(pop %||% .getExProjectPop(pathProject))
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Must not specify both marker and chnl")
  }
  .assertStringVector(pop)
  if (!is.null(.data)) {
    .data <- .asStimGatingSet(.data, pop)
  }
  exList <- purrr::map(pop, function(popCurr) {
    indCurrVec <- as.character(ind %||% .getExProjectInd(pathProject, popCurr))
    .assertStringVector(indCurrVec)
    purrr::map(indCurrVec, function(indCurr) {
      isMarker <- !is.null(marker)
      chnl <- if (isMarker) {
        stimgateMetaReadMarkerLab(pathProject)[as.character(marker)]
      } else {
        as.character(
          chnl %||% .getExProjectChnl(pathProject, popCurr, indCurr)
        )
      }
      .assertStringVector(chnl)
      ex <- .dataGetExInit(
        if (is.null(.data)) NULL else .data[[as.integer(indCurr)]],
        popCurr,
        chnl,
        indCurr,
        pathProject
      )
      ex <- .dataGetExExcMin(ex, excMin)
      ex <- .dataGetExCytPos(
        ex = ex,
        chnlGate = chnlGate,
        markerGate = markerGate,
        pop = popCurr,
        ind = indCurr,
        combnExc = combnExc,
        gateTypeCytPos = gateTypeCytPos,
        mult = mult,
        pathProject = pathProject
      )
      ex <- .dataGetExBias(
        ex,
        ind = indCurr,
        pathProject = pathProject,
        bias = bias
      )
      ex <- .dataGetExRenamed(ex, isMarker, pathProject)
      transChnlFinal <- if (isMarker) transMarker else transChnl
      ex <- .dataGetExTrans(ex, transFn, transChnlFinal)
      attr(ex, "chnl") <- paste0(chnl, collapse = "&*&")
      ex <- .dataGetExMeta(ex, popCurr, indCurr)
      ex
    }) |>
      stats::setNames(as.character(indCurrVec))
  }) |>
    stats::setNames(as.character(pop))
  exDf <- exList |> purrr::map_df(function(x) x |> dplyr::bind_rows())
  probGMinList <- purrr::map(
    exList,
    function(exIndList) {
      purrr::map(
        exIndList,
        function(ex) {
          list(attr(ex, "probGMin") %||% 1) |>
            stats::setNames(attr(ex, "chnl"))
        }
      ) |>
        stats::setNames(names(exIndList))
    }
  ) |>
    stats::setNames(names(exList))
  attr(exDf, "probGMin") <- probGMinList
  nCellPosDf <- if (length(exList) == 0L) {
    tibble::tibble(
      pop = character(0),
      ind = character(0),
      nCellPos = integer(0)
    )
  } else {
    purrr::map_df(names(exList), function(popCurr) {
      purrr::map_df(names(exList[[popCurr]]), function(indCurr) {
        tibble::tibble(
          pop = as.character(popCurr),
          ind = as.character(indCurr),
          nCellPos = as.integer(nrow(exList[[popCurr]][[indCurr]]))
        )
      })
    })
  }
  attr(exDf, "nCellPos") <- nCellPosDf
  exDf
}

.dataGetExInit <- function(.data, pop, chnl, ind, pathProject) {
  chnlCut <- chnl[[1]]
  extraChnl <- setdiff(chnl, chnlCut)
  extraChnl <- if (length(extraChnl) == 0L) NULL else extraChnl
  .getEx(
    .data,
    pop,
    chnlCut,
    ind,
    NULL,
    NULL,
    extraChnl,
    pathProject,
    FALSE
  )
}

#' @keywords internal
.getExProjectPop <- function(pathProject) {
  .assertString(pathProject)
  pathDir <- file.path(pathProject, "sampleData")
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  popVec <- .gateGetDirs(pathDir, "pop_")
  .assertStringVector(popVec)
  popVec
}

#' @keywords internal
.getExProjectInd <- function(pathProject, pop = NULL) {
  pop <- pop %||% .getExProjectPop(pathProject)
  pop <- pop[[1]]
  .assertString(pop)
  pathDir <- file.path(pathProject, "sampleData", paste0("pop_", pop))
  if (!dir.exists(pathDir)) {
    return(character(0))
  }
  indVec <- .gateGetDirs(pathDir, "ind_")
  .assertStringVector(indVec)
  indVec
}

#' @keywords internal
.getExProjectChnl <- function(pathProject, pop = NULL, ind = NULL) {
  pop <- pop %||% .getExProjectPop(pathProject)
  pop <- pop[[1]]
  .assertString(pop)
  ind <- ind %||% .getExProjectInd(pathProject, pop)
  ind <- ind[[1]]
  .assertString(ind)
  pathChnlDir <- .getExChnlPathDir(ind, pop, pathProject)
  .assertString(pathChnlDir)
  if (!dir.exists(pathChnlDir)) {
    return(character(0))
  }
  fnVec <- list.files(pathChnlDir)
  fnVec <- fnVec[grepl("^chnl_.*\\.rds$", fnVec)]
  chnlVec <- unique(sub("^chnl_(.*)\\.rds$", "\\1", fnVec))
  .assertStringVector(chnlVec)
  chnlVec
}

#' @keywords internal
.dataGetExBias <- function(ex, ind, pathProject, bias) {
  if (!bias) {
    return(ex)
  }
  indBatchList <- stimgateMetaReadBatchList(pathProject)
  indUns <- .getIndUns(ind, indBatchList)
  # only apply bias to unstim
  isUns <- ind == indUns
  if (!isUns) {
    return(ex)
  }
  # apply bias
  chnlList <- stimgateMetaReadSettingsChnls(pathProject)
  chnlLab <- stimgateMetaReadChnlLab(pathProject)
  for (chnl in colnames(ex)) {
    # Completed settings are marker-keyed; older projects may use channels.
    settings <- chnlList[[chnlLab[chnl]]] %||% chnlList[[chnl]]
    bias <- settings[["biasUns"]] %||% 0
    ex[[chnl]] <- ex[[chnl]] + bias
  }
  ex
}

#' @keywords internal
.dataGetExExcMin <- function(ex, excMin) {
  if (!excMin) {
    attr(ex, "probGMin") <- NULL
    return(ex)
  }
  nCellInit <- nrow(ex)
  cnVec <- setdiff(colnames(ex), c("pop", "ind"))
  minVec <- vapply(
    cnVec,
    function(x) min(ex[[x]], na.rm = TRUE),
    numeric(1)
  ) |>
    stats::setNames(cnVec)
  for (cn in cnVec) {
    incVec <- ex[[cn]] > minVec[[cn]]
    ex <- ex[incVec, ]
  }
  attr(ex, "probGMin") <- nrow(ex) / nCellInit
  ex
}

#' @keywords internal
.dataGetExCytPos <- function(
  ex,
  chnlGate,
  markerGate,
  pop,
  ind,
  combnExc = NULL,
  gateTypeCytPos = "cyt",
  mult = FALSE,
  pathProject
) {
  if (is.null(chnlGate) && is.null(markerGate)) {
    return(ex)
  }
  if (!is.null(chnlGate) && !is.null(markerGate)) {
    stop("Must not specify both chnlGate and markerGate")
  }
  chnlGate <- chnlGate %||% stimgateMetaReadMarkerLab(pathProject)[markerGate]
  gateTblInd <- .gateGetGateTblAll(pop, chnlGate, pathProject) |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .dataGetExCytPosInc(
    ex,
    gateTblInd,
    mult,
    chnlGate,
    gateTypeCytPos
  )

  if (nrow(ex) == 0L) {
    message("No stimulation-positive cells.")
    return(ex)
  }

  ex <- .dataGetExCytPosExc(
    ex,
    combnExc,
    gateTblInd,
    chnlGate,
    gateTypeCytPos
  )

  if (nrow(ex) == 0L) {
    message(
      "No stimulation-positive cells after excluding specified cytokine combinations."
    ) # nolint
    return(ex)
  }

  ex
}

#' @keywords internal
.dataGetExCytPosInc <- function(
  ex,
  gateTblInd,
  mult,
  chnl,
  gateTypeCytPos
) {
  incVec <- rep(FALSE, nrow(ex))

  if (!mult) {
    incVec <- .getPosInd(
      # nolint
      ex = ex,
      gateTbl = gateTblInd,
      chnl = chnl,
      chnlAlt = NULL,
      gateTypeCytPos = gateTypeCytPos
    )
  } else {
    incVec <- .getPosIndMult(
      # nolint
      ex = ex,
      gateTbl = gateTblInd,
      chnl = chnl,
      gateTypeCytPos = gateTypeCytPos
    )
  }
  ex[incVec, , drop = FALSE]
}

#' @keywords internal
.dataGetExCytPosExc <- function(
  ex,
  combnExc,
  gateTblInd,
  chnlGate,
  gateTypeCytPos
) {
  if (is.null(combnExc)) {
    return(ex)
  }
  for (chnlPos in combnExc) {
    if (nrow(ex) == 0) {
      break
    }
    excVec <- .getPosIndCytCombn(
      # nolint
      ex = ex,
      gateTbl = gateTblInd,
      chnlPos = chnlPos,
      chnlNeg = setdiff(chnlGate, chnlPos),
      gateTypeCytPos = gateTypeCytPos
    )
    ex <- ex[!excVec, , drop = FALSE]
  }
  ex
}

#' @keywords internal
.dataGetExRenamed <- function(ex, isMarker, pathProject) {
  # if user specified markers, then give them back a table
  # with column names as markers
  if (!isMarker) {
    return(ex)
  }
  colnames(ex) <- stimgateMetaReadChnlLab(pathProject)[
    colnames(ex)
  ]
  ex
}

#' @keywords internal
.dataGetExTrans <- function(ex, transFn, transChnl) {
  # transform
  if (is.null(transFn)) {
    return(ex)
  }
  if (is.null(transChnl)) {
    ex <- transFn(ex)
  } else {
    for (nm in transChnl) {
      ex[, nm] <- transFn(ex[, nm])
    }
  }
  ex
}

#' @keywords internal
.dataGetExMeta <- function(ex, pop, ind) {
  metaDf <- tibble::tibble(
    pop = rep(pop, nrow(ex)),
    ind = rep(ind, nrow(ex))
  )
  attrList <- attributes(ex)
  attrVecNmOrig <- names(attrList)
  attrVecNmAdd <- c(
    "isUns",
    "probGMin",
    "chnl"
  )
  attrVecNmAdd <- intersect(attrVecNmAdd, attrVecNmOrig)
  ex <- tibble::as_tibble(cbind(metaDf, ex))
  for (i in seq_along(attrVecNmAdd)) {
    attr(ex, attrVecNmAdd[i]) <- attrList[[attrVecNmAdd[i]]]
  }
  ex
}
