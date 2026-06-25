strDetectAny <- function(string, pattern) {
  vapply(
    pattern,
    function(patternCurr) grepl(patternCurr, string),
    logical(1)
  ) |>
    any()
}

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
      .data = .data[[i]],
      pop = pop,
      chnlCut = chnlCut,
      extraChnl = extraChnl,
      ind = i,
      # specify corresponding unstim
      indUns = indBatch[length(indBatch)],
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
  # collect all the channels we need
  # get expression information as a tibble
  # get .data
  allSaved <- .getExCheckChnlSaved(
    chnl = c(chnlCut, extraChnl),
    ind = ind,
    pop = pop,
    pathProject = pathProject
  )
  ex <- if (allSaved) {
    .getExOld(
      pop = pop,
      chnl = c(chnlCut, extraChnl),
      ind = ind,
      pathProject = pathProject
    )
  } else {
    .getExNew(
      .data = .data,
      chnl = c(chnlCut, extraChnl),
      ind = ind,
      pop = pop,
      pathProject = pathProject,
      save = TRUE
    )
  }
  .getExAddAttributes(
    ex = ex,
    ind = ind,
    indUns = indUns,
    batch = batch,
    chnlCut = chnlCut,
    pop = pop,
    addAttributes = addAttributes
  )
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
.getExOld <- function(pop, chnl, ind, pathProject) {
  # get expression information as a tibble
  # get .data
  ex <- tibble::tibble(
    V1 = .getExOldChnlReadInd(chnl[[1]], ind, pop, pathProject)
  )
  names(ex) <- chnl[[1]]
  chnlRemaining <- chnl[-1]
  for (i in seq_along(chnlRemaining)) {
    chnlCurr <- chnlRemaining[[i]]
    ex[[chnlCurr]] <-
      .getExOldChnlReadInd(chnlCurr, ind, pop, pathProject)
  }
  ex
}

#' @keywords internal
.getExOldChnlReadInd <- function(chnl, ind, pop, pathProject) {
  pathChnl <- .getExChnlPath(chnl, ind, pop, pathProject)
  readRDS(pathChnl)
}

#' @keywords internal
.getExNew <- function(.data, pop, chnl, ind, pathProject, save) {
  # get expression information as a tibble
  # get .data
  fr <- flowWorkspace::gh_pop_get_data(.data, y = pop)
  ex <- flowCore::exprs(fr)[, chnl, drop = FALSE] |>
    tibble::as_tibble()
  .getExNewChnlSave(
    ex = ex,
    ind = ind,
    pop = pop,
    pathProject = pathProject,
    save = save
  )
  ex
}

#' @keywords internal
.getExNewChnlSave <- function(ex, ind, pop, pathProject, save) {
  if (!save) {
    return(invisible(FALSE))
  }
  for (chnlCurr in colnames(ex)) {
    .getExNewChnlSaveInd(
      ex = ex,
      chnl = chnlCurr,
      ind = ind,
      pop = pop,
      pathProject = pathProject
    )
  }
  invisible(TRUE)
}

#' @keywords internal
.getExNewChnlSaveInd <- function(ex, chnl, ind, pop, pathProject) {
  pathChnl <- .getExChnlPath(chnl, ind, pop, pathProject)
  if (file.exists(pathChnl)) {
    return(invisible(FALSE))
  }
  if (!dir.exists(dirname(pathChnl))) {
    dir.create(dirname(pathChnl), recursive = TRUE)
  }
  saveRDS(ex[[chnl]], pathChnl)
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

#' @keywords internal
.getExAddAttributes <- function(
  ex,
  ind,
  indUns,
  batch,
  chnlCut,
  pop,
  addAttributes
) {
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

.getInd <- function(ex) {
  attr(ex, "ind")
}

#' @keywords internal
.getCut <- function(ex) {
  ex[[attr(ex, "chnlCut")]]
}

#' @keywords internal
.getBatchEx <- function(ex) {
  attr(ex, "batch")
}

#' @keywords internal
.getIndUns <- function(ind, indBatchList) {
  hasInd <- vapply(
    indBatchList,
    function(x) ind %in% x,
    logical(1)
  )
  if (sum(hasInd) > 1L) {
    # this is an unstim, as it appears
    # in more than one batch
    return(ind)
  }
  hasInd <- which(hasInd)
  indBatch <- indBatchList[hasInd] |>
    unlist()
  indBatch[[length(indBatch)]]
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

#' @title Read saved expression data from project
#' @description Read channel expression vectors saved under a project's
#'   sampleData directory and return them as a tibble with sample metadata
#'   columns.
#' @param pathProject character Path to project.
#' @param .data GatingSet or NULL GatingSet object to extract expression data
#'   from. Default is NULL.
#' @param pop character or NULL Population name(s). Default is detected from
#'   project sampleData.
#' @param ind character or NULL Index/indices of samples. Default is detected
#'   from project sampleData.
#' @param chnl character or NULL Channel name(s) to return. Default is
#'   detected from project sampleData.
#' @param marker character or NULL Marker name(s) to return. Cannot be
#'   specified with `chnl`. Default is NULL.
#' @param bias logical Whether to add bias to unstimulated sample used in the
#'   gating. Default is `FALSE`.
#' @param excMin logical Whether to exclude cells with the minimum
#'   expression for any channels. Default is FALSE.
#' @param combnExc list or NULL Combinations of channels to exclude. Default
#'   is NULL.
#' @param chnlGate character or NULL Channel name(s) to use for gating.
#'   Cannot be specified with `marker_gate`. Default is NULL.
#' @param markerGate character or NULL Marker name(s) to use for gating.
#'   Cannot be specified with `chnl_gate`. Default is NULL.
#' @param gateTypeCytPos character Gate type to use for cytokine-positive
#'   cells. Default is "cyt".
#' @param gateTypeSinglePos character Gate type to use for single-positive
#'   cells. Default is "single".
#' @param mult logical Whether to return only multi-functional cells (positive
#'   for multiple markers). Default is FALSE.
#' @param gateUnsMethod character Method for gating unstimulated cells.
#'   Default is "min".
#' @param transFn function or NULL Transformation function to apply to
#'   expression values. Default is NULL.
#' @param transChnl character or NULL Channel name(s) to transform when using
#'   channel names. Default is NULL (transforms all channels).
#' @param transMarker character or NULL Marker name(s) to transform when
#'   using marker names. Default is NULL (transforms all markers).
#' @return A tibble with columns `pop`, `ind` and one column per requested
#'   channel. Rows correspond to cells.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "sampleData", "POP1", "ind_1"),
#'   recursive = TRUE
#' )
#' saveRDS(
#'   c(1, 2, 3),
#'   file.path(tmp, "sampleData", "POP1", "ind_1", "chnl_BC1.rds")
#' )
#' saveRDS(
#'   c(4, 5, 6),
#'   file.path(tmp, "sampleData", "POP1", "ind_1", "chnl_BC2.rds")
#' )
#' getStimExpr(tmp)
#' getStimExpr(tmp, chnl = "BC1")
#' }
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
  gateTypeSinglePos = "single",
  mult = FALSE,
  gateUnsMethod = "min",
  transFn = NULL,
  transChnl = NULL,
  transMarker = NULL
) {
  .assertString(pathProject)
  pop <- pop %|c|% .getExProjectPop(pathProject)
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Must not specify both marker and chnl")
  }
  .assertStringVector(pop)
  exList <- purrr::map(pop, function(popCurr) {
    ind <- ind %|c|% .getExProjectInd(pathProject, popCurr)
    .assertStringVector(ind)
    purrr::map(ind, function(indCurr) {
      chnl <- if (!is.null(marker)) {
        isMarker <- TRUE
        marker <- as.character(marker)
        stimgateMetaReadMarkerLab(pathProject)[marker]
      } else {
        isMarker <- FALSE
        chnl %|c|%
          .getExProjectChnl(pathProject, popCurr, indCurr)
      }
      .assertStringVector(chnl)
      ex <- .dataGetExInit(
        .data,
        popCurr,
        chnl,
        indCurr,
        pathProject
      )
      ex <- .dataGetExExcMin(ex, excMin, pop, chnl, indCurr)
      ex <- .dataGetExCytPos(
        ex = ex,
        chnlGate = chnlGate,
        markerGate = markerGate,
        pop = popCurr,
        ind = indCurr,
        combnExc = combnExc,
        gateTypeCytPos = gateTypeCytPos,
        gateTypeSinglePos = gateTypeSinglePos,
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
      stats::setNames(ind)
  }) |>
    stats::setNames(pop)
  exDf <- exList |> purrr::map_df(function(x) x |> dplyr::bind_rows())
  probGMinList <- purrr::map(
    exList,
    function(exIndList) {
      purrr::map(
        exIndList,
        function(ex) {
          probGMin <- attr(ex, "probGMin") %||% 1
          if (is.null(probGMin)) {
            return(1.0)
          }
          chnl <- attr(ex, "chnl")
          list(probGMin) |> stats::setNames(chnl)
        }
      ) |>
        stats::setNames(names(exIndList))
    }
  ) |>
    stats::setNames(names(exList))
  attr(exDf, "probGMin") <- probGMinList
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
  popVec <- list.dirs(
    file.path(pathProject, "sampleData"),
    recursive = FALSE
  ) |>
    basename() |>
    sub("^pop_(.*)$", "\\1", x = _)
  .assertStringVector(popVec)
  popVec
}

#' @keywords internal
.getExProjectInd <- function(pathProject, pop = NULL) {
  pop <- pop %||% .getExProjectPop(pathProject)
  pop <- pop[[1]]
  .assertString(pop)
  indVec <- list.dirs(
    file.path(pathProject, "sampleData", paste0("pop_", pop)),
    recursive = FALSE
  ) |>
    basename() |>
    sub("^ind_", "", x = _)
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
  pathChnlDir <- file.path(
    pathProject,
    "sampleData",
    paste0("pop_", pop),
    paste0("ind_", ind)
  )
  .assertString(pathChnlDir)
  chnlVec <- list.files(pathChnlDir) |>
    sub("^chnl_(.*)\\.rds$", "\\1", x = _)
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
  for (chnl in colnames(ex)) {
    bias <- chnlList[[chnl]][["biasUns"]]
    ex[[chnl]] <- ex[[chnl]] + bias
  }
  ex
}

#' @keywords internal
.dataGetExExcMin <- function(
  ex,
  excMin,
  pop = NULL,
  chnl = NULL,
  ind = NULL
) {
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
  gateTypeSinglePos = "single",
  mult = FALSE,
  pathProject
) {
  if (is.null(chnlGate) && is.null(markerGate)) {
    return(ex)
  }
  if (!is.null(chnlGate) && !is.null(markerGate)) {
    stop("Must not specify both chnlGate and markerGate")
  }
  cnVec <- colnames(ex)
  chnlGate <- if (!is.null(markerGate)) {
    isMarker <- TRUE
    stimgateMetaReadMarkerLab(pathProject)[markerGate]
  } else {
    isMarker <- FALSE
    chnlGate %||%
      .getExProjectChnl(pathProject, pop, ind)
  }
  gateTblInd <- .gateGetGateTblAll(NULL, pop, chnlGate, pathProject) |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .dataGetExCytPosInc(
    ex,
    gateTblInd,
    mult,
    chnlGate,
    gateTypeCytPos,
    gateTypeSinglePos
  )

  if (nrow(ex) == 0L) {
    message("No stimulation-positive cells.")
    return(.dataGetExZeroTbl(cnVec))
  }

  ex <- .dataGetExCytPosExc(
    ex,
    combnExc,
    gateTblInd,
    chnlGate,
    gateTypeCytPos,
    gateTypeSinglePos
  )

  if (nrow(ex) == 0L) {
    message(
      "No stimulation-positive cells after excluding specified cytokine combinations."
    ) # nolint
    return(.dataGetExZeroTbl(cnVec))
  }

  ex
}

.dataGetExZeroTbl <- function(cn) {
  outDf <- matrix(rep(NA_real_, length(cn)), ncol = length(cn))
  colnames(outDf) <- cn
  tibble::as_tibble(outDf)
}

#' @keywords internal
.dataGetExCytPosInc <- function(
  ex,
  gateTblInd,
  mult,
  chnl,
  gateTypeCytPos,
  gateTypeSinglePos
) {
  incVec <- rep(FALSE, nrow(ex))

  if (!mult) {
    incVec <- .getPosInd(
      # nolint
      ex = ex,
      gateTbl = gateTblInd,
      chnl = chnl,
      chnlAlt = NULL,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos
    )
  } else {
    incVec <- .getPosIndMult(
      # nolint
      ex = ex,
      gateTbl = gateTblInd,
      chnl = chnl,
      chnlAlt = NULL,
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
  gateTypeCytPos,
  gateTypeSinglePos
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
      chnlAlt = NULL,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos
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
    pop = pop,
    ind = ind
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
