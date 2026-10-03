#' @keywords internal
.gateCytPos <- function(
  chnlSettings,
  indBatchList,
  .data,
  gateName = NULL,
  calcCytPos = TRUE,
  stage,
  pathProject
) {
  .debug("-------------") # nolint
  .debug("getting cytokine-positive gates") # nolint
  .debug("-------------") # nolint

  # prep
  # -------------------------------

  # vector of chanls
  chnlVec <- purrr::map_chr(chnlSettings, "chnlCut")

  # chnlLab
  chnlLabVec <- .getLabs(.data = .data[[1]], chnlCut = chnlVec)

  # get max bwMin for densities from chnlSettings elements
  bwMin <- stats::quantile(purrr::map_dbl(chnlSettings, "bwMin"), 0.8)

  # get original gates
  gateTbl <- .getCytPosGatesGateTblGet(
    chnlVec = chnlVec,
    pop = chnlSettings[[1]]$popGate,
    pathProject = pathProject,
    chnlLab = chnlLabVec
  )

  # keep cytokine-positive gates (gateCyt)
  # as the original gates, if not actually
  # gating on cytokine-positive-only cells
  if (!calcCytPos) {
    .debug("Returning original gates as cyt+ gates") # nolint
    return(gateTbl |> dplyr::mutate(gateCyt = gate)) # nolint
  }

  # get cyt+ gates for each of the different gate types
  gnVec <- unique(gateTbl$gateName)
  if (any(grepl("Clust$", gnVec))) {
    gnVec <- gnVec[grepl("Clust$", gnVec)]
  } else if (any(grepl("Adj$", gnVec))) {
    gnVec <- gnVec[grepl("Adj$", gnVec)]
  }

  purrr::map_df(gnVec, function(gn) {
    gateTblGn <- gateTbl |> dplyr::filter(gateName == gn)
    .debug(
      "Getting cyt+ gates for gateName: ",
      gateTblGn$gateName[[1]]
    ) # nolint
    indVec <- unlist(indBatchList)

    cpTblCyt <- purrr::map_df(indVec, function(ind) {
      indUns <- .getIndUns(ind, indBatchList)
      batch <- .getBatch(ind, indBatchList)
      .getCytPosGatesInd(
        ind = ind,
        .data = .data,
        indUns = indUns,
        gateTblGn = gateTblGn,
        chnlVec = chnlVec,
        chnlLabVec = chnlLabVec,
        popGate = chnlSettings[[1]]$popGate,
        bwMin = bwMin,
        stage = stage,
        pathProject = pathProject,
        batch = batch
      )
    })

    if (nrow(cpTblCyt) == 0L) {
      return(dplyr::mutate(gateTblGn, gateCyt = NA_real_))
    }

    # join gateCyt onto gateTbl
    gateTblGn |>
      dplyr::left_join(
        cpTblCyt,
        by = c("batch", "ind", "chnl", "marker")
      )
  })
}

#' @keywords internal
.getCytPosGatesInd <- function(
  ind,
  .data,
  indUns,
  gateTblGn,
  chnlVec,
  chnlLabVec,
  popGate,
  bwMin,
  stage,
  batch,
  pathProject
) {
  .debug("Getting cyt+ gates for ind: ", ind) # nolint

  # return if ind in batch is the unstim ind
  if (ind == indUns) {
    return(NULL)
  }

  # get expression dataframe
  ex <- .getEx(
    .data = .data[[ind]],
    pop = popGate,
    chnlCut = chnlVec,
    ind = ind,
    indUns = indUns,
    batch = batch,
    pathProject = pathProject
  )

  # gates
  # -----------------
  gateTblInd <- gateTblGn |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  basePos <- .getCytPosBasePos(
    ex = ex,
    gateTblInd = gateTblInd
  )

  # ==============
  # Calculate cyt+ cutpoints
  # ==============
  cpVecCytPos <- purrr::map_dbl(
    seq_along(chnlVec),
    function(i) {
      .getCpPosGatesChnl(
        chnlCurr = chnlVec[[i]],
        ex = ex,
        gateTblInd = gateTblInd,
        basePos = basePos,
        bwMin = bwMin,
        ind = ind,
        stage = stage,
        pathProject = pathProject
      )
    }
  )

  tibble::tibble(
    batch = batch,
    ind = attr(ex, "ind"),
    chnl = chnlVec,
    marker = chnlLabVec[chnlVec],
    gateCyt = cpVecCytPos
  )
}

#' @keywords internal
.getCpPosGatesChnl <- function(
  chnlCurr,
  ex,
  gateTblInd,
  basePos,
  bwMin,
  ind,
  stage,
  pathProject
) {
  .debug("chnlCurr: ", chnlCurr) # nolint
  cpOrig <- gateTblInd$gate[gateTblInd$chnl == chnlCurr]
  if (length(cpOrig) == 0L || is.na(cpOrig)) {
    return(NA_real_)
  }

  # subset only cells pos for at least one other cyt
  # --------------
  posCurr <- basePos$pos[[chnlCurr]]

  incVec <- basePos$nPos - as.integer(posCurr) > 0L
  .intSaveNm(
    paste0(chnlCurr, "_incVec"),
    incVec,
    ind,
    stage,
    pathProject
  )

  .intSaveNm(
    paste0(chnlCurr, "_cpOrig"),
    cpOrig,
    ind,
    stage,
    pathProject
  )

  # The full stimulated marginal distribution supplies the safe lower
  # boundary. The conditional distribution among cells positive for at least
  # one other cytokine supplies any candidate antimode.
  shapeRef <- .getCytPosMarginalReference(
    ex = ex,
    chnl = chnlCurr,
    bwMin = bwMin
  )
  .intSaveNm(
    paste0(chnlCurr, "_shapeReference"),
    shapeRef,
    ind,
    stage,
    pathProject
  )

  cpTaut <- .getCpPosTautString(
    ex = ex,
    inc = incVec,
    chnl = chnlCurr,
    cpOrig = cpOrig,
    peakX = shapeRef$peakX,
    windowWidth = shapeRef$windowWidth,
    lower = shapeRef$lowerX
  )
  .intSaveNm(
    paste0(chnlCurr, "_cpTaut"),
    cpTaut$threshold,
    ind,
    stage,
    pathProject
  )

  # =====================
  # Final cutpoint
  # =====================
  cpCytPos <- if (is.finite(cpTaut$threshold)) {
    min(cpOrig, cpTaut$threshold)
  } else {
    cpOrig
  }

  .intSaveNm(
    paste0(chnlCurr, "_tautStringInfo"),
    cpTaut,
    ind,
    stage,
    pathProject
  )
  .intSaveNm(
    paste0(chnlCurr, "_cpCytPosFinal"),
    cpCytPos,
    ind,
    stage,
    pathProject
  )

  cpCytPos
}
