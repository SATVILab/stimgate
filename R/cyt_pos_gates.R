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
  chnlVec <- .getCytPosGatesChnlVecFromChnlList(
    chnlSettings = chnlSettings
  )

  # chnlLab
  chnlLabVec <- .getLabs(.data = .data[[1]], chnlCut = chnlVec)

  # get max bwMin for densities from chnlSettings elements
  bwMin <- .gateCytPosMaxBwMin(chnlSettings)

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
    force(gateTblGn)
    .getCytPosGatesGateName(
      gateTblGn = gateTblGn,
      .data = .data,
      indBatchList = indBatchList,
      chnlVec = chnlVec,
      chnlLabVec = chnlLabVec,
      popGate = chnlSettings[[1]]$popGate,
      bwMin = bwMin,
      calcCytPos = calcCytPos,
      stage = stage,
      pathProject = pathProject
    )
  })
}

#' @keywords internal
.getCytPosGatesGateName <- function(
  gateTblGn,
  .data,
  indBatchList,
  chnlVec,
  chnlLabVec,
  popGate,
  bwMin,
  calcCytPos,
  stage,
  pathProject
) {
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
      popGate = popGate,
      bwMin = bwMin,
      calcCytPos = calcCytPos,
      stage = stage,
      pathProject = pathProject,
      batch = batch
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  # join gateCyt onto gateTbl
  gateTblGn |>
    dplyr::left_join(
      cpTblCyt |>
        dplyr::select(batch, ind, chnl, marker, gateCyt), # nolint
      by = c("batch", "ind", "chnl", "marker")
    )
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
  calcCytPos,
  stage,
  batch,
  pathProject
) {
  .debug("Getting cyt+ gates for ind: ", ind) # nolint

  # return if ind in batch is the last one, as that is the unstim ind
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

  exUns <- .getEx(
    .data = .data[[indUns]],
    pop = popGate, # nolint
    chnlCut = chnlVec,
    ind = indUns, # BUG FIX: Cache under unstim ind, not stim ind
    indUns = indUns,
    batch = batch,
    pathProject = pathProject
  )

  # gates
  # -----------------
  gateTblInd <- gateTblGn |>
    dplyr::filter(.data$ind == .env$ind) # nolint

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
  bwMin,
  ind,
  stage,
  pathProject
) {
  .debug("chnlCurr: ", chnlCurr) # nolint
  if (is.na(gateTblInd$gate[gateTblInd$chnl == chnlCurr])) {
    return(NA)
  }

  # subset only cells pos for at least one other cyt
  # --------------
  nonNaChnlVec <- gateTblInd |>
    dplyr::filter(!is.na(gate)) |> # nolint
    dplyr::pull("chnl")
  incVec <- .getPosInd(
    ex = ex,
    gateTbl = gateTblInd,
    chnl = setdiff(nonNaChnlVec, chnlCurr),
    chnlAlt = setdiff(nonNaChnlVec, chnlCurr),
    gateTypeCytPos = "base",
    gateTypeSinglePos = "base"
  )
  .intSaveNm(
    paste0(chnlCurr, "_incVec"),
    incVec,
    ind,
    stage,
    pathProject
  )

  # get original cutpoint
  cpOrig <- gateTblInd |>
    dplyr::filter(chnl == chnlCurr) |> # nolint
    dplyr::pull("gate")

  .intSaveNm(
    paste0(chnlCurr, "_cpOrig"),
    cpOrig,
    ind,
    stage,
    pathProject
  )

  # =====================
  # Cutpoint - based on cyt+ cells
  # =====================
  cpPos <- .getCpPos(
    ex = ex,
    inc = incVec,
    chnl = chnlCurr,
    bwMin = bwMin,
    trustNoOrHighAm = FALSE,
    minCell = 10,
    cpOrig = cpOrig,
    nLoop = 5
  )
  .intSaveNm(
    paste0(chnlCurr, "_cpPos"),
    cpPos,
    ind,
    stage,
    pathProject
  )

  # BUG FIX: Actually calculate cpNeg here instead of checking for its existence
  cpNeg <- .getCpNeg(
    ex = ex,
    inc = incVec,
    chnl = chnlCurr,
    bwMin = bwMin,
    minCell = 10 
  )

  # =====================
  # Final cutpoint
  # =====================
  cpCytPos <- min(cpPos, cpOrig, na.rm = TRUE)
  if (!is.na(cpNeg)) {
    cpCytPos <- min(cpCytPos, cpNeg, na.rm = TRUE)
  }
  
  .intSaveNm(
    paste0(chnlCurr, "_cpCytPosFinal"),
    cpCytPos,
    ind,
    stage,
    pathProject
  )

  cpCytPos
}

#' @keywords internal
.gateCytPosMaxBwMin <- function(chnl, quant = 0.8) {
  chnl |>
    purrr::map_dbl(~ .x$bwMin) |>
    quantile(quant)
}
