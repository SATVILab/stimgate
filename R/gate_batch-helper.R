#' @keywords internal
.gateBatchAll <- function(
    indBatch,
    batch,
    exList,
    .data,
    chnlSettings,
    stage,
    pathProject) {
  .debug("chnlSettings$gateTbl is NULL") # nolint
  .debug(
    "gating ",
    paste0(indBatch, collapse = "-") # nolint
  )

  # create bare list
  gateList <- .getCpUnsLoc(
    exList = exList,
    .data = .data,
    chnlSettings = chnlSettings,
    stage = stage,
    pathProject = pathProject
  )

  if (!is.null(chnlSettings$tolCtrl)) {
    for (tol in chnlSettings$tolCtrl) {
      .debug("getting tg-based cutpoint as a control") # nolint
      gateList[[paste0("tgCtrl_", tol)]] <- .getCpTg(
        exList = exList,
        chnlSettings = chnlSettings,
        tgType = "tolCtrl",
        stage = stage,
        pathProject = pathProject
      )
    }
  }

  .gateBatchTbl(gateList, attr(exList[[1]], "batch")) # nolint
}

#' @keywords internal
.gateBatchTbl <- function(gateList, batch) {
  purrr::map_df(seq_along(gateList), function(i) {
    .gateBatchTblAlongType(gateList, batch, i)
  })
}

#' @keywords internal
.gateBatchTblAlongType <- function(gateList, batch, i) {
  .debug("gate list index", i) # nolint
  cpList <- .gateBatchTblCp(gateList[[i]])
  gateType <- .gateBatchTblType(gateList, i)
  purrr::map_df(seq_along(cpList), function(j) {
    .gateBatchTblAlongCombn(cpList, gateType, batch, j)
  })
}

#' @keywords internal
.gateBatchTblAlongCombn <- function(cpList, gateType, batch, j) {
  .debug("gate list sub-index", j) # nolint
  gateCombn <- .gateBatchTblCombn(cpList, j)
  gateVec <- cpList[[j]]
  n <- length(gateVec)
  locGenerated <- attr(gateVec, "locGenerated") %||% rep(FALSE, n)
  locGeneratedDirect <- attr(gateVec, "locGeneratedDirect") %||% rep(FALSE, n)
  locSource <- attr(gateVec, "locSource") %||% rep(NA_character_, n)
  locReason <- attr(gateVec, "locReason") %||% rep(NA_character_, n)

  tibble::tibble(
    gateName = .gateBatchTblName(gateType, gateCombn),
    gateType = gateType,
    gateCombn = gateCombn,
    batch = batch,
    ind = .gateBatchTblInd(cpList, j),
    gate = .gateBatchTblGate(cpList, j),
    gateUse = .gateBatchTblUse(.env$gateType),
    locGenerated = locGenerated %in% TRUE,
    locGeneratedDirect = locGeneratedDirect %in% TRUE,
    locSource = as.character(locSource),
    locReason = as.character(locReason)
  )
}


#' @keywords internal
.gateBatchTblName <- function(gateType, gateCombn) {
  paste0(gateType, "_", gateCombn)
}

#' @keywords internal
.gateBatchTblCp <- function(gateListElem) {
  if (!"cp" %in% names(gateListElem)) {
    return(gateListElem)
  }
  gateListElem[["cp"]]
}

#' @keywords internal
.gateBatchTblType <- function(gateList, i) {
  names(gateList)[i]
}
#' @keywords internal
.gateBatchTblCombn <- function(cpList, j) {
  names(cpList)[[j]]
}

#' @keywords internal
.gateBatchTblInd <- function(cpList, j) {
  as.character(names(cpList[[j]]))
}
#' @keywords internal
.gateBatchTblGate <- function(cpList, j) {
  cpList[[j]]
}

#' @keywords internal
.gateBatchTblUse <- function(gateType) {
  if (grepl("tgCtrl_", gateType)) {
    return("ctrl")
  }
  "gate"
}
