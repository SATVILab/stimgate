#' @keywords internal
.gateBatchAll <- function(
  indBatch,
  batch,
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
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

  if (!is.null(chnlSettings$tolClust)) {
    .debug("getting tolerance gate") # nolint
    gateList[["tgClust"]] <- .getCpTg(
      exList = exList,
      chnlSettings = chnlSettings,
      tgType = "tolClust",
      stage = stage,
      pathProject = pathProject
    )
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
  tibble::tibble(
    gateName = .gateBatchTblName(gateType, gateCombn),
    gateType = gateType,
    gateCombn = gateCombn,
    batch = batch,
    ind = .gateBatchTblInd(cpList, j),
    gate = .gateBatchTblGate(cpList, j),
    gateUse = .gateBatchTblUse(.env$gateType)
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
  if (grepl("tgClust", gateType)) "tgClust" else "gate"
}

#' @keywords internal
.gateBatchSingle <- function(
  indBatch,
  batch,
  exList,
  .data,
  chnlSettings,
  calcCytPosGates,
  stage,
  pathProject
) {
  .debug("chnlSettings$gateTbl is not NULL") # nolint
  .debug(paste0("Gating batch ", batch))

  # =================================
  # get pre-adj and -clust gates for each gate type
  # =================================
  gateTbl <- .gateBatchSingleTblFormat(chnlSettings$gateTbl)
  gateNameVec <- unique(gateTbl$gateName)

  # get single-pos gates for each of the gate types already done
  purrr::map_df(gateNameVec, function(gateNameCurr) {
    .debug("getting single-pos gate", gateNameCurr) # nolint
    gateTblGnMarker <- gateTbl |>
      dplyr::filter(gateName == gateNameCurr, chnl == chnlSettings$chnlCut) # nolint

    gateNameTblRow <- gateTblGnMarker[1, , drop = FALSE]
    gateMethod <- gateNameTblRow$gateMethod
    adjInd <- stringr::str_detect(gateNameCurr, "Adj")
    clustInd <- stringr::str_detect(gateNameCurr, "Clust")

    # filter to yield cells negative for all cytokine combinations
    # except possibly this cytokine single-positive
    exListNegButSinglePosCurr <- purrr::map(
      seq_along(exList),
      function(i) {
        if (i == indInBatchUns) {
          return(exList[[i]])
        }

        gateTblGnInd <- gateTbl |>
          dplyr::filter(
            ind == attr(exList[[i]], "ind"), # nolint
            gateName == gateNameCurr # nolint
          )

        posIndVecButSinglePosCurr <-
          .getPosIndButSinglePosForOneCyt(
            ex = exList[[i]],
            gateTbl = gateTblGnInd,
            chnlSingleExc = chnlSettings$chnlCut,
            chnl = NULL,
            gateTypeCytPos = if (calcCytPosGates) "cyt" else "base",
            gateTypeSinglePos = "base"
          )
        exList[[i]][!posIndVecButSinglePosCurr, , drop = FALSE]
      }
    ) |>
      stats::setNames(names(exList))

    # Temporarily overlay single execution scope overrides into local settings clone
    chnlSettingsExec <- chnlSettings
    chnlSettingsExec$gateCombn <- gateNameTblRow$gateCombn

    gateList <- switch(
      gateMethod,
      "tg" = .getCpTg(
        exList = exListNegButSinglePosCurr,
        chnlSettings = chnlSettingsExec,
        tgType = "tg",
        stage = stage,
        pathProject = pathProject
      ),
      "loc" = {
        chnlSettingsExec$gateNameCurr <- gateNameCurr
        chnlSettingsExec$exUns <- exList[[length(exList)]]
        .getCpUnsLoc(
          exList = exListNegButSinglePosCurr,
          .data = .data,
          chnlSettings = chnlSettingsExec,
          stage = stage,
          pathProject = pathProject
        )[[1]][[1]]
      }
    )

    if (names(gateList)[[1]] == "cp") {
      gateList <- gateList[["cp"]]
    }

    gateTblOut <- switch(
      gateMethod,
      "tg" = ,
      "loc" = ,
      "uns" = tibble::tibble(
        gateName = gateNameCurr,
        ind = names(gateList[[1]]),
        gate = gateList[[1]],
        gateUse = "gate"
      )
    )

    if (adjInd) {
      chnlSettingsAdj <- chnlSettings
      chnlSettingsAdj$gateCombn <- "no"
      chnlSettingsAdj$tol <- tolCtrl

      gateList <- .getCpTg(
        exList = exListNegButSinglePosCurr,
        chnlSettings = chnlSettingsAdj,
        tgType = "Adj",
        stage = stage,
        pathProject = pathProject
      )

      if (names(gateList)[[1]] == "cp") {
        gateList <- gateList[["cp"]]
      }

      gateTblAdj <- tibble::tibble(
        gateName = gateNameCurr,
        ind = names(gateList[[1]]),
        gate = gateList[[1]],
        gateUse = "ctrl"
      )
      gateTblOut <- gateTblOut |> dplyr::bind_rows(gateTblAdj)
    }

    if (clustInd) {
      chnlSettingsClust <- chnlSettings
      chnlSettingsClust$gateCombn <- "no"
      chnlSettingsClust$tol <- tolClustSingle

      gateList <- .getCpTg(
        exList = exListNegButSinglePosCurr,
        chnlSettings = chnlSettingsClust,
        tgType = "Clust",
        stage = stage,
        pathProject = pathProject
      )

      if (names(gateList)[[1]] == "cp") {
        gateList <- gateList[["cp"]]
      }

      gateTblClust <- tibble::tibble(
        gateName = gateNameCurr,
        ind = names(gateList[[1]]),
        gate = gateList[[1]],
        gateUse = "tgClust"
      )
      gateTblOutFinal <- gateTblOut |>
        dplyr::bind_rows(gateTblClust)
    } else {
      gateTblOutFinal <- gateTblOut
    }

    gateTblOutFinal
  }) |>
    dplyr::mutate(
      gateType = purrr::map_chr(
        gateName, # nolint
        function(gn) stringr::str_split(gn, "_")[[1]][1]
      ),
      gateCombn = gateName |>
        stringr::str_remove("Adj") |>
        stringr::str_remove("Clust") |>
        stringr::str_remove(gateType) |> # nolint
        stringr::str_remove("_")
    ) |>
    dplyr::mutate(ind = as.numeric(ind)) # nolint
}

#' @keywords internal
.gateBatchSingleTblFormat <- function(gateTbl) {
  gateTbl |>
    dplyr::mutate(
      gateMethod = .gateBatchSingleTblFormatMethod(gateName) # nolint
    ) |>
    dplyr::mutate(
      gateCombn = .gateBatchSingleTblFormatCombn(gateName) # nolint
    ) |>
    dplyr::mutate(
      clust = purrr::map_chr(
        gateName,
        function(x) ifelse(stringr::str_detect(x, "Clust"), "Clust", "")
      )
    ) |>
    dplyr::mutate(
      adj = purrr::map_chr(
        gateName,
        function(x) ifelse(stringr::str_detect(x, "Adj"), "Adj", "")
      )
    )
}

#' @keywords internal
.gateBatchSingleTblFormatMethod <- function(gateName) {
  gateMethod <- ifelse(
    stringr::str_detect(gateName, "loc"),
    "loc",
    NA # nolint
  )
  gateMethod <- ifelse(
    stringr::str_detect(gateName, "tg"),
    "tg",
    gateMethod # nolint
  )
  ifelse(
    stringr::str_detect(gateName, "uns"),
    "uns",
    gateMethod
  )
}

#' @keywords internal
.gateBatchSingleTblFormatCombn <- function(gateName) {
  purrr::map_chr(gateName, function(x) {
    if (stringr::str_detect(x, "No")) {
      return("no")
    }
    if (stringr::str_detect(x, "Min")) {
      return("min")
    }
    if (stringr::str_detect(x, "Prejoin")) {
      return("prejoin")
    }
  })
}
