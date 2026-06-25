#' @keywords internal
.getStatsChnlLabGet <- function(chnlLab, .data, chnl) {
  if (is.null(chnlLab)) {
    chnlLab <- .getLabs(
      .data = .data[[1]],
      chnlCut = chnl
    )
  }
  chnlLab
}

#' @keywords internal
.getStatsGateTblGet <- function(
  gateTbl,
  chnlLab,
  pathProject,
  popGate,
  gateName = NULL,
  tolClust = NULL
) {
  if (!is.null(gateTbl)) {
    return(gateTbl)
  }
  purrr::map_df(
    names(chnlLab),
    function(chnlCurr) {
      gateTblCurr <- .gatesGetPathAll(
        pathProject = pathProject,
        pop = popGate,
        chnlCut = chnlCurr,
        init = FALSE
      ) |>
        readRDS()
      if (!is.null(gateName)) {
        gateTblCurr <- gateTblCurr |>
          dplyr::filter(gateName == .env$gateName) # nolint
      }
      if (!is.null(tolClust)) {
        if (tolClust) {
          gateTblCurr <- gateTblCurr |>
            dplyr::filter(grepl("Clust$", gateName))
        }
      }

      gateTblCurr |>
        dplyr::mutate(
          chnl = chnlCurr,
          marker = chnlLab[chnlCurr]
        ) |>
        dplyr::select(
          chnl,
          marker,
          gateName, # nolint
          batch,
          ind,
          gate,
          gateCyt,
          gateSingle # nolint
        )
    }
  )
}

#' @keywords internal
.getStatsChnlGet <- function(chnl, gateTbl) {
  if (!is.null(chnl)) {
    return(chnl)
  }
  unique(gateTbl$chnl)
}

#' @keywords internal
.getStatsGateNameGet <- function(gateName, gateTbl) {
  if (!is.null(gateName)) {
    return(gateName)
  }
  unique(gateTbl$gateName)
}

#' @keywords internal
.getStatsCombnMatListGet <- function(nChnl, nPos) {
  purrr::map(
    seq_len(nChnl),
    function(nPos) gtools::combinations(n = nChnl, r = nPos)
  ) |>
    stats::setNames(seq_len(nChnl))
}

#' @keywords internal
.getStatsCytCombnVecListGet <- function(combnMatList, chnl) {
  purrr::map(
    names(combnMatList),
    function(nPosNm) {
      combnMat <- combnMatList[[nPosNm]]
      purrr::map_chr(seq_len(nrow(combnMat)), function(i) {
        chnlPos <- chnl[combnMat[i, , drop = TRUE]]
        purrr::map_chr(chnl, function(chnlCurr) {
          if (chnlCurr %in% chnlPos) {
            return(paste0(chnlCurr, "~+~"))
          }
          paste0(chnlCurr, "~-~")
        }) |>
          paste0(collapse = "")
      })
    }
  ) |>
    stats::setNames(names(combnMatList))
}

#' @keywords internal
.getStatsGateTblSave <- function(
  gateTbl,
  pathProject,
  popGate,
  chnlLab,
  chnl,
  save
) {
  if (!save) {
    return(invisible(FALSE))
  }
  if (!"chnl" %in% colnames(gateTbl)) {
    gateTbl <- gateTbl |>
      dplyr::mutate(
        chnl = chnl[[1]],
        marker = chnlLab[chnl[[1]]]
      )
  }

  gateTbl <- gateTbl |>
    dplyr::select(
      gateName,
      chnl,
      marker,
      ind,
      dplyr::everything() # nolint
    )
  gateTbl[, "ind"] <- as.character(gateTbl[["ind"]])

  gateTbl <- gateTbl |>
    dplyr::arrange(gateName, chnl, marker, ind) # nolint
  pathSaveRds <- .gatesGetPathAll(
    pathProject = pathProject,
    pop = popGate,
    chnlCut = chnl[1],
    init = FALSE
  )
  pathSaveCsv <- sub("\\.rds$", ".csv", pathSaveRds)
  if (file.exists(pathSaveRds)) {
    file.remove(pathSaveRds)
  }
  if (file.exists(pathSaveCsv)) {
    file.remove(pathSaveCsv)
  }
  if (!dir.exists(dirname(pathSaveCsv))) {
    dir.create(dirname(pathSaveCsv), recursive = TRUE)
  }
  utils::write.csv(gateTbl, pathSaveCsv, row.names = FALSE)
  saveRDS(gateTbl, pathSaveRds)
}

#' @keywords internal
.statsSave <- function(save, statTbl, pathProject) {
  if (!save) {
    return(invisible(statTbl))
  }
  if (!dir.exists(pathProject)) {
    dir.create(pathProject, recursive = TRUE)
  }
  if ("ind" %in% colnames(statTbl)) {
    statTbl[, "ind"] <- as.character(statTbl[["ind"]])
  }
  if ("batch" %in% colnames(statTbl)) {
    statTbl[, "batch"] <- as.character(statTbl[["batch"]])
  }
  fnRds <- "gateStats.rds"
  fnCsv <- "gateStats.csv"
  pathSaveFnRds <- file.path(pathProject, fnRds)
  pathSaveFnCsv <- file.path(pathProject, fnCsv)
  if (file.exists(pathSaveFnRds)) {
    file.remove(pathSaveFnRds)
  }
  if (file.exists(pathSaveFnCsv)) {
    file.remove(pathSaveFnCsv)
  }
  utils::write.csv(statTbl, pathSaveFnCsv, row.names = FALSE)
  saveRDS(statTbl, pathSaveFnRds)
  invisible(pathProject)
}
