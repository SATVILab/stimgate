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
          dplyr::filter(gateName %in% .env$gateName) # nolint
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
          gateCyt # nolint
        )
    }
  )
}

#' @keywords internal
.getStatsCombnMatListGet <- function(nChnl) {
  if (nChnl > 30L) {
    stop("Combination statistics support at most 30 channels.")
  }
  purrr::map(
    seq_len(nChnl),
    function(nPos) t(utils::combn(nChnl, nPos))
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
        paste0(chnl, ifelse(chnl %in% chnlPos, "~+~", "~-~"), collapse = "")
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
  uniqueChnls <- unique(gateTbl$chnl)
  for (chnlCurr in uniqueChnls) {
    gateTblCurr <- gateTbl |>
      dplyr::filter(chnl == chnlCurr)
    pathSaveRds <- .gatesGetPathAll(
      pathProject = pathProject,
      pop = popGate,
      chnlCut = chnlCurr,
      init = FALSE
    )
    pathSaveCsv <- sub("\\.rds$", ".csv", pathSaveRds)
    if (!dir.exists(dirname(pathSaveCsv))) {
      dir.create(dirname(pathSaveCsv), recursive = TRUE)
    }
    utils::write.csv(gateTblCurr, pathSaveCsv, row.names = FALSE)
    saveRDS(gateTblCurr, pathSaveRds)
  }
}

#' @keywords internal
.statsSave <- function(save, statTbl, pathProject) {
  if (!save) {
    return(invisible(statTbl))
  }
  if (!dir.exists(pathProject)) {
    dir.create(pathProject, recursive = TRUE)
  }
  fnRds <- "gateStats.rds"
  fnCsv <- "gateStats.csv"
  pathSaveFnRds <- file.path(pathProject, fnRds)
  pathSaveFnCsv <- file.path(pathProject, fnCsv)
  utils::write.csv(statTbl, pathSaveFnCsv, row.names = FALSE)
  saveRDS(statTbl, pathSaveFnRds)
  invisible(pathProject)
}
