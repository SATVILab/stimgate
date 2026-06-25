#' @keywords internal
.getIncVec <- function(chnlCurr, chnlVec, ex, gateTblInd) {
  incVec <- rep(FALSE, nrow(ex))

  for (chnlAlt in setdiff(chnlVec, chnlCurr)) {
    cp <- gateTblInd |>
      dplyr::filter(.data$chnl == chnlAlt) |> # nolint
      dplyr::pull("gate")
    incVec <- incVec | ex[[chnlAlt]] > cp
  }

  incVec
}

#' @keywords internal
.getCpNeg <- function(
  ex,
  inc,
  chnl,
  bwMin,
  minCell = 1e3,
  maxPeakRatio = 1e3
) {
  # get cytokine-negative cells
  exNeg <- ex[!inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]

  # if too few, then return
  if (length(exNeg) <= minCell) {
    return(NA)
  }

  # calculate density
  densNeg <- density(exNeg, bw = "SJ")
  if (densNeg$bw < bwMin) {
    densNeg <- density(exNeg, bw = bwMin)
  }

  # calculate mode furthest to the right
  densLen <- length(densNeg$x)

  modeIndVec <- (2:(densLen - 1))[
    (densNeg$y[-1] > densNeg$y[-densLen])[-(densLen - 1)] &
      (densNeg$y[-densLen] > densNeg$y[-1])[-1]
  ]
  modeVec <- densNeg$x[modeIndVec]
  modeVecY <- densNeg$y[modeIndVec]
  modeIndVec <- modeIndVec[modeVecY > 0.01 * max(modeVecY)]
  modeVec <- densNeg$x[modeIndVec]

  # get all points that are to the
  # right of the right-most mode that
  # also are much smaller than peak
  lowNegDensPtsVec <- densNeg$x[
    densNeg$y < max(densNeg$y) / maxPeakRatio &
      densNeg$x >= max(modeVec)
  ]

  # return NA if no such no low neg dens points found, other return left-most such point
  ifelse(length(lowNegDensPtsVec) == 0, NA, min(lowNegDensPtsVec))
}


#' @keywords internal
.getCpPos <- function(
  ex,
  inc,
  chnl,
  bwMin,
  trustNoOrHighAm = FALSE,
  minCell = 10,
  cpOrig,
  nLoop = 5
) {
  cpPos <- .getCpPosInd(
    ex = ex,
    inc = inc,
    chnl = chnl,
    bwMin = bwMin,
    adjust = 1,
    trustNoOrHighAm = FALSE,
    minCell = 10,
    cpOrig = cpOrig
  )

  k <- 1
  while (is.na(cpPos) && k <= nLoop) {
    if (is.na(cpPos)) {
      cpPos <- .getCpPosInd(
        ex = ex,
        inc = inc,
        chnl = chnl,
        bwMin = bwMin,
        adjust = 0.5^k,
        trustNoOrHighAm = FALSE,
        minCell = 10,
        cpOrig = cpOrig
      )
    }
    k <- k + 1
  }

  cpPos
}


#' @keywords internal
.getCpPosInd <- function(
  ex,
  inc,
  chnl,
  bwMin,
  adjust = 1,
  trustNoOrHighAm = FALSE,
  minCell = 10,
  cpOrig
) {
  # .data
  exPos <- ex[inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]
  if (length(exPos) < minCell) {
    return(NA)
  }
  exNeg <- ex[!inc & ex[[chnl]] > min(ex[[chnl]]), ][[chnl]]

  # bw
  bwDens <- density(exPos, bw = "SJ")$bw
  bwSd <- sd(exNeg)
  bwFinal <- max(bwDens, bwSd, bwMin)

  # calculate density
  # ------------------
  densPos <- density(exPos, bw = bwFinal, adjust = adjust)
  densNeg <- density(
    exNeg,
    bw = bwFinal,
    adjust = adjust,
    from = min(densPos$x),
    to = max(densPos$x)
  )

  # calculate modes and antimodes
  # -----------------------------
  densLen <- length(densPos$y)

  # antimodes
  amIndVec <- (2:(densLen - 1))[
    (densPos$y[-1] < densPos$y[-densLen])[-(densLen - 1)] &
      (densPos$y[-densLen] < densPos$y[-1])[-1]
  ]
  amVec <- densPos$x[amIndVec]
  amVecHeight <- densPos$y[amIndVec]

  # modes
  modeIndVec <- (2:(densLen - 1))[
    (densPos$y[-1] > densPos$y[-densLen])[-(densLen - 1)] &
      (densPos$y[-densLen] > densPos$y[-1])[-1]
  ]
  modeVec <- densPos$x[modeIndVec]
  modeVecHeight <- densPos$y[modeIndVec]
  modeIndVec <- modeIndVec[modeVecHeight > 0.01 * max(modeVecHeight)]
  modeVec <- densPos$x[modeIndVec]
  modeVecHeight <- densPos$y[modeIndVec]

  # calculate lowest mode for cyt-neg cells
  # --------------------------
  modeIndVecNeg <- (2:(densLen - 1))[
    (densNeg$y[-1] > densNeg$y[-densLen])[-(densLen - 1)] &
      (densNeg$y[-densLen] > densNeg$y[-1])[-1]
  ]
  modeVecHeightNeg <- densNeg$y[modeIndVecNeg]
  modeIndVecNeg <- modeIndVecNeg[
    modeVecHeightNeg > 0.01 * max(modeVecHeightNeg)
  ]
  highestModeNeg <- max(densNeg$x[modeIndVecNeg])

  # calculate cpShape
  # -------------------
  if (length(amVec) == 0) {
    if (!trustNoOrHighAm) {
      return(NA)
    }
    if (trustNoOrHighAm) {
      cpShape <- ifelse(cpOrig < max(modeVec), highestModeNeg, NA)
      return(cpShape)
    }
  }

  # nearest mode to cpOrig
  modeVecAboveCpOrig <- modeVec[modeVec > cpOrig]

  if (length(modeVecAboveCpOrig) == 0) {
    return(NA)
  }
  modeAboveCpOrigMin <- min(modeVecAboveCpOrig)

  amVecMoreThanCpOrig <- amVec[amVec > cpOrig]
  amRightMinInd <- ifelse(
    length(amVecMoreThanCpOrig) > 0,
    which(amVec == min(amVecMoreThanCpOrig)),
    NA
  )
  amRightMin <- amVec[amRightMinInd]

  if (!is.na(amRightMin[1])) {
    if (amRightMin < modeAboveCpOrigMin) {
      return(amRightMin)
    }
  }

  # get left-most antimode of antimodes less than cpOrig
  amVecLessThanCpOrig <- amVec[amVec < cpOrig]
  maxLeftAmInd <- ifelse(
    length(amVecLessThanCpOrig) > 0,
    which(amVec == max(amVecLessThanCpOrig)),
    NA
  )
  maxLeftAm <- amVec[maxLeftAmInd]
  maxLeftAmHeight <- amVecHeight[maxLeftAmInd]

  if (length(maxLeftAm) == 0 || all(is.na(maxLeftAm))) {
    cpShape <- ifelse(trustNoOrHighAm, highestModeNeg, NA)
    return(cpShape)
  }

  maxModeLessThanCpOrig <- max(modeVec[modeVec < cpOrig])
  maxModeLessThanCpOrigHeight <- modeVecHeight[
    which(modeVec == maxModeLessThanCpOrig)
  ]
  if (maxLeftAmHeight > 0.5 * maxModeLessThanCpOrigHeight) {
    if (trustNoOrHighAm) {
      minModeAboveCpOrigInd <- which(
        modeVec == min(modeVec[modeVec > cpOrig])
      )
      minModeAboveCpOrigHeight <- modeVecHeight[
        minModeAboveCpOrigInd
      ]
      minModeAboveCpOrig <- modeVec[minModeAboveCpOrigInd]

      propMoveToRight <- maxModeLessThanCpOrigHeight /
        (maxModeLessThanCpOrigHeight + minModeAboveCpOrigHeight)
      cpShape <- max(
        maxLeftAm,
        maxModeLessThanCpOrig +
          propMoveToRight *
            (minModeAboveCpOrig - maxModeLessThanCpOrig)
      )
      return(cpShape)
    } else {
      return(NA)
    }
  }

  maxLeftAm
}


#' @keywords internal
.getCytPosGatesChnlVecFromChnlList <- function(chnlSettings) {
  purrr::map_chr(chnlSettings, function(x) x$chnlCut)
}

#' @keywords internal
.getCytPosGatesGateTblGet <- function(
  chnlVec,
  pop,
  pathProject,
  chnlLab
) {
  .debug("Getting gateTbl") # nolint
  purrr::map_df(chnlVec, function(chnlCurr) {
    .gatesGetPathAll(
      pathProject = pathProject,
      pop = pop,
      chnlCut = chnlCurr,
      init = TRUE
    ) |>
      readRDS() |>
      dplyr::mutate(chnl = chnlCurr, marker = chnlLab[chnlCurr]) |>
      dplyr::select(chnl, marker, gateName, batch, ind, gate) # nolint
  })
}
