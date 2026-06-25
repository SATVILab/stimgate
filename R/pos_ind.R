# Get logical indicator for cytokine-positive cells
# Returns a logical vector indicating whether each cell is positive
# for any of the specified channels using a single threshold type
#' @keywords internal
.getPosIndSimple <- function(ex, gateTbl, chnl = NULL, gateType) {
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }
  posVec <- rep(FALSE, nrow(ex))
  for (chnlCurr in chnl) {
    gateTblChnlInd <- which(gateTbl$chnl == chnlCurr)
    gateTblColName <- switch(
      gateType,
      "base" = "gate",
      "cyt" = "gate_cyt",
      "single" = "gate_single",
      stop(paste0(
        "gateType ",
        gateType,
        " not recognised in .getPosIndSimple"
      ))
    )
    posVec <- posVec |
      ex[[chnlCurr]] > gateTbl[[gateTblColName]][[gateTblChnlInd]]
  }
  posVec
}

# Identify cells that are positive for at least two cytokines.
# Get logical indicator for multi-functional cells
# Identifies cells positive for multiple cytokines with customizable requirements
#' @keywords internal
.getPosIndMult <- function(
  ex,
  gateTbl,
  chnl = NULL,
  chnlAlt = NULL,
  gateTypeCytPos
) {
  # must specify types of gates to use for cyt+ cells
  if (!gateTypeCytPos %in% c("base", "cyt")) {
    stop(paste0(
      "gateTypeCytPos value of ",
      ifelse(missing(gateTypeCytPos), "blank", gateTypeCytPos),
      ' not either "cyt" or "base" in function .getPosIndMult.'
    ))
  }

  # set defaults
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }
  if (is.null(chnlAlt)) {
    chnlAlt <- chnl
  }

  # if not using cyt+ gates, simply count any cell as multifunctional
  # if it is above base threshold for at least two cyts
  if (gateTypeCytPos == "base") {
    posVecCytPosMult <- rep(0, nrow(ex))
    for (chnlCurr in chnl) {
      posVecCytPosCurr <- .getPosIndSimple(
        ex = ex,
        gateTbl = gateTbl,
        chnl = chnlCurr,
        gateType = "base"
      )
      posVecCytPosMult <- posVecCytPosMult + posVecCytPosCurr
    }
    return(posVecCytPosMult >= 2)
  }

  # if using cyt+ gates, then for each cytokine in chnl, determine
  # which cells are positive for all other cytokines in chnlAlt and
  # then determine which of these are also positive for the current chnl.
  if (gateTypeCytPos == "cyt") {
    posVecCytPosMult <- rep(FALSE, nrow(ex))
    for (chnlCurr in chnl) {
      posVecCurr <- .getPosIndSimple(
        ex = ex,
        gateTbl = gateTbl,
        chnl = chnlCurr,
        gateType = "base"
      )
      posVecAlt <- .getPosIndSimple(
        ex = ex,
        gateTbl = gateTbl,
        chnl = setdiff(c(chnlAlt, chnl), chnlCurr),
        gateType = "cyt"
      )

      # positive for chnlCurr using base threshold and any other using
      # cyt+ thresholds
      posVecCytPosMultCurr <- posVecAlt & posVecCurr
      posVecCytPosMult <- posVecCytPosMult | posVecCytPosMultCurr

      for (chnlAltCurr in setdiff(c(chnl, chnlAlt), chnlCurr)) {
        posVecCurrCyt <- .getPosIndSimple(
          ex = ex,
          gateTbl = gateTbl,
          chnl = chnlCurr,
          gateType = "cyt"
        )

        posVecAltBase <- .getPosIndSimple(
          ex = ex,
          gateTbl = gateTbl,
          chnl = chnlAltCurr,
          gateType = "base"
        )

        posVecCytPosMultCurr <- posVecAltBase & posVecCurrCyt
        posVecCytPosMult <- posVecCytPosMult | posVecCytPosMultCurr
      }
    }
  }

  posVecCytPosMult
}

# Identify cells positive for all cytokines except one
# Finds cells positive for every cytokine except one specified channel
#' @keywords internal
.getPosIndButSinglePosForOneCyt <- function(
  ex,
  gateTbl,
  chnlSingleExc,
  chnl = NULL,
  gateTypeCytPos,
  gateTypeSinglePos
) {
  .getPosIndButSinglePosForOneCytCheck(
    gateTypeSinglePos = gateTypeSinglePos
  )
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }
  chnl <- c(chnlSingleExc, chnl) |>
    unique()

  # cells positive for any cytokine except
  # current cytokine using base or single threshold
  posVecSingleIndAnyCytButCurr <- .getPosIndSimple(
    ex = ex,
    gateTbl = gateTbl,
    chnl = setdiff(chnl, chnlSingleExc),
    gateType = gateTypeSinglePos
  )

  # above specifies all cells that are polyfunctional,
  # if no adjusted thresholds are used
  if (gateTypeCytPos == "base" && gateTypeSinglePos == "base") {
    return(posVecSingleIndAnyCytButCurr)
  }

  # cells positive for any two cytokines
  posVecMultiCyt <- .getPosIndMult(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    chnlAlt = chnl,
    gateTypeCytPos = gateTypeCytPos
  )

  posVecSingleIndAnyCytButCurr | posVecMultiCyt
}

#' @keywords internal
.getPosIndButSinglePosForOneCytCheck <- function(
  gateTypeSinglePos
) {
  if (missing(gateTypeSinglePos)) {
    stop("gateTypeSinglePos missing")
  }
  if (!gateTypeSinglePos %in% c("base", "single")) {
    stop(paste0(
      "gateTypeSinglePos value of ",
      gateTypeSinglePos,
      ' not either "single" or "base" in function .getPosIndButSinglePosForOneCytCheck.'
    ))
  }
  invisible(TRUE)
}

# Identify cells that express at least one cytokine
# Returns a logical vector indicating cytokine-positive cells using flexible thresholds
#' @keywords internal
.getPosInd <- function(
  ex,
  gateTbl,
  chnl,
  chnlAlt = NULL,
  gateTypeCytPos,
  gateTypeSinglePos
) {
  # must specify types of gates to use for single+ cells
  if (!gateTypeSinglePos %in% c("base", "single")) {
    stop(paste0(
      "gateTypeSinglePos value of ",
      ifelse(missing(gateTypeSinglePos), "blank", gateTypeSinglePos),
      " not either 'single' or 'base' in function .getPosInd"
    ))
  }

  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }
  if (is.null(chnlAlt)) {
    chnlAlt <- setdiff(unique(gateTbl$chnl), chnl)
  }
  chnlAlt <- setdiff(chnlAlt, chnl)

  # if only base thresholds are used, then it's sufficient to
  # look only for cells that are positive for the required channels
  if (gateTypeCytPos == "base" && gateTypeSinglePos == "base") {
    posIndVec <- .getPosIndSimple(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      gateType = "base"
    )
    return(posIndVec)
  }

  # cells positive for any cytokine that is required
  posVecSingle <- .getPosIndSimple(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    gateType = gateTypeSinglePos
  )

  # cells positive for at least two cytokines, for some cytokines that are required
  posVecMulti <- .getPosIndMult(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    chnlAlt = chnlAlt,
    gateTypeCytPos = gateTypeCytPos
  )

  # cells positive for either one of the require cytokines (possibly on its own)
  # or any other cytokine together with the required cytokine
  posIndVec <- posVecSingle | posVecMulti

  posIndVec
}

#' @keywords internal
.getPosIndCytCombn <- function(
  ex,
  gateTbl,
  chnlPos,
  chnlNeg,
  chnlAlt,
  gateTypeCytPos,
  gateTypeSinglePos
) {
  # must specify types of gates to use for single+ cells
  if (!gateTypeSinglePos %in% c("base", "single")) {
    stop(paste0(
      "gateTypeSinglePos value of ",
      ifelse(missing(gateTypeSinglePos), "blank", gateTypeSinglePos),
      " not either 'single' or 'base' in function .getPosIndCytCombn"
    ))
  }
  chnl <- unique(c(chnlPos, chnlNeg, chnlAlt))
  chnlPosIndVecPos <- rep(TRUE, nrow(ex))

  # get all cells that are positive for any marker
  for (chnlCurr in chnlPos) {
    chnlPosIndVecPosCurr <- .getPosInd(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnlCurr,
      chnlAlt = setdiff(chnl, chnlCurr),
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos
    )
    chnlPosIndVecPos <- chnlPosIndVecPos & chnlPosIndVecPosCurr
  }

  if (length(chnlNeg) > 0) {
    # get all cells that are negative for any marker that they should be negative
    chnlNegIndVecPos <- .getPosInd(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnlNeg,
      chnlAlt = chnl,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos
    )

    # all cells positive only for positive markers are those that
    # are positive for those markers and negative for all other markers
    chnlPosIndVecPos <- chnlPosIndVecPos & !chnlNegIndVecPos
  }

  chnlPosIndVecPos
}
