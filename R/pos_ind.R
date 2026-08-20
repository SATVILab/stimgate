# Get cached per-cell base and cytokine-positive threshold comparisons
#' @keywords internal
.getPosIndCache <- function(
  ex,
  gateTbl,
  chnl = NULL,
  posCache = NULL
) {
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  chnl <- unique(as.character(chnl))

  if (is.null(posCache)) {
    posCache <- list(
      base = list(),
      cyt = list(),
      n = nrow(ex)
    )
  }

  if (!identical(posCache$n, nrow(ex))) {
    stop("posCache does not correspond to the supplied expression table.")
  }

  if (is.null(posCache$base)) {
    posCache$base <- list()
  }

  if (is.null(posCache$cyt)) {
    posCache$cyt <- list()
  }

  hasGateCyt <- "gateCyt" %in% colnames(gateTbl)

  for (chnlCurr in chnl) {
    gateTblChnlInd <- which(
      as.character(gateTbl$chnl) == chnlCurr
    )

    if (is.null(posCache$base[[chnlCurr]])) {
      posCache$base[[chnlCurr]] <-
        ex[[chnlCurr]] > gateTbl$gate[[gateTblChnlInd]]
    }

    if (
      hasGateCyt &&
        is.null(posCache$cyt[[chnlCurr]])
    ) {
      posCache$cyt[[chnlCurr]] <-
        ex[[chnlCurr]] > gateTbl$gateCyt[[gateTblChnlInd]]
    }
  }

  posCache
}


# Get one cached logical threshold-comparison vector
#' @keywords internal
.getPosIndCacheGet <- function(
  posCache,
  chnl,
  gateType
) {
  gateList <- switch(
    gateType,
    "base" = posCache$base,
    "cyt" = posCache$cyt,
    stop(
      paste0(
        "gateType ",
        gateType,
        " not recognised."
      )
    )
  )

  out <- gateList[[as.character(chnl)]]

  if (is.null(out)) {
    stop(
      paste0(
        "No ",
        gateType,
        " positivity vector available for channel ",
        chnl,
        "."
      )
    )
  }

  out
}


# Logical OR across cached positivity vectors while preserving NA semantics
#' @keywords internal
.getPosIndCacheAny <- function(
  posCache,
  chnl,
  gateType
) {
  out <- rep(FALSE, posCache$n)

  for (chnlCurr in chnl) {
    out <- out |
      .getPosIndCacheGet(
        posCache = posCache,
        chnl = chnlCurr,
        gateType = gateType
      )
  }

  out
}


# Count TRUE and NA values across cached positivity vectors
#' @keywords internal
.getPosIndCacheCount <- function(
  posCache,
  chnl,
  gateType
) {
  nTrue <- integer(posCache$n)
  nNa <- integer(posCache$n)

  for (chnlCurr in chnl) {
    posCurr <- .getPosIndCacheGet(
      posCache = posCache,
      chnl = chnlCurr,
      gateType = gateType
    )

    nTrue <- nTrue +
      as.integer(
        !is.na(posCurr) &
          posCurr
      )

    nNa <- nNa +
      as.integer(
        is.na(posCurr)
      )
  }

  list(
    nTrue = nTrue,
    nNa = nNa
  )
}


# Convert TRUE/NA counts to the result of a logical OR
#' @keywords internal
.getPosIndCacheAnyFromCount <- function(count) {
  out <- rep(
    FALSE,
    length(count$nTrue)
  )

  out[count$nNa > 0L] <- NA
  out[count$nTrue > 0L] <- TRUE

  out
}


# Get positivity for at least one channel other than the current channel
#' @keywords internal
.getPosIndCacheAnyExcept <- function(
  posCache,
  count,
  chnlCurr,
  gateType
) {
  posCurr <- .getPosIndCacheGet(
    posCache = posCache,
    chnl = chnlCurr,
    gateType = gateType
  )

  countOther <- list(
    nTrue = count$nTrue -
      as.integer(
        !is.na(posCurr) &
          posCurr
      ),
    nNa = count$nNa -
      as.integer(
        is.na(posCurr)
      )
  )

  .getPosIndCacheAnyFromCount(
    countOther
  )
}


# Get logical indicator for cytokine-positive cells
# Returns a logical vector indicating whether each cell is positive
# for any of the specified channels using a single threshold type
#' @keywords internal
.getPosIndSimple <- function(
  ex,
  gateTbl,
  chnl = NULL,
  gateType,
  posCache = NULL
) {
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    posCache = posCache
  )

  .getPosIndCacheAny(
    posCache = posCache,
    chnl = chnl,
    gateType = gateType
  )
}


# Identify cells that are positive for at least two cytokines
#' @keywords internal
.getPosIndMult <- function(
  ex,
  gateTbl,
  chnl = NULL,
  chnlAlt = NULL,
  gateTypeCytPos,
  posCache = NULL
) {
  if (!gateTypeCytPos %in% c("base", "cyt")) {
    stop(
      paste0(
        "gateTypeCytPos value of ",
        ifelse(
          missing(gateTypeCytPos),
          "blank",
          gateTypeCytPos
        ),
        ' not either "cyt" or "base" in function .getPosIndMult.'
      )
    )
  }

  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  if (is.null(chnlAlt)) {
    chnlAlt <- chnl
  }

  if (gateTypeCytPos == "base") {
    posCache <- .getPosIndCache(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      posCache = posCache
    )

    count <- .getPosIndCacheCount(
      posCache = posCache,
      chnl = chnl,
      gateType = "base"
    )

    # Preserve the previous arithmetic behaviour:
    # any NA among the contributing channels produces NA.
    out <- count$nTrue >= 2L
    out[count$nNa > 0L] <- NA

    return(out)
  }

  # In cyt mode a cell is multifunctional when, for at least one
  # requested cytokine, either:
  #
  # 1. that cytokine clears its base threshold and another cytokine
  #    clears its cyt+ threshold, or
  # 2. that cytokine clears its cyt+ threshold and another cytokine
  #    clears its base threshold.
  contextChnl <- unique(
    c(
      as.character(chnlAlt),
      as.character(chnl)
    )
  )

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = contextChnl,
    posCache = posCache
  )

  baseCount <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = contextChnl,
    gateType = "base"
  )

  cytCount <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = contextChnl,
    gateType = "cyt"
  )

  posVecCytPosMult <- rep(
    FALSE,
    nrow(ex)
  )

  for (chnlCurr in chnl) {
    posCurrBase <- .getPosIndCacheGet(
      posCache = posCache,
      chnl = chnlCurr,
      gateType = "base"
    )

    posCurrCyt <- .getPosIndCacheGet(
      posCache = posCache,
      chnl = chnlCurr,
      gateType = "cyt"
    )

    posOtherBase <- .getPosIndCacheAnyExcept(
      posCache = posCache,
      count = baseCount,
      chnlCurr = chnlCurr,
      gateType = "base"
    )

    posOtherCyt <- .getPosIndCacheAnyExcept(
      posCache = posCache,
      count = cytCount,
      chnlCurr = chnlCurr,
      gateType = "cyt"
    )

    posVecCytPosMult <-
      posVecCytPosMult |
      (posCurrBase &
        posOtherCyt) |
      (posCurrCyt &
        posOtherBase)
  }

  posVecCytPosMult
}


# Get context-dependent positivity separately for every supplied cytokine
#' @keywords internal
.getPosIndByChnl <- function(
  ex,
  gateTbl,
  chnl = NULL,
  gateTypeCytPos,
  posCache = NULL
) {
  if (!gateTypeCytPos %in% c("base", "cyt")) {
    stop(
      paste0(
        "gateTypeCytPos value of ",
        gateTypeCytPos,
        ' not either "cyt" or "base".'
      )
    )
  }

  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  chnl <- unique(
    as.character(chnl)
  )

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    posCache = posCache
  )

  if (gateTypeCytPos == "base") {
    return(
      posCache$base[chnl]
    )
  }

  # The current cyt+ rule for one cytokine simplifies exactly to:
  #
  # base-positive for that cytokine
  # OR
  # cyt+-positive for that cytokine AND base-positive for
  # at least one other cytokine.
  baseCount <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = chnl,
    gateType = "base"
  )

  out <- stats::setNames(
    vector(
      "list",
      length(chnl)
    ),
    chnl
  )

  for (chnlCurr in chnl) {
    posCurrBase <- .getPosIndCacheGet(
      posCache = posCache,
      chnl = chnlCurr,
      gateType = "base"
    )

    posCurrCyt <- .getPosIndCacheGet(
      posCache = posCache,
      chnl = chnlCurr,
      gateType = "cyt"
    )

    posOtherBase <- .getPosIndCacheAnyExcept(
      posCache = posCache,
      count = baseCount,
      chnlCurr = chnlCurr,
      gateType = "base"
    )

    out[[chnlCurr]] <-
      posCurrBase |
      (posCurrCyt &
        posOtherBase)
  }

  out
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
  posCache = NULL
) {
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  chnl <- c(
    chnlSingleExc,
    chnl
  ) |>
    unique()

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    posCache = posCache
  )

  # Cells positive for any cytokine except the current cytokine
  # using its base threshold.
  posVecSingleIndAnyCytButCurr <-
    .getPosIndCacheAny(
      posCache = posCache,
      chnl = setdiff(
        chnl,
        chnlSingleExc
      ),
      gateType = "base"
    )

  if (gateTypeCytPos == "base") {
    return(
      posVecSingleIndAnyCytButCurr
    )
  }

  posVecMultiCyt <- .getPosIndMult(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    chnlAlt = chnl,
    gateTypeCytPos = gateTypeCytPos,
    posCache = posCache
  )

  posVecSingleIndAnyCytButCurr |
    posVecMultiCyt
}

# Identify cells to exclude separately for each cytokine
#' @keywords internal
.getPosIndButSinglePosByChnl <- function(
  ex,
  gateTbl,
  chnl = NULL,
  gateTypeCytPos,
  posCache = NULL
) {
  if (!gateTypeCytPos %in% c("base", "cyt")) {
    stop(
      paste0(
        "gateTypeCytPos value of ",
        gateTypeCytPos,
        ' not either "cyt" or "base".'
      )
    )
  }

  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  chnl <- unique(as.character(chnl))

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = chnl,
    posCache = posCache
  )

  baseCount <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = chnl,
    gateType = "base"
  )

  posVecMultiCyt <- if (gateTypeCytPos == "cyt") {
    .getPosIndMult(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      chnlAlt = chnl,
      gateTypeCytPos = "cyt",
      posCache = posCache
    )
  } else {
    NULL
  }

  out <- stats::setNames(
    vector(
      "list",
      length(chnl)
    ),
    chnl
  )

  for (chnlCurr in chnl) {
    otherBasePos <- .getPosIndCacheAnyExcept(
      posCache = posCache,
      count = baseCount,
      chnlCurr = chnlCurr,
      gateType = "base"
    )

    out[[chnlCurr]] <- if (gateTypeCytPos == "base") {
      otherBasePos
    } else {
      otherBasePos | posVecMultiCyt
    }
  }

  out
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
  posCache = NULL
) {
  if (is.null(chnl)) {
    chnl <- unique(gateTbl$chnl)
  }

  if (is.null(chnlAlt)) {
    chnlAlt <- setdiff(
      unique(gateTbl$chnl),
      chnl
    )
  }

  chnlAlt <- setdiff(
    chnlAlt,
    chnl
  )

  if (gateTypeCytPos == "base") {
    posCache <- .getPosIndCache(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      posCache = posCache
    )

    return(
      .getPosIndCacheAny(
        posCache = posCache,
        chnl = chnl,
        gateType = "base"
      )
    )
  }

  if (gateTypeCytPos != "cyt") {
    stop(
      paste0(
        "gateTypeCytPos value of ",
        gateTypeCytPos,
        ' not either "cyt" or "base".'
      )
    )
  }

  contextChnl <- unique(
    c(
      as.character(chnl),
      as.character(chnlAlt)
    )
  )

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = contextChnl,
    posCache = posCache
  )

  posByChnl <- .getPosIndByChnl(
    ex = ex,
    gateTbl = gateTbl,
    chnl = contextChnl,
    gateTypeCytPos = "cyt",
    posCache = posCache
  )

  posIndVec <- rep(
    FALSE,
    nrow(ex)
  )

  for (chnlCurr in chnl) {
    posIndVec <-
      posIndVec |
      posByChnl[[as.character(chnlCurr)]]
  }

  posIndVec
}


# Get cell membership for one exact cytokine combination
#' @keywords internal
.getPosIndCytCombn <- function(
  ex,
  gateTbl,
  chnlPos,
  chnlNeg,
  chnlAlt,
  gateTypeCytPos,
  posCache = NULL,
  posByChnl = NULL
) {
  chnl <- unique(
    c(
      chnlPos,
      chnlNeg,
      chnlAlt
    )
  )

  if (is.null(posByChnl)) {
    posCache <- .getPosIndCache(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      posCache = posCache
    )

    posByChnl <- .getPosIndByChnl(
      ex = ex,
      gateTbl = gateTbl,
      chnl = chnl,
      gateTypeCytPos = gateTypeCytPos,
      posCache = posCache
    )
  }

  chnlPosIndVecPos <- rep(
    TRUE,
    nrow(ex)
  )

  for (chnlCurr in chnlPos) {
    chnlPosIndVecPos <-
      chnlPosIndVecPos &
      posByChnl[[as.character(chnlCurr)]]
  }

  if (length(chnlNeg) > 0L) {
    chnlNegIndVecPos <- rep(
      FALSE,
      nrow(ex)
    )

    for (chnlCurr in chnlNeg) {
      chnlNegIndVecPos <-
        chnlNegIndVecPos |
        posByChnl[[as.character(chnlCurr)]]
    }

    chnlPosIndVecPos <-
      chnlPosIndVecPos &
      !chnlNegIndVecPos
  }

  chnlPosIndVecPos
}
