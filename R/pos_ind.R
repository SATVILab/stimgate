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

  hasGateCyt <- "gateCyt" %in% colnames(gateTbl)

  for (chnlCurr in chnl) {
    gateTblChnlInd <- which(
      as.character(gateTbl$chnl) == chnlCurr
    )

    # Unstimulated samples have saved expression but no stimulation gates.
    if (length(gateTblChnlInd) == 0L) {
      posCache$base[[chnlCurr]] <- rep(FALSE, nrow(ex))
      if (hasGateCyt) {
        posCache$cyt[[chnlCurr]] <- rep(FALSE, nrow(ex))
      }
      next
    }

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
  gateList <- switch(gateType,
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


# Get positivity for at least one channel other than the current channel,
# with the NA semantics of a logical OR across those channels
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

  out <- count$nTrue - as.integer(!is.na(posCurr) & posCurr) > 0L
  out[!out & count$nNa - as.integer(is.na(posCurr)) > 0L] <- NA
  out
}


# Identify cells that are positive for at least two cytokines
#' @keywords internal
.getPosIndMult <- function(
  ex,
  gateTbl,
  chnl = NULL,
  gateTypeCytPos,
  posCache = NULL
) {
  gateTypeCytPos <- match.arg(gateTypeCytPos, c("base", "cyt"))

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

  if (gateTypeCytPos == "base") {
    # Any NA among the contributing channels produces NA.
    out <- baseCount$nTrue >= 2L
    out[baseCount$nNa > 0L] <- NA
    return(out)
  }

  # In cyt mode a cell is multifunctional when, for at least one
  # requested cytokine, either:
  #
  # 1. that cytokine clears its base threshold and another cytokine
  #    clears its cyt+ threshold, or
  # 2. that cytokine clears its cyt+ threshold and another cytokine
  #    clears its base threshold.
  cytCount <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = chnl,
    gateType = "cyt"
  )

  posList <- lapply(chnl, function(chnlCurr) {
    (posCache$base[[chnlCurr]] &
      .getPosIndCacheAnyExcept(posCache, cytCount, chnlCurr, "cyt")) |
      (.getPosIndCacheGet(posCache, chnlCurr, "cyt") &
        .getPosIndCacheAnyExcept(posCache, baseCount, chnlCurr, "base"))
  })

  Reduce("|", posList, rep(FALSE, nrow(ex)))
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
  gateTypeCytPos <- match.arg(gateTypeCytPos, c("base", "cyt"))

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

  if (gateTypeCytPos == "base") {
    return(posCache$base[chnl])
  }

  # The cyt+ rule for one cytokine is:
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

  lapply(chnl, function(chnlCurr) {
    posCache$base[[chnlCurr]] |
      (.getPosIndCacheGet(posCache, chnlCurr, "cyt") &
        .getPosIndCacheAnyExcept(posCache, baseCount, chnlCurr, "base"))
  }) |>
    stats::setNames(chnl)
}


# Identify cells positive for any other cytokine than one specified channel
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

  .getPosIndButSinglePosByChnl(
    ex = ex,
    gateTbl = gateTbl,
    chnl = unique(c(chnlSingleExc, chnl)),
    gateTypeCytPos = gateTypeCytPos,
    posCache = posCache
  )[[as.character(chnlSingleExc)]]
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
  gateTypeCytPos <- match.arg(gateTypeCytPos, c("base", "cyt"))

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
      gateTypeCytPos = "cyt",
      posCache = posCache
    )
  } else {
    FALSE
  }

  lapply(chnl, function(chnlCurr) {
    .getPosIndCacheAnyExcept(posCache, baseCount, chnlCurr, "base") |
      posVecMultiCyt
  }) |>
    stats::setNames(chnl)
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
    chnlAlt <- unique(gateTbl$chnl)
  }

  # chnlAlt supplies context for the cyt+ rule but is not itself tested.
  posByChnl <- .getPosIndByChnl(
    ex = ex,
    gateTbl = gateTbl,
    chnl = unique(c(as.character(chnl), as.character(chnlAlt))),
    gateTypeCytPos = gateTypeCytPos,
    posCache = posCache
  )

  Reduce("|", posByChnl[as.character(chnl)], rep(FALSE, nrow(ex)))
}


# Get cell membership for one exact cytokine combination
#' @keywords internal
.getPosIndCytCombn <- function(
  ex,
  gateTbl,
  chnlPos,
  chnlNeg,
  gateTypeCytPos,
  posCache = NULL,
  posByChnl = NULL
) {
  if (is.null(posByChnl)) {
    posByChnl <- .getPosIndByChnl(
      ex = ex,
      gateTbl = gateTbl,
      chnl = c(chnlPos, chnlNeg),
      gateTypeCytPos = gateTypeCytPos,
      posCache = posCache
    )
  }

  posAll <- Reduce(
    "&",
    posByChnl[as.character(chnlPos)],
    rep(TRUE, nrow(ex))
  )
  posAny <- Reduce(
    "|",
    posByChnl[as.character(chnlNeg)],
    rep(FALSE, nrow(ex))
  )

  posAll & !posAny
}
