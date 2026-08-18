# Prepare expression data list with bias and noise adjustments
# Gets measurements for samples and applies bias and noise modifications
# Returns a list where each element is a numeric vector
#' @keywords internal
.prepareExListWithBiasAndNoise <- function(
  exList,
  ind,
  excMin,
  bias = 0,
  noiseSd = NULL
) {
  purrr::map(ind, function(indCurr) {
    cutTbl <- exList[[as.character(indCurr)]]
    attrList <- attributes(cutTbl)
    if (excMin) {
      nRowInit <- nrow(cutTbl)
      cutTbl <- cutTbl[
        .getCut(cutTbl) > min(.getCut(cutTbl)),
      ] # nolint
      nRowFin <- nrow(cutTbl)
      attr(cutTbl, "probGMin") <- nRowFin / nRowInit
    }
    cutTbl[[attr(cutTbl, "chnlCut")]] <- .getCut(cutTbl) + bias # nolint
    if (!is.null(noiseSd)) {
      cutTbl <- cutTbl[[attr(cutTbl, "chnlCut")]] +
        stats::rnorm(nrow(cutTbl), sd = noiseSd) # nolint
    }
    cutTbl |>
      .prepareExListWithBiasAndNoiseAddAttr(attrList)
  }) |>
    stats::setNames(as.character(ind))
}

.prepareExListWithBiasAndNoiseAddAttr <- function(ex, attrList) {
  attrVecNmOrig <- names(attrList)
  attrVecNmAdd <- c(
    "ind",
    "indUns",
    "isUns",
    "chnlCut",
    "batch",
    "popGate",
    "probGMin"
  )
  attrVecNmAdd <- intersect(attrVecNmAdd, attrVecNmOrig)
  for (i in seq_along(attrVecNmAdd)) {
    attr(ex, attrVecNmAdd[i]) <- attrList[[attrVecNmAdd[i]]]
  }
  ex
}


# Linear interpolation function
# Interpolates a value from input and output vectors
#' @keywords internal
.interp <- function(val, x, y) {
  xLow <- x[x <= val] |> max()
  if (xLow == val) {
    return(y[which(x == xLow)])
  }
  xHigh <- x[x >= val] |> min()
  xLowInd <- which(x == xLow)
  xHighInd <- which(x == xHigh)

  y[xLowInd] +
    (val - xLow) * (y[xHighInd] - y[xLowInd]) / (xHigh - xLow)
}

#' @keywords internal
.getCpTg <- function(
  exList,
  chnlSettings,
  tgType,
  stage,
  pathProject
) {
  # get native StimGate tailgate cutpoint
  .debug("Getting tg cutpoint")

  # Extract explicit parameters cleanly from chnlSettings at local execution scope
  chnlCut <- chnlSettings$chnlCut
  gateCombn <- chnlSettings$gateCombn
  excMin <- chnlSettings$excMin
  minCell <- chnlSettings$minCell
  cpMin <- chnlSettings$cpMin
  bw <- chnlSettings$bw

  stageChnl <- file.path(stage, chnlCut)
  cpList <- list()

  .getNativeTailgateCp <- function(ex) {
    ex <- ex[!is.na(.getCut(ex)), , drop = FALSE]
    if (excMin) {
      ex <- ex[.getCut(ex) > min(.getCut(ex)), , drop = FALSE]
    }
    if (nrow(ex) < max(minCell, 5L)) {
      return(
        if (nrow(ex) == 0L) {
          NA_real_
        } else {
          max(
            cpMin,
            max(.getCut(ex)) +
              (max(.getCut(ex)) - min(.getCut(ex))) / 5
          )
        }
      )
    }

    bwEst <- try(
      suppressWarnings(ks::hpi(.getCut(ex), deriv.order = 1)),
      silent = TRUE
    )
    if (inherits(bwEst, "try-error") || !is.finite(bwEst)) {
      bwEst <- max(bw, 1e-8)
    }
    bwTg <- max(bw, bwEst)
    densityObj <- try(
      suppressWarnings(stats::density(
        x = .getCut(ex),
        bw = bwTg,
        n = 512L,
        from = min(.getCut(ex)),
        to = max(.getCut(ex))
      )),
      silent = TRUE
    )
    if (inherits(densityObj, "try-error") || is.null(densityObj)) {
      return(
        max(
          cpMin,
          max(.getCut(ex)) +
            (max(.getCut(ex)) - min(.getCut(ex))) / 5
        )
      )
    }

    tailgate <- .getCpTailgate(
      density = data.frame(x = densityObj$x, y = densityObj$y),
      peakX = densityObj$x[which.max(densityObj$y)],
      fraction = 1 / 200
    )
    cp <- suppressWarnings(as.numeric(tailgate$lowerBoundX)[1L])
    if (!is.finite(cp)) {
      cp <- max(
        cpMin,
        max(.getCut(ex)) +
          (max(.getCut(ex)) - min(.getCut(ex))) / 5
      )
    }

    cp
  }

  if ("prejoin" %in% gateCombn) {
    .debug("prejoin")
    indGate <- names(exList)[-1]
    ex <- dplyr::bind_rows(exList[indGate])
    ex <- ex[!is.na(.getCut(ex)), ]
    if (excMin) {
      ex <- ex[.getCut(ex) > min(.getCut(ex)), ]
    }
    if (nrow(ex) < max(minCell, 5)) {
      .intSaveNm(
        file.path(tgType, "cpTgPrejoinTooFewToCutReliably"),
        cpMin,
        paste0(indGate, collapse = ","),
        stageChnl,
        pathProject
      )
      cpVec <- stats::setNames(rep(NA, length(indGate)), indGate)
    } else {
      cp <- .getNativeTailgateCp(ex)
      .intSaveNm(
        file.path(tgType, "cpTgPrejoinInit"),
        cp,
        paste0(indGate, collapse = ","),
        stageChnl,
        pathProject
      )
      if (is.na(cp) || length(.getCut(ex)) < minCell) {
        .intSaveNm(
          file.path(tgType, "cpTgPrejoinIsNaOrTooFewToCutReliably"),
          cp,
          paste0(indGate, collapse = ","),
          stageChnl,
          pathProject
        )
        cp <- max(
          cpMin,
          max(.getCut(ex)) +
            (max(.getCut(ex)) - min(.getCut(ex))) / 5
        )
      }
      .intSaveNm(
        file.path(tgType, "cpTgPrejoinFinal"),
        cp,
        paste0(indGate, collapse = ","),
        stageChnl,
        pathProject
      )
      cpVec <- stats::setNames(rep(cp, length(indGate)), indGate)
    }

    cpList <- cpList |> append(list(prejoin = cpVec))
  }

  # get cutpoint if group method is not prejoin
  nonPrejoinCombnVec <- setdiff(gateCombn, "prejoin")

  if (length(nonPrejoinCombnVec) > 0) {
    .debug("non-prejoin")
    indGate <- names(exList)[-1]
    cpTgVec <- purrr::map_dbl(indGate, function(ind) {
      .debug("ind", ind)
      ex <- exList[[as.character(ind)]]
      cp <- .getNativeTailgateCp(ex)
      .intSaveNm(
        file.path(tgType, "cpTgIndInit"),
        cp,
        ind,
        stageChnl,
        pathProject
      )
      cp
    }) |>
      stats::setNames(indGate)

    .debug("combining thresholds") # nolint

    cpTgList <- .combineCp(
      cp = cpTgVec,
      gateCombn = gateCombn
    ) |>
      stats::setNames(gateCombn)

    .intSaveNm(
      file.path(tgType, "cpTgList"),
      cpTgList,
      names(cpTgList),
      stageChnl,
      pathProject
    )

    cpList <- cpList |> append(cpTgList)
  }
  .debug("Done tg cutpoint")

  cpList
}


# Get mid-probability cut
#' @keywords internal
.getCpPwmid <- function(highIndTbl, cpScp) {
  # get table to model - all values above changepoint
  modTbl <- highIndTbl |>
    dplyr::filter(chnlCut >= cpScp)

  # check if too little .data to model with here
  if (nrow(modTbl) < 5 || length(unique(highIndTbl$high)) == 1) {
    message(
      "cpPwmid set to cpScp + 5 due to too few obs above cpSc or too few postive obs"
    )
    return(cpScp + 5)
  }

  # model the probability of positivity above this point
  if (nrow(modTbl) < 40) {
    fitPw <- stats::glm(high ~ chnlCut, family = stats::binomial, .data = modTbl)
  } else {
    fitPw <- stats::glm(
      high ~ splines::ns(chnlCut, df = 3),
      family = stats::binomial,
      .data = modTbl
    )
  }

  # get predictions over a range of cut values
  predTbl <- tibble::tibble(
    chnlCut = seq(min(modTbl$chnlCut), max(modTbl$chnlCut))
  )
  predVec <- stats::predict(fitPw, predTbl, type = "response")
  predTbl <- predTbl |> dplyr::mutate(pred = predVec)

  # find prediction in between middle and max
  minVal <- mean(
    highIndTbl |> dplyr::filter(chnlCut < cpScp) |> dplyr::pull("high")
  )
  maxVal <- max(predTbl$pred)
  midProb <- mean(c(maxVal, minVal))

  # find the cutpoint that minimises the difference objective function
  optFunc <- function(x) {
    predTbl <- tibble::tibble(chnlCut = x)
    predVec <- stats::predict(fitPw, predTbl, type = "response")
    (predVec - midProb)^2
  }

  # initial search parameter
  initPar <- mean(c(cpScp, max(highIndTbl$chnlCut)))
  # optimisation
  stats::optim(
    initPar,
    fn = optFunc,
    method = "Brent",
    lower = cpScp,
    upper = max(highIndTbl$chnlCut)
  )$par
}


# Get axis labels from annotated data frame
#' @keywords internal
.getLabs <- function(.data, chnlCut, high = NULL) {
  force(.data)
  adfData <- flowWorkspace::gh_pop_get_data(.data) |>
    flowCore::parameters() |>
    flowCore::pData()

  if (!is.null(high)) {
    cutLab <- adfData[["desc"]][[which(adfData$name == chnlCut)]] |>
      stats::setNames(chnlCut)
    return(cutLab)
  }

  purrr::map_chr(chnlCut, function(cutCurr) {
    adfData[["desc"]][[which(adfData$name == cutCurr)]]
  }) |>
    stats::setNames(chnlCut)
}


#' @keywords internal
.combineCp <- function(cp, gateCombn) {
  purrr::map(gateCombn, function(gateCombnCurr) {
    if (all(purrr::map_lgl(cp, is.na))) {
      return(stats::setNames(cp, names(cp)))
    }
    if (is.null(gateCombnCurr) || gateCombnCurr %in% c("no", "prejoin")) {
      return(cp)
    }
    if (gateCombnCurr == "min") {
      return(stats::setNames(
        rep(
          min(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gateCombnCurr == "mean") {
      return(stats::setNames(
        rep(
          mean(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gateCombnCurr == "trim20") {
      return(stats::setNames(
        rep(
          mean(cp, trim = 0.2, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gateCombnCurr == "median") {
      return(stats::setNames(
        rep(
          stats::median(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
    if (gateCombnCurr == "max") {
      return(stats::setNames(
        rep(
          max(cp, na.rm = TRUE),
          length(cp)
        ),
        names(cp)
      ))
    }
  }) |>
    stats::setNames(gateCombn)
}
