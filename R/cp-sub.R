# Prepare expression data list with bias adjustment
# Optionally drops minimum-expression cells and shifts the cut channel by bias
# Returns a list of expression tables
#' @keywords internal
.prepareExListWithBiasAndNoise <- function(
  exList,
  ind,
  excMin,
  bias = 0
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
