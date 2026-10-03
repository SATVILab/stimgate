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
  attrsToKeep <- c(
    "ind",
    "indUns",
    "isUns",
    "chnlCut",
    "batch",
    "popGate",
    "probGMin"
  )
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
    for (nm in intersect(attrsToKeep, names(attrList))) {
      attr(cutTbl, nm) <- attrList[[nm]]
    }
    cutTbl
  }) |>
    stats::setNames(as.character(ind))
}


# Get axis labels from annotated data frame
#' @keywords internal
.getLabs <- function(.data, chnlCut, high = NULL) {
  force(.data)
  adfData <- flowWorkspace::gh_pop_get_data(.data) |>
    flowCore::parameters() |>
    flowCore::pData()

  descMap <- stats::setNames(as.character(adfData[["desc"]]), adfData[["name"]])
  descMap[chnlCut]
}


#' @keywords internal
.combineCp <- function(cp, gateCombn) {
  purrr::map(gateCombn, function(gateCombnCurr) {
    if (all(is.na(cp))) {
      return(stats::setNames(cp, names(cp)))
    }
    if (is.null(gateCombnCurr) || gateCombnCurr %in% c("no", "prejoin")) {
      return(cp)
    }
    val <- switch(gateCombnCurr,
      min = min(cp, na.rm = TRUE),
      mean = mean(cp, na.rm = TRUE),
      trim20 = mean(cp, trim = 0.2, na.rm = TRUE),
      median = stats::median(cp, na.rm = TRUE),
      max = max(cp, na.rm = TRUE),
      NULL
    )
    if (!is.null(val)) {
      stats::setNames(rep(val, length(cp)), names(cp))
    }
  }) |>
    stats::setNames(gateCombn)
}
