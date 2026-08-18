#' @rdname chnlLab
#'
#' @title Get markers and channels
#'
#' @description From a cytometry object (e.g. flowFrame or flowSet),
#' either get a character vector of markers
#' or channels (getChnl and getMarker),
#' or get a named vector that converts
#' between channel names and marker names (e.g. chnlToMarker).
#'
#' @param data object of class flowFrame, flowSet. Channel and corresponding
#' marker names are drawn from here.
#'
#' @details
#' Note that chnlLab is equivalent to chnlToMarker,
#' and markerLab is equivalent to markerToChnl.
#'
#' @return A named character vector.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("flowCore", quietly = TRUE)) {
#'   exprs <- matrix(
#'     seq_len(8),
#'     ncol = 2,
#'     dimnames = list(NULL, c("FSC-A", "FL1-H"))
#'   )
#'   ff <- flowCore::flowFrame(exprs)
#'
#'   # Get channel to marker mapping
#'   chnlLab(ff)
#' }
#' }
#'
#' @aliases markerLab, chnlToMarker, markerToChnl, getMarker, getChnl
#' @export
chnlLab <- function(data) {
  adf <- switch(
    class(data)[1],
    "GatingSet" = {
      gh <- data[[flowWorkspace::sampleNames(data)[1]]]
      fr <- flowWorkspace::gh_pop_get_data(gh)
      flowCore::parameters(fr)@data
    },
    "GatingHierarchy" = {
      fr <- flowWorkspace::gh_pop_get_data(data)
      flowCore::parameters(fr)@data
    },
    "flowFrame" = flowCore::parameters(data)@data,
    "flowSet" = flowCore::parameters(data[[1]])@data,
    "cytoframe" = flowCore::parameters(data)@data,
    "cytoset" = flowCore::parameters(data[[1]])@data,
    stop("classOfDataNotRecognised")
  )

  labVec <- stats::setNames(adf$desc, adf$name)
  for (i in seq_along(labVec)) {
    if (is.na(labVec[i])) {
      labVec[i] <- names(labVec)[i]
    }
  }

  labVec
}
