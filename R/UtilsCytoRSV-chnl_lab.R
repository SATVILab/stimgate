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
#' # Create example flowFrame-like data structure
#' data(GvHD, package = "flowCore")
#' fs <- GvHD[1:2]
#'
#' # Get channel to marker mapping
#' chnlLab(fs[[1]])
#' }
#'
#' @aliases markerLab, chnlToMarker, markerToChnl, getMarker, getChnl
#' @export
chnlLab <- function(data) {
  adf <- switch(
    class(data)[1],
    "flowFrame" = flowCore::parameters(data)@data,
    "flowSet" = flowCore::parameters(data[[1]])@data,
    "cytoframe" = flowCore::parameters(data)@data,
    "cytoset" = flowCore::parameters(data[[1]])@data,
    stop("class of data not recognised")
  )

  labVec <- setNames(adf$desc, adf$name)
  for (i in seq_along(labVec)) {
    if (is.na(labVec[i])) {
      labVec[i] <- names(labVec)[i]
    }
  }

  labVec
}
