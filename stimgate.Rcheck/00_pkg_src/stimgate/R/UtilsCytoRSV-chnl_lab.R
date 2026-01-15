#' @rdname chnl_lab
#'
#' @title Get markers and channels
#'
#' @description From a cytometry object (e.g. flowFrame or flowSet),
#' either get a character vector of markers
#' or channels (get_chnl and get_marker),
#' or get a named vector that converts
#' between channel names and marker names (e.g. chnl_to_marker).
#'
#' @param data object of class flowFrame, flowSet. Channel and corresponding
#' marker names are drawn from here.
#'
#' @details
#' Note that chnl_lab is equivalent to chnl_to_marker,
#' and marker_lab is equivalent to marker_to_chnl.
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
#' chnl_lab(fs[[1]])
#' }
#'
#' @aliases marker_lab, chnl_to_marker, marker_to_chnl, get_marker, get_chnl
#' @export
chnl_lab <- function(data) {
  adf <- switch(class(data)[1],
    "flowFrame" = flowCore::parameters(data)@data,
    "flowSet" = flowCore::parameters(data[[1]])@data,
    "cytoframe" = flowCore::parameters(data)@data,
    "cytoset" = flowCore::parameters(data[[1]])@data,
    stop("class of data not recognised")
  )

  lab_vec <- setNames(adf$desc, adf$name)
  for (i in seq_along(lab_vec)) {
    if (is.na(lab_vec[i])) {
      lab_vec[i] <- names(lab_vec)[i]
    }
  }

  lab_vec
}
