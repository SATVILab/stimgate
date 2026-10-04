#' @rdname chnlLab
#' @title Get channel-to-marker labels
#' @description Map channel names to marker labels in a cytometry object.
#'   Channels without a marker label use their channel name.
#' @param data flowFrame, flowSet, GatingSet, GatingHierarchy, cytoframe or
#'   cytoset Cytometry data. For sets, labels come from the first sample.
#' @return A character vector of marker labels, named by channel.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' chnlLab(gs)
#' @export
chnlLab <- function(data) {
  adf <- switch(class(data)[1],
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
  isNa <- is.na(labVec)
  labVec[isNa] <- names(labVec)[isNa]

  labVec
}
