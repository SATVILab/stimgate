#' @title Generate a batch list of sample indices
#' @description Groups sample rows by batch/donor identifiers, screens out samples
#'   falling below a minimum cell count threshold, and structures the output so that
#'   the unstimulated control index is always positioned as the first element of each batch.
#' @param fnTblInfo data.frame. Sample metadata containing annotations.
#' @param colGrp character vector. One or more column names used to define batches/groups.
#' @param colStim character. Column name containing stimulation identifiers.
#' @param unsChr character. Name/string specifying the unstimulated control sample.
#' @param colNCell character. Column name containing the cell count for each sample.
#' @param minCell numeric. Minimum number of cells required to retain a sample.
#' @return A named list where each element contains a numeric vector of sample
#'   indices representing a batch, with the unstimulated control index at the beginning.
#' @examples
#' fnTblInfo <- data.frame(
#'   donor = c("d1", "d1", "d2", "d2", "d3", "d3"),
#'   stim = c("stim", "uns", "uns", "stim", "uns", "stim"),
#'   nCell = c(5000, 4000, 6000, 5500, 3000, 50)
#' )
#' # Donor d3 is dropped because its stimulated sample has too few cells
#' getBatchList(
#'   fnTblInfo,
#'   colGrp = "donor",
#'   colStim = "stim",
#'   unsChr = "uns",
#'   colNCell = "nCell",
#'   minCell = 100
#' )
#' @export
getBatchList <- function(
  fnTblInfo,
  colGrp,
  colStim,
  unsChr,
  colNCell,
  minCell
) {
  grpVec <- do.call(paste, c(fnTblInfo[colGrp], list(sep = "_")))
  grpVecUnique <- unique(grpVec)

  outList <- stats::setNames(lapply(grpVecUnique, function(grp) {
    selVecInd <- which(grpVec == grp)
    selVecInd <- selVecInd[fnTblInfo[[colNCell]][selVecInd] >= minCell]
    if (length(selVecInd) <= 1L) {
      return(NULL)
    }

    isUns <- fnTblInfo[[colStim]][selVecInd] == unsChr
    if (!any(isUns)) {
      return(NULL)
    }

    c(selVecInd[isUns], selVecInd[!isUns])
  }), grpVecUnique)

  outList[!vapply(outList, is.null, logical(1))]
}
