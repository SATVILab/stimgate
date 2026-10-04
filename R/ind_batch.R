#' @title Group samples into batches
#' @description Group metadata rows by donor or batch, keeping the unstimulated
#'   control first. Drop samples below `minCell` and groups with fewer than two
#'   remaining samples or no control.
#' @param fnTblInfo data.frame Sample metadata, one row per sample.
#' @param colGrp character vector Column names defining batches.
#' @param colStim character Column containing stimulation labels.
#' @param unsChr character Label identifying unstimulated controls.
#' @param colNCell character Column containing sample cell counts.
#' @param minCell numeric Minimum cell count to retain a sample.
#' @return A named list of integer row indices per batch. Names join group values
#'   with underscores; control indices precede stimulated indices.
#' @examples
#' samples <- data.frame(
#'   donor = c("d1", "d1", "d2", "d2"),
#'   stim = c("stim", "uns", "uns", "stim"),
#'   nCell = c(5000, 4000, 3000, 50)
#' )
#' # d2 is dropped: only its control meets the cell-count limit
#' getBatchList(samples, "donor", "stim", "uns", "nCell", minCell = 100)
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
