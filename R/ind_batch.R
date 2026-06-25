#' @title Generate a batch list of sample indices
#' @description Groups sample rows by batch/donor identifiers, screens out samples 
#'   falling below a minimum cell count threshold, and structures the output so that 
#'   the unstimulated control index is always positioned as the final element of each batch.
#' @param fnTblInfo data.frame. Sample metadata containing annotations.
#' @param colGrp character vector. One or more column names used to define batches/groups.
#' @param colStim character. Column name containing stimulation identifiers.
#' @param unsChr character. Name/string specifying the unstimulated control sample.
#' @param colNCell character. Column name containing the cell count for each sample.
#' @param minCell numeric. Minimum number of cells required to retain a sample.
#' @return A named list where each element contains a numeric vector of sample 
#'   indices representing a batch, with the unstimulated control index at the end.
#' @export
getBatchList <- function(
  fnTblInfo,
  colGrp,
  colStim,
  unsChr,
  colNCell,
  minCell
) {
  fnTblInfo[["rowNumber"]] <- seq_len(nrow(fnTblInfo))
  
  # Construct group vector combining all grouping columns
  grpVec <- fnTblInfo[[colGrp[[1]]]]
  for (i in seq_along(colGrp)[-1]) {
    grpVec <- paste0(grpVec, "_", fnTblInfo[[colGrp[[i]]]])
  }
  
  grpVecUnique <- unique(grpVec)
  
  outList <- lapply(grpVecUnique, function(grp) {
    selVecInd <- which(grpVec == grp)
    selVecStim <- fnTblInfo[[colStim]][selVecInd]
    
    # Early return if no unstim present in the initial group slice
    if (!unsChr %in% selVecStim) {
      return(NULL)
    }
    
    # Filter by minimum cell count
    selVecNCell <- fnTblInfo[[colNCell]][selVecInd]
    selVecInd <- selVecInd[selVecNCell >= minCell]
    
    # Batch must have at least one stimulated and one unstimulated sample
    if (length(selVecInd) <= 1L) {
      return(NULL)
    }
    
    selVecStimFinal <- fnTblInfo[[colStim]][selVecInd]
    if (!unsChr %in% selVecStimFinal) {
      return(NULL)
    }

    fnTblSel <- fnTblInfo[selVecInd, ]
    rowNumberUns <- fnTblSel[["rowNumber"]][
      fnTblSel[[colStim]] == unsChr
    ]
    
    # Order so that unstim indices always come last
    c(
      setdiff(fnTblSel[["rowNumber"]], rowNumberUns),
      rowNumberUns
    )
  }) |>
    stats::setNames(grpVecUnique)
  
  # Remove batches that were excluded during screening
  outList[!vapply(outList, is.null, logical(1))]
}
