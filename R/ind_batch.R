#' @title Generate a batch list of sample indices
#' @description Groups sample rows by batch/donor identifiers, screens out samples 
#'   falling below a minimum cell count threshold, and structures the output so that 
#'   the unstimulated control index is always positioned as the final element of each batch.
#' @param fn_tbl_info data.frame. Sample metadata containing annotations.
#' @param col_grp character vector. One or more column names used to define batches/groups.
#' @param col_stim character. Column name containing stimulation identifiers.
#' @param uns_chr character. Name/string specifying the unstimulated control sample.
#' @param col_n_cell character. Column name containing the cell count for each sample.
#' @param min_cell numeric. Minimum number of cells required to retain a sample.
#' @return A named list where each element contains a numeric vector of sample 
#'   indices representing a batch, with the unstimulated control index at the end.
#' @export
get_batch_list <- function(
  fn_tbl_info,
  col_grp,
  col_stim,
  uns_chr,
  col_n_cell,
  min_cell
) {
  fn_tbl_info[["row_number"]] <- seq_len(nrow(fn_tbl_info))
  
  # Construct group vector combining all grouping columns
  grp_vec <- fn_tbl_info[[col_grp[[1]]]]
  for (i in seq_along(col_grp)[-1]) {
    grp_vec <- paste0(grp_vec, "_", fn_tbl_info[[col_grp[[i]]]])
  }
  
  grp_vec_unique <- unique(grp_vec)
  
  out_list <- lapply(grp_vec_unique, function(grp) {
    sel_vec_ind <- which(grp_vec == grp)
    sel_vec_stim <- fn_tbl_info[[col_stim]][sel_vec_ind]
    
    # Early return if no unstim present in the initial group slice
    if (!uns_chr %in% sel_vec_stim) {
      return(NULL)
    }
    
    # Filter by minimum cell count
    sel_vec_n_cell <- fn_tbl_info[[col_n_cell]][sel_vec_ind]
    sel_vec_ind <- sel_vec_ind[sel_vec_n_cell >= min_cell]
    
    # Batch must have at least one stimulated and one unstimulated sample
    if (length(sel_vec_ind) <= 1L) {
      return(NULL)
    }
    
    sel_vec_stim_final <- fn_tbl_info[[col_stim]][sel_vec_ind]
    if (!uns_chr %in% sel_vec_stim_final) {
      return(NULL)
    }

    fn_tbl_sel <- fn_tbl_info[sel_vec_ind, ]
    row_number_uns <- fn_tbl_sel[["row_number"]][
      fn_tbl_sel[[col_stim]] == uns_chr
    ]
    
    # Order so that unstim indices always come last
    c(
      setdiff(fn_tbl_sel[["row_number"]], row_number_uns),
      row_number_uns
    )
  }) |>
    stats::setNames(grp_vec_unique)
  
  # Remove batches that were excluded during screening
  out_list[!vapply(out_list, is.null, logical(1))]
}
