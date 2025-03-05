get_batch_list <- function(fn_tbl_info,
                           col_grp,
                           col_stim,
                           uns_chr,
                           col_n_cell,
                           min_cell) {
  fn_tbl_info[["row_number"]] <- seq_len(nrow(fn_tbl_info))
  grp_vec <- fn_tbl_info[[col_grp[[1]]]]
  for (i in seq_along(col_grp)[-1]) {
    grp_vec <- paste0(grp_vec, "_", fn_tbl_info[[col_grp[[i]]]])
  }
  grp_vec_unique <- unique(grp_vec)
  out_list <- lapply(grp_vec_unique, function(grp) {
    sel_vec_ind <- which(grp_vec == grp)
    sel_vec_stim <- fn_tbl_info[[col_stim]][sel_vec_ind]
    if (!uns_chr %in% sel_vec_stim) {
      return(NULL)
    }
    sel_vec_n_cell <- fn_tbl_info[[col_n_cell]][sel_vec_ind]
    sel_vec_ind <- sel_vec_ind[sel_vec_n_cell >= min_cell]
    if (length(sel_vec_ind) %in% c(0, 1)) {
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
    c(
      setdiff(fn_tbl_sel[["row_number"]], row_number_uns),
      row_number_uns
    )
  })
  out_list[!vapply(out_list, is.null, logical(1))]
}
