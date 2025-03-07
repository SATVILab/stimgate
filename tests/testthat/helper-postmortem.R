
get_fn_tbl_info_postmortem <- function(gs) {
  
  fn_vec_orig <- NULL
  n_cell_vec <- NULL
  for (i in seq_along(gs)) {
    fr <- flowWorkspace::gh_pop_get_data(gs[[i]])
    fn_vec_orig[[i]] <- basename(flowCore::keyword(fr, "FILENAME")[["FILENAME"]])
    n_cell_vec[[i]] <- nrow(flowCore::exprs(fr))
  }
  fn_vec_orig <- fn_vec_orig |> unlist()
  fn_vec <- tolower(fn_vec_orig)
  fn_vec <- gsub("_processed.*|-processed.*", "", fn_vec)
  fn_vec <- gsub(" ", "-", fn_vec)
  pid_vec <- gsub("^pm_", "", fn_vec) |>
    gsub("-.*|_.*", "", x = _)
  pid_vec <- paste0("pid_", pid_vec)
  # remove pid
  fn_vec <- gsub(
    "^pm_\\d+-+|^pm-\\d+-+|^pm-\\d+_+|^pm_\\d+_+", "", fn_vec
  )
  fn_vec <- gsub("-1$|_1$", "", fn_vec)
  fn_vec <- gsub("_", "-", fn_vec)
  location_vec <- gsub("-.*", "", fn_vec)
  stim_vec <- gsub(".*-", "", fn_vec)
  fn_tbl_info <- tibble::tibble(
    fn = fn_vec_orig,
    pid = pid_vec,
    location = location_vec,
    stim = stim_vec,
    n_cell_pop = n_cell_vec |> unlist()
  )
  fn_tbl_info
}

get_batch_list_postmortem <- function(fn_tbl_info,
                                      col_grp,
                                      col_stim,
                                      uns_chr,
                                      col_n_cell,
                                      min_cell_uns,
                                      min_cell_stim) {
  fn_tbl_info[["row_number"]] <- seq_len(nrow(fn_tbl_info))
  grp_vec <- fn_tbl_info[[col_grp[[1]]]]
  for (i in seq_along(col_grp)[-1]) {
    grp_vec <- paste0(grp_vec, "_", fn_tbl_info[[col_grp[[i]]]])
  }
  grp_vec_unique <- unique(grp_vec)
  out_list <- lapply(grp_vec_unique, function(grp) {
    sel_vec_ind <- which(grp_vec == grp)
    sel_vec_ind_uns <- sel_vec_ind[
      which(fn_tbl_info[[col_stim]][sel_vec_ind] == uns_chr)
      ]
    # check only one unstim sample
    if (length(sel_vec_ind_uns) > 1L) {
      stop(paste0("Multiple unstim samples for group ", grp))
    }
    # no unstim, ignore
    if (length(sel_vec_ind_uns) == 0L) {
      return(NULL)
    }
    # check sufficient unstim cells
    if (fn_tbl_info[[col_n_cell]][sel_vec_ind_uns] < min_cell_uns) {
      return(NULL)
    }
    sel_vec_ind_stim <- setdiff(sel_vec_ind, sel_vec_ind_uns)
    # remove stim samples with insufficient cells
    sel_vec_ind_stim <- sel_vec_ind_stim[
      fn_tbl_info[[col_n_cell]][sel_vec_ind_stim] >= min_cell_stim
    ]
    # return NULL if there is no stim sample
    # with sufficient cells
    if (length(sel_vec_ind_stim) == 0L) {
      return(NULL)
    }
    # return with stim indices first
    c(sel_vec_ind_stim, sel_vec_ind_uns)
  })
  # remove stim-and-unstim sets with no stim
  # with sufficient cells and/or no unstim with sufficient
  # cells
  out_list[!vapply(out_list, is.null, logical(1))]
}

.get_ind_uns <- function(ind, ind_batch_list) {
  has_ind <- vapply(
    ind_batch_list, function(x) ind %in% x, logical(1)
  )
  if (sum(has_ind) > 1L) {
    # this is an unstim, as it appears
    # in more than one batch
    return(ind)
  }
  ind_batch_list_match <- ind_batch_list[[has_ind]]
  ind_batch_list_match[[length(ind_batch_list_match)]]
}
