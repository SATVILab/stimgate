.acs_assert_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " not found at: ", path)
  }
  invisible(path)
}

.acs_assert_dir <- function(path, description) {
  if (!dir.exists(path)) {
    stop(description, " not found at: ", path)
  }
  invisible(path)
}

.acs_get_expr <- function(gs, ind, chnl, pop = "root") {
  gh <- gs[[as.integer(ind)]]
  ff <- flowWorkspace::gh_pop_get_data(gh, pop)
  ex <- flowCore::exprs(ff)
  if (!chnl %in% colnames(ex)) {
    stop("Channel '", chnl, "' is not present in sample index ", ind, ".")
  }
  as.numeric(ex[, chnl])
}

.acs_estimate_from_threshold <- function(x_stim, x_uns, threshold) {
  threshold <- suppressWarnings(as.numeric(threshold))[1]
  n_stim <- length(x_stim)
  n_uns <- length(x_uns)

  if (!is.finite(threshold)) {
    return(tibble::tibble(
      threshold = NA_real_,
      nCellStim = n_stim,
      nCellUns = n_uns,
      nPosStim = NA_integer_,
      nPosUns = NA_integer_,
      propStim = NA_real_,
      propUns = NA_real_,
      propRespEst = NA_real_
    ))
  }

  n_pos_stim <- sum(x_stim > threshold, na.rm = TRUE)
  n_pos_uns <- sum(x_uns > threshold, na.rm = TRUE)
  prop_stim <- n_pos_stim / n_stim
  prop_uns <- n_pos_uns / n_uns

  tibble::tibble(
    threshold = threshold,
    nCellStim = n_stim,
    nCellUns = n_uns,
    nPosStim = n_pos_stim,
    nPosUns = n_pos_uns,
    propStim = prop_stim,
    propUns = prop_uns,
    propRespEst = prop_stim - prop_uns
  )
}

.acs_stim_pair_tbl <- function(batch_list, sample_metadata) {
  purrr::imap_dfr(batch_list, function(ind_batch, batch_nm) {
    ind_uns <- ind_batch[[1]]
    ind_stim <- ind_batch[-1]

    tibble::tibble(
      batch = batch_nm,
      indUns = ind_uns,
      indStim = ind_stim,
      sampleUns = sample_metadata$sample_name[ind_uns],
      sampleStim = sample_metadata$sample_name[ind_stim],
      conditionStim = sample_metadata$condition[ind_stim]
    )
  })
}
