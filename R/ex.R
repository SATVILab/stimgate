str_detect_any <- function(string, pattern) {
  vapply(
    pattern,
    function(pattern_curr) grepl(pattern_curr, string),
    logical(1)
  ) |>
    any()
}

.get_ex_list <- function(.data,
                         ind_batch,
                         batch,
                         pop,
                         chnl_cut,
                         extra_chnl = NULL,
                         path_project) {
  is_path_given <- is.character(path_project) && nzchar(path_project)
  if (!is_path_given) {
    stop("path_project must be a non-empty character string.")
  }
  # get expression .data for each batch
  lapply(ind_batch, function(i) {
    .get_ex(
      .data = .data[[i]],
      pop = pop,
      chnl_cut = chnl_cut,
      extra_chnl = extra_chnl,
      ind = i,
      # specify corresponding unstim
      ind_uns = ind_batch[length(ind_batch)],
      batch = batch,
      path_project = path_project
    )
  }) |>
    stats::setNames(as.character(ind_batch))
}


.get_ex <- function(.data,
                    pop,
                    chnl_cut,
                    ind,
                    ind_uns,
                    batch,
                    extra_chnl = NULL,
                    path_project) {
  # collect all the channels we need
  # get expression information as a tibble
  # get .data
  all_saved <- .get_ex_check_chnl_saved(
    chnl = c(chnl_cut, extra_chnl),
    ind = ind,
    pop = pop,
    path_project = path_project
  )
  ex <- if (all_saved) {
    .get_ex_old(
      .data = .data,
      pop = pop,
      chnl = c(chnl_cut, extra_chnl),
      ind = ind,
      path_project = path_project
    )
  } else {
    .get_ex_new(
      .data = .data,
      chnl = c(chnl_cut, extra_chnl),
      ind = ind,
      pop = pop,
      path_project = path_project
    )
  }
  .get_ex_add_attributes(
    ex = ex,
    ind = ind,
    ind_uns = ind_uns,
    batch = batch,
    chnl_cut = chnl_cut,
    pop = pop
  )
}

.get_ex_check_chnl_saved <- function(chnl,
                                     ind,
                                     pop,
                                     path_project) {
  path_chnl_dir <- .get_ex_chnl_path_dir(ind, pop, path_project)
  if (!dir.exists(path_chnl_dir)) {
    return(FALSE)
  }
  fn_vec <- list.files(path_chnl_dir)
  req_vec <- paste0("chnl_", chnl, ".rds")
  all(req_vec %in% fn_vec)
}

.get_ex_old <- function(.data,
                        pop,
                        chnl,
                        ind,
                        path_project) {
  # get expression information as a tibble
  # get .data
  ex <- tibble::tibble(
    V1 = .get_ex_old_chnl_read_ind(chnl[[1]], ind, pop, path_project)
  )
  names(ex) <- chnl[[1]]
  chnl_remaining <- chnl[-1]
  for (i in seq_along(chnl_remaining)) {
    chnl_curr <- chnl_remaining[[i]]
    ex[[chnl_curr]] <-
      .get_ex_old_chnl_read_ind(chnl_curr, ind, pop, path_project)
  }
  ex
}

.get_ex_old_chnl_read_ind <- function(chnl,
                                      ind,
                                      pop,
                                      path_project) {
  path_chnl <- .get_ex_chnl_path(chnl, ind, pop, path_project)
  readRDS(path_chnl)
}

.get_ex_new <- function(.data,
                        pop,
                        chnl,
                        ind,
                        path_project) {
  # get expression information as a tibble
  # get .data
  fr <- flowWorkspace::gh_pop_get_data(.data, y = pop)
  ex <- flowCore::exprs(fr)[, chnl, drop = FALSE] |>
    tibble::as_tibble()
  .get_ex_new_chnl_save(
    ex = ex,
    ind = ind,
    pop = pop,
    path_project = path_project
  )
  ex
}

.get_ex_new_chnl_save <- function(ex,
                                  ind,
                                  pop,
                                  path_project) {
  for (chnl_curr in colnames(ex)) {
    .get_ex_new_chnl_save_ind(
      ex = ex,
      chnl = chnl_curr,
      ind = ind,
      pop = pop,
      path_project = path_project
    )
  }
  invisible(TRUE)
}

.get_ex_new_chnl_save_ind <- function(ex,
                                      chnl,
                                      ind,
                                      pop,
                                      path_project) {
  path_chnl <- .get_ex_chnl_path(chnl, ind, pop, path_project)
  if (file.exists(path_chnl)) {
    return(invisible(FALSE))
  }
  if (!dir.exists(dirname(path_chnl))) {
    dir.create(dirname(path_chnl), recursive = TRUE)
  }
  saveRDS(ex[[chnl]], path_chnl)
}

.get_ex_chnl_path <- function(chnl, ind, pop, path_project) {
  file.path(
    .get_ex_chnl_path_dir(ind, pop, path_project),
    paste0("chnl_", chnl, ".rds")
  )
}
.get_ex_chnl_path_dir <- function(ind, pop, path_project) {
  file.path(
    path_project,
    "sample_data",
    pop,
    paste0("ind_", ind)
  )
}

.get_ex_add_attributes <- function(ex,
                                   ind,
                                   ind_uns,
                                   batch,
                                   chnl_cut,
                                   pop) {
  attr(ex, "ind") <- ind |> as.character()
  attr(ex, "ind_uns") <- ind_uns |> as.character()
  attr(ex, "is_uns") <- ind == ind_uns
  attr(ex, "chnl_cut") <- chnl_cut
  attr(ex, "batch") <- batch
  attr(ex, "pop") <- pop

  ex
}

.get_cut <- function(ex) {
  ex[[attr(ex, "chnl_cut")]]
}

.get_batch_ex <- function(ex) {
  attr(ex, "batch")
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
  has_ind <- which(has_ind)
  ind_batch <- ind_batch_list[has_ind] |>
    unlist()
  ind_batch[[length(ind_batch)]]
}

.get_batch <- function(ind, ind_batch_list) {
  has_ind <- vapply(
    ind_batch_list, function(x) ind %in% x, logical(1)
  )
  names(ind_batch_list)[has_ind]
}
