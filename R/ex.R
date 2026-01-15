str_detect_any <- function(string, pattern) {
  vapply(
    pattern,
    function(pattern_curr) grepl(pattern_curr, string),
    logical(1)
  ) |>
    any()
}

#' @keywords internal
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


#' @keywords internal
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
      path_project = path_project,
      save = TRUE
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

#' @keywords internal
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

#' @keywords internal
.get_ex_old <- function(pop,
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

#' @keywords internal
.get_ex_old_chnl_read_ind <- function(chnl,
                                      ind,
                                      pop,
                                      path_project) {
  path_chnl <- .get_ex_chnl_path(chnl, ind, pop, path_project)
  readRDS(path_chnl)
}

#' @keywords internal
.get_ex_new <- function(.data,
                        pop,
                        chnl,
                        ind,
                        path_project,
                        save) {
  # get expression information as a tibble
  # get .data
  fr <- flowWorkspace::gh_pop_get_data(.data, y = pop)
  ex <- flowCore::exprs(fr)[, chnl, drop = FALSE] |>
    tibble::as_tibble()
  .get_ex_new_chnl_save(
    ex = ex,
    ind = ind,
    pop = pop,
    path_project = path_project,
    save = save
  )
  ex
}

#' @keywords internal
.get_ex_new_chnl_save <- function(ex,
                                  ind,
                                  pop,
                                  path_project,
                                  save) {
  if (!save) {
    return(invisible(FALSE))
  }
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

#' @keywords internal
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

#' @keywords internal
.get_ex_chnl_path <- function(chnl, ind, pop, path_project) {
  file.path(
    .get_ex_chnl_path_dir(ind, pop, path_project),
    paste0("chnl_", chnl, ".rds")
  )
}
#' @keywords internal
.get_ex_chnl_path_dir <- function(ind, pop, path_project) {
  file.path(
    path_project,
    "sample_data",
    paste0("pop_", pop),
    paste0("ind_", ind)
  )
}

#' @keywords internal
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

#' @keywords internal
.get_cut <- function(ex) {
  ex[[attr(ex, "chnl_cut")]]
}

#' @keywords internal
.get_batch_ex <- function(ex) {
  attr(ex, "batch")
}

#' @keywords internal
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

#' @keywords internal
.get_batch <- function(ind, ind_batch_list) {
  has_ind <- vapply(
    ind_batch_list, function(x) ind %in% x, logical(1)
  )
  names(ind_batch_list)[has_ind]
}

#' @title Read saved expression data from project
#' @description Read channel expression vectors saved under a project's sample_data
#'   directory and return them as a tibble with sample metadata columns.
#' @param path_project character Path to project.
#' @param pop character or NULL Population name(s). Default is detected from project sample_data.
#' @param ind character or NULL Index/indices of samples. Default is detected from project sample_data.
#' @param chnl character or NULL Channel name(s) to return. Default is detected from project sample_data.
#' @param bias logical. Whether to add bias to unstimulated sample used in the gating. Default is `FALSE`.
#' @param exc_min logical. Whether to exclude cells with the minimum expression for any channels.
#' @return A tibble with columns `pop`, `ind` and one column per requested channel. Rows correspond to cells.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "sample_data", "POP1", "ind_1"), recursive = TRUE)
#' saveRDS(c(1, 2, 3), file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC1.rds"))
#' saveRDS(c(4, 5, 6), file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC2.rds"))
#' stimgate_data_get_ex(tmp)
#' stimgate_data_get_ex(tmp, chnl = "BC1")
#' }
#' @export
stimgate_data_get_ex <- function(path_project,
                                 pop = NULL,
                                 ind = NULL,
                                 chnl = NULL,
                                 bias = FALSE,
                                 exc_min = FALSE,
                                 combn_exc = NULL,
                                 cyt_pos = FALSE,
                                 gate_type_cyt_pos = "cyt",
                                 gate_type_single_pos = "single",
                                 mult = FALSE,
                                 gate_uns_method = "min") {
  .assert_string(path_project)
  pop <- pop %||% .get_ex_project_pop(path_project)
  .assert_string_vector(pop)
  purrr::map_df(pop, function(pop_curr) {
    ind <- ind %||% .get_ex_project_ind(path_project, pop_curr)
    .assert_string_vector(ind)
    purrr::map_df(ind, function(ind_curr) {
      chnl <- chnl %||%
        .get_ex_project_chnl(path_project, pop_curr, ind_curr)
      .assert_string_vector(chnl)
      ex <- .get_ex_old(
        pop = pop_curr,
        chnl = chnl,
        ind = ind_curr,
        path_project = path_project
      )
      ex <- .data_get_ex_exc_min(ex, exc_min)
      ex <- .data_get_ex_cyt_pos(
        ex,
        cyt_pos = cyt_pos,
        chnl = chnl,
        combn_exc = combn_exc,
        gate_type_cyt_pos = gate_type_cyt_pos,
        gate_type_single_pos = gate_type_single_pos,
        mult = mult
      )
      ex <- .data_get_ex_bias(
        ex,
        ind = ind_curr,
        path_project = path_project,
        bias = bias
      )
      meta_df <- tibble::tibble(
        pop = pop_curr,
        ind = ind_curr
      )
      tibble::as_tibble(cbind(meta_df, ex))
    })
  })
}

#' @keywords internal
.get_ex_project_pop <- function(path_project) {
  .assert_string(path_project)
  pop_vec <- list.dirs(
    file.path(path_project, "sample_data"),
    recursive = FALSE
  ) |>
    basename() |>
    sub("^pop_(.*)$", "\\1", x = _)
  .assert_string_vector(pop_vec)
  pop_vec
}

#' @keywords internal
.get_ex_project_ind <- function(path_project, pop = NULL) {
  pop <- pop %||% .get_ex_project_pop(path_project)
  pop <- pop[[1]]
  .assert_string(pop)
  ind_vec <- list.dirs(
    file.path(path_project, "sample_data", paste0("pop_", pop)),
    recursive = FALSE
  ) |>
    basename() |>
    sub("^ind_", "", x = _)
  .assert_string_vector(ind_vec)
  ind_vec
}

#' @keywords internal
.get_ex_project_chnl <- function(path_project, pop = NULL, ind = NULL) {
  pop <- pop %||% .get_ex_project_pop(path_project)
  pop <- pop[[1]]
  .assert_string(pop)
  ind <- ind %||% .get_ex_project_ind(path_project, pop)
  ind <- ind[[1]]
  .assert_string(ind)
  path_chnl_dir <- file.path(
    path_project,
    "sample_data",
    paste0("pop_", pop),
    paste0("ind_", ind)
  )
  .assert_string(path_chnl_dir)
  chnl_vec <- list.files(path_chnl_dir) |>
    sub("^chnl_(.*)\\.rds$", "\\1", x = _)
  .assert_string_vector(chnl_vec)
  chnl_vec
}

#' @keywords internal
.data_get_ex_bias <- function(ex,
                              ind,
                              path_project,
                              bias) {
  if (!bias) {
    return(ex)
  }
  ind_batch_list <- stimgate_meta_read_batch_list(path_project)
  ind_uns <- .get_ind_uns(ind, ind_batch_list)
  # only apply bias to unstim
  is_uns <- ind == ind_uns
  if (!is_uns) {
    return(ex)
  }
  # apply bias
  chnl_list <- stimgate_meta_read_settings_chnls(path_project)
  for (chnl in colnames(ex)) {
    bias <- chnl_list[[chnl]][["bias_uns"]]
    ex[[chnl]] <- ex[[chnl]] + bias
  }
  ex
}

#' @keywords internal
.data_get_ex_exc_min <- function(ex, exc_min) {
  if (!exc_min) {
    return(ex)
  }
  cn_vec <- setdiff(colnames(ex), c("pop", "ind"))
  min_vec <- vapply(
    cn_vec,
    function(x) min(ex[[x]], na.rm = TRUE),
    numeric(1)
  ) |>
    stats::setNames(cn_vec)
  for (cn in cn_vec) {
    inc_vec <- ex[[cn]] > min_vec[[cn]]
    ex <- ex[inc_vec, ]
  }
  ex
}

#' @keywords internal
.data_get_ex_cyt_pos <- function(ex,
                                 cyt_pos,
                                 chnl = NULL,
                                 combn_exc = NULL,
                                 gate_type_cyt_pos = "cyt",
                                 gate_type_single_pos = "single",
                                 mult = FALSE) {
  if (!cyt_pos) {
    return(ex)
  }
  cn_vec <- setdiff(colnames(ex), c("pop", "ind"))
  exc_vec <- vapply(
    cn_vec,
    function(x) {
      is_cyt <- str_detect_any(
        x,
        pattern = gate_type_cyt_pos
      )
      is_single <- str_detect_any(
        x,
        pattern = gate_type_single_pos
      )
      if (is_cyt) {
        return(TRUE)
      }
      if (is_single) {
        return(FALSE)
      }
      stop(
        "Channel name '", x,
        "' does not contain either '",
        gate_type_cyt_pos,
        "' or '",
        gate_type_single_pos,
        "' to determine its gate type."
      )
    },
    logical(1)
  ) |>
    stats::setNames(cn_vec)
  exc_mat <- vapply(
    combn_exc,
    function(x) exc_vec[[x]],
    logical(n = 1)
  ) |>
    t() |>
    matrix(nrow = length(combn_exc), ncol = length(cn_vec)) |>
    stats::setDimnames(
      list(
        NULL,
        cn_vec
      )
    )
  if (mult) {
    inc_vec <- apply(exc_mat, 2, all)
  } else {
    inc_vec <- apply(exc_mat, 2, any)
  }
  ex[, inc_vec, drop = FALSE]
}
