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
                    path_project,
                    add_attributes = TRUE) {
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
    pop = pop,
    add_attributes = add_attributes
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
                                   pop,
                                   add_attributes) {
  if (!add_attributes) {
    return(ex)
  }
  attr(ex, "ind") <- ind |> as.character()
  attr(ex, "ind_uns") <- ind_uns |> as.character()
  attr(ex, "is_uns") <- ind == ind_uns
  attr(ex, "chnl_cut") <- chnl_cut
  attr(ex, "batch") <- batch
  attr(ex, "pop") <- pop

  ex
}

.get_ind <- function(ex) {
  attr(ex, "ind")
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
#' @description Read channel expression vectors saved under a project's
#'   sample_data directory and return them as a tibble with sample metadata
#'   columns.
#' @param path_project character Path to project.
#' @param .data GatingSet or NULL GatingSet object to extract expression data
#'   from. Default is NULL.
#' @param pop character or NULL Population name(s). Default is detected from
#'   project sample_data.
#' @param ind character or NULL Index/indices of samples. Default is detected
#'   from project sample_data.
#' @param chnl character or NULL Channel name(s) to return. Default is
#'   detected from project sample_data.
#' @param marker character or NULL Marker name(s) to return. Cannot be
#'   specified with `chnl`. Default is NULL.
#' @param bias logical Whether to add bias to unstimulated sample used in the
#'   gating. Default is `FALSE`.
#' @param exc_min logical Whether to exclude cells with the minimum
#'   expression for any channels. Default is FALSE.
#' @param combn_exc list or NULL Combinations of channels to exclude. Default
#'   is NULL.
#' @param chnl_gate character or NULL Channel name(s) to use for gating.
#'   Cannot be specified with `marker_gate`. Default is NULL.
#' @param marker_gate character or NULL Marker name(s) to use for gating.
#'   Cannot be specified with `chnl_gate`. Default is NULL.
#' @param gate_type_cyt_pos character Gate type to use for cytokine-positive
#'   cells. Default is "cyt".
#' @param gate_type_single_pos character Gate type to use for single-positive
#'   cells. Default is "single".
#' @param mult logical Whether to return only multi-functional cells (positive
#'   for multiple markers). Default is FALSE.
#' @param gate_uns_method character Method for gating unstimulated cells.
#'   Default is "min".
#' @param trans_fn function or NULL Transformation function to apply to
#'   expression values. Default is NULL.
#' @param trans_chnl character or NULL Channel name(s) to transform when using
#'   channel names. Default is NULL (transforms all channels).
#' @param trans_marker character or NULL Marker name(s) to transform when
#'   using marker names. Default is NULL (transforms all markers).
#' @return A tibble with columns `pop`, `ind` and one column per requested
#'   channel. Rows correspond to cells.
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' dir.create(file.path(tmp, "sample_data", "POP1", "ind_1"),
#'   recursive = TRUE
#' )
#' saveRDS(
#'   c(1, 2, 3),
#'   file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC1.rds")
#' )
#' saveRDS(
#'   c(4, 5, 6),
#'   file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC2.rds")
#' )
#' stimgate_data_get_ex(tmp)
#' stimgate_data_get_ex(tmp, chnl = "BC1")
#' }
#' @export
stimgate_data_get_ex <- function(path_project,
                                 .data = NULL,
                                 pop = NULL,
                                 ind = NULL,
                                 chnl = NULL,
                                 marker = NULL,
                                 bias = FALSE,
                                 exc_min = FALSE,
                                 combn_exc = NULL,
                                 chnl_gate = NULL,
                                 marker_gate = NULL,
                                 gate_type_cyt_pos = "cyt",
                                 gate_type_single_pos = "single",
                                 mult = FALSE,
                                 gate_uns_method = "min",
                                 trans_fn = NULL,
                                 trans_chnl = NULL,
                                 trans_marker = NULL) {
  .assert_string(path_project)
  pop <- pop %|c|% .get_ex_project_pop(path_project)
  if (!is.null(chnl) && !is.null(marker)) {
    stop("Must not specify both marker and chnl")
  }
  .assert_string_vector(pop)
  purrr::map_df(pop, function(pop_curr) {
    ind <- ind %|c|% .get_ex_project_ind(path_project, pop_curr)
    .assert_string_vector(ind)
    purrr::map_df(ind, function(ind_curr) {
      chnl <- if (!is.null(marker)) {
        is_marker <- TRUE
        marker <- as.character(marker)
        stimgate_meta_read_marker_lab(path_project)[marker]
      } else {
        is_marker <- FALSE
        chnl %|c|%
          .get_ex_project_chnl(path_project, pop_curr, ind_curr)
      }
      .assert_string_vector(chnl)
      ex <- .data_get_ex_init(
        .data, pop_curr, chnl, ind_curr, path_project
      )
      ex <- .data_get_ex_exc_min(ex, exc_min)
      ex <- .data_get_ex_cyt_pos(
        ex = ex,
        chnl_gate = chnl_gate,
        marker_gate = marker_gate,
        pop = pop_curr,
        ind = ind_curr,
        combn_exc = combn_exc,
        gate_type_cyt_pos = gate_type_cyt_pos,
        gate_type_single_pos = gate_type_single_pos,
        mult = mult,
        path_project = path_project
      )
      ex <- .data_get_ex_bias(
        ex,
        ind = ind_curr,
        path_project = path_project,
        bias = bias
      )
      ex <- .data_get_ex_renamed(ex, is_marker, path_project)
      trans_chnl_final <- if (is_marker) trans_marker else trans_chnl
      ex <- .data_get_ex_trans(ex, trans_fn, trans_chnl_final)
      .data_get_ex_meta(ex, pop_curr, ind_curr)
    })
  })
}

.data_get_ex_init <- function(.data,
                              pop,
                              chnl,
                              ind,
                              path_project) {
  chnl_cut <- chnl[[1]]
  extra_chnl <- setdiff(chnl, chnl_cut)
  extra_chnl <- if (length(extra_chnl) == 0L) NULL else extra_chnl
  .get_ex(
    .data, pop, chnl_cut, ind, NULL, NULL,
    extra_chnl, path_project, FALSE
  )
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
                                 chnl_gate,
                                 marker_gate,
                                 pop,
                                 ind,
                                 combn_exc = NULL,
                                 gate_type_cyt_pos = "cyt",
                                 gate_type_single_pos = "single",
                                 mult = FALSE,
                                 path_project) {
  if (is.null(chnl_gate) && is.null(marker_gate)) {
    return(ex)
  }
  if (!is.null(chnl_gate) && !is.null(marker_gate)) {
    stop("Must not specify both chnl_gate and marker_gate")
  }
  cn_vec <- colnames(ex)
  chnl_gate <- if (!is.null(marker_gate)) {
    is_marker <- TRUE
    stimgate_meta_read_marker_lab(path_project)[marker_gate]
  } else {
    is_marker <- FALSE
    chnl_gate %||%
      .get_ex_project_chnl(path_project, pop, ind)
  }
  gate_tbl_ind <- .gate_get_gate_tbl_all(NULL, pop, chnl_gate, path_project) |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .data_get_ex_cyt_pos_inc(
    ex, gate_tbl_ind, mult, chnl_gate,
    gate_type_cyt_pos, gate_type_single_pos
  )

  if (nrow(ex) == 0L) {
    message("No stimulation-positive cells.")
    return(.data_get_ex_zero_tbl(cn_vec))
  }

  ex <- .data_get_ex_cyt_pos_exc(
    ex, combn_exc, gate_tbl_ind, chnl,
    gate_type_cyt_pos, gate_type_single_pos
  )

  if (nrow(ex) == 0L) {
    message("No stimulation-positive cells after excluding specified cytokine combinations.") # nolint
    return(.data_get_ex_zero_tbl(cn_vec))
  }

  ex
}

.data_get_ex_zero_tbl <- function(cn) {
  out_df <- matrix(rep(NA_real_, length(cn)), ncol = length(cn))
  colnames(out_df) <- cn
  tibble::as_tibble(out_df)
}

#' @keywords internal
.data_get_ex_cyt_pos_inc <- function(ex,
                                     gate_tbl_ind,
                                     mult,
                                     chnl,
                                     gate_type_cyt_pos,
                                     gate_type_single_pos) {
  inc_vec <- rep(FALSE, nrow(ex))

  if (!mult) {
    inc_vec <- .get_pos_ind( # nolint
      ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
      gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )
  } else {
    inc_vec <- .get_pos_ind_mult( # nolint
      ex = ex, gate_tbl = gate_tbl_ind, chnl = chnl, chnl_alt = NULL,
      gate_type_cyt_pos = gate_type_cyt_pos
    )
  }
  ex[inc_vec, , drop = FALSE]
}

#' @keywords internal
.data_get_ex_cyt_pos_exc <- function(ex,
                                     combn_exc,
                                     gate_tbl_ind,
                                     chnl_gate,
                                     gate_type_cyt_pos,
                                     gate_type_single_pos) {
  if (is.null(combn_exc)) {
    return(ex)
  }
  for (chnl_pos in combn_exc) {
    if (nrow(ex) == 0) break
    exc_vec <- .get_pos_ind_cyt_combn( # nolint
      ex = ex, gate_tbl = gate_tbl_ind,
      chnl_pos = chnl_pos, chnl_neg = setdiff(chnl_gate, chnl_pos),
      chnl_alt = NULL, gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )
    ex <- ex[!exc_vec, , drop = FALSE]
  }
  ex
}

#' @keywords internal
.data_get_ex_renamed <- function(ex,
                                 is_marker,
                                 path_project) {
  # if user specified markers, then give them back a table
  # with column names as markers
  if (!is_marker) {
    return(ex)
  }
  colnames(ex) <- stimgate_meta_read_chnl_lab(path_project)[
    colnames(ex)
  ]
  ex
}

#' @keywords internal
.data_get_ex_trans <- function(ex,
                               trans_fn,
                               trans_chnl) {
  # transform
  if (is.null(trans_fn)) {
    return(ex)
  }
  if (is.null(trans_chnl)) {
    ex <- trans_fn(ex)
  } else {
    for (nm in trans_chnl) {
      ex[, nm] <- trans_fn(ex[, nm])
    }
  }
  ex
}

#' @keywords internal
.data_get_ex_meta <- function(ex, pop, ind) {
  meta_df <- tibble::tibble(
    pop = pop,
    ind = ind
  )
  tibble::as_tibble(cbind(meta_df, ex))
}
