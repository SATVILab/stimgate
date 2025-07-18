#' Detect if any pattern matches a string
#'
#' @description Helper function to check if any of multiple patterns match a given string
#'
#' @param string character. String to search in
#' @param pattern character vector. Patterns to search for
#'
#' @return logical. TRUE if any pattern matches, FALSE otherwise
#'
#' @keywords internal
str_detect_any <- function(string, pattern) {
  vapply(
    pattern,
    function(pattern_curr) grepl(pattern_curr, string),
    logical(1)
  ) |>
    any()
}

#' Get list of expression matrices for specified samples
#'
#' @description Internal helper function to extract expression matrices from GatingSet
#' for specified samples and populations
#'
#' @param .data GatingSet. Source data
#' @param ind_batch numeric vector. Indices of samples in batch
#' @param batch character. Batch identifier  
#' @param pop character. Population to extract
#' @param chnl_cut character vector. Channels to include
#' @param extra_chnl character vector. Additional channels to include (optional)
#'
#' @return list of expression matrices for each sample
#'
#' @keywords internal
.get_ex_list <- function(.data,
                         ind_batch,
                         batch,
                         pop,
                         chnl_cut,
                         extra_chnl = NULL) {
  # get expression .data for each batch
  lapply(ind_batch, function(i) {
    .get_ex(
      .data = .data[[i]],
      pop = pop,
      chnl_cut,
      extra_chnl = extra_chnl,
      ind = i,
      # specify corresponding unstim
      ind_uns = ind_batch[length(ind_batch)],
      batch = batch
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
                    extra_chnl = NULL) {
  # get expression information as a tibble
  # get .data
  fr <- flowWorkspace::gh_pop_get_data(.data, y = pop)
  ex <- flowCore::exprs(fr)
  ex <- ex[, c(chnl_cut, extra_chnl), drop = FALSE] |> tibble::as_tibble()

  attr(ex, "ind") <- ind |> as.character()
  attr(ex, "ind_uns") <- ind_uns |> as.character()
  attr(ex, "is_uns") <- ind == ind_uns
  attr(ex, "chnl_cut") <- chnl_cut
  attr(ex, "batch") <- batch

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
