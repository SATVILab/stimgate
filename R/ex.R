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
                         cut,
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