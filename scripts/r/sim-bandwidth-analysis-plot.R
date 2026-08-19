format_bw_lab <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(x))
  lab <- format(signif(x, digits), scientific = FALSE, trim = TRUE)
  lab <- sub("\\.?0+$", "", lab)
  lab[is.na(x)] <- NA_character_
  lab
}

format_bw_file <- function(x) {
  x <- format_bw_lab(x)
  x <- gsub("\\.", "p", x)
  gsub("[^A-Za-z0-9]+", "_", x)
}

safe_file_lab <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

.file_safe <- safe_file_lab

make_bw_colour_values <- function(bw_vec, base_col_vec = NULL) {
  if (is.null(base_col_vec) || length(base_col_vec) == 0L) {
    base_col_vec <- c(
      "#e0ecf4",
      "#bfd3e6",
      "#9ebcda",
      "#8c96c6",
      "#8c6bb1",
      "#88419d",
      "#810f7c"
    )
  }

  bw_num <- sort(unique(as.numeric(bw_vec)))
  bw_lab <- format_bw_lab(bw_num)
  n_bw <- length(bw_num)

  if (n_bw <= length(base_col_vec)) {
    col_vec <- base_col_vec[seq(length(base_col_vec) - n_bw + 1L, length(base_col_vec))]
  } else {
    col_vec <- grDevices::colorRampPalette(base_col_vec)(n_bw)
  }

  stats::setNames(col_vec, bw_lab)
}

make_bw_linetype_scale <- function(bw_vec) {
  bw_num <- sort(unique(as.numeric(bw_vec)))
  bw_lab <- format_bw_lab(bw_num)
  n_bw <- length(bw_num)

  if (n_bw == 1L) {
    lty <- "solid"
  } else {
    dash_len <- round(seq(8, 2, length.out = n_bw - 1L))
    gap_len <- 4L

    lty <- c(
      "solid",
      paste0(
        as.character(as.hexmode(dash_len)),
        as.character(as.hexmode(gap_len))
      )
    )
  }

  list(
    levels = bw_lab,
    values = stats::setNames(lty, bw_lab)
  )
}

add_bw_labs <- function(.data) {
  if (!is.data.frame(.data)) {
    return(.data)
  }

  .data |>
    dplyr::mutate(
      bw_core_lab = format_bw_lab(.data$bw_core),
      bw_extra_lab = format_bw_lab(.data$bw_extra)
    )
}
