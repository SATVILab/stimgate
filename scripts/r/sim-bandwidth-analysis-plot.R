format_bw_lab <- function(x, digits = 4L) {
  x_num <- suppressWarnings(as.numeric(x))
  if (length(x_num) == 0L) {
    return(character(0))
  }

  lab <- vapply(x_num, function(x_curr) {
    if (length(x_curr) == 0L || is.na(x_curr) || !is.finite(x_curr)) {
      return("NA")
    }
    format(signif(x_curr, digits), scientific = FALSE, trim = TRUE)
  }, FUN.VALUE = character(1))

  sub("\\.?0+$", "", lab)
}

make_bw_colour_values <- function(bw_vec, base_col_vec = NULL) {
  bw_lab <- format_bw_lab(bw_vec)
  bw_unique <- unique(bw_lab)

  if (is.null(base_col_vec) || length(base_col_vec) == 0L) {
    base_col_vec <- grDevices::hcl.colors(length(bw_unique), palette = "Dark 2")
  }

  if (length(base_col_vec) < length(bw_unique)) {
    base_col_vec <- rep(base_col_vec, length.out = length(bw_unique))
  }

  names(base_col_vec) <- bw_unique
  out <- base_col_vec[match(bw_lab, names(base_col_vec))]
  stats::setNames(out, bw_lab)
}

make_bw_linetype_scale <- function(bw_vec, lty_vec = NULL) {
  bw_lab <- format_bw_lab(bw_vec)
  bw_unique <- unique(bw_lab)

  if (is.null(lty_vec) || length(lty_vec) == 0L) {
    lty_vec <- c(
      "solid",
      "dashed",
      "dotted",
      "dotdash",
      "longdash",
      "twodash",
      "1F",
      "F1",
      "4C88C488"
    )
  }

  if (length(lty_vec) < length(bw_unique)) {
    lty_vec <- rep(lty_vec, length.out = length(bw_unique))
  }

  names(lty_vec) <- bw_unique
  out <- lty_vec[match(bw_lab, names(lty_vec))]
  stats::setNames(out, bw_lab)
}

add_bw_labs <- function(.data, bw_cols = c("bw", "bw_core", "bw_extra")) {
  if (!is.data.frame(.data)) {
    return(.data)
  }

  for (bw_col in bw_cols) {
    if (bw_col %in% names(.data)) {
      .data[[paste0(bw_col, "_lab")]] <- format_bw_lab(.data[[bw_col]])
    }
  }

  .data
}
