# Shared figure style for the analysis QMDs: transformation labels and order,
# number formats, theme, figure sizes and section headings for plot loops.
# Source after analysis-runtime.R.

# Transformations always appear in this order, with these labels.
.analysis_trans_labels <- c(gaussian = "Gaussian", skew = "Skew", gamma = "Gamma")

# Factor of transformation labels in the standard order. Unknown values keep
# their own label after the standard ones.
.analysis_trans_factor <- function(x) {
  x <- as.character(x)
  lab <- ifelse(
    x %in% names(.analysis_trans_labels),
    unname(.analysis_trans_labels[x]),
    x
  )
  extra <- setdiff(unique(lab), .analysis_trans_labels)
  factor(lab, levels = c(.analysis_trans_labels, extra))
}

# Standard-order transformation names (lower case) present in `x`.
.analysis_trans_order <- function(x) {
  x <- unique(as.character(x))
  c(
    intersect(names(.analysis_trans_labels), x),
    setdiff(x, names(.analysis_trans_labels))
  )
}

# Plain decimal with only the digits needed, never scientific notation:
# 0.0025 -> "0.0025", 1e5 -> "100,000".
.analysis_label_number <- function(x, digits = 3) {
  out <- formatC(signif(x, digits), format = "fg", digits = digits, big.mark = ",")
  out <- trimws(out)
  # formatC pads with zeros after the decimal point; integers keep their zeros.
  decimal <- grepl(".", out, fixed = TRUE)
  out[decimal] <- sub("\\.?0+$", "", out[decimal])
  out[is.na(x)] <- NA_character_
  out
}

# Proportion as a percentage with only the digits needed:
# 0.0001 -> "0.01%", 0.002 -> "0.2%", 0.2 -> "20%".
.analysis_label_percent <- function(x, digits = 3) {
  out <- paste0(.analysis_label_number(100 * x, digits), "%")
  out[is.na(x)] <- NA_character_
  out
}

# Facet labeller for response probabilities and other proportions.
.analysis_labeller_percent <- function(prefix = "") {
  function(x) paste0(prefix, .analysis_label_percent(as.numeric(x)))
}

# Base text size (pt) and figure width (cm, an A4 text block). Saving every
# figure at this width keeps text the same size across figures.
.analysis_fig_font_size <- 10
.analysis_fig_width <- 16
.analysis_fig_max_height <- 22

# Shared theme: white background, light grid, facet strips with a black
# border and no fill, legend underneath.
.analysis_theme <- function(base_size = .analysis_fig_font_size, grid = "xy") {
  cowplot::theme_cowplot(font_size = base_size) +
    cowplot::background_grid(major = grid) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center"
    )
}

# Save at the standard width so text matches across figures. Heights above a
# page are capped unless `allow_tall` is TRUE for figures that cannot fit.
.analysis_save_fig <- function(
    plot,
    path,
    height = 12,
    width = .analysis_fig_width,
    allow_tall = FALSE) {
  if (!isTRUE(allow_tall)) {
    height <- min(height, .analysis_fig_max_height)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path, plot = plot, width = width, height = height,
    units = "cm", limitsize = FALSE
  )
  invisible(path)
}

# Markdown heading for one level of a plot loop. Use in chunks with
# `results: asis`; `level` counts `#`s, so nested loops go one level deeper.
.analysis_heading <- function(text, level) {
  cat("\n\n", strrep("#", level), " ", text, "\n\n", sep = "")
  invisible(NULL)
}

# Print a plot inside a `results: asis` loop, separated from the next heading.
.analysis_print_fig <- function(plot) {
  print(plot)
  cat("\n\n")
  invisible(plot)
}

# Colour roles, kept distinct so a colour means one thing across the analyses:
# - methods (QMDs 7-10): raspberry StimGate, slate-blue Tailgate, saffron
#   F-beta; distinct in hue and lightness, so colour-blind and greyscale safe;
# - over/under direction: ColorBrewer BrBG teal and brown (see
#   `.simBandwidthSignedErrorColours`);
# - error statistic (median, upper percentile, maximum): blues and lavender;
# - bandwidth: a sequential purple ramp (`make_bw_colour_values()`).
.analysis_method_colours <- c(
  stimgate = "#C0395A", tailgate = "#3D5A80", fbeta = "#E9A23B"
)
.analysis_method_labels <- c(
  stimgate = "StimGate", tailgate = "Tailgate", fbeta = "F-beta"
)
.analysis_stat_colours <- c(median = "#0072B2", upper = "#56B4E9", max = "#8C8DBA")

# Colour scale for methods; `aesthetics = "fill"` for fills.
.analysis_scale_method <- function(aesthetics = "colour", ...) {
  ggplot2::scale_colour_manual(
    values = .analysis_method_colours,
    labels = .analysis_method_labels,
    aesthetics = aesthetics,
    ...
  )
}
