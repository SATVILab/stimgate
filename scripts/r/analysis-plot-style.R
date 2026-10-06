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

# Construct variants from one completed plot. Removing marked MC layers leaves
# all points, axes, non-MC intervals and ratio coordinates unchanged.
.analysis_mcse_plot_variants <- function(plot, mcse_mode = NULL) {
  if (is.null(mcse_mode)) return(list(original = plot))
  mode <- .analysis_mcse_mode(mcse_mode)
  variants <- mode
  stats::setNames(lapply(variants, function(version) {
    out <- plot
    if (version == "off") {
      out$layers <- lapply(plot$layers, function(layer) {
        if (!isTRUE(attr(layer, "analysis_mcse"))) return(layer)
        # Keep interval scale training (including facet ranges and bar widths),
        # but draw nothing. Copy the layer instead of mutating the original.
        blank <- ggplot2::ggproto(NULL, ggplot2::GeomBlank,
          setup_data = layer$geom$setup_data, setup_params = layer$geom$setup_params,
          required_aes = layer$geom$required_aes, default_aes = layer$geom$default_aes)
        ggplot2::ggproto(NULL, layer, geom = blank, show.legend = FALSE)
      })
    }
    out
  }), variants)
}

.analysis_mcse_fig_path <- function(path, version) {
  if (version == "original") path else file.path(dirname(path), paste0("mcse_", version), basename(path))
}

# Save at the standard width so text matches across figures. Heights above a
# page are capped unless `allow_tall` is TRUE for figures that cannot fit.
.analysis_save_fig <- function(
    plot,
    path,
    height = 12,
    width = .analysis_fig_width,
    allow_tall = FALSE,
    mcse_mode = NULL) {
  if (!isTRUE(allow_tall)) {
    height <- min(height, .analysis_fig_max_height)
  }
  variants <- .analysis_mcse_plot_variants(plot, mcse_mode)
  paths <- vapply(names(variants), function(version) .analysis_mcse_fig_path(path, version), character(1))
  for (version in names(variants)) {
    dir.create(dirname(paths[[version]]), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = paths[[version]], plot = variants[[version]], width = width, height = height,
      units = "cm", limitsize = FALSE
    )
  }
  invisible(unname(paths))
}

# Markdown heading for one level of a plot loop. Use in chunks with
# `results: asis`; `level` counts `#`s, so nested loops go one level deeper.
.analysis_heading <- function(text, level) {
  cat("\n\n", strrep("#", level), " ", text, "\n\n", sep = "")
  invisible(NULL)
}

# Print a plot inside a `results: asis` loop, separated from the next heading.
.analysis_print_fig <- function(plot, mcse_mode = NULL) {
  variants <- .analysis_mcse_plot_variants(plot, mcse_mode)
  for (version in names(variants)) {
    if (version != "original") {
      cat("\n\n**Monte Carlo error bars: ", version, "**\n\n", sep = "")
    }
    print(variants[[version]])
    cat("\n\n")
  }
  invisible(plot)
}

# Print before opening a save device so asis headings keep their own figure.
.analysis_print_save_fig <- function(plot, path, ..., mcse_mode = NULL) {
  .analysis_print_fig(plot, mcse_mode = mcse_mode)
  .analysis_save_fig(plot, path, ..., mcse_mode = mcse_mode)
}

# Colour roles, kept distinct so a colour means one thing across the analyses:
# - methods (QMDs 7-10): Okabe-Ito blue StimGate, bluish-green Tailgate,
#   orange F-beta, for colour-blind accessibility;
# - over/under direction: ColorBrewer BrBG teal and brown (see
#   `.simBandwidthSignedErrorColours`);
# - error statistic (median, upper percentile, maximum): blues and lavender;
# - bandwidth: a sequential purple ramp (`make_bw_colour_values()`).
.analysis_method_colours <- c(
  stimgate = "#0072B2", tailgate = "#009E73", fbeta = "#E69F00"
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
