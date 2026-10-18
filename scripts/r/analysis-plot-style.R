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

# Actual parent-cell counts, compact labels; never alter the simulation grid.
.analysis_label_cells <- function(x) {
  x <- as.numeric(as.character(x))
  ifelse(x >= 1000, paste0(.analysis_label_number(x / 1000), "k"),
    .analysis_label_number(x)
  )
}

.analysis_cell_x <- function(x, discrete = length(unique(x[!is.na(x)])) == 2L) {
  if (discrete) factor(x, levels = sort(unique(x))) else x
}

.analysis_scale_cells <- function(x, discrete = is.factor(x), log = TRUE) {
  if (discrete) {
    return(ggplot2::scale_x_discrete(labels = .analysis_label_cells))
  }
  if (log) ggplot2::scale_x_log10(breaks = sort(unique(x)), labels = .analysis_label_cells) else
    ggplot2::scale_x_continuous(breaks = sort(unique(x)), labels = .analysis_label_cells)
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

# Shared theme: white background, light grid, facet strips with a soft grey
# border and no fill, legend underneath.
# Thesis figures (Manuscript22RodoThesis, stanford-thesis Typst template: A4,
# 1.5 in left and 1 in other margins, so a 146.5 x 246.2 mm text area and
# 12 pt body text) are drawn at their printed size and inserted at full width.
# Save one to the projr `output-thesis` folder only when the thesis uses it.
.analysis_thesis_width <- 146.5

# Facet strips: white background with a soft grey border.
.analysis_strip_theme <- function() {
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white", colour = "grey65"),
    text = ggplot2::element_text(colour = "#333333"),
    axis.text = ggplot2::element_text(colour = "#333333"),
    axis.ticks = ggplot2::element_line(colour = "#333333")
  )
}

# Text sizes (pt) for a figure printed at the thesis text width, with white
# facet strips.
.analysis_thesis_theme <- function() {
  .analysis_strip_theme() + ggplot2::theme(
    text = ggplot2::element_text(size = 8, colour = "#333333"),
    axis.text = ggplot2::element_text(size = 7.5),
    axis.title = ggplot2::element_text(size = 8.5),
    strip.text = ggplot2::element_text(size = 7.5),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(2, 4, 2, 2)
  )
}

# Save as vector PDF and 300 dpi PNG at the printed size (mm).
.analysis_save_thesis_fig <- function(plot, path_no_ext, height,
                                      width = .analysis_thesis_width) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    ggplot2::ggsave(
      paste0(path_no_ext, ".", ext), plot, width = width, height = height,
      units = "mm", dpi = 300,
      device = if (ext == "pdf") grDevices::cairo_pdf else "png"
    )
  }
  invisible(paste0(path_no_ext, c(".pdf", ".png")))
}

.analysis_theme <- function(base_size = .analysis_fig_font_size, grid = "xy") {
  cowplot::theme_cowplot(font_size = base_size) +
    cowplot::background_grid(major = grid) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = NA, colour = "grey65"),
      text = ggplot2::element_text(colour = "#333333"),
      axis.text = ggplot2::element_text(colour = "#333333"),
      axis.line = ggplot2::element_line(colour = "#333333"),
      axis.ticks = ggplot2::element_line(colour = "#333333"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center"
    )
}

# Train every facet in data space without setting limits or clipping intervals.
# Omitting facet variables broadcasts the two invisible endpoints to every panel,
# including free-y panels, without training their x scales or aesthetic legends.
.analysis_y_floor <- function(range = c(0, 0.1)) {
  ggplot2::geom_blank(
    data = data.frame(.axis_floor = range),
    mapping = ggplot2::aes(y = .axis_floor), inherit.aes = FALSE
  )
}

# Reserve room for strips, axes and legends, plus 5 cm per facet row.
.analysis_facet_height <- function(plot) {
  layout <- ggplot2::ggplot_build(plot)$layout$layout
  5 + 5 * max(layout$ROW)
}

# Construct variants from one completed plot. Removing marked MC layers leaves
# all points, axes, non-MC intervals and ratio coordinates unchanged.
.analysis_mcse_plot_variants <- function(plot, mcse_mode = NULL) {
  if (is.null(mcse_mode)) {
    return(list(original = plot))
  }
  mode <- .analysis_mcse_mode(mcse_mode)
  variants <- mode
  stats::setNames(lapply(variants, function(version) {
    out <- plot
    if (version == "off") {
      out$layers <- lapply(plot$layers, function(layer) {
        if (!isTRUE(attr(layer, "analysis_mcse"))) {
          return(layer)
        }
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
.analysis_print_fig <- function(plot, mcse_mode = NULL, width = NULL, height = NULL) {
  variants <- .analysis_mcse_plot_variants(plot, mcse_mode)
  for (version in names(variants)) {
    if (version != "original") {
      cat("\n\n**Monte Carlo error bars: ", version, "**\n\n", sep = "")
    }
    if (!is.null(height) && isTRUE(getOption("knitr.in.progress"))) {
      # A loop can contain different facet counts. A PNG device per plot makes
      # HTML use the same dimensions as the saved figure, independent of the
      # chunk's fixed fig.width/fig.height. Embed it as a data URI: Quarto only
      # embeds figures knitr recorded, and drops other files in its figure
      # folder, while a knit_print() of include_graphics() in a
      # `results: asis` chunk prints only the path.
      path <- tempfile("analysis-figure-", fileext = ".png")
      on.exit(unlink(path), add = TRUE)
      (function() {
        grDevices::png(path, width = width, height = height, units = "cm", res = 150)
        on.exit(grDevices::dev.off())
        print(variants[[version]])
      })()
      cat("![](", knitr::image_uri(path), "){width=100%}", sep = "")
    } else {
      print(variants[[version]])
    }
    cat("\n\n")
  }
  invisible(plot)
}

# Print before opening a save device so asis headings keep their own figure.
.analysis_print_save_fig <- function(plot, path, ..., mcse_mode = NULL,
                                     fit_panels = FALSE) {
  sizes <- if (isTRUE(fit_panels)) {
    list(
    width = .analysis_fig_width, height = .analysis_facet_height(plot), allow_tall = TRUE
  ) } else {
    list()
  }
  if (isTRUE(fit_panels)) {
    .analysis_print_fig(plot, mcse_mode = mcse_mode,
      width = sizes$width, height = sizes$height)
  } else {
    .analysis_print_fig(plot, mcse_mode = mcse_mode)
  }
  if (exists(".simBandwidthDisplayNote", mode = "function")) .simBandwidthDisplayNote(.analysis_mcse_plot_variants(plot, mcse_mode)[[1L]], path)
  args <- list(...)
  args[names(sizes)] <- sizes
  do.call(.analysis_save_fig, c(list(plot = plot, path = path, mcse_mode = mcse_mode), args))
}

# Colour roles, kept distinct so a colour means one thing across the analyses:
# - methods (QMDs 7-10): Okabe-Ito orange StimGate, bluish-green Tailgate,
#   blue F-beta, reddish-purple Tailgate at default settings, for colour-blind
#   accessibility;
# - over/under direction: ColorBrewer BrBG teal and brown (see
#   `.simBandwidthSignedErrorColours`);
# - error statistic (median, upper percentile, maximum): blues and lavender;
# - bandwidth: a sequential purple ramp (`make_bw_colour_values()`).
.analysis_method_colours <- c(
  stimgate = "#E69F00", fbeta = "#0072B2", tailgate = "#009E73",
  tailgate_default = "#CC79A7"
)
.analysis_method_labels <- c(
  stimgate = "StimGate", fbeta = "F-beta", tailgate = "Tailgate",
  tailgate_default = "Tailgate (default settings)"
)
.analysis_method_shapes <- c(stimgate = 16, fbeta = 15, tailgate = 17, tailgate_default = 2)
.analysis_method_linetypes <- c(
  stimgate = "solid", fbeta = "42", tailgate = "22", tailgate_default = "13"
)
.analysis_stat_colours <- c(median = "#0072B2", upper = "#56B4E9", max = "#8C8DBA")

# 10th-90th percentile bands behind median lines (thesis figures of Analyses
# 2a, 7 and 8). The band is the series colour at `.analysis_band_alpha` (0.18):
# at 0.10-0.12 the orange StimGate band almost vanishes on white, while at 0.25
# three overlapping bands turn a muddy grey-green that hides which method is
# where; 0.18 keeps each hue identifiable with the medians drawn on top. Band
# boundary lines are left out ("none"), dotted ("dotted") or thin solid lines
# ("thin") whose darkness lies between the full-colour median line and the
# pale fill.
.analysis_band_alpha <- 0.18
.analysis_band_median_linewidth <- 0.6
.analysis_band_line_styles <- list(
  dotted = list(linetype = "11", linewidth = 0.45, alpha = 0.9),
  thin = list(linetype = "solid", linewidth = 0.25, alpha = 0.55)
)

# Ribbon and optional boundary lines for `lo`-`hi` bands. The ribbon fill is
# mapped to `fill` (add the matching fill scale); boundary lines take the
# plot's colour and group mappings. A ribbon needs two x values, so series
# with a single x value (identified by the `series` columns, which must
# include the facet variables) get a pale vertical bar of the same alpha,
# dodged by `single_dodge` (x units; dodge their median points the same way).
# Rows of `d` whose series (`series` columns) has a single x value.
.analysis_band_single <- function(d, x, series) {
  key <- interaction(d[series], drop = TRUE)
  n_x <- tapply(d[[x]], key, function(v) length(unique(v)))
  as.vector(n_x[as.character(key)] == 1L)
}

.analysis_band_layers <- function(lines = c("none", "dotted", "thin"), fill,
                                  lo = "lo", hi = "hi", x = "x", series = NULL,
                                  single_dodge = 0.18,
                                  alpha = .analysis_band_alpha) {
  lines <- match.arg(lines)
  out <- list(ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]], fill = .data[[fill]]),
    alpha = alpha, colour = NA, show.legend = FALSE
  ))
  if (length(series)) {
    out <- c(out, list(ggplot2::geom_linerange(
      data = function(d) d[.analysis_band_single(d, x, series), , drop = FALSE],
      ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]], colour = .data[[fill]]),
      alpha = 2 * alpha, linewidth = 1.6, show.legend = FALSE,
      position = ggplot2::position_dodge(width = single_dodge)
    )))
  }
  if (identical(lines, "none")) {
    return(out)
  }
  style <- .analysis_band_line_styles[[lines]]
  c(out, lapply(c(lo, hi), function(y) {
    ggplot2::geom_line(ggplot2::aes(y = .data[[y]]), linetype = style$linetype,
      linewidth = style$linewidth, alpha = style$alpha)
  }))
}

.analysis_method_order <- names(.analysis_method_labels)
.analysis_method_factor <- function(x) {
  x <- as.character(x)
  factor(x, levels = c(.analysis_method_order, setdiff(unique(x), .analysis_method_order)))
}

.analysis_method_data <- function(data) {
  if ("method" %in% names(data)) data$method <- .analysis_method_factor(data$method)
  data
}

# Scale limits reorder only the observed methods, including after filtering.
.analysis_method_limits <- function(x) {
  c(intersect(.analysis_method_order, x), setdiff(x, .analysis_method_order))
}

# Shared method encodings; matching names and labels merge their legends.
.analysis_scale_method <- function(aesthetics = "colour", ...) {
  colour_aesthetics <- intersect(aesthetics, c("colour", "color", "fill"))
  scales <- list()
  if (length(colour_aesthetics)) {
    scales <- c(scales, list(ggplot2::scale_colour_manual(
    values = .analysis_method_colours, labels = .analysis_method_labels,
    name = "Method", aesthetics = colour_aesthetics, limits = .analysis_method_limits, ...
  )))
  }
  if ("shape" %in% aesthetics) {
    scales <- c(scales, list(ggplot2::scale_shape_manual(
    values = .analysis_method_shapes, labels = .analysis_method_labels, name = "Method", limits = .analysis_method_limits, ...
  )))
  }
  if ("linetype" %in% aesthetics) {
    scales <- c(scales, list(ggplot2::scale_linetype_manual(
    values = .analysis_method_linetypes, labels = .analysis_method_labels, name = "Method", limits = .analysis_method_limits, ...
  )))
  }
  if (length(scales) == 1L) scales[[1]] else scales
}
