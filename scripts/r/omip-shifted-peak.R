# Analyses 14b and 15b: OMIP-111 and OMIP-016 with StimGate's shifted-peak
# rule (`stimControl(locShiftedPeakRef = TRUE)`), compared with the saved
# default results of Analyses 14 and 15. Source after analysis-runtime.R and
# analysis-plot-style.R.

# Analysis 14's settings with only the shifted-peak rule switched on.
.omipShiftedPeakSettings <- function(settings) {
  settings$semantics <- paste0(settings$semantics, "-shifted-peak")
  settings$control$locShiftedPeakRef <- TRUE
  settings
}

# Whether the rule applied to each stimulated tube's own local-FDR gate, read
# from the final `loc_minClust` gate table of each StimGate project. `labels`
# has one row of identifying columns per project.
.omipShiftedPeakFired <- function(paths, labels, gateName = "loc_minClust") {
  if (length(paths) != nrow(labels)) stop("One label row is needed per project.")
  dplyr::bind_rows(lapply(seq_along(paths), function(i) {
    gates <- stimgate::getStimGates(paths[[i]])
    gates <- gates[gates$gateName == gateName, , drop = FALSE]
    if (!"locShiftedPeakRef" %in% names(gates)) {
      stop("StimGate project lacks the shifted-peak flag: ", paths[[i]])
    }
    data.frame(
      labels[rep(i, nrow(gates)), , drop = FALSE],
      ind = as.character(gates$ind), chnl = as.character(gates$chnl),
      marker = as.character(gates$marker),
      shiftedPeakRef = gates$locShiftedPeakRef %in% TRUE,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }))
}

# One row per OMIP-111 outcome with the default (Analysis 14) and
# shifted-peak results side by side. Tube indices are positions within each
# strain's sample map, as in .omip111Run().
.omip111ShiftedPeakPaired <- function(default, shifted, fired, samples) {
  key <- c("strain", "mouse", "population", "marker", "method")
  if (nrow(default) != nrow(shifted) || anyDuplicated(default[key]) ||
    anyDuplicated(shifted[key])) {
    stop("Default and shifted-peak OMIP-111 cohorts differ.")
  }
  tubes <- dplyr::bind_rows(lapply(split(samples, samples$strain), function(x) {
    data.frame(strain = x$strain, ind = as.character(seq_len(nrow(x))), sampleStim = x$sample)
  }))
  fired <- dplyr::left_join(fired, tubes, by = c("strain", "ind"))
  defaultCols <- default[c(key, "threshold", "propBs", "errorPp")]
  names(defaultCols)[-seq_along(key)] <- paste0(names(defaultCols)[-seq_along(key)], "Default")
  out <- shifted |>
    dplyr::left_join(defaultCols, by = key) |>
    dplyr::left_join(
      fired[c("sampleStim", "population", "marker", "shiftedPeakRef")],
      by = c("sampleStim", "population", "marker")
    )
  stimgate <- out$method == "StimGate"
  if (anyNA(out$shiftedPeakRef[stimgate])) stop("Missing shifted-peak flags.")
  out$shiftedPeakRef[!stimgate] <- FALSE
  out$unit <- out$mouse
  out
}

# Comparators do not depend on StimGate, so their results must not change.
.omipShiftedPeakComparatorsUnchanged <- function(paired) {
  other <- paired[paired$method != "StimGate", , drop = FALSE]
  same <- function(a, b) all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
  same(other$threshold, other$thresholdDefault) && same(other$propBs, other$propBsDefault)
}

# Mean StimGate errors (percentage points) with and without the rule.
.omipShiftedPeakErrorSummary <- function(paired, groups) {
  paired[paired$method %in% c("StimGate", "stimgate"), , drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
    dplyr::summarise(
      n = dplyr::n(), rule_applied = sum(.data$shiftedPeakRef),
      mean_error_default_pp = mean(.data$errorPpDefault),
      mean_error_shifted_pp = mean(.data$errorPp),
      mean_abs_error_default_pp = mean(abs(.data$errorPpDefault)),
      mean_abs_error_shifted_pp = mean(abs(.data$errorPp)),
      changed = sum(.data$threshold != .data$thresholdDefault, na.rm = TRUE),
      .groups = "drop"
    )
}

# StimGate signed errors per marker, default against shifted peak, one line
# per unit (mouse or tube); filled points where the rule applied.
.omipShiftedPeakErrorPlot <- function(paired, facet) {
  sg <- paired[paired$method %in% c("StimGate", "stimgate"), , drop = FALSE]
  long <- dplyr::bind_rows(
    data.frame(sg[c(facet, "marker", "unit", "shiftedPeakRef")],
      setting = "default", errorPp = sg$errorPpDefault
    ),
    data.frame(sg[c(facet, "marker", "unit", "shiftedPeakRef")],
      setting = "shifted peak", errorPp = sg$errorPp
    )
  )
  long$setting <- factor(long$setting, levels = c("default", "shifted peak"))
  long$applied <- ifelse(long$shiftedPeakRef, "rule applied", "rule not applied")
  ggplot2::ggplot(long, ggplot2::aes(.data$setting, .data$errorPp)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_line(ggplot2::aes(group = .data$unit), colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$applied, colour = .data$setting),
      size = 2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(c(facet, "marker"))),
      scales = "free_y", ncol = 5
    ) +
    .analysis_y_floor(c(-10, 10)) +
    ggplot2::scale_shape_manual(values = c(`rule applied` = 16, `rule not applied` = 1)) +
    ggplot2::scale_colour_manual(values = c(default = "#0072B2", `shifted peak` = "#D55E00")) +
    ggplot2::labs(
      x = NULL, y = "StimGate − reference (percentage points)",
      shape = NULL, colour = "StimGate setting"
    ) +
    .analysis_theme()
}
