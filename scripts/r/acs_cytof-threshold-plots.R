# Analysis 9b: distributions of final base gates, with real stimulated cells.
# Thresholds stay on the cached asinh scale; conditional gates are separate.
.acsCytofThresholdTable <- function(gates, comparators, sampleMap) {
  channels <- .acsCytofChannelMap()
  stimMap <- sampleMap[sampleMap$stim != "uns", c("ind", "SampleID", "stim")]
  stimMap$ind <- as.character(stimMap$ind)
  tables <- c(list(stimgate = gates |>
    dplyr::filter(.data$gateName == "loc_minClust") |>
    dplyr::mutate(
      threshold = .data$gate,
      thresholdOrigin = .data$locSource,
      thresholdFallbackUsed = !(.data$locGenerated %in% TRUE)
    )), comparators)
  dplyr::bind_rows(lapply(names(tables), function(method) {
    d <- tables[[method]]
    d$ind <- as.character(d$ind)
    d <- d[d$ind %in% stimMap$ind & d$chnl %in% names(channels), ]
    expected <- expand.grid(ind = stimMap$ind, chnl = names(channels))
    if (anyDuplicated(d[c("ind", "chnl")]) ||
      nrow(dplyr::anti_join(expected, d, by = c("ind", "chnl")))) {
      stop("Missing or duplicate ACS thresholds for ", method, ".")
    }
    d |>
      dplyr::select(dplyr::all_of(c(
        "ind", "chnl", "threshold", "thresholdOrigin", "thresholdFallbackUsed"
      )), dplyr::any_of(c(
        "locGenerated", "locGeneratedDirect", "locSource", "locReason"
      ))) |>
      dplyr::left_join(stimMap, by = "ind") |>
      dplyr::mutate(method = method, cyt = unname(channels[.data$chnl]))
  }))
}

.acsCytofReadThresholds <- function(
    paths, sampleMap, locThresholdMethod, tailgateSelection) {
  manifest <- .acsCytofReadStimGateManifest(paths$stimgate, locThresholdMethod)
  objects <- stats::setNames(lapply(c("fbeta", "tailgate"), function(method) {
    .acsCytofReadComparatorCache(paths[[method]], method,
      nSample = nrow(sampleMap), tailgateSelection = tailgateSelection
    )
  }), c("fbeta", "tailgate"))
  .acsCytofValidateManifests(c(list(manifest), lapply(objects, `[[`, "manifest")))
  if (!identical(manifest$context$preprocessing$sampleMap, sampleMap)) {
    stop("ACS thresholds do not match the loaded GatingSet sample map.")
  }
  .acsCytofThresholdTable(
    stimgate::getStimGates(paths$stimgate),
    lapply(objects, `[[`, "thresholds"), sampleMap
  )
}

# Sample rows once per tube, retaining zeros and all six markers together.
# Equal caps keep memory bounded and ensure every nonempty tube contributes.
.acsCytofThresholdCells <- function(gs, sampleMap, nCell = 2000L, seed = 917L) {
  channels <- .acsCytofChannelMap()
  rows <- sampleMap[sampleMap$stim != "uns", ]
  withr::with_seed(seed, dplyr::bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    gh <- gs[[as.integer(rows$ind[[i]])]]
    ex <- flowCore::exprs(flowWorkspace::gh_pop_get_data(gh, "root"))
    ex <- ex[, names(channels), drop = FALSE]
    index <- sample.int(nrow(ex), min(nCell, nrow(ex)))
    tibble::as_tibble(ex[index, , drop = FALSE]) |>
      tidyr::pivot_longer(
        dplyr::everything(), names_to = "chnl", values_to = "expression"
      ) |>
      dplyr::mutate(
        ind = as.character(rows$ind[[i]]), SampleID = rows$SampleID[[i]],
        stim = rows$stim[[i]], cyt = unname(channels[.data$chnl])
      )
  })))
}

# One population/stimulation/cytokine per figure, including finite fallbacks.
.acsCytofThresholdPlot <- function(thresholds, cells, title, bins = 20L) {
  methods <- c("stimgate", "fbeta", "tailgate")
  labels <- .analysis_method_labels[methods]
  labels[["fbeta"]] <- "F-beta (zeros removed)"
  coverage <- vapply(methods, function(m) {
    d <- thresholds[thresholds$method == m, ]
    paste0(
      labels[[m]], ": ",
      sum(is.finite(d$threshold)), "/", nrow(d),
      " finite, ", sum(d$thresholdFallbackUsed %in% TRUE), " fallback"
    )
  }, character(1))
  d <- thresholds[is.finite(thresholds$threshold), ]
  x <- cells$expression[is.finite(cells$expression)]
  # Always retain the whole expression distribution, even if every gate failed.
  span <- range(c(d$threshold, x), finite = TRUE)
  if (!all(is.finite(span))) {
    stop("No finite ACS thresholds or expression for ", title, ".")
  }
  if (diff(span) == 0) span <- span + c(-0.5, 0.5)
  breaks <- seq(span[[1]], span[[2]], length.out = bins + 1L)
  steps <- dplyr::bind_rows(lapply(methods, function(m) {
    values <- d$threshold[d$method == m]
    if (!length(values)) {
      return(NULL)
    }
    h <- graphics::hist(values, breaks = breaks, plot = FALSE)
    data.frame(
      method = m, x = c(breaks, tail(breaks, 1)), y = c(0, h$density, 0)
    )
  }))
  peak <- if (nrow(steps)) max(steps$y) else 1
  density <- data.frame(x = numeric(), y = numeric())
  if (length(x) > 1L) {
    den <- stats::density(x, from = span[[1]], to = span[[2]])
    density <- data.frame(x = den$x, y = den$y / max(den$y) * peak)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_area(
      data = density, ggplot2::aes(x = .data$x, y = .data$y),
      fill = "grey88", colour = "grey60", linewidth = 0.3
    )
  if (nrow(steps)) {
    steps$method <- .analysis_method_factor(steps$method)
    p <- p + ggplot2::geom_step(
      data = steps,
      ggplot2::aes(
        x = .data$x, y = .data$y,
        colour = .data$method, linetype = .data$method
      ),
      direction = "hv", linewidth = 0.6
    )
  }
  p +
    ggplot2::scale_colour_manual(
      values = .analysis_method_colours[methods], labels = labels,
      limits = .analysis_method_limits, name = "Method"
    ) +
    ggplot2::scale_linetype_manual(
      values = .analysis_method_linetypes[methods], labels = labels,
      limits = .analysis_method_limits, name = "Method"
    ) +
    ggplot2::labs(
      x = "Threshold / expression (asinh scale)", y = "Threshold density",
      caption = paste0(
        paste(coverage, collapse = "\n"), "\n",
        "Grey: sampled stimulated expression (", length(x), " cells, ",
        length(unique(cells$ind)), " tubes); ",
        "peak scaled to tallest threshold bin.\n",
        "Up to 2,000 cells per tube; exact zeros retained. ",
        "Finite fallback gates included."
      )
    ) +
    .analysis_theme(grid = "none") + ggplot2::theme(legend.position = "bottom")
}
