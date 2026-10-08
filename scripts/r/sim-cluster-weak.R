# Analysis 12B: evaluate threshold sharing on unchanged, truth-labelled cells.
# Source analysis-runtime.R, analysis-plot-style.R and sim-compare-freq_bs.R first.

.simClusterWeakCells <- function(settings) {
  counts <- unlist(settings[c("n_strong", "n_weak", "n_cell")])
  probs <- unlist(settings[c("strong_prob", "weak_prob")])
  if (length(counts) != 3L || any(!is.finite(counts)) ||
      any(counts < 1 | counts != floor(counts)) ||
      length(probs) != 2L || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1) || probs[[2]] >= probs[[1]]) {
    stop("Use positive integer sample/cell counts and 0 < weak_prob < strong_prob < 1.")
  }
  if (!is.finite(settings$mean_pos) || settings$mean_pos <= 0 ||
      !is.finite(settings$variance) || settings$variance <= 0) {
    stop("mean_pos and variance must be finite and positive.")
  }
  response <- c(rep("Strong", settings$n_strong), rep("Weak", settings$n_weak))
  probability <- ifelse(response == "Strong", settings$strong_prob, settings$weak_prob)
  # Draw all sample seeds before gating; algorithm RNG cannot change the cells.
  seeds <- sample.int(.Machine$integer.max, length(response), replace = TRUE)
  experiments <- lapply(seq_along(response), function(i) {
    .analysis_with_seed(seeds[[i]], simcyto::simCytExperiment(
      nSample = 1L, nMarker = 1L, nCondition = 2L, nCluster = 2L,
      nCellByCondition = rep(settings$n_cell, 2L),
      transformationFunc = simcyto::simCytTransformIdentity(),
      mixtureType = "gaussianOnly",
      meanExprMat = matrix(c(0, settings$mean_pos), ncol = 1L),
      clusterLabelVec = c("gn", "gp"), probVecUns = c(1, 0),
      probExact = TRUE,
      probResponseVecByStimCondition = list(c(-probability[[i]], probability[[i]])),
      samplePerturbationSd = 0, conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      covEvMin = settings$variance, covEvMax = settings$variance
    ))
  })
  frames <- unlist(lapply(experiments, `[[`, "flowFrameList"), recursive = FALSE)
  labels <- unlist(lapply(experiments, `[[`, "labelsList"), recursive = FALSE)
  names(frames) <- paste0("tube", seq_along(frames))
  gs <- flowWorkspace::GatingSet(methods::as(frames, "flowSet"))
  metadata <- tibble::tibble(
    sample = paste0(tolower(response), "-", stats::ave(seq_along(response), response, FUN = seq_along)),
    ind = as.character(2L * seq_along(response)), response = response,
    response_probability = probability, sample_seed = seeds,
    scenario_seed = settings$seed, n_cell_per_tube = settings$n_cell,
    negative_mean = 0, responder_mean = settings$mean_pos,
    component_variance = settings$variance,
    bandwidth = settings$bw, unstimulated_bias = settings$bias_uns
  )
  cells <- dplyr::bind_rows(lapply(seq_along(frames), function(i) {
    x <- flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]], "root"))[, "F1"]
    if (length(labels[[i]]) != length(x) || any(!labels[[i]] %in% c("gn", "gp"))) {
      stop("Simulation labels do not match the GatingSet cells.")
    }
    tibble::tibble(
      ind = as.character(if (i %% 2L) i + 1L else i),
      condition = if (i %% 2L) "Unstimulated" else "Stimulated",
      cell = seq_along(x), expression = x, label = as.character(labels[[i]])
    )
  })) |>
    dplyr::left_join(metadata, by = "ind")
  list(gs = gs, cells = cells, settings_table = metadata,
       batches = lapply(seq_along(response), function(i) c(2L * i - 1L, 2L * i)))
}

.simClusterWeakScore <- function(cells, gates) {
  dplyr::bind_rows(lapply(seq_len(nrow(gates)), function(i) {
    row <- gates[i, , drop = FALSE]
    tube <- cells[cells$ind == row$ind & cells$condition == "Stimulated", ]
    if (!nrow(tube)) stop("No cells for gate sample: ", row$ind)
    counts <- .simCompareConfusionCounts(tube$expression, tube$label, row$threshold)
    tp <- counts$nTruePos
    fp <- counts$nFalsePos
    n_truth <- sum(tube$label == "gp")
    n_negative <- sum(tube$label == "gn")
    dplyr::bind_cols(row, tibble::as_tibble(counts), tibble::tibble(
      n_responders = n_truth, n_negative = n_negative,
      nPosStim = tp + fp,
      sensitivity = if (n_truth > 0) tp / n_truth else NA_real_,
      fdp = if (!is.na(tp + fp) && tp + fp > 0) fp / (tp + fp) else NA_real_,
      false_positive_rate = if (n_negative > 0) fp / n_negative else NA_real_
    ))
  }))
}

.simClusterWeakChanges <- function(scores) {
  before <- scores[scores$stage == "Before", ]
  after <- scores[scores$stage == "After", ]
  if (anyDuplicated(before$ind) || anyDuplicated(after$ind) ||
      !setequal(before$ind, after$ind)) stop("Expected one before/after pair per sample.")
  dplyr::inner_join(
    before |> dplyr::select("ind", "sample", "response", "threshold", "nTruePos", "nFalsePos", "sensitivity", "fdp"),
    after |> dplyr::select("ind", "threshold", "nTruePos", "nFalsePos", "sensitivity", "fdp"),
    by = "ind", suffix = c("_before", "_after")
  ) |>
    dplyr::mutate(
      threshold_change = .data$threshold_after - .data$threshold_before,
      responders_recovered = .data$nTruePos_after - .data$nTruePos_before,
      additional_false_positives = .data$nFalsePos_after - .data$nFalsePos_before,
      sensitivity_change = .data$sensitivity_after - .data$sensitivity_before,
      fdp_change = .data$fdp_after - .data$fdp_before
    )
}

# Cached-result semantics. v2: `settings$loc_threshold_method` is passed to
# stimControl(locThresholdMethod = ) and recorded in the gate rows; v1 caches
# used probability-sum matching without recording it.
.simClusterWeakSemantics <- "cluster-weak-v5"

# The local-FDR threshold method StimGate resolved and saved for every channel.
.simClusterWeakThresholdMethod <- function(path_project, requested) {
  saved <- stimgate::stimgateMetaReadSettingsChnls(path_project)
  method <- unique(vapply(saved, function(x) {
    as.character(x$locThresholdMethod %||% NA_character_)
  }, character(1L)))
  if (!identical(method, requested)) {
    stop("StimGate saved locThresholdMethod '",
      paste(method, collapse = "', '"), "' but '", requested,
      "' was requested.")
  }
  method
}

.simClusterWeakRun <- function(settings) {
  .analysis_require_packages(c("simcyto", "flowCore", "flowWorkspace"))
  loc_threshold_method <- settings$loc_threshold_method
  if (!is.character(loc_threshold_method) ||
      length(loc_threshold_method) != 1L) {
    stop("settings$loc_threshold_method must be set explicitly.")
  }
  .analysis_with_seed(settings$seed, {
    data <- .simClusterWeakCells(settings)
    path <- tempfile("cluster-weak-")
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    old <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_)
    on.exit({
      if (is.na(old)) Sys.unsetenv("STIMGATE_INTERMEDIATE") else
        Sys.setenv(STIMGATE_INTERMEDIATE = old)
    }, add = TRUE)
    Sys.setenv(STIMGATE_INTERMEDIATE = "all")
    stimgate::gateStim(
      pathProject = path, .data = data$gs, batchList = data$batches,
      popGate = "root", marker = "MarkerF1", bw = settings$bw,
      biasUns = settings$bias_uns,
      control = stimgate::stimControl(
        clusterGates = TRUE, calcCytPosGates = FALSE,
        bwCluster = settings$bw, locEnforceShapeThreshold = FALSE,
        locThresholdMethod = loc_threshold_method
      )
    )
    method <- .simClusterWeakThresholdMethod(path, loc_threshold_method)
    details <- stimgate::getStimGatesDetailed(path, chnl = "F1")
    allocations <- details |>
      dplyr::filter(.data$detailObject == "locClusterQuantileTbl", .data$detailPathStage == "init")
    if (nrow(allocations) != nrow(data$settings_table) || anyDuplicated(allocations$ind) ||
        !setequal(as.character(allocations$ind), data$settings_table$ind)) {
      stop("Expected one initial cluster diagnostic row for every stimulated sample.")
    }
    applied <- stimgate::getStimGates(path, chnl = "F1")
    # getStimGates() names the clustered gate "loc_minClust" and the original "loc_min".
    final <- applied |>
      dplyr::filter(.data$gateName == "loc_minClust")
    original <- applied |>
      dplyr::filter(.data$gateName == "loc_min")
    check <- dplyr::left_join(
      allocations, final |> dplyr::select("ind", final_gate = "gate"), by = "ind"
    )
    if (nrow(final) != nrow(allocations) ||
        !isTRUE(all.equal(check$cpJoinTgOrig, check$final_gate, check.attributes = FALSE))) {
      stop("Cluster diagnostic gates disagree with the final applied clustered gates.")
    }
    original_check <- dplyr::left_join(
      allocations, original |> dplyr::select("ind", original_gate = "gate"), by = "ind"
    )
    if (nrow(original) != nrow(allocations) ||
        !isTRUE(all.equal(original_check$cpOrigQuantMin, original_check$original_gate, check.attributes = FALSE))) {
      stop("Pre-cluster diagnostic gates disagree with the applied original gates.")
    }
    gates <- dplyr::bind_rows(
      allocations |> dplyr::transmute(ind = as.character(.data$ind), stage = "Before", threshold = .data$cpOrigQuantMin),
      allocations |> dplyr::transmute(ind = as.character(.data$ind), stage = "After", threshold = .data$cpJoinTgOrig)
    ) |> dplyr::left_join(data$settings_table, by = "ind") |>
      dplyr::mutate(locThresholdMethod = .env$method)
    scores <- .simClusterWeakScore(data$cells, gates)
    # Both gate versions are applied by the package and appear in its statistics.
    stats <- stimgate::getStimStats(path) |>
      dplyr::filter(.data$gateName %in% c(original$gateName, final$gateName), grepl("~\\+~", .data$cytCombn)) |>
      dplyr::transmute(ind = as.character(.data$ind),
        stage = ifelse(.data$gateName %in% original$gateName, "Before", "After"),
        package_nPosStim = .data$countStim)
    scored <- dplyr::left_join(scores, stats, by = c("ind", "stage"))
    if (nrow(scored) != nrow(scores) || anyNA(scored$package_nPosStim) ||
        !isTRUE(all.equal(as.numeric(scored$nPosStim), as.numeric(scored$package_nPosStim)))) {
      stop("Truth-based positive counts disagree with package statistics.")
    }
    list(semantics = .simClusterWeakSemantics, settings = settings,
         settings_table = data$settings_table, cells = data$cells,
         allocations = allocations, scores = scored,
         changes = .simClusterWeakChanges(scored))
  })
}

.simClusterWeakValidate <- function(result, settings) {
  if (!identical(result$semantics, .simClusterWeakSemantics) || !identical(result$settings, settings)) {
    stop("Weak-response cache settings changed; rerun simulations with the same profile.")
  }
  if (!"locThresholdMethod" %in% names(result$scores) ||
      !identical(unique(result$scores$locThresholdMethod), settings$loc_threshold_method)) {
    stop("Weak-response cache does not record the requested locThresholdMethod.")
  }
  expected <- 2L * (settings$n_strong + settings$n_weak)
  if (nrow(result$scores) != expected || anyDuplicated(result$scores[, c("ind", "stage")]) ||
      nrow(result$allocations) != expected / 2L ||
      !setequal(result$scores$ind, result$settings_table$ind)) {
    stop("Incomplete cached sample/gate pairs.")
  }
  counts <- .simClusterWeakScore(result$cells, result$scores[, c("ind", "stage", "threshold")])
  totals <- counts$n_responders + counts$n_negative
  if (any(totals != settings$n_cell)) stop("Incomplete cached stimulated cells.")
  for (nm in c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg", "nPosStim")) {
    if (!identical(counts[[nm]], result$scores[[nm]])) stop("Invalid cached truth counts: ", nm)
  }
  if (!isTRUE(all.equal(as.numeric(counts$nPosStim), as.numeric(result$scores$package_nPosStim)))) {
    stop("Cached counts disagree with package statistics.")
  }
  if (!isTRUE(all.equal(.simClusterWeakChanges(result$scores), result$changes))) {
    stop("Cached paired changes disagree with scored gates.")
  }
  invisible(result)
}

.simClusterWeakDistributions <- function(result, tail = FALSE) {
  # Histogram heights are fractions of the whole tube, including truth subgroups.
  # This keeps a rare response's size visible instead of normalising it to one.
  cells <- result$cells |>
    dplyr::mutate(component = dplyr::case_when(
      .data$condition == "Unstimulated" ~ "Unstimulated",
      .data$label == "gp" ~ "Stimulated responders",
      TRUE ~ "Stimulated negatives"
    ))
  breaks <- seq(min(cells$expression) - 0.1, max(cells$expression) + 0.1, length.out = 180L)
  histograms <- cells |>
    dplyr::group_by(.data$sample, .data$condition, .data$component) |>
    dplyr::group_modify(function(.x, .y) {
      h <- graphics::hist(.x$expression, breaks = breaks, plot = FALSE)
      tibble::tibble(expression = h$mids, fraction = h$counts / result$settings$n_cell)
    }) |>
    dplyr::ungroup()
  gates <- result$scores |>
    dplyr::mutate(stage = factor(.data$stage, levels = c("Before", "After")))
  p <- ggplot2::ggplot(histograms, ggplot2::aes(.data$expression, .data$fraction, colour = .data$component)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data = gates,
      ggplot2::aes(xintercept = .data$threshold, linetype = .data$stage, group = .data$stage),
      inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(sample), ncol = 2) +
    ggplot2::scale_colour_manual(values = c("Unstimulated" = "#0072B2",
      "Stimulated negatives" = "#999999", "Stimulated responders" = "#D55E00")) +
    ggplot2::scale_linetype_manual(values = c(Before = "dashed", After = "solid")) +
    ggplot2::labs(x = "Marker expression", y = "Fraction of tube per bin", colour = NULL, linetype = "Gate") +
    .analysis_theme()
  if (tail) {
    upper <- max(cells$expression)
    p <- p + ggplot2::coord_cartesian(
      xlim = c(result$settings$mean_pos / 2, upper),
      ylim = c(0, max(histograms$fraction[histograms$expression >= result$settings$mean_pos / 2]))
    )
  }
  p
}

.simClusterWeakOutcomes <- function(result) {
  long <- result$scores |>
    tidyr::pivot_longer(c("threshold", "sensitivity", "fdp", "false_positive_rate"),
      names_to = "metric", values_to = "value") |>
    dplyr::mutate(stage = factor(.data$stage, levels = c("Before", "After")),
      metric = factor(.data$metric, levels = c("threshold", "sensitivity", "fdp", "false_positive_rate"),
        labels = c("Gate", "Sensitivity", "False discovery proportion", "False positive rate")))
  ggplot2::ggplot(long, ggplot2::aes(.data$stage, .data$value, group = .data$sample, colour = .data$response)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", ncol = 2) +
    ggplot2::scale_colour_manual(values = c(Strong = "#0072B2", Weak = "#D55E00")) +
    ggplot2::labs(x = NULL, y = NULL, colour = "Response") + .analysis_theme()
}
