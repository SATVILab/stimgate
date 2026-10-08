# Analysis 12A: lab location-shift demonstration. Source analysis-runtime.R and
# analysis-plot-style.R first. Plot builders have no filesystem side effects.
.simClusterLabData <- function(n_sample_lab = 6L, n_cell = 2000L, lab_shift = 4) {
  stopifnot(
    length(n_sample_lab) == 1L, is.finite(n_sample_lab),
    n_sample_lab >= 2L, n_sample_lab == as.integer(n_sample_lab),
    length(n_cell) == 1L, is.finite(n_cell), n_cell >= 20L,
    n_cell == as.integer(n_cell),
    length(lab_shift) == 1L, is.finite(lab_shift), lab_shift >= 0
  )
  n_sample <- 2L * n_sample_lab
  experiment <- simcyto::simCytExperiment(
    nSample = n_sample, nMarker = 1L, nCondition = 2L, nCluster = 2L,
    nCellByCondition = c(n_cell, n_cell),
    transformationFunc = simcyto::simCytTransformIdentity(),
    mixtureType = "gaussianOnly", meanExprMat = matrix(c(10, 20), ncol = 1L),
    clusterLabelVec = c("gn", "gp"), probVecUns = c(0.98, 0.02),
    probExact = TRUE, probResponseVecByStimCondition = list(c(-0.2, 0.2)),
    samplePerturbationSd = 0.05, conditionPerturbationSd = 0,
    clusterPerturbationSd = 0, covEvMin = 1, covEvMax = 1
  )
  samples <- tibble::tibble(
    sample = sprintf("sample-%02d", seq_len(n_sample)),
    lab = rep(c("Lab A", "Lab B"), each = n_sample_lab),
    shift = rep(c(0, lab_shift), each = n_sample_lab),
    ind_uns = as.character(2L * seq_len(n_sample) - 1L),
    ind = as.character(2L * seq_len(n_sample))
  )
  matrices <- lapply(seq_along(experiment$flowFrameList), function(i) {
    sample_i <- (i + 1L) %/% 2L
    x <- flowCore::exprs(experiment$flowFrameList[[i]])[, 1L]
    matrix(x + samples$shift[[sample_i]], ncol = 1L,
      dimnames = list(NULL, "Marker"))
  })
  names(matrices) <- paste0(rep(samples$sample, each = 2L), c("-uns", "-stim"))
  batch_list <- stats::setNames(lapply(seq_len(n_sample), function(i) {
    c(2L * i - 1L, 2L * i)
  }), samples$sample)
  expression <- dplyr::bind_rows(lapply(seq_along(matrices), function(i) {
    sample_i <- (i + 1L) %/% 2L
    tibble::tibble(
      sample = samples$sample[[sample_i]], lab = samples$lab[[sample_i]],
      condition = if (i %% 2L == 1L) "Unstimulated" else "Stimulated",
      expression = matrices[[i]][, 1L]
    )
  }))
  list(matrices = matrices, batch_list = batch_list,
    samples = samples, expression = expression)
}

.simClusterLabAllocations <- function(samples, details) {
  # The quantile table carries both original and adjusted thresholds and grp.
  # Do not read all cluster_final rows: they also include final frequency rows.
  cluster_rows <- tibble::tibble(
    ind = character(), cluster = character(), threshold_original = double(),
    threshold_adjusted = double(), direct_threshold = logical(),
    cluster_action = character(), cluster_reason = character()
  )
  if (nrow(details) > 0L) {
    selected <- details |>
      dplyr::filter(.data$detailObject == "locClusterQuantileTbl",
        .data$detailPathStage == "init")
    if (nrow(selected) > 0L) {
      cluster_rows <- selected |>
        dplyr::transmute(
          ind = as.character(.data$ind), cluster = as.character(.data$grp),
          threshold_original = .data$cpOrigQuantMin,
          threshold_adjusted = .data$cpJoinTgOrig,
          direct_threshold = .data$locGeneratedDirect,
          cluster_action = .data$locClusterAction,
          cluster_reason = .data$locClusterReason
        )
      if (anyDuplicated(cluster_rows$ind)) {
        stop("Expected one initial cluster row per stimulated sample.")
      }
      if (any(!cluster_rows$ind %in% samples$ind)) {
        stop("Cluster details contain an unknown stimulated sample.")
      }
    }
  }
  samples |>
    dplyr::left_join(cluster_rows, by = "ind")
}

.simClusterLabSummary <- function(allocations) {
  assigned <- allocations |>
    dplyr::filter(!is.na(.data$cluster))
  composition <- assigned |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(n_labs = dplyr::n_distinct(.data$lab), .groups = "drop")
  tibble::tibble(
    n_samples = nrow(allocations), n_assigned = nrow(assigned),
    n_unassigned = nrow(allocations) - nrow(assigned),
    n_clusters = nrow(composition),
    n_mixed_lab_clusters = sum(composition$n_labs > 1L),
    labs_separated = if (nrow(assigned) != nrow(allocations)) NA else
      nrow(composition) > 0L && all(composition$n_labs == 1L)
  )
}

# The local-FDR threshold method StimGate resolved and saved for every channel.
# It must equal the requested method, so output rows record what was used.
.simClusterThresholdMethod <- function(path_project, requested) {
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

# `loc_threshold_method` is passed to stimControl(locThresholdMethod = ).
.simClusterLabRun <- function(seed = 558L, n_sample_lab = 6L,
    n_cell = 2000L, lab_shift = 4, loc_threshold_method = "region") {
  if (!is.character(loc_threshold_method) ||
      length(loc_threshold_method) != 1L) {
    stop("loc_threshold_method must be one character value.")
  }
  path_project <- tempfile("stimgate-cluster-lab-")
  on.exit(unlink(path_project, recursive = TRUE), add = TRUE)
  old_intermediate <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_)
  on.exit({
    if (is.na(old_intermediate)) Sys.unsetenv("STIMGATE_INTERMEDIATE") else
      Sys.setenv(STIMGATE_INTERMEDIATE = old_intermediate)
  }, add = TRUE)
  .analysis_with_seed(seed, {
    data <- .simClusterLabData(n_sample_lab, n_cell, lab_shift)
    Sys.setenv(STIMGATE_INTERMEDIATE = "all")
    stimgate::gateStim(
      pathProject = path_project, .data = data$matrices,
      batchList = data$batch_list, chnl = "Marker", bw = 0.3, biasUns = 0.3,
      control = stimgate::stimControl(
        clusterGates = TRUE, calcCytPosGates = FALSE, bwCluster = 0.3,
        bwMin = "none", bwMax = "none",
        locThresholdMethod = loc_threshold_method
      )
    )
    method <- .simClusterThresholdMethod(path_project, loc_threshold_method)
    details <- stimgate::getStimGatesDetailed(path_project, chnl = "Marker")
    allocations <- .simClusterLabAllocations(data$samples, details) |>
      dplyr::mutate(locThresholdMethod = .env$method)
    final_gates <- stimgate::getStimGates(path_project, chnl = "Marker") |>
      dplyr::mutate(locThresholdMethod = .env$method)
    list(expression = data$expression, allocations = allocations,
      summary = .simClusterLabSummary(allocations), details = details,
      final_gates = final_gates)
  })
}

.simClusterLabDistributionPlot <- function(expression) {
  ggplot2::ggplot(expression,
    ggplot2::aes(x = .data$expression, colour = .data$lab, group = .data$sample)) +
    ggplot2::geom_density(bw = 0.3, linewidth = 0.4, alpha = 0.65) +
    ggplot2::facet_wrap(~condition, ncol = 1L) +
    ggplot2::scale_colour_manual(values = c("Lab A" = "#0072B2", "Lab B" = "#D55E00")) +
    ggplot2::labs(x = "Marker expression", y = "Density", colour = "Lab") +
    .analysis_theme()
}

.simClusterLabAllocationPlot <- function(allocations) {
  plot_data <- allocations |>
    dplyr::mutate(cluster_display = dplyr::coalesce(.data$cluster, "Unassigned"))
  ggplot2::ggplot(plot_data,
    ggplot2::aes(x = .data$cluster_display, y = .data$sample, colour = .data$lab)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c("Lab A" = "#0072B2", "Lab B" = "#D55E00")) +
    ggplot2::labs(x = "Assigned cluster (arbitrary label)", y = "Sample", colour = "Lab") +
    .analysis_theme()
}
