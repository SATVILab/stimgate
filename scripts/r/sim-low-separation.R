# Low-separation demonstration of the cytokine-positive gates (Analysis 11).
#
# Two markers (TNF, IFNg) are simulated with simcyto on the Gaussian
# transformation. Each scenario is one dataset gated once with gateStim(), and
# the ordinary and cytokine-positive gates are then applied to the same cells.
# Truth comes from the simulated component labels.

.simLowSepMarkers <- c(tnf = "TNF", ifng = "IFNg")
.simLowSepChnl <- c(tnf = "F1", ifng = "F2")
.simLowSepCombn <- c("TNF-IFNg-", "TNF+IFNg-", "TNF-IFNg+", "TNF+IFNg+")
.simLowSepGateTypes <- c(ordinary = "Ordinary gates", cytpos = "Cytokine-positive gates")
.simLowSepGateColours <- c(ordinary = "#7F7F7F", cytpos = "#D55E00")

# Quantities whose background-subtracted frequencies are compared. Marginal
# frequencies count every cell positive for a cytokine; exclusive frequencies
# count one combination only.
.simLowSepQuantities <- tibble::tribble(
  ~quantity, ~quantity_type, ~quantity_lab,
  "TNF+", "marginal", "TNF+ (any)",
  "IFNg+", "marginal", "IFN\u03b3+ (any)",
  "TNF+IFNg+", "exclusive", "TNF+ IFN\u03b3+",
  "TNF+IFNg-", "exclusive", "TNF+ IFN\u03b3- only",
  "TNF-IFNg+", "exclusive", "TNF- IFN\u03b3+ only"
)

# Fixed scientific settings shared by every scenario.
.simLowSepMainSettings <- function() {
  list(
    transformation = "gaussian",
    cov_ev = 1.5,
    background_relative_to_response = 0.2,
    n_cell_uns_relative_to_stim = 1,
    # Share of the responding cells in each positive combination.
    response_split = c("TNF+IFNg-" = 0.25, "TNF-IFNg+" = 0.15, "TNF+IFNg+" = 0.6),
    prob_exact = TRUE,
    calc_cyt_pos_gates = TRUE,
    cluster_gates = TRUE
  )
}

# Full scenario grid. IDs and seeds are assigned here, before any dev/quick
# filtering, so a scenario keeps its data whatever subset is run.
.simLowSepGrid <- function(simulation_seed) {
  separation <- tibble::tibble(
    separation = c("low", "lower"),
    mean_pos = c(4.5, 3.5)
  )
  response <- tibble::tibble(
    response_level = c("lower", "higher"),
    prob_response = c(0.01, 0.05)
  )
  grid <- tidyr::expand_grid(
    separation,
    response,
    n_cell = c(1e5, 5e3)
  )
  seeds <- .analysis_with_seed(simulation_seed, {
    sample.int(.Machine$integer.max, nrow(grid))
  })
  dplyr::mutate(
    grid,
    sim_id = seq_len(nrow(grid)),
    sim_seed = as.integer(seeds),
    .before = 1L
  )
}

.simLowSepScenarioLab <- function(grid) {
  paste0(
    grid$separation, " separation (mean ", grid$mean_pos, "), ",
    grid$response_level, " response (", 100 * grid$prob_response, "%)"
  )
}

# Cell count without scientific notation, for figure file names.
.simLowSepCellKey <- function(n_cell) {
  format(n_cell, scientific = FALSE, trim = TRUE)
}

.simLowSepCellLab <- function(n_cell) {
  paste0(format(n_cell, big.mark = ",", scientific = FALSE, trim = TRUE), " cells")
}

# Simulate one dataset: `n_sample` unstimulated/stimulated pairs.
.simLowSepSimulate <- function(row, n_sample, settings = .simLowSepMainSettings()) {
  if (!identical(settings$transformation, "gaussian")) {
    stop("This demonstration uses the Gaussian transformation only.")
  }
  split <- settings$response_split[.simLowSepCombn[-1L]]
  if (!isTRUE(all.equal(sum(split), 1))) {
    stop("response_split must sum to one.")
  }
  prob_resp <- row$prob_response * split
  prob_uns <- prob_resp * settings$background_relative_to_response
  n_cell_uns <- round(row$n_cell * settings$n_cell_uns_relative_to_stim)
  res <- simcyto::simCytExperiment(
    nSample = n_sample,
    nMarker = 2L,
    nCondition = 2L,
    nCluster = 4L,
    nCellByCondition = c(n_cell_uns, row$n_cell),
    transformationFunc = simcyto::simCytTransformGaussian(),
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(
      c(0, 0, row$mean_pos, 0, 0, row$mean_pos, row$mean_pos, row$mean_pos),
      ncol = 2L, byrow = TRUE
    ),
    clusterLabelVec = .simLowSepCombn,
    probVecUns = c(1 - sum(prob_uns), prob_uns),
    probExact = settings$prob_exact,
    probResponseVecByStimCondition = list(c(
      -sum(prob_resp - prob_uns), prob_resp - prob_uns
    )),
    covEvMin = settings$cov_ev,
    covEvMax = settings$cov_ev
  )
  res$flowFrameList <- lapply(res$flowFrameList, function(ff) {
    flowCore::markernames(ff) <- stats::setNames(
      unname(.simLowSepMarkers), unname(.simLowSepChnl)
    )
    ff
  })
  res
}

# Combination calls for one tube. The cytokine-positive rule matches the
# package: positive by the ordinary gate, or above the cytokine-positive gate
# while ordinarily positive for the other cytokine.
.simLowSepCall <- function(x_tnf, x_ifng, gate, gate_type) {
  pos_tnf <- x_tnf > gate[["tnf"]]
  pos_ifng <- x_ifng > gate[["ifng"]]
  if (identical(gate_type, "cytpos")) {
    pos_tnf_cyt <- pos_tnf | (x_tnf > gate[["tnf_cyt"]] & pos_ifng)
    pos_ifng <- pos_ifng | (x_ifng > gate[["ifng_cyt"]] & pos_tnf)
    pos_tnf <- pos_tnf_cyt
  } else if (!identical(gate_type, "ordinary")) {
    stop("gate_type must be \"ordinary\" or \"cytpos\".")
  }
  .simLowSepCombn[1L + pos_tnf + 2L * pos_ifng]
}

# Gate values for one stimulated sample, from getStimGates().
.simLowSepSampleGate <- function(gate_tbl, ind) {
  rows <- gate_tbl[as.character(gate_tbl$ind) == as.character(ind), ]
  out <- vapply(names(.simLowSepChnl), function(m) {
    r <- rows[rows$chnl == .simLowSepChnl[[m]], ]
    if (nrow(r) != 1L) {
      stop("Expected one gate for channel ", .simLowSepChnl[[m]], " in sample ", ind, ".")
    }
    c(r$gate, r$gateCyt)
  }, numeric(2))
  c(
    tnf = out[[1, "tnf"]], ifng = out[[1, "ifng"]],
    tnf_cyt = out[[2, "tnf"]], ifng_cyt = out[[2, "ifng"]]
  )
}

# Confusion counts (truth by call) for every tube and gate type.
.simLowSepCounts <- function(expr_list, labels_list, gate_tbl, n_sample) {
  purrr::map_dfr(seq_len(n_sample), function(s) {
    ind_uns <- 2L * s - 1L
    ind_stim <- 2L * s
    gate <- .simLowSepSampleGate(gate_tbl, ind_stim)
    purrr::map_dfr(c(uns = ind_uns, stim = ind_stim), function(ind) {
      ex <- expr_list[[ind]]
      truth <- factor(labels_list[[ind]], levels = .simLowSepCombn)
      if (anyNA(truth)) stop("Unexpected simulated labels in tube ", ind, ".")
      purrr::map_dfr(names(.simLowSepGateTypes), function(gt) {
        call <- factor(
          .simLowSepCall(ex[, .simLowSepChnl[["tnf"]]], ex[, .simLowSepChnl[["ifng"]]], gate, gt),
          levels = .simLowSepCombn
        )
        tibble::as_tibble(as.data.frame(table(truth = truth, call = call), stringsAsFactors = FALSE)) |>
          dplyr::rename(n = "Freq") |>
          dplyr::mutate(gate_type = gt, .before = 1L)
      })
    }, .id = "tube") |>
      dplyr::mutate(sample = s, .before = 1L)
  }) |>
    dplyr::mutate(n = as.integer(.data$n))
}

# The package's saved combination counts use the cytokine-positive rule; the
# counts recomputed here must reproduce them exactly.
.simLowSepValidateCounts <- function(counts, stats_tbl) {
  combn_code <- c(
    "TNF+IFNg-" = "F1~+~F2~-~", "TNF-IFNg+" = "F1~-~F2~+~",
    "TNF+IFNg+" = "F1~+~F2~+~", "TNF-IFNg-" = "F1~-~F2~-~"
  )
  mine <- counts |>
    dplyr::filter(.data$gate_type == "cytpos") |>
    dplyr::group_by(.data$sample, .data$tube, .data$call) |>
    dplyr::summarise(n = sum(.data$n), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "tube", values_from = "n") |>
    dplyr::mutate(
      ind = as.character(2L * .data$sample),
      cytCombn = unname(combn_code[.data$call])
    )
  joined <- dplyr::inner_join(
    mine,
    dplyr::mutate(stats_tbl, ind = as.character(.data$ind)),
    by = c("ind", "cytCombn")
  )
  if (nrow(joined) != nrow(mine) || nrow(stats_tbl) != nrow(mine)) {
    stop("Package combination statistics do not cover every sample and combination.")
  }
  bad <- joined$stim != joined$countStim | joined$uns != joined$countUns
  if (any(bad)) {
    stop(
      "Recomputed cytokine-positive combination counts differ from getStimStats() in ",
      sum(bad), " rows."
    )
  }
  invisible(TRUE)
}

# Simulate and gate one scenario. Returns gates, confusion counts and the
# stimulated cells of the first sample (prespecified for the figures).
.simLowSepRunScenario <- function(row, n_sample, settings = .simLowSepMainSettings(),
                                  path_project = NULL) {
  if (nrow(row) != 1L) stop("row must have exactly one scenario.")
  .analysis_with_seed(row$sim_seed, {
    sim <- .simLowSepSimulate(row, n_sample, settings)
    gs <- flowWorkspace::GatingSet(methods::as(sim$flowFrameList, "flowSet"))
    path_project <- path_project %||% tempfile(paste0("lowsep-", row$sim_id, "-"))
    on.exit(unlink(path_project, recursive = TRUE), add = TRUE)
    batch_list <- lapply(seq_len(n_sample), function(s) c(2L * s - 1L, 2L * s))
    suppressMessages(gateStim(
      pathProject = path_project, .data = gs, batchList = batch_list,
      marker = unname(.simLowSepMarkers),
      control = stimControl(
        calcCytPosGates = settings$calc_cyt_pos_gates,
        clusterGates = settings$cluster_gates
      )
    ))
    gate_tbl <- getStimGates(path_project)
    # Classify the single-precision expression that StimGate gated.
    expr_list <- lapply(seq_along(gs), function(i) {
      flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]], "root"))
    })
    counts <- .simLowSepCounts(expr_list, sim$labelsList, gate_tbl, n_sample)
    .simLowSepValidateCounts(counts, getStimStats(path_project))
    gates <- purrr::map_dfr(seq_len(n_sample), function(s) {
      g <- .simLowSepSampleGate(gate_tbl, 2L * s)
      rows <- gate_tbl[as.character(gate_tbl$ind) == as.character(2L * s), ]
      tibble::tibble(
        sample = s,
        marker = unname(.simLowSepMarkers),
        gate = unname(g[c("tnf", "ifng")]),
        gate_cyt = unname(g[c("tnf_cyt", "ifng_cyt")]),
        gate_name = rows$gateName[match(.simLowSepChnl, rows$chnl)],
        loc_source = if ("locSource" %in% names(rows)) rows$locSource[match(.simLowSepChnl, rows$chnl)] else NA_character_
      )
    })
    cells <- tibble::tibble(
      tnf = expr_list[[2L]][, .simLowSepChnl[["tnf"]]],
      ifng = expr_list[[2L]][, .simLowSepChnl[["ifng"]]],
      truth = sim$labelsList[[2L]]
    )
    scen <- dplyr::select(row, dplyr::all_of(c(
      "sim_id", "separation", "mean_pos", "response_level", "prob_response", "n_cell"
    )))
    list(
      gates = dplyr::bind_cols(scen, gates),
      counts = dplyr::bind_cols(scen, counts),
      cells = dplyr::bind_cols(scen, cells)
    )
  })
}

.simLowSepRunGrid <- function(grid, n_sample, settings = .simLowSepMainSettings()) {
  out <- lapply(seq_len(nrow(grid)), function(i) {
    message("Low-separation scenario ", grid$sim_id[[i]], " (", i, " of ", nrow(grid), ")")
    .simLowSepRunScenario(grid[i, ], n_sample, settings)
  })
  list(
    gates = purrr::map_dfr(out, "gates"),
    counts = purrr::map_dfr(out, "counts"),
    cells = purrr::map_dfr(out, "cells")
  )
}

# Cache ---------------------------------------------------------------------

.simLowSepCacheSettings <- function(grid, n_sample, settings, simulation_seed, profile) {
  list(
    analysis_semantics_version = "sim-low-separation-v1",
    grid = as.data.frame(grid),
    n_sample = as.integer(n_sample),
    settings = settings,
    simulation_seed = as.integer(simulation_seed),
    profile = profile
  )
}

.simLowSepWriteCache <- function(results, settings, path) {
  .write_rds_atomic(list(
    settings = settings, results = results,
    run_id = .sanitize_run_id(Sys.getenv("ANALYSIS_RUN_ID", unset = ""))
  ), path)
}

.simLowSepReadCache <- function(path, settings, analysis_key) {
  qmd <- "analysis/11-sim-low-separation-cyt-pos.qmd"
  cached <- .analysis_read_rds(path, analysis_key, qmd)
  if (!is.list(cached) || !all(c("settings", "results") %in% names(cached))) {
    .analysis_cache_error(analysis_key, paste0("Cached results at ", path, " are malformed."), qmd)
  }
  .analysis_check_expected_run(cached, analysis_key, qmd)
  if (!isTRUE(all.equal(cached$settings, settings))) {
    .analysis_cache_error(analysis_key, paste0(
      "Cached results at ", path, " were made with different settings ",
      "(scenario grid, sample count, simulation settings or seed)."
    ), qmd)
  }
  cached$results
}

# Summaries -----------------------------------------------------------------

# Per stimulated sample and marker: true positives, cells called positive,
# true positives recovered and negative-cell contamination.
.simLowSepRecovery <- function(counts) {
  stim <- dplyr::filter(counts, .data$tube == "stim")
  scen_cols <- c("sim_id", "separation", "mean_pos", "response_level", "prob_response", "n_cell")
  purrr::map_dfr(c(TNF = "TNF\\+", IFNg = "IFNg\\+"), function(pattern) {
    stim |>
      dplyr::mutate(
        truth_pos = grepl(pattern, .data$truth),
        call_pos = grepl(pattern, .data$call)
      ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(scen_cols, "sample", "gate_type")))) |>
      dplyr::summarise(
        n_cell_tube = sum(.data$n),
        n_true_pos = sum(.data$n[.data$truth_pos]),
        n_called_pos = sum(.data$n[.data$call_pos]),
        n_recovered = sum(.data$n[.data$truth_pos & .data$call_pos]),
        n_contaminating = sum(.data$n[!.data$truth_pos & .data$call_pos]),
        .groups = "drop"
      )
  }, .id = "marker") |>
    dplyr::mutate(
      sensitivity = dplyr::if_else(.data$n_true_pos > 0, .data$n_recovered / .data$n_true_pos, NA_real_),
      fdp = dplyr::if_else(.data$n_called_pos > 0, .data$n_contaminating / .data$n_called_pos, NA_real_)
    )
}

# Background-subtracted marginal and exclusive frequencies, estimated and
# true, per stimulated sample and gate type.
.simLowSepFrequencies <- function(counts) {
  scen_cols <- c("sim_id", "separation", "mean_pos", "response_level", "prob_response", "n_cell")
  group_cols <- c(scen_cols, "sample", "gate_type", "tube")
  totals <- counts |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(n_tube = sum(.data$n), .groups = "drop")
  purrr::map_dfr(.simLowSepQuantities$quantity, function(q) {
    hit <- function(x) {
      if (q == "TNF+") grepl("TNF\\+", x) else if (q == "IFNg+") grepl("IFNg\\+", x) else x == q
    }
    counts |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        n_est = sum(.data$n[hit(.data$call)]),
        n_truth = sum(.data$n[hit(.data$truth)]),
        .groups = "drop"
      ) |>
      dplyr::mutate(quantity = q)
  }) |>
    dplyr::left_join(totals, by = group_cols) |>
    dplyr::mutate(prop_est = .data$n_est / .data$n_tube, prop_truth = .data$n_truth / .data$n_tube) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(c(scen_cols, "sample", "gate_type", "quantity")),
      names_from = "tube",
      values_from = c("n_est", "n_truth", "prop_est", "prop_truth")
    ) |>
    dplyr::mutate(
      prop_bs_est = .data$prop_est_stim - .data$prop_est_uns,
      prop_bs_truth = .data$prop_truth_stim - .data$prop_truth_uns,
      error_bs = .data$prop_bs_est - .data$prop_bs_truth
    ) |>
    dplyr::left_join(.simLowSepQuantities, by = "quantity")
}

# Mean over the samples of one dataset, with how often the cytokine-positive
# gates moved the estimate closer to the truth.
.simLowSepFrequencySummary <- function(freq) {
  scen_cols <- c("sim_id", "separation", "mean_pos", "response_level", "prob_response", "n_cell")
  paired <- freq |>
    dplyr::select(dplyr::all_of(c(scen_cols, "sample", "quantity", "gate_type", "error_bs"))) |>
    tidyr::pivot_wider(names_from = "gate_type", values_from = "error_bs") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen_cols, "quantity")))) |>
    dplyr::summarise(
      n_closer = sum(abs(.data$cytpos) < abs(.data$ordinary)),
      n_further = sum(abs(.data$cytpos) > abs(.data$ordinary)),
      n_unchanged = sum(abs(.data$cytpos) == abs(.data$ordinary)),
      .groups = "drop"
    )
  freq |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen_cols, "quantity", "quantity_type", "gate_type")))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      mean_prop_bs_truth = mean(.data$prop_bs_truth),
      mean_prop_bs_est = mean(.data$prop_bs_est),
      mean_error_bs = mean(.data$error_bs),
      mean_abs_error_bs = mean(abs(.data$error_bs)),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "gate_type",
      values_from = c("mean_prop_bs_est", "mean_error_bs", "mean_abs_error_bs")
    ) |>
    dplyr::left_join(paired, by = c(scen_cols, "quantity"))
}

# Plots (return ggplot objects; no files are written here) -------------------

.simLowSepFacetLabels <- function(tbl) {
  ordered <- dplyr::arrange(tbl, dplyr::desc(.data$mean_pos), .data$prob_response)
  tbl$scenario_lab <- factor(
    .simLowSepScenarioLab(tbl),
    levels = unique(.simLowSepScenarioLab(ordered))
  )
  tbl
}

# Hex plot of TNF against IFNg in the first sample's stimulated tube. Ordinary
# gates span the plot; each cytokine-positive gate is drawn only where the
# other cytokine is ordinarily positive, the only cells it applies to.
.simLowSepPlotHex <- function(cells, gates, bins = 70) {
  cells <- .simLowSepFacetLabels(cells)
  g <- gates |>
    dplyr::filter(.data$sample == 1L) |>
    dplyr::select(dplyr::all_of(c("sim_id", "marker", "gate", "gate_cyt"))) |>
    tidyr::pivot_wider(names_from = "marker", values_from = c("gate", "gate_cyt"))
  lim <- cells |>
    dplyr::group_by(.data$sim_id) |>
    dplyr::summarise(
      x_max = max(.data$tnf), y_max = max(.data$ifng),
      x_min = min(.data$tnf), y_min = min(.data$ifng), .groups = "drop"
    )
  g <- dplyr::left_join(g, lim, by = "sim_id") |>
    dplyr::left_join(dplyr::distinct(cells, .data$sim_id, .data$scenario_lab), by = "sim_id")
  # Regions recovered only by the cytokine-positive gates (empty if no move).
  rect <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab,
      xmin = .data$gate_TNF, xmax = .data$x_max,
      ymin = pmin(.data$gate_cyt_IFNg, .data$gate_IFNg), ymax = .data$gate_IFNg
    ),
    dplyr::transmute(g, .data$scenario_lab,
      xmin = pmin(.data$gate_cyt_TNF, .data$gate_TNF), xmax = .data$gate_TNF,
      ymin = .data$gate_IFNg, ymax = .data$y_max
    )
  ) |>
    dplyr::filter(.data$xmax > .data$xmin, .data$ymax > .data$ymin)
  ordinary <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_TNF, xend = .data$gate_TNF,
      y = .data$y_min, yend = .data$y_max, line = "TNF"),
    dplyr::transmute(g, .data$scenario_lab, x = .data$x_min, xend = .data$x_max,
      y = .data$gate_IFNg, yend = .data$gate_IFNg, line = "IFNg")
  )
  cytpos <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_cyt_TNF, xend = .data$gate_cyt_TNF,
      y = .data$gate_IFNg, yend = .data$y_max, line = "TNF"),
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_TNF, xend = .data$x_max,
      y = .data$gate_cyt_IFNg, yend = .data$gate_cyt_IFNg, line = "IFNg")
  )
  # Distinct groups keep coincident lines (unmoved gates) from being merged.
  ordinary$group <- paste("ordinary", seq_len(nrow(ordinary)))
  cytpos$group <- paste("cytpos", seq_len(nrow(cytpos)))
  ggplot(cells, aes(x = .data$tnf, y = .data$ifng)) +
    geom_hex(bins = bins) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_rect(
      data = rect, inherit.aes = FALSE,
      aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = .simLowSepGateColours[["cytpos"]], alpha = 0.18
    ) +
    geom_segment(
      data = ordinary, inherit.aes = FALSE,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend,
        group = .data$group, colour = "ordinary", linetype = "ordinary"),
      linewidth = 0.6
    ) +
    geom_segment(
      data = cytpos, inherit.aes = FALSE,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend,
        group = .data$group, colour = "cytpos", linetype = "cytpos"),
      linewidth = 0.8
    ) +
    scale_colour_manual(values = .simLowSepGateColours, labels = .simLowSepGateTypes, name = NULL) +
    scale_linetype_manual(values = c(ordinary = "22", cytpos = "solid"), labels = .simLowSepGateTypes, name = NULL) +
    facet_wrap(~ scenario_lab, ncol = 2, labeller = ggplot2::label_wrap_gen(width = 36)) +
    coord_equal() +
    labs(x = "TNF", y = "IFN\u03b3") +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# IFNg density among all stimulated cells and among TNF+ cells (ordinary TNF
# gate), with the ordinary and cytokine-positive IFNg gates.
.simLowSepPlotConditionalDensity <- function(cells, gates, n_grid = 512L) {
  g <- gates |>
    dplyr::filter(.data$sample == 1L) |>
    dplyr::select(dplyr::all_of(c("sim_id", "marker", "gate", "gate_cyt"))) |>
    tidyr::pivot_wider(names_from = "marker", values_from = c("gate", "gate_cyt"))
  dens <- cells |>
    dplyr::left_join(g, by = "sim_id") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(
      "sim_id", "separation", "mean_pos", "response_level", "prob_response", "n_cell"
    )))) |>
    dplyr::group_modify(function(d, key) {
      sub <- list(all = d$ifng, tnf_pos = d$ifng[d$tnf > d$gate_TNF[[1]]])
      rng <- range(d$ifng)
      purrr::map_dfr(names(sub), function(nm) {
        x <- sub[[nm]]
        if (length(x) < 2L) {
          return(tibble::tibble(ifng = numeric(), density = numeric(), cells = character(), n = integer()))
        }
        de <- stats::density(x, n = n_grid, from = rng[[1]], to = rng[[2]])
        tibble::tibble(ifng = de$x, density = de$y, cells = nm, n = length(x))
      })
    }) |>
    dplyr::ungroup() |>
    .simLowSepFacetLabels()
  cell_labs <- c(all = "All stimulated cells", tnf_pos = "TNF+ cells (ordinary TNF gate)")
  g <- g |>
    dplyr::left_join(dplyr::distinct(dens, .data$sim_id, .data$scenario_lab), by = "sim_id") |>
    dplyr::mutate(
      shift = .data$gate_IFNg - .data$gate_cyt_IFNg,
      shift_lab = dplyr::if_else(
        .data$shift > 0,
        paste0("Lowered by ", formatC(.data$shift, format = "f", digits = 2), "\namong TNF+ cells"),
        "Gate not lowered"
      )
    )
  lines <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_IFNg, gate_type = "ordinary"),
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_cyt_IFNg, gate_type = "cytpos")
  ) |>
    dplyr::mutate(group = seq_len(dplyr::n()))
  arrows <- dplyr::filter(g, .data$shift > 0)
  ggplot(dens, aes(x = .data$ifng, y = .data$density)) +
    geom_line(aes(linetype = .data$cells), linewidth = 0.6) +
    scale_linetype_manual(values = c(all = "solid", tnf_pos = "42"), labels = cell_labs, name = NULL) +
    geom_vline(
      data = lines,
      aes(xintercept = .data$x, colour = .data$gate_type, group = .data$group),
      linewidth = 0.7
    ) +
    geom_segment(
      data = arrows, inherit.aes = FALSE,
      aes(x = .data$gate_IFNg, xend = .data$gate_cyt_IFNg, y = 0, yend = 0),
      colour = .simLowSepGateColours[["cytpos"]], linewidth = 0.8,
      arrow = grid::arrow(length = grid::unit(0.2, "cm"))
    ) +
    geom_text(
      data = g, inherit.aes = FALSE,
      aes(x = Inf, y = Inf, label = .data$shift_lab),
      hjust = 1.02, vjust = 1.3, size = 2.7
    ) +
    scale_colour_manual(values = .simLowSepGateColours, labels = c(
      ordinary = "Ordinary IFN\u03b3 gate", cytpos = "Cytokine-positive IFN\u03b3 gate"
    ), name = NULL) +
    scale_y_sqrt() +
    facet_wrap(~ scenario_lab, ncol = 2, scales = "free_y", labeller = ggplot2::label_wrap_gen(width = 36)) +
    labs(x = "IFN\u03b3", y = "Density (square-root scale)") +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# How far each sample's cytokine-positive gate moved below its ordinary gate.
.simLowSepPlotGateShift <- function(gates) {
  gates |>
    .simLowSepFacetLabels() |>
    dplyr::mutate(
      shift = .data$gate - .data$gate_cyt,
      marker = factor(.data$marker, levels = unname(.simLowSepMarkers), labels = c("TNF", "IFN\u03b3")),
      cell_lab = factor(.simLowSepCellLab(.data$n_cell), levels = .simLowSepCellLab(sort(unique(.data$n_cell))))
    ) |>
    ggplot(aes(x = .data$cell_lab, y = .data$shift)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(
      position = position_jitter(width = 0.15, height = 0, seed = 1L),
      colour = .simLowSepGateColours[["cytpos"]], alpha = 0.7
    ) +
    facet_grid(marker ~ scenario_lab, labeller = ggplot2::labeller(scenario_lab = ggplot2::label_wrap_gen(width = 18))) +
    labs(x = NULL, y = "Ordinary gate minus cytokine-positive gate") +
    .analysis_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Paired per-sample recovery of true positives and contamination by negative
# cells, before and after the cytokine-positive gates.
.simLowSepPlotRecovery <- function(recovery) {
  recovery |>
    .simLowSepFacetLabels() |>
    tidyr::pivot_longer(
      c("sensitivity", "fdp"),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      panel = factor(
        paste0(
          ifelse(.data$marker == "TNF", "TNF", "IFN\u03b3"), ": ",
          ifelse(.data$metric == "sensitivity", "true positives recovered", "called positives that are negative")
        ),
        levels = c(
          "TNF: true positives recovered", "IFN\u03b3: true positives recovered",
          "TNF: called positives that are negative", "IFN\u03b3: called positives that are negative"
        )
      ),
      gate_type = factor(.data$gate_type, levels = names(.simLowSepGateTypes))
    ) |>
    ggplot(aes(x = .data$gate_type, y = .data$value)) +
    geom_line(aes(group = .data$sample), colour = "grey70") +
    geom_point(aes(colour = .data$gate_type)) +
    scale_colour_manual(values = .simLowSepGateColours, labels = .simLowSepGateTypes, name = NULL) +
    scale_x_discrete(labels = c(ordinary = "Ordinary", cytpos = "Cyt+")) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    .simLowSepYFloor(c(0, 0.1)) +
    facet_wrap(~ panel + scenario_lab, ncol = 4, scales = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 24, multi_line = FALSE)) +
    labs(x = NULL, y = NULL) +
    .analysis_theme() +
    theme(legend.position = "bottom")
}

# Signed error of background-subtracted frequencies (estimate minus truth),
# paired by sample, for ordinary and cytokine-positive gates.
.simLowSepPlotFrequencyError <- function(freq) {
  freq |>
    .simLowSepFacetLabels() |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = factor(.data$gate_type, levels = names(.simLowSepGateTypes)),
      error_pp = 100 * .data$error_bs
    ) |>
    ggplot(aes(x = .data$gate_type, y = .data$error_pp)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(aes(group = .data$sample), colour = "grey70") +
    geom_point(aes(colour = .data$gate_type)) +
    scale_colour_manual(values = .simLowSepGateColours, labels = .simLowSepGateTypes, name = NULL) +
    scale_x_discrete(labels = c(ordinary = "Ordinary", cytpos = "Cyt+")) +
    facet_grid(scenario_lab ~ quantity_lab, labeller = ggplot2::labeller(
      scenario_lab = ggplot2::label_wrap_gen(width = 18),
      quantity_lab = ggplot2::label_wrap_gen(width = 12)
    )) +
    labs(x = NULL, y = "Estimated minus true frequency (percentage points)") +
    .analysis_theme() +
    theme(legend.position = "bottom")
}

# `.analysis_y_floor()` where available; otherwise a blank-point floor with
# the same effect (trains the y range without censoring data).
.simLowSepYFloor <- function(range) {
  if (exists(".analysis_y_floor", mode = "function")) {
    return(.analysis_y_floor(range))
  }
  geom_blank(data = data.frame(.y_floor = range), aes(y = .data$.y_floor), inherit.aes = FALSE)
}

# Tables --------------------------------------------------------------------

.simLowSepTableDir <- function(fig_key, path_root = NULL, create = TRUE) {
  .analysis_project_dir("output", c("table", fig_key), path_root, create)
}

.simLowSepWriteTable <- function(tbl, name, fig_key, path_root = NULL) {
  path <- file.path(.simLowSepTableDir(fig_key, path_root), paste0(name, ".csv"))
  utils::write.csv(tbl, path, row.names = FALSE)
  invisible(path)
}
