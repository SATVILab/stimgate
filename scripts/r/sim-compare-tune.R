# Tuning the F-beta and Tailgate comparators (Analysis 6).
#
# Paired stimulated/unstimulated tubes are simulated as in Analysis 7 (one
# marker, two clusters, `simcyto`). Each tube is gated by Tailgate over a grid
# of derivative tolerances and biases and by F-beta over a grid of beta values.
# Only the comparators run: StimGate is not part of this analysis. Every gate
# is scored against the simulated labels with the package's `x > gate` rule.
#
# Requires sim-misc.R, sim-bandwidth.R and sim-compare-freq_bs.R.

# Scenario grid -------------------------------------------------------------

# Every transformation and separation of Analysis 7 at one response
# probability and two cell counts. IDs and seeds are assigned here, before any
# dev/quick filtering. The seeds differ from Analysis 7's, so the settings are
# chosen on data that Analysis 7 does not use.
.simTuneGrid <- function(simulation_seed, prob_response = 0.002,
                         n_cell = c(5e3, 1e5)) {
  grid <- tidyr::expand_grid(
    .simMiscGetMeanPosTbl(),
    prob_response = prob_response,
    n_cell = n_cell
  ) |>
    dplyr::mutate(transformation = factor(
      .data$transformation, levels = c("gaussian", "skew", "gamma")
    )) |>
    dplyr::arrange(.data$transformation, .data$mean_pos_setting, .data$n_cell) |>
    dplyr::mutate(transformation = as.character(.data$transformation))
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

# Fixed simulation settings (those of Analysis 7) and the comparator grids.
# Tolerances are fractions of the largest absolute first derivative of the
# density, as cytoUtils' `auto_tol` uses (1e-2), in half steps on a log10
# scale; 1e-2 itself is on the grid. The absolute tolerance 0.01 is
# cytoUtils' own default, with `auto_tol = FALSE`.
.simTuneSettings <- function() {
  list(
    sim = list(
      n_cluster = 2L,
      cov_ev = 1.5,
      prob_exact = TRUE,
      background_relative_to_response = 0.2,
      n_cell_uns_relative_to_stim = 1,
      sample_perturbation_sd = 0,
      condition_perturbation_sd = 0,
      cluster_perturbation_sd = 0
    ),
    tailgate = list(
      x = "stim",
      adjust = 1,
      method = "firstDeriv",
      tol_rel = 10^seq(-6, -1, by = 0.5),
      tol_abs = 0.01,
      bias = c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4, 0.5)
    ),
    fbeta = list(
      beta = c(0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.5, 2),
      theta = 2,
      width = 10
    ),
    fallback_margin = 0.05
  )
}

# Reference settings shown beside the selected ones.
.simTuneReferenceSettings <- function() {
  tibble::tribble(
    ~setting_ref, ~method, ~tol_type, ~tol, ~bias, ~beta, ~setting_ref_lab,
    "fbeta_published", "fbeta", NA, NA, NA, 0.8,
    "F-beta, published (beta 0.8)",
    "tailgate_cytoutils", "tailgate", "absolute", 0.01, 0, NA,
    "Tailgate, cytoUtils default (tol 0.01, no bias)",
    "tailgate_default", "tailgate", "relative", 0.01, 0, NA,
    "Tailgate, Analysis 7 default (auto tol 1%, no bias)",
    "tailgate_current", "tailgate", "relative", 0.01, 0.1, NA,
    "Tailgate, Analysis 7 current (auto tol 1%, bias 0.1)"
  )
}

.simTuneSettingCols <- c("method", "tol_type", "tol", "bias", "beta")
.simTuneScenarioCols <- c(
  "sim_id", "transformation", "mean_pos_setting", "mean_pos",
  "prob_response", "n_cell"
)

.simTuneSettingLab <- function(method, tol_type, tol, bias, beta) {
  dplyr::case_when(
    method == "fbeta" ~ paste0("F-beta, beta ", beta),
    tol_type == "absolute" ~ paste0("Tailgate, tol ", tol, ", bias ", bias),
    TRUE ~ paste0("Tailgate, tol ", .simTuneTolText(tol), " of max slope, bias ", bias)
  )
}

# A relative tolerance as a power of ten, e.g. 1e-2 or 3.2e-5.
.simTuneTolText <- function(tol) {
  out <- formatC(tol, format = "e", digits = 1)
  sub("[.]0e", "e", sub("e([+-])0?", "e\\1", out))
}

.simTuneScenarioLab <- function(tbl) {
  paste0(
    tbl$transformation, ", ", tbl$mean_pos_setting, " separation (",
    tbl$mean_pos, ")"
  )
}

.simTuneCellLab <- function(n_cell) {
  paste0(format(n_cell, big.mark = ",", scientific = FALSE, trim = TRUE), " cells")
}

# Simulation ----------------------------------------------------------------

# One dataset of `n_sample` unstimulated/stimulated pairs, with the same
# simulator call as `.simCompareFreqBs()`.
.simTuneSimulate <- function(row, n_sample, settings = .simTuneSettings()) {
  s <- settings$sim
  n_cell_uns <- round(row$n_cell * s$n_cell_uns_relative_to_stim)
  prob_uns <- row$prob_response * s$background_relative_to_response
  .simCompareSimCytExperiment(
    nSample = n_sample,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = s$n_cluster,
    nCellByCondition = c(n_cell_uns, row$n_cell),
    transformationFunc = .simCompareGetTrans(row$transformation),
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(c(0, row$mean_pos), byrow = TRUE, ncol = 1),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(1 - prob_uns, prob_uns),
    probExact = s$prob_exact,
    probResponseVecByStimCondition = list(c(-row$prob_response, row$prob_response)),
    conditionPerturbationSd = s$condition_perturbation_sd,
    clusterPerturbationSd = s$cluster_perturbation_sd,
    samplePerturbationSd = s$sample_perturbation_sd,
    covEvMin = s$cov_ev,
    covEvMax = s$cov_ev
  )
}

# Gates ---------------------------------------------------------------------

# Tailgate cutpoints (before bias) for every tolerance. The bandwidth is
# estimated once with ks::hpi(), as `.simCompareTailgateThreshold()` does, and
# a relative tolerance is that fraction of the largest absolute first
# derivative, so a fraction of 0.01 reproduces `autoTol = TRUE` exactly.
.simTuneTailgateCuts <- function(x, settings = .simTuneSettings()) {
  tg <- settings$tailgate
  tols <- tibble::tibble(
    tol_type = c(rep("relative", length(tg$tol_rel)), "absolute"),
    tol = c(tg$tol_rel, tg$tol_abs)
  )
  setup_error <- NA_character_
  setup <- tryCatch({
    bw <- suppressWarnings(ks::hpi(x, deriv.order = 1L))
    deriv <- cytoUtils:::.deriv_density(
      x = x, deriv = 1, bandwidth = bw, adjust = tg$adjust
    )
    list(bw = bw, max_deriv = max(abs(deriv$y)))
  }, error = function(e) {
    setup_error <<- conditionMessage(e)
    NULL
  })
  purrr::pmap_dfr(tols, function(tol_type, tol) {
    if (is.null(setup)) {
      return(tibble::tibble(
        tol_type = tol_type, tol = tol, tol_value = NA_real_,
        cut = NA_real_, error = setup_error
      ))
    }
    tol_value <- if (tol_type == "relative") tol * setup$max_deriv else tol
    err <- NA_character_
    cut <- tryCatch(
      .simCompareTailgateThreshold(
        x = x, adjust = tg$adjust, bandwidth = setup$bw, method = tg$method,
        tol = tol_value, autoTol = FALSE, bias = 0
      )$threshold,
      error = function(e) {
        err <<- conditionMessage(e)
        NA_real_
      }
    )
    tibble::tibble(
      tol_type = tol_type, tol = tol, tol_value = tol_value,
      cut = cut, error = err
    )
  })
}

# Tailgate gates for every tolerance and bias.
.simTuneTailgateGates <- function(x, settings = .simTuneSettings()) {
  tidyr::expand_grid(
    .simTuneTailgateCuts(x, settings),
    bias = settings$tailgate$bias
  ) |>
    dplyr::mutate(
      method = "tailgate", beta = NA_real_, threshold = .data$cut + .data$bias
    )
}

.simTuneFbetaGates <- function(x_uns, x_stim, fbeta_env,
                               settings = .simTuneSettings()) {
  fb <- settings$fbeta
  purrr::map_dfr(fb$beta, function(beta) {
    err <- NA_character_
    threshold <- tryCatch(
      .simCompareFbetaThreshold(
        xUns = x_uns, xStim = x_stim, fbetaEnv = fbeta_env,
        beta = beta, theta = fb$theta, width = fb$width
      )$threshold,
      error = function(e) {
        err <<- conditionMessage(e)
        NA_real_
      }
    )
    tibble::tibble(
      method = "fbeta", tol_type = NA_character_, tol = NA_real_,
      tol_value = NA_real_, bias = NA_real_, beta = beta,
      cut = threshold, threshold = threshold, error = err
    )
  })
}

# Score each gate as Analysis 7 does: a comparator that finds no cutpoint
# without an error gets the high-value fallback gate; errors stay missing.
.simTuneScoreGates <- function(gates, x_stim, x_uns, labels_stim,
                               settings = .simTuneSettings()) {
  purrr::pmap_dfr(gates, function(...) {
    g <- list(...)
    est <- .simCompareEstimateFromThreshold(
      xStim = x_stim, xUns = x_uns, threshold = g$threshold,
      fallbackHighValue = is.na(g$error),
      fallbackMargin = settings$fallback_margin,
      labelsStim = labels_stim
    )
    tibble::as_tibble(c(g[setdiff(names(g), "threshold")], est))
  })
}

# Histogram counts of one sample, for drawing its gates without keeping cells.
.simTuneHist <- function(x_stim, x_uns, labels_stim, n_bins = 120L) {
  rng <- range(c(x_stim, x_uns))
  breaks <- seq(rng[[1]], rng[[2]], length.out = n_bins + 1L)
  count <- function(x) graphics::hist(x, breaks = breaks, plot = FALSE)$counts
  gp <- labels_stim %in% "gp"
  tibble::tibble(
    x_lo = utils::head(breaks, -1L),
    x_hi = breaks[-1L],
    stim_negative = count(x_stim[!gp]),
    stim_positive = count(x_stim[gp]),
    unstim = count(x_uns)
  )
}

# Run -----------------------------------------------------------------------

# Dataset seeds of one scenario, drawn from its `sim_seed` as
# `.simCompareFreqBs()` draws its iteration seeds.
.simTuneIterSeeds <- function(sim_seed, n_iter) {
  .analysis_with_seed(sim_seed, {
    sample.int(.Machine$integer.max, n_iter, replace = TRUE)
  })
}

# Prefix rows with the scenario columns of `row`.
.simTuneWithScenario <- function(row, tbl) {
  if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
  scen <- dplyr::select(row, dplyr::all_of(.simTuneScenarioCols))
  dplyr::bind_cols(scen[rep(1L, nrow(tbl)), ], tbl)
}

# One dataset of one scenario, seeded by its own `iter_seed`, so the result
# does not depend on which process runs it or in what order. Histograms are
# kept for the first `n_hist_sample` samples of the first dataset.
.simTuneRunDataset <- function(row, iter, iter_seed, n_sample, fbeta_env,
                               settings = .simTuneSettings(),
                               n_hist_sample = 4L) {
  if (nrow(row) != 1L) stop("row must have exactly one scenario.")
  .analysis_with_seed(iter_seed, {
    sim <- .simTuneSimulate(row, n_sample, settings)
    truth <- .simCompareTruthTable(sim$labelsList, n_sample, 2L)
    per_sample <- lapply(seq_len(n_sample), function(s) {
      x_uns <- as.numeric(flowCore::exprs(sim$flowFrameList[[2L * s - 1L]])[, "F1"])
      x_stim <- as.numeric(flowCore::exprs(sim$flowFrameList[[2L * s]])[, "F1"])
      labels <- sim$labelsList[[2L * s]]
      gates <- dplyr::bind_rows(
        .simTuneTailgateGates(x_stim, settings),
        .simTuneFbetaGates(x_uns, x_stim, fbeta_env, settings)
      )
      scores <- .simTuneScoreGates(gates, x_stim, x_uns, labels, settings) |>
        dplyr::mutate(iter = iter, sample = as.character(s), .before = 1L)
      hist <- if (iter == 1L && s <= n_hist_sample) {
        dplyr::mutate(.simTuneHist(x_stim, x_uns, labels),
          iter = iter, sample = as.character(s), .before = 1L)
      }
      list(scores = scores, hist = hist)
    })
    scores <- purrr::map_dfr(per_sample, "scores") |>
      dplyr::left_join(
        dplyr::select(truth, "sample", "propRespTruth", "propStimTruth"),
        by = "sample"
      )
    list(
      scores = .simTuneWithScenario(row, scores),
      hist = .simTuneWithScenario(row, purrr::map_dfr(per_sample, "hist"))
    )
  })
}

# Every dataset of one scenario, in this process.
.simTuneRunScenario <- function(row, n_sample, n_iter, fbeta_env,
                                settings = .simTuneSettings(),
                                n_hist_sample = 4L) {
  if (nrow(row) != 1L) stop("row must have exactly one scenario.")
  iter_seeds <- .simTuneIterSeeds(row$sim_seed, n_iter)
  out <- lapply(seq_len(n_iter), function(iter) {
    .simTuneRunDataset(row, iter, iter_seeds[[iter]], n_sample, fbeta_env,
      settings, n_hist_sample)
  })
  list(scores = purrr::map_dfr(out, "scores"), hist = purrr::map_dfr(out, "hist"))
}

# Helper files a worker sources, in dependency order.
.simTuneHelperFiles <- c(
  "analysis-runtime.R", "sim-misc.R", "sim-bandwidth.R",
  "sim-compare-freq_bs.R", "sim-compare-tune.R"
)

# Every dataset of every scenario. With `workers > 1`, datasets run in
# parallel in multisession workers, largest scenarios first. Each worker loads
# the current checkout, sources the helpers and creates its own F-beta Python
# environment, which is never sent between processes. Seeds are fixed per
# dataset, so serial and parallel runs give the same results.
.simTuneRunGrid <- function(grid, n_sample, n_iter, path_fbeta,
                            settings = .simTuneSettings(), n_hist_sample = 4L,
                            workers = 1L, root_dir = NULL) {
  tasks <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    tibble::tibble(
      row_index = i, n_cell = grid$n_cell[[i]], iter = seq_len(n_iter),
      iter_seed = .simTuneIterSeeds(grid$sim_seed[[i]], n_iter)
    )
  }) |>
    dplyr::arrange(dplyr::desc(.data$n_cell), .data$row_index, .data$iter)
  workers <- max(1L, as.integer(workers))

  if (workers == 1L) {
    # The reticulate environment is created in this process and kept local.
    fbeta_env <- .simCompareFbetaEnvironment(pathFbeta = path_fbeta)
    out <- lapply(seq_len(nrow(tasks)), function(k) {
      t <- tasks[k, ]
      message("Comparator tuning: scenario ", grid$sim_id[[t$row_index]],
        ", dataset ", t$iter, " (", k, " of ", nrow(tasks), ")")
      .simTuneRunDataset(grid[t$row_index, ], t$iter, t$iter_seed, n_sample,
        fbeta_env, settings, n_hist_sample)
    })
  } else {
    if (is.null(root_dir)) stop("root_dir is required for parallel runs.")
    root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
    path_fbeta <- normalizePath(path_fbeta, winslash = "/", mustWork = TRUE)
    helper_paths <- file.path(root_dir, "scripts", "r", .simTuneHelperFiles)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    message("Comparator tuning: ", nrow(tasks), " datasets on ", workers, " workers")
    out <- furrr::future_map(seq_len(nrow(tasks)), function(k) {
      t <- tasks[k, ]
      env <- new.env(parent = globalenv())
      for (f in helper_paths) source(f, local = env)
      env$.simCompareEnsureCurrentCheckout(root_dir)
      fbeta_env <- env$.simCompareFbetaEnvironment(pathFbeta = path_fbeta)
      env$.simTuneRunDataset(grid[t$row_index, ], t$iter, t$iter_seed, n_sample,
        fbeta_env, settings, n_hist_sample)
    }, .options = furrr::furrr_options(
      seed = TRUE, scheduling = Inf, globals = c(
        "tasks", "grid", "helper_paths", "root_dir", "path_fbeta",
        "n_sample", "settings", "n_hist_sample"
      )
    ))
  }
  order_tbl <- function(x) {
    dplyr::arrange(x, .data$sim_id, .data$iter, as.integer(.data$sample))
  }
  list(
    scores = order_tbl(purrr::map_dfr(out, "scores")),
    hist = order_tbl(purrr::map_dfr(out, "hist"))
  )
}

# Cache ---------------------------------------------------------------------

.simTuneCacheSettings <- function(grid, n_sample, n_iter, settings,
                                  simulation_seed, profile) {
  list(
    analysis_semantics_version = "sim-tune-comparators-v1",
    grid = as.data.frame(grid),
    n_sample = as.integer(n_sample),
    n_iter = as.integer(n_iter),
    settings = settings,
    simulation_seed = as.integer(simulation_seed),
    profile = profile
  )
}

.simTuneWriteCache <- function(results, settings, path) {
  .write_rds_atomic(list(
    settings = settings, results = results,
    run_id = .sanitize_run_id(Sys.getenv("ANALYSIS_RUN_ID", unset = ""))
  ), path)
}

.simTuneReadCache <- function(path, settings, analysis_key) {
  qmd <- "analysis/6-sim-tune-comparators.qmd"
  cached <- .analysis_read_rds(path, analysis_key, qmd)
  if (!is.list(cached) || !all(c("settings", "results") %in% names(cached))) {
    .analysis_cache_error(analysis_key, paste0("Cached results at ", path, " are malformed."), qmd)
  }
  .analysis_check_expected_run(cached, analysis_key, qmd)
  if (!isTRUE(all.equal(cached$settings, settings))) {
    .analysis_cache_error(analysis_key, paste0(
      "Cached results at ", path, " were made with different settings ",
      "(scenario grid, dataset size, comparator grids or seed)."
    ), qmd)
  }
  cached$results
}

# Summaries -----------------------------------------------------------------

.simTuneQuantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  unname(stats::quantile(x, p, names = FALSE))
}

# F1 per stimulated tube. Errors leave the counts and F1 missing; F1 is
# undefined only when a tube has no genuine positives and nothing is selected.
.simTuneMetrics <- function(scores) {
  failed <- !is.na(scores$error)
  scores |>
    dplyr::mutate(
      dplyr::across(
        c("nTruePos", "nFalsePos", "nFalseNeg", "nTrueNeg", "propRespEst"),
        ~ dplyr::if_else(failed, NA, .x)
      ),
      f1_den = 2 * .data$nTruePos + .data$nFalsePos + .data$nFalseNeg,
      f1 = dplyr::if_else(.data$f1_den > 0, 2 * .data$nTruePos / .data$f1_den, NA_real_),
      est_ratio = .data$propRespEst / .data$propRespTruth,
      setting_lab = .simTuneSettingLab(
        .data$method, .data$tol_type, .data$tol, .data$bias, .data$beta
      )
    ) |>
    dplyr::select(-"f1_den")
}

# Tube-level distributions pooled over the datasets of each scenario.
.simTuneSummary <- function(metrics) {
  metrics |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(
      .simTuneScenarioCols, .simTuneSettingCols, "setting_lab"
    )))) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_datasets = dplyr::n_distinct(.data$iter),
      n_error = sum(!is.na(.data$error)),
      n_fallback = sum(.data$thresholdFallbackUsed %in% TRUE),
      n_f1_finite = sum(is.finite(.data$f1)),
      f1_p10 = .simTuneQuantile(.data$f1, 0.1),
      f1_median = .simTuneQuantile(.data$f1, 0.5),
      f1_p90 = .simTuneQuantile(.data$f1, 0.9),
      n_est_finite = sum(is.finite(.data$propRespEst)),
      truth = mean(.data$propRespTruth),
      est_p10 = .simTuneQuantile(.data$propRespEst, 0.1),
      est_median = .simTuneQuantile(.data$propRespEst, 0.5),
      est_p90 = .simTuneQuantile(.data$propRespEst, 0.9),
      gate_median = .simTuneQuantile(.data$threshold, 0.5),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # How far the central 80% of estimates strays from the truth: the larger
      # of the relative errors of the 10th and 90th percentiles.
      est_band_error = pmax(
        abs(.data$est_p10 / .data$truth - 1), abs(.data$est_p90 / .data$truth - 1)
      ),
      truth_in_band = .data$est_p10 <= .data$truth & .data$truth <= .data$est_p90
    )
}

# Equal-weight averages over scenarios for each setting, and the selection:
# keep settings whose mean median F1 is within `f1_tolerance` of the best for
# that method, then take the smallest mean estimate-band error, then the
# highest mean 10th-percentile F1.
.simTuneRank <- function(summary, f1_tolerance = 0.02) {
  summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(.simTuneSettingCols, "setting_lab")))) |>
    dplyr::summarise(
      n_scenarios = dplyr::n(),
      n_scenarios_f1 = sum(is.finite(.data$f1_median)),
      n_error = sum(.data$n_error),
      n_fallback = sum(.data$n_fallback),
      mean_f1_median = mean(.data$f1_median),
      min_f1_median = min(.data$f1_median),
      mean_f1_p10 = mean(.data$f1_p10),
      mean_f1_p90 = mean(.data$f1_p90),
      mean_est_band_error = mean(.data$est_band_error),
      max_est_band_error = max(.data$est_band_error),
      n_truth_in_band = sum(.data$truth_in_band %in% TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$method) |>
    dplyr::mutate(
      f1_eligible = is.finite(.data$mean_f1_median) &
        .data$mean_f1_median >= max(.data$mean_f1_median, na.rm = TRUE) - f1_tolerance
    ) |>
    dplyr::arrange(
      .data$method, dplyr::desc(.data$f1_eligible), .data$mean_est_band_error,
      dplyr::desc(.data$mean_f1_p10), .by_group = TRUE
    ) |>
    dplyr::mutate(
      selected = dplyr::row_number() == 1L & .data$f1_eligible,
      rank = dplyr::row_number()
    ) |>
    dplyr::ungroup()
}

# The selected setting of each method beside the reference settings.
.simTuneCompareSettings <- function(rank) {
  ref <- .simTuneReferenceSettings()
  sel <- rank |>
    dplyr::filter(.data$selected) |>
    dplyr::transmute(
      setting_ref = paste0(.data$method, "_selected"),
      method = .data$method, tol_type = .data$tol_type, tol = .data$tol,
      bias = .data$bias, beta = .data$beta,
      setting_ref_lab = paste0(
        ifelse(.data$method == "fbeta", "F-beta", "Tailgate"), ", selected (",
        sub("^[^,]*, ", "", .data$setting_lab), ")"
      )
    )
  out <- dplyr::bind_rows(ref, sel)
  levels <- c(
    "fbeta_published", "fbeta_selected", "tailgate_cytoutils",
    "tailgate_default", "tailgate_current", "tailgate_selected"
  )
  out$setting_ref <- factor(out$setting_ref, levels = levels)
  dplyr::arrange(out, .data$setting_ref)
}

.simTuneSettingColours <- c(
  fbeta_published = "#56B4E9", fbeta_selected = "#0072B2",
  tailgate_cytoutils = "#999999", tailgate_default = "#CC79A7",
  tailgate_current = "#009E73", tailgate_selected = "#D55E00"
)

.simTuneSettingLinetypes <- c(
  fbeta_published = "solid", fbeta_selected = "22",
  tailgate_cytoutils = "solid", tailgate_default = "solid",
  tailgate_current = "solid", tailgate_selected = "22"
)

# Join tube or summary rows to the compared settings. A selected setting can
# equal a reference one, in which case both rows are kept.
.simTuneJoinSettings <- function(tbl, compare) {
  key <- function(x) {
    paste(x$method, x$tol_type, signif(x$tol, 6), signif(x$bias, 6), signif(x$beta, 6))
  }
  compare$.key <- key(compare)
  tbl$.key <- key(tbl)
  dplyr::inner_join(
    tbl,
    dplyr::select(compare, ".key", "setting_ref", "setting_ref_lab"),
    by = ".key", relationship = "many-to-many"
  ) |>
    dplyr::select(-".key")
}

# Paired comparison within datasets: the median F1 and the median estimate of
# each dataset under the selected setting minus under a reference setting.
.simTuneDatasetDifferences <- function(metrics, compare) {
  per_dataset <- .simTuneJoinSettings(metrics, compare) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(.simTuneScenarioCols, "iter", "method", "setting_ref")))) |>
    dplyr::summarise(
      f1_median = .simTuneQuantile(.data$f1, 0.5),
      est_median = .simTuneQuantile(.data$propRespEst, 0.5),
      .groups = "drop"
    )
  selected <- dplyr::filter(per_dataset, grepl("_selected$", .data$setting_ref)) |>
    dplyr::select(-"setting_ref") |>
    dplyr::rename(f1_selected = "f1_median", est_selected = "est_median")
  per_dataset |>
    dplyr::filter(!grepl("_selected$", .data$setting_ref)) |>
    dplyr::inner_join(selected, by = c(.simTuneScenarioCols, "iter", "method")) |>
    dplyr::mutate(
      f1_diff = .data$f1_selected - .data$f1_median,
      est_diff = .data$est_selected - .data$est_median
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(.simTuneScenarioCols, "method", "setting_ref")))) |>
    dplyr::summarise(
      n_datasets = sum(is.finite(.data$f1_diff)),
      f1_diff_mean = mean(.data$f1_diff, na.rm = TRUE),
      f1_diff_se = stats::sd(.data$f1_diff, na.rm = TRUE) / sqrt(.data$n_datasets),
      est_diff_mean = mean(.data$est_diff, na.rm = TRUE),
      est_diff_se = stats::sd(.data$est_diff, na.rm = TRUE) / sqrt(.data$n_datasets),
      .groups = "drop"
    )
}

# Plots ---------------------------------------------------------------------

.simTuneFacet <- function() {
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$scenario_lab), cols = ggplot2::vars(.data$cell_lab),
    labeller = ggplot2::label_wrap_gen(width = 18)
  )
}

.simTuneAddLabs <- function(tbl) {
  tbl$scenario_lab <- factor(.simTuneScenarioLab(tbl), levels = unique(.simTuneScenarioLab(
    dplyr::arrange(tbl, factor(.data$transformation, c("gaussian", "skew", "gamma")), .data$mean_pos_setting)
  )))
  tbl$cell_lab <- factor(.simTuneCellLab(tbl$n_cell), levels = .simTuneCellLab(sort(unique(tbl$n_cell))))
  tbl
}

.simTuneTolLab <- function(tol_type, tol) {
  ifelse(tol_type == "absolute", paste0("abs ", tol), .simTuneTolText(tol))
}

# Tailgate tolerance-by-bias heat map of one summary statistic.
.simTunePlotTailgateHeat <- function(summary, stat, fill_lab, selected = NULL,
                                     direction = 1) {
  tbl <- summary |>
    dplyr::filter(.data$method == "tailgate") |>
    .simTuneAddLabs() |>
    dplyr::mutate(
      tol_lab = factor(
        .simTuneTolLab(.data$tol_type, .data$tol),
        levels = unique(.simTuneTolLab(.data$tol_type, .data$tol)[order(.data$tol_type == "absolute", .data$tol)])
      ),
      bias_lab = factor(.data$bias),
      value = .data[[stat]]
    )
  p <- ggplot(tbl, aes(x = .data$tol_lab, y = .data$bias_lab, fill = .data$value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", .data$value)), size = 2) +
    scale_fill_viridis_c(name = fill_lab, direction = direction, na.value = "grey85") +
    .simTuneFacet() +
    labs(x = "Derivative tolerance (fraction of the steepest slope, or absolute)", y = "Bias") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(selected)) {
    sel <- selected |>
      dplyr::mutate(
        tol_lab = factor(.simTuneTolLab(.data$tol_type, .data$tol), levels = levels(tbl$tol_lab)),
        bias_lab = factor(.data$bias, levels = levels(tbl$bias_lab))
      )
    p <- p + geom_tile(
      data = sel, aes(x = .data$tol_lab, y = .data$bias_lab),
      inherit.aes = FALSE, fill = NA, colour = "red", linewidth = 0.8
    )
  }
  p
}

# Tailgate statistic against bias, one line per tolerance, with the 10th-90th
# percentile ribbon only for the tolerances in `tol_show`.
.simTunePlotTailgateLines <- function(summary, which = c("f1", "est")) {
  which <- match.arg(which)
  tbl <- summary |>
    dplyr::filter(.data$method == "tailgate") |>
    .simTuneAddLabs() |>
    dplyr::mutate(tol_lab = factor(
      .simTuneTolLab(.data$tol_type, .data$tol),
      levels = unique(.simTuneTolLab(.data$tol_type, .data$tol)[order(.data$tol_type == "absolute", .data$tol)])
    ))
  cols <- if (which == "f1") c("f1_p10", "f1_median", "f1_p90") else c("est_p10", "est_median", "est_p90")
  tbl$lo <- tbl[[cols[1]]]
  tbl$mid <- tbl[[cols[2]]]
  tbl$hi <- tbl[[cols[3]]]
  if (which == "est") {
    tbl <- dplyr::mutate(tbl, dplyr::across(c("lo", "mid", "hi"), ~ .x / .data$truth))
  }
  p <- ggplot(tbl, aes(x = .data$bias, y = .data$mid, colour = .data$tol_lab, fill = .data$tol_lab)) +
    geom_ribbon(aes(ymin = .data$lo, ymax = .data$hi), alpha = 0.08, colour = NA) +
    geom_line() +
    geom_point(size = 0.8) +
    scale_colour_viridis_d(name = "Tolerance", end = 0.9, aesthetics = c("colour", "fill")) +
    .simTuneFacet() +
    labs(x = "Bias") +
    theme_bw()
  if (which == "f1") {
    p + labs(y = "F1 (median; band: 10th-90th percentile)") + coord_cartesian(ylim = c(0, 1))
  } else {
    p + geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Estimate / truth (median; band: 10th-90th percentile)") +
      .analysis_y_floor(c(0, 2))
  }
}

# Tailgate statistic against log10 relative tolerance, one line per bias,
# with 10th-90th percentile ribbons.
.simTunePlotTailgateTol <- function(summary, which = c("f1", "est")) {
  which <- match.arg(which)
  tbl <- summary |>
    dplyr::filter(.data$method == "tailgate", .data$tol_type == "relative") |>
    .simTuneAddLabs() |>
    dplyr::mutate(bias_lab = factor(.data$bias))
  cols <- if (which == "f1") c("f1_p10", "f1_median", "f1_p90") else c("est_p10", "est_median", "est_p90")
  tbl$lo <- tbl[[cols[1]]]
  tbl$mid <- tbl[[cols[2]]]
  tbl$hi <- tbl[[cols[3]]]
  if (which == "est") {
    tbl <- dplyr::mutate(tbl, dplyr::across(c("lo", "mid", "hi"), ~ .x / .data$truth))
  }
  p <- ggplot(tbl, aes(x = .data$tol, y = .data$mid, colour = .data$bias_lab, fill = .data$bias_lab)) +
    geom_vline(xintercept = 0.01, linetype = "dotted") +
    geom_ribbon(aes(ymin = .data$lo, ymax = .data$hi), alpha = 0.06, colour = NA) +
    geom_line() +
    geom_point(size = 0.8) +
    scale_x_log10(breaks = 10^(-6:-1), labels = .simTuneTolText) +
    scale_colour_viridis_d(name = "Bias", end = 0.9, aesthetics = c("colour", "fill")) +
    .simTuneFacet() +
    labs(x = "Tolerance (fraction of the steepest slope, log10 scale)") +
    theme_bw()
  if (which == "f1") {
    p + labs(y = "F1 (median; band: 10th-90th percentile)") + coord_cartesian(ylim = c(0, 1))
  } else {
    p + geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Estimate / truth (median; band: 10th-90th percentile)") +
      .analysis_y_floor(c(0, 2))
  }
}

# F-beta statistic against beta: median with 10th-90th percentile bars.
.simTunePlotFbeta <- function(summary, which = c("f1", "est")) {
  which <- match.arg(which)
  tbl <- summary |> dplyr::filter(.data$method == "fbeta") |> .simTuneAddLabs()
  cols <- if (which == "f1") c("f1_p10", "f1_median", "f1_p90") else c("est_p10", "est_median", "est_p90")
  tbl$lo <- tbl[[cols[1]]]
  tbl$mid <- tbl[[cols[2]]]
  tbl$hi <- tbl[[cols[3]]]
  if (which == "est") {
    tbl <- dplyr::mutate(tbl, dplyr::across(c("lo", "mid", "hi"), ~ .x / .data$truth))
  }
  p <- ggplot(tbl, aes(x = .data$beta, y = .data$mid)) +
    geom_vline(xintercept = 0.8, colour = .simTuneSettingColours[["fbeta_published"]], linetype = "dotted") +
    geom_errorbar(aes(ymin = .data$lo, ymax = .data$hi), width = 0, colour = "#0072B2") +
    geom_line(colour = "#0072B2") +
    geom_point(colour = "#0072B2") +
    scale_x_log10(breaks = sort(unique(tbl$beta))) +
    .simTuneFacet() +
    labs(x = "beta (log scale; dotted: published 0.8)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (which == "f1") {
    p + labs(y = "F1 (median; bars: 10th-90th percentile)") + coord_cartesian(ylim = c(0, 1))
  } else {
    p + geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Estimate / truth (median; bars: 10th-90th percentile)") +
      .analysis_y_floor(c(0, 2))
  }
}

# Selected against reference settings: median and 10th-90th percentiles.
.simTunePlotCompare <- function(summary, compare, which = c("f1", "est")) {
  which <- match.arg(which)
  tbl <- .simTuneJoinSettings(summary, compare) |> .simTuneAddLabs()
  cols <- if (which == "f1") c("f1_p10", "f1_median", "f1_p90") else c("est_p10", "est_median", "est_p90")
  tbl$lo <- tbl[[cols[1]]]
  tbl$mid <- tbl[[cols[2]]]
  tbl$hi <- tbl[[cols[3]]]
  if (which == "est") {
    tbl <- dplyr::mutate(tbl, dplyr::across(c("lo", "mid", "hi"), ~ .x / .data$truth))
  }
  labs_vec <- stats::setNames(compare$setting_ref_lab, as.character(compare$setting_ref))
  p <- ggplot(tbl, aes(x = .data$setting_ref, y = .data$mid, colour = .data$setting_ref)) +
    geom_pointrange(aes(ymin = .data$lo, ymax = .data$hi), size = 0.3) +
    scale_colour_manual(values = .simTuneSettingColours, labels = labs_vec, name = NULL) +
    scale_x_discrete(labels = NULL) +
    .simTuneFacet() +
    labs(x = NULL) +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "vertical",
      axis.ticks.x = element_blank())
  if (which == "f1") {
    p + labs(y = "F1 (median and 10th-90th percentile)") + coord_cartesian(ylim = c(0, 1))
  } else {
    p + geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Estimate / truth (median and 10th-90th percentile)") +
      .analysis_y_floor(c(0, 2))
  }
}

# Gates of every tube under the compared settings.
.simTunePlotGateSpread <- function(metrics, compare) {
  tbl <- .simTuneJoinSettings(metrics, compare) |> .simTuneAddLabs()
  labs_vec <- stats::setNames(compare$setting_ref_lab, as.character(compare$setting_ref))
  ggplot(tbl, aes(x = .data$setting_ref, y = .data$threshold, colour = .data$setting_ref)) +
    geom_boxplot(outlier.size = 0.4, fill = NA) +
    scale_colour_manual(values = .simTuneSettingColours, labels = labs_vec, name = NULL) +
    scale_x_discrete(labels = NULL) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$scenario_lab, .data$cell_lab), scales = "free_y", ncol = 2,
      labeller = ggplot2::label_wrap_gen(width = 30, multi_line = FALSE)
    ) +
    labs(x = NULL, y = "Gate") +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "vertical", axis.ticks.x = element_blank())
}

# Gates of the compared settings drawn on the samples kept as histograms
# (one scenario and cell count). Stimulated cells are stacked by simulated
# label; the unstimulated tube is a step outline. Counts use a pseudo-log
# scale so the few responding cells stay visible.
.simTunePlotGates <- function(hist, metrics, compare) {
  gates <- .simTuneJoinSettings(
    dplyr::semi_join(metrics, hist, by = c("sim_id", "iter", "sample")), compare
  )
  labs_vec <- stats::setNames(compare$setting_ref_lab, as.character(compare$setting_ref))
  # Stack responders on top of the negative cells, so the bar height is the
  # stimulated tube's count.
  bars <- dplyr::bind_rows(
    dplyr::mutate(hist, label = "Stimulated, simulated negative",
      ymin = 0, ymax = .data$stim_negative),
    dplyr::mutate(hist, label = "Stimulated, simulated responder",
      ymin = .data$stim_negative, ymax = .data$stim_negative + .data$stim_positive)
  ) |>
    dplyr::filter(.data$ymax > .data$ymin) |>
    dplyr::mutate(label = factor(.data$label, c(
      "Stimulated, simulated responder", "Stimulated, simulated negative"
    )))
  uns <- dplyr::mutate(hist, x = (.data$x_lo + .data$x_hi) / 2)
  sample_lab <- function(x) paste("Sample", x)
  ggplot() +
    geom_rect(
      data = bars,
      aes(xmin = .data$x_lo, xmax = .data$x_hi, ymin = .data$ymin, ymax = .data$ymax, fill = .data$label),
      alpha = 0.7
    ) +
    geom_step(data = uns, aes(x = .data$x_lo, y = .data$unstim), colour = "grey30", linewidth = 0.3) +
    geom_vline(
      data = gates,
      aes(xintercept = .data$threshold, colour = .data$setting_ref,
        linetype = .data$setting_ref, group = .data$setting_ref),
      linewidth = 0.6, alpha = 0.9
    ) +
    # Selected settings are dashed, so they stay visible where they give the
    # same gate as a reference setting.
    scale_linetype_manual(values = .simTuneSettingLinetypes, labels = labs_vec, name = NULL) +
    scale_fill_manual(values = c("#E69F00", "grey70"), name = NULL) +
    scale_colour_manual(values = .simTuneSettingColours, labels = labs_vec, name = NULL) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1e3, 1e4)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$sample), ncol = 2,
      labeller = ggplot2::as_labeller(sample_lab)) +
    labs(x = "Expression", y = "Cells (pseudo-log scale; outline: unstimulated)") +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "vertical")
}
