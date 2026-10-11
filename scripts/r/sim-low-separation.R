# Low-separation simulations of the cytokine-positive gates (Analyses 11a
# and 11b).
#
# Two markers (TNF, IFNg) are simulated with simcyto on the Gaussian
# transformation. Each simulation row is one dataset (`iter`) of `n_sample`
# unstimulated/stimulated pairs, gated once with gateStim(); four rules are
# then applied to the same cells: the ordinary gates, the package's antimode
# cytokine-positive gates, and the co-expression gates of
# `coexpression-gates.R` with and without the double-positive screen. Truth
# comes from the simulated component labels. Analysis 11b adds a
# stimulated-tube batch mismatch to the all-negative component.
#
# Rows run through the shared runner of sim-bandwidth-analysis-run.R
# (`.simBandwidthRunRow()` etc.), so source analysis-runtime.R, sim-misc.R,
# sim-bandwidth.R, sim-bandwidth-analysis-io.R, sim-bandwidth-analysis-run.R
# and coexpression-gates.R first.

.simLowSepMarkers <- c(tnf = "TNF", ifng = "IFNg")
.simLowSepChnl <- c(tnf = "F1", ifng = "F2")
.simLowSepCombn <- c("TNF-IFNg-", "TNF+IFNg-", "TNF-IFNg+", "TNF+IFNg+")

# The four rules compared, with one label, colour and short label each.
.simLowSepGateTypes <- c(
  ordinary = "Ordinary gates",
  antimode = "Antimode gates",
  coex_screened = "Co-expression gates (screened)",
  coex_unscreened = "Co-expression gates (unscreened)"
)
.simLowSepGateColours <- c(
  ordinary = "#7F7F7F", antimode = "#D55E00",
  coex_screened = "#0072B2", coex_unscreened = "#56B4E9"
)
.simLowSepGateShort <- c(
  ordinary = "Ordinary", antimode = "Antimode",
  coex_screened = "Co-exp.\nscreened", coex_unscreened = "Co-exp.\nunscreened"
)
# Co-expression rules and their `.coexLowerGate(screen = )` value.
.simLowSepCoexScreen <- c(coex_screened = TRUE, coex_unscreened = FALSE)

# Scenario columns carried into every result row. 11a has no mismatch
# (mismatch_type "none"); in 11b, `base_scenario_id` identifies the
# biological scenario whose data are shared across mismatch settings.
.simLowSepScenCols <- c(
  "scenario_id", "base_scenario_id", "family", "separation", "mean_pos", "mean_pos_ifng",
  "response_level", "prob_response", "response_split", "background_relative",
  "background_split", "n_cell", "mismatch_type", "mismatch_label", "stim_mean_shift",
  "stim_sd_multiplier"
)
# Columns identifying one simulation row (one dataset of one scenario).
.simLowSepRowCols <- c("sim_id", "sim_seed", "iter")

# Largest number of datasets per scenario. Fixed, so a row keeps its sim_id
# and seed whatever the run size (final runs use all of them).
.simLowSepIterMax <- 20L

# Share of the responding cells in each positive combination (TNF+IFNg-,
# TNF-IFNg+, TNF+IFNg+), and of the unstimulated tube's background cells.
.simLowSepSplits <- list(
  default = c("TNF+IFNg-" = 0.25, "TNF-IFNg+" = 0.15, "TNF+IFNg+" = 0.6),
  single = c("TNF+IFNg-" = 0.5, "TNF-IFNg+" = 0.5, "TNF+IFNg+" = 0),
  double = c("TNF+IFNg-" = 0, "TNF-IFNg+" = 0, "TNF+IFNg+" = 1),
  weak_double = c("TNF+IFNg-" = 0.47, "TNF-IFNg+" = 0.47, "TNF+IFNg+" = 0.06)
)

# Scenario families, in display order.
.simLowSepFamilies <- c(
  main = "Main grid",
  high_separation = "High separation",
  no_coexpression = "Responders single-positive only",
  weak_double = "Weak double-positive response",
  double_only = "Responders double-positive only",
  weak_ifng = "Weak IFNγ, stronger TNF",
  null = "No response (stimulated = unstimulated)",
  null_background_coexpression = "No response, co-expressing background",
  background_coexpression = "Co-expressing background"
)

# Quantities whose background-subtracted frequencies are compared. Marginal
# frequencies count every cell positive for a cytokine; exclusive frequencies
# count one combination only.
.simLowSepQuantities <- tibble::tribble(
  ~quantity, ~quantity_type, ~quantity_lab,
  "TNF+", "marginal", "TNF+ (any)",
  "IFNg+", "marginal", "IFNγ+ (any)",
  "TNF+IFNg+", "exclusive", "TNF+ IFNγ+",
  "TNF+IFNg-", "exclusive", "TNF+ IFNγ- only",
  "TNF-IFNg+", "exclusive", "TNF- IFNγ+ only"
)

# Response tests, in display order.
.simLowSepTests <- c(
  fisher = "Fisher's exact test, each sample",
  wilcoxon = "Wilcoxon signed-rank test, each dataset",
  t = "Paired t-test, each dataset"
)

# Fixed scientific settings shared by every scenario. `loc_threshold_method`
# is passed to stimControl(locThresholdMethod = ). The antimode rule is run
# by gateStim(); both co-expression rules are applied afterwards to the
# ordinary gates with `coexpr`.
.simLowSepMainSettings <- function(loc_threshold_method = "region", n_sample = 12L) {
  list(
    transformation = "gaussian",
    cov_ev = 1.5,
    coexpr = list(nBin = 20L, rMin = 3.5, frac = 0.75, zMin = 2),
    n_cell_uns_relative_to_stim = 1,
    # Shares of the responding (and background) cells in each positive
    # combination, chosen per scenario by name.
    splits = .simLowSepSplits,
    prob_exact = TRUE,
    calc_cyt_pos_gates = TRUE,
    cyt_pos_method = "antimode",
    cluster_gates = TRUE,
    loc_threshold_method = loc_threshold_method,
    n_sample = as.integer(n_sample)
  )
}

# Datasets per scenario and samples per dataset. Final runs use 20 datasets
# of 12 samples and draft runs 5; quick and dev runs use 2 datasets of 3
# samples (smoke checks only).
.simLowSepSizes <- function(sim_size, quick = FALSE, dev = FALSE) {
  if (isTRUE(quick) || isTRUE(dev)) {
    return(list(n_dataset = 2L, n_sample = 3L))
  }
  if (!sim_size %in% c("final", "draft")) {
    stop("sim_size must be \"final\" or \"draft\".")
  }
  list(n_dataset = if (identical(sim_size, "final")) .simLowSepIterMax else 5L, n_sample = 12L)
}

# The local-FDR threshold method StimGate resolved and saved for every marker.
# It must equal the requested method, so output rows record what was used.
.simLowSepThresholdMethod <- function(path_project, requested) {
  saved <- stimgateMetaReadSettingsChnls(path_project)
  method <- unique(vapply(saved, function(x) {
    as.character(x$locThresholdMethod %||% NA_character_)
  }, character(1L)))
  if (!identical(method, requested)) {
    stop(
      "StimGate saved locThresholdMethod '", paste(method, collapse = "', '"),
      "' but '", requested, "' was requested."
    )
  }
  method
}

# Grids ---------------------------------------------------------------------

# Columns of a grid without batch mismatch (Analysis 11a).
.simLowSepNoMismatch <- function(grid) {
  dplyr::mutate(grid,
    base_scenario_id = .data$scenario_id,
    mismatch_type = "none", mismatch_label = "No mismatch",
    stim_mean_shift = 0, stim_sd_multiplier = 1
  )
}

# Full scenario grid of Analysis 11a. IDs and seeds are assigned here, before
# any dev/quick filtering, so a scenario keeps its data whatever subset is
# run. Each family block is appended after the earlier ones with seeds from
# its own stream, so earlier scenarios keep their IDs and seeds: the main grid
# (stream `simulation_seed`), the first families (`simulation_seed + 1`) and
# the weak double-positive and co-expressing no-response families
# (`simulation_seed + 2`).
.simLowSepGrid <- function(simulation_seed) {
  main <- tidyr::expand_grid(
    tibble::tibble(separation = c("low", "lower"), mean_pos = c(4.5, 3.5)),
    tibble::tibble(response_level = c("lower", "higher"), prob_response = c(0.01, 0.05)),
    n_cell = c(1e5, 5e3)
  )
  seeds_main <- .analysis_with_seed(simulation_seed, {
    sample.int(.Machine$integer.max, nrow(main))
  })
  main <- dplyr::mutate(main, family = "main", mean_pos_ifng = .data$mean_pos,
    response_split = "default", background_relative = 0.2, background_split = "response",
    scenario_seed = as.integer(seeds_main))
  resp <- tibble::tibble(response_level = c("lower", "higher"), prob_response = c(0.01, 0.05))
  sep_null <- tibble::tibble(separation = c("low", "lower"), mean_pos = c(4.5, 3.5))
  cells <- c(1e5, 5e3)
  extra <- dplyr::bind_rows(
    tidyr::expand_grid(family = "high_separation", separation = "high", mean_pos = 6,
      mean_pos_ifng = 6, resp, response_split = "default", background_relative = 0.2,
      background_split = "response", n_cell = cells),
    tidyr::expand_grid(family = "no_coexpression", separation = "low", mean_pos = 4.5,
      mean_pos_ifng = 4.5, resp, response_split = "single", background_relative = 0.2,
      background_split = "response", n_cell = cells),
    tidyr::expand_grid(family = "double_only", separation = "low", mean_pos = 4.5,
      mean_pos_ifng = 4.5, resp, response_split = "double", background_relative = 0.2,
      background_split = "response", n_cell = cells),
    tidyr::expand_grid(family = "weak_ifng", separation = "weak IFNg", mean_pos = 4.5,
      mean_pos_ifng = 3, resp, response_split = "default", background_relative = 0.2,
      background_split = "response", n_cell = cells),
    tidyr::expand_grid(family = "null", sep_null,
      response_level = "none", prob_response = 0.01, response_split = "default",
      background_relative = 1, background_split = "response", n_cell = cells) |>
      dplyr::mutate(mean_pos_ifng = .data$mean_pos),
    tidyr::expand_grid(family = "background_coexpression", separation = "low", mean_pos = 4.5,
      mean_pos_ifng = 4.5, resp, response_split = "default", background_relative = 0.5,
      background_split = "double", n_cell = cells)
  )
  seeds_extra <- .analysis_with_seed(simulation_seed + 1L, {
    sample.int(.Machine$integer.max, nrow(extra))
  })
  extra$scenario_seed <- as.integer(seeds_extra)
  # Appended families: a very weak double-positive response (responders
  # mostly single-positive) and no response with an all double-positive
  # background shared by both tubes.
  added <- dplyr::bind_rows(
    tidyr::expand_grid(family = "weak_double", separation = "low", mean_pos = 4.5,
      mean_pos_ifng = 4.5, resp, response_split = "weak_double", background_relative = 0.2,
      background_split = "response", n_cell = cells),
    tidyr::expand_grid(family = "null_background_coexpression", sep_null,
      response_level = "none", prob_response = 0.01, response_split = "double",
      background_relative = 1, background_split = "response", n_cell = cells) |>
      dplyr::mutate(mean_pos_ifng = .data$mean_pos)
  )
  seeds_added <- .analysis_with_seed(simulation_seed + 2L, {
    sample.int(.Machine$integer.max, nrow(added))
  })
  added$scenario_seed <- as.integer(seeds_added)
  grid <- dplyr::bind_rows(main, extra, added)
  grid$scenario_id <- seq_len(nrow(grid))
  grid <- .simLowSepNoMismatch(grid)
  dplyr::select(grid, dplyr::all_of(c("scenario_id", "scenario_seed", .simLowSepScenCols[-1])))
}

# Seeds of a scenario's datasets, drawn before any method runs: the first
# dataset uses the scenario seed itself, later ones are drawn from it.
.simLowSepIterSeeds <- function(scenario_seed, n_iter = .simLowSepIterMax) {
  extra <- .analysis_with_seed(scenario_seed, {
    sample.int(.Machine$integer.max, n_iter - 1L)
  })
  as.integer(c(scenario_seed, extra))
}

# One row per scenario and dataset. `sim_id` comes from the scenario ID and
# dataset number on the full grid (`.simLowSepIterMax` datasets), and each
# dataset's seed from the scenario's data seed (`scenario_seed`), so rows
# sharing a scenario seed (11b mismatch settings) share their data.
.simLowSepRowGrid <- function(scenario_grid, n_iter = .simLowSepIterMax) {
  rows <- lapply(seq_len(nrow(scenario_grid)), function(i) {
    scen <- scenario_grid[i, , drop = FALSE]
    seeds <- .simLowSepIterSeeds(scen$scenario_seed[[1]], n_iter)
    dplyr::bind_cols(
      scen[rep(1L, n_iter), , drop = FALSE],
      tibble::tibble(
        iter = seq_len(n_iter),
        sim_id = as.integer((scen$scenario_id[[1]] - 1L) * n_iter + seq_len(n_iter)),
        sim_seed = seeds
      )
    )
  })
  out <- dplyr::bind_rows(rows)
  dplyr::select(out, dplyr::all_of(c("sim_id", "sim_seed", "iter")), dplyr::everything())
}

# Batch mismatch settings of Analysis 11b, applied to the stimulated tube's
# all-negative (TNF- IFNg-) component only: Analysis 2b's Gaussian mean
# shifts (0 is the no-mismatch baseline) and 10% SD inflation.
.simLowSepBatchMismatch <- function() {
  dplyr::bind_rows(
    tibble::tibble(
      mismatch_type = "mean_shift", stim_mean_shift = c(0, 0.025, 0.05, 0.2),
      stim_sd_multiplier = 1
    ),
    tibble::tibble(mismatch_type = "sd_inflation", stim_mean_shift = 0, stim_sd_multiplier = 1.1)
  ) |>
    dplyr::mutate(mismatch_label = dplyr::case_when(
      .data$mismatch_type == "sd_inflation" ~ "SD inflation 10%",
      .data$stim_mean_shift == 0 ~ "No mismatch",
      TRUE ~ paste0("Mean shift ", .data$stim_mean_shift)
    ))
}

# Baseline scenarios of Analysis 11b: a subset of Analysis 11a's families
# (main low separation, single-positive only, weak double-positive,
# double-positive only, and the two no-response families at low separation),
# every response level and both cell counts.
.simLowSepBatchBaseline <- function(grid_11a) {
  dplyr::filter(grid_11a,
    .data$family %in% c("main", "no_coexpression", "weak_double", "double_only",
      "null", "null_background_coexpression"),
    .data$separation == "low"
  )
}

# Full scenario grid of Analysis 11b: each baseline scenario (with Analysis
# 11a's scenario ID as `base_scenario_id` and its data seed) crossed with
# every mismatch setting. Mismatch settings of one baseline share its seed
# (common random numbers): the mismatch is deterministic, so the
# no-mismatch rows reproduce Analysis 11a's datasets.
.simLowSepBatchGrid <- function(simulation_seed) {
  base <- .simLowSepBatchBaseline(.simLowSepGrid(simulation_seed)) |>
    dplyr::select(-dplyr::all_of(c("mismatch_type", "mismatch_label", "stim_mean_shift",
      "stim_sd_multiplier", "base_scenario_id"))) |>
    dplyr::rename(base_scenario_id = "scenario_id")
  grid <- tidyr::expand_grid(base, .simLowSepBatchMismatch())
  grid$scenario_id <- seq_len(nrow(grid))
  dplyr::select(grid, dplyr::all_of(c("scenario_id", "scenario_seed", .simLowSepScenCols[-1])))
}

# Labels ---------------------------------------------------------------------

# Scenario label. With `mismatch = TRUE` the 11b mismatch is appended, and
# with `cells = TRUE` the cell count.
.simLowSepScenarioLab <- function(grid, mismatch = TRUE, cells = FALSE) {
  resp <- ifelse(grid$response_level == "none", "no response (background 1%)",
    paste0(grid$response_level, " response (", 100 * grid$prob_response, "%)"))
  means <- ifelse(grid$mean_pos == grid$mean_pos_ifng,
    paste0(grid$separation, " separation (mean ", grid$mean_pos, ")"),
    paste0("TNF mean ", grid$mean_pos, ", IFNγ mean ", grid$mean_pos_ifng))
  # Families other than the main grid share separations and responses, so
  # they are named in the label.
  fam <- ifelse(grid$family == "main", "", paste0(.simLowSepFamilies[grid$family], ": "))
  out <- paste0(fam, means, ", ", resp)
  if (isTRUE(cells)) {
    out <- paste0(out, ", ", .simLowSepCellLab(grid$n_cell))
  }
  if (isTRUE(mismatch) && "mismatch_type" %in% names(grid)) {
    mm <- grid$mismatch_type != "none"
    out[mm] <- paste0(out[mm], ", stimulated TNF- IFNγ- cells: ", tolower(grid$mismatch_label[mm]))
  }
  out
}

# Cell count without scientific notation, for figure file names.
.simLowSepCellKey <- function(n_cell) {
  format(n_cell, scientific = FALSE, trim = TRUE)
}

.simLowSepCellLab <- function(n_cell) {
  paste0(format(n_cell, big.mark = ",", scientific = FALSE, trim = TRUE), " cells")
}

.simLowSepFacetLabels <- function(tbl, mismatch = TRUE, cells = FALSE) {
  ord <- order(
    match(tbl$family, names(.simLowSepFamilies)), -tbl$mean_pos, -tbl$mean_pos_ifng,
    tbl$prob_response, -tbl$n_cell,
    if ("stim_sd_multiplier" %in% names(tbl)) tbl$stim_sd_multiplier else 0,
    if ("stim_mean_shift" %in% names(tbl)) tbl$stim_mean_shift else 0
  )
  lab <- .simLowSepScenarioLab(tbl, mismatch = mismatch, cells = cells)
  tbl$scenario_lab <- factor(lab, levels = unique(lab[ord]))
  tbl
}

# Expected net (stimulated minus unstimulated) frequency of each quantity in
# each scenario, from the simulation settings. Quantities with no expected
# net response are those where a positive test is a false positive (the
# batch mismatch of 11b is not a response).
.simLowSepExpectedNet <- function(grid, splits = .simLowSepSplits) {
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    r <- grid[i, , drop = FALSE]
    split <- splits[[r$response_split[[1]]]][.simLowSepCombn[-1L]]
    split_bg <- if (identical(r$background_split[[1]], "response")) {
      split
    } else {
      splits[[r$background_split[[1]]]][.simLowSepCombn[-1L]]
    }
    net <- r$prob_response[[1]] * (split - r$background_relative[[1]] * split_bg)
    vals <- c(
      "TNF+" = net[["TNF+IFNg-"]] + net[["TNF+IFNg+"]],
      "IFNg+" = net[["TNF-IFNg+"]] + net[["TNF+IFNg+"]],
      net
    )
    tibble::tibble(scenario_id = r$scenario_id[[1]], quantity = names(vals), net_expected = unname(vals))
  }) |>
    dplyr::filter(.data$quantity %in% .simLowSepQuantities$quantity) |>
    dplyr::mutate(response_truth = dplyr::if_else(abs(.data$net_expected) > 1e-12, "response", "none"))
}

# Simulation and gating ------------------------------------------------------

# Simulate one dataset: `n_sample` unstimulated/stimulated pairs. A batch
# mismatch in `row` (stim_mean_shift, stim_sd_multiplier) is applied by
# simcyto to the stimulated tube's TNF- IFNg- component on both markers; it
# is deterministic, so it uses no random numbers.
.simLowSepSimulate <- function(row, n_sample, settings = .simLowSepMainSettings()) {
  if (!identical(settings$transformation, "gaussian")) {
    stop("This demonstration uses the Gaussian transformation only.")
  }
  split <- settings$splits[[row$response_split]][.simLowSepCombn[-1L]]
  split_bg <- if (identical(row$background_split, "response")) {
    split
  } else {
    settings$splits[[row$background_split]][.simLowSepCombn[-1L]]
  }
  if (!isTRUE(all.equal(sum(split), 1)) || !isTRUE(all.equal(sum(split_bg), 1))) {
    stop("Response and background splits must each sum to one.")
  }
  prob_resp <- row$prob_response * split
  prob_uns <- row$prob_response * row$background_relative * split_bg
  n_cell_uns <- round(row$n_cell * settings$n_cell_uns_relative_to_stim)
  args <- list(
    nSample = n_sample,
    nMarker = 2L,
    nCondition = 2L,
    nCluster = 4L,
    nCellByCondition = c(n_cell_uns, row$n_cell),
    transformationFunc = simcyto::simCytTransformGaussian(),
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(
      c(0, 0, row$mean_pos, 0, 0, row$mean_pos_ifng, row$mean_pos, row$mean_pos_ifng),
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
  args <- c(args, .simLowSepMismatchArgs(row))
  res <- do.call(simcyto::simCytExperiment, args)
  res$flowFrameList <- lapply(res$flowFrameList, function(ff) {
    flowCore::markernames(ff) <- stats::setNames(
      unname(.simLowSepMarkers), unname(.simLowSepChnl)
    )
    ff
  })
  res
}

# simcyto arguments for the row's batch mismatch: none without one (so 11a
# and the 11b baseline simulate exactly as before), otherwise the shift or
# SD multiplier restricted to the TNF- IFNg- component. Fails when simcyto
# cannot restrict the mismatch to one component.
.simLowSepMismatchArgs <- function(row) {
  shift <- if ("stim_mean_shift" %in% names(row)) row$stim_mean_shift[[1]] else 0
  sd_mult <- if ("stim_sd_multiplier" %in% names(row)) row$stim_sd_multiplier[[1]] else 1
  out <- list()
  if (!identical(as.numeric(shift), 0) || !identical(as.numeric(sd_mult), 1)) {
    supported <- c("stimMeanShiftClusters", "stimSdMultiplierClusters") %in%
      names(formals(simcyto::simCytExperiment))
    if (!all(supported)) {
      stop("The installed simcyto cannot apply a mismatch to one component ",
        "(stimMeanShiftClusters / stimSdMultiplierClusters); update simcyto.")
    }
  }
  if (!identical(as.numeric(shift), 0)) {
    out$stimMeanShift <- shift
    out$stimMeanShiftClusters <- .simLowSepCombn[[1]]
  }
  if (!identical(as.numeric(sd_mult), 1)) {
    out$stimSdMultiplier <- sd_mult
    out$stimSdMultiplierClusters <- .simLowSepCombn[[1]]
  }
  out
}

# Combination calls for one tube. The antimode rule matches the package:
# positive by the ordinary gate, or above the cytokine-positive gate while
# ordinarily positive for the other cytokine. The co-expression rules
# (`.coexPositive()`) use the lowered gates `low` (`.coexLowerGates()`) of
# that rule.
.simLowSepCall <- function(x_tnf, x_ifng, gate, gate_type, low = NULL) {
  pos_tnf <- x_tnf > gate[["tnf"]]
  pos_ifng <- x_ifng > gate[["ifng"]]
  if (gate_type %in% names(.simLowSepCoexScreen)) {
    pos <- .coexPositive(
      data.frame(TNF = x_tnf, IFNg = x_ifng),
      c(TNF = gate[["tnf"]], IFNg = gate[["ifng"]]), low
    )
    pos_tnf <- pos[, "TNF"]
    pos_ifng <- pos[, "IFNg"]
  } else if (identical(gate_type, "antimode")) {
    pos_tnf_cyt <- pos_tnf | (x_tnf > gate[["tnf_cyt"]] & pos_ifng)
    pos_ifng <- pos_ifng | (x_ifng > gate[["ifng_cyt"]] & pos_tnf)
    pos_tnf <- pos_tnf_cyt
  } else if (!identical(gate_type, "ordinary")) {
    stop("gate_type must be one of: ", paste0("\"", names(.simLowSepGateTypes), "\"", collapse = ", "), ".")
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

# Co-expression gates for one stimulated sample, for both rules: every
# ordered pair lowered from the ordinary gates with its own and its
# unstimulated tube's cells, with (`coex_screened`) and without
# (`coex_unscreened`) the double-positive screen.
.simLowSepCoexGates <- function(expr_list, gate, s, settings) {
  tube <- function(ind) {
    ex <- expr_list[[ind]]
    data.frame(TNF = ex[, .simLowSepChnl[["tnf"]]], IFNg = ex[, .simLowSepChnl[["ifng"]]])
  }
  dat <- list(stim = tube(2L * s), uns = tube(2L * s - 1L),
    gate = c(TNF = gate[["tnf"]], IFNg = gate[["ifng"]]))
  cs <- settings$coexpr
  purrr::map_dfr(names(.simLowSepCoexScreen), function(gt) {
    low <- .coexLowerGates(dat, nBin = cs$nBin, rMin = cs$rMin, frac = cs$frac, zMin = cs$zMin,
      screen = .simLowSepCoexScreen[[gt]])
    dplyr::mutate(tibble::as_tibble(low), gate_type = gt, .before = 1L)
  })
}

# Confusion counts (truth by call) for every tube and rule. `low_list` holds
# each sample's co-expression gates (`.simLowSepCoexGates()`).
.simLowSepCounts <- function(expr_list, labels_list, gate_tbl, n_sample, low_list) {
  purrr::map_dfr(seq_len(n_sample), function(s) {
    ind_uns <- 2L * s - 1L
    ind_stim <- 2L * s
    gate <- .simLowSepSampleGate(gate_tbl, ind_stim)
    purrr::map_dfr(c(uns = ind_uns, stim = ind_stim), function(ind) {
      ex <- expr_list[[ind]]
      truth <- factor(labels_list[[ind]], levels = .simLowSepCombn)
      if (anyNA(truth)) stop("Unexpected simulated labels in tube ", ind, ".")
      purrr::map_dfr(names(.simLowSepGateTypes), function(gt) {
        low <- if (gt %in% names(.simLowSepCoexScreen)) {
          low_list[[s]][low_list[[s]]$gate_type == gt, , drop = FALSE]
        }
        call <- factor(
          .simLowSepCall(ex[, .simLowSepChnl[["tnf"]]], ex[, .simLowSepChnl[["ifng"]]], gate, gt, low),
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

# The package's saved combination counts use the antimode rule; the counts
# recomputed here must reproduce them exactly.
.simLowSepValidateCounts <- function(counts, stats_tbl) {
  combn_code <- c(
    "TNF+IFNg-" = "F1~+~F2~-~", "TNF-IFNg+" = "F1~-~F2~+~",
    "TNF+IFNg+" = "F1~+~F2~+~", "TNF-IFNg-" = "F1~-~F2~-~"
  )
  mine <- counts |>
    dplyr::filter(.data$gate_type == "antimode") |>
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
      "Recomputed antimode combination counts differ from getStimStats() in ",
      sum(bad), " rows."
    )
  }
  invisible(TRUE)
}

# Scenario function for `.simBandwidthRunRow()`: simulate and gate one
# dataset (it does not seed; the runner seeds with `row$sim_seed`). gateStim()
# runs the antimode rule (`cytPosMethod = "antimode"`, set explicitly), whose
# recomputed counts must equal getStimStats(); the co-expression rules are
# applied to the ordinary gates afterwards. Returns one row with list columns:
# gates, co-expression gates, confusion counts and, for the first dataset,
# the first sample's stimulated cells (chosen in advance for the figures).
.simLowSepScenario <- function(row, settings) {
  if (nrow(row) != 1L) stop("row must have exactly one scenario.")
  loc_threshold_method <- settings$loc_threshold_method
  if (!is.character(loc_threshold_method) || length(loc_threshold_method) != 1L) {
    stop("settings$loc_threshold_method must be set explicitly.")
  }
  if (!isTRUE(settings$calc_cyt_pos_gates) || !identical(settings$cyt_pos_method, "antimode")) {
    stop("settings must request the antimode gates (calc_cyt_pos_gates = TRUE, cyt_pos_method = ",
      "\"antimode\"): the co-expression rules are applied afterwards.")
  }
  n_sample <- as.integer(settings$n_sample)
  if (length(n_sample) != 1L || is.na(n_sample) || n_sample < 1L) {
    stop("settings$n_sample must be one positive integer.")
  }
  sim <- .simLowSepSimulate(row, n_sample, settings)
  gs <- flowWorkspace::GatingSet(methods::as(sim$flowFrameList, "flowSet"))
  id <- if ("sim_id" %in% names(row)) row$sim_id[[1]] else row$scenario_id[[1]]
  path_project <- tempfile(paste0("lowsep-", id, "-"))
  on.exit(unlink(path_project, recursive = TRUE), add = TRUE)
  batch_list <- lapply(seq_len(n_sample), function(s) c(2L * s - 1L, 2L * s))
  suppressMessages(gateStim(
    pathProject = path_project, .data = gs, batchList = batch_list,
    marker = unname(.simLowSepMarkers),
    calcCytPosGates = settings$calc_cyt_pos_gates,
    control = stimControl(
      cytPosMethod = settings$cyt_pos_method,
      clusterGates = settings$cluster_gates,
      locThresholdMethod = loc_threshold_method
    )
  ))
  method <- .simLowSepThresholdMethod(path_project, loc_threshold_method)
  gate_tbl <- getStimGates(path_project)
  # Classify the single-precision expression that StimGate gated.
  expr_list <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]], "root"))
  })
  low_list <- lapply(seq_len(n_sample), function(s) {
    .simLowSepCoexGates(expr_list, .simLowSepSampleGate(gate_tbl, 2L * s), s, settings)
  })
  counts <- .simLowSepCounts(expr_list, sim$labelsList, gate_tbl, n_sample, low_list)
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
      loc_source = if ("locSource" %in% names(rows)) rows$locSource[match(.simLowSepChnl, rows$chnl)] else NA_character_,
      locThresholdMethod = method
    )
  })
  coex <- purrr::map_dfr(seq_len(n_sample), function(s) {
    dplyr::mutate(low_list[[s]], sample = s, .before = 1L)
  })
  first_dataset <- !"iter" %in% names(row) || isTRUE(row$iter[[1]] == 1L)
  cells <- if (first_dataset) {
    tibble::tibble(
      sample = 1L,
      tnf = expr_list[[2L]][, .simLowSepChnl[["tnf"]]],
      ifng = expr_list[[2L]][, .simLowSepChnl[["ifng"]]],
      truth = sim$labelsList[[2L]]
    )
  } else {
    tibble::tibble(sample = integer(), tnf = numeric(), ifng = numeric(), truth = character())
  }
  tibble::tibble(
    n_sample_run = n_sample,
    gates = list(gates),
    coex_gates = list(coex),
    counts = list(counts),
    cells = list(cells)
  )
}

# Validation of collated rows (after the shared checks of sim_id, seeds and
# runtime errors): every dataset has the expected samples, both
# co-expression rules for both ordered pairs, a full confusion table per
# tube and rule, and the requested local-FDR threshold method.
.simLowSepValidator <- function(n_sample, loc_threshold_method) {
  function(tbl) {
    need <- c("n_sample_run", "gates", "coex_gates", "counts", "cells")
    if (!all(need %in% names(tbl))) {
      return("missing low-separation result columns")
    }
    n_rows <- function(col) vapply(tbl[[col]], NROW, integer(1))
    n_combn <- length(.simLowSepCombn)^2
    bad <- !(tbl$n_sample_run %in% n_sample) |
      n_rows("gates") != 2L * n_sample |
      n_rows("coex_gates") != 2L * length(.simLowSepCoexScreen) * n_sample |
      n_rows("counts") != n_sample * 2L * length(.simLowSepGateTypes) * n_combn
    problems <- character()
    if (any(bad)) {
      problems <- paste0("unexpected result sizes for sim_id: ",
        paste(sort(as.integer(tbl$sim_id[bad])), collapse = ", "))
    }
    method_ok <- vapply(tbl$gates, function(g) {
      !is.null(g) && "locThresholdMethod" %in% names(g) &&
        all(g$locThresholdMethod %in% loc_threshold_method)
    }, logical(1))
    if (!all(method_ok)) {
      problems <- c(problems, paste0(
        "gates without the requested locThresholdMethod for sim_id: ",
        paste(sort(as.integer(tbl$sim_id[!method_ok])), collapse = ", ")
      ))
    }
    problems
  }
}

# Unnest the per-dataset list columns into one table each, keyed by the
# scenario and row columns. Cells keep only their scenario and dataset keys
# (join `scenarios` for labels), since they are large.
.simLowSepCollate <- function(tbl) {
  tbl <- dplyr::arrange(tbl, .data$sim_id)
  keys <- intersect(c(.simLowSepScenCols, .simLowSepRowCols), names(tbl))
  unnest_col <- function(col, cols = keys) {
    out <- lapply(seq_len(nrow(tbl)), function(i) {
      x <- tbl[[col]][[i]]
      if (is.null(x) || nrow(x) == 0L) return(NULL)
      dplyr::bind_cols(tbl[rep(i, nrow(x)), cols, drop = FALSE], tibble::as_tibble(x))
    })
    dplyr::bind_rows(out)
  }
  list(
    scenarios = dplyr::distinct(tbl[, intersect(.simLowSepScenCols, names(tbl)), drop = FALSE]),
    gates = unnest_col("gates"),
    coex_gates = unnest_col("coex_gates"),
    counts = unnest_col("counts"),
    cells = unnest_col("cells", intersect(c("scenario_id", "sim_id", "iter"), keys))
  )
}

# Summaries -----------------------------------------------------------------

# Scenario and row key columns present in `tbl`.
.simLowSepKeys <- function(tbl) {
  intersect(c(.simLowSepScenCols, .simLowSepRowCols), names(tbl))
}

.simLowSepScen <- function(tbl) {
  intersect(.simLowSepScenCols, names(tbl))
}

# Per tube, sample and marker: true positives, cells called positive, true
# positives recovered and negative-cell contamination (truth by marker:
# a cell is truly TNF+ if its simulated combination contains TNF+).
.simLowSepRecovery <- function(counts) {
  keys <- .simLowSepKeys(counts)
  purrr::map_dfr(c(TNF = "TNF\\+", IFNg = "IFNg\\+"), function(pattern) {
    counts |>
      dplyr::mutate(
        truth_pos = grepl(pattern, .data$truth),
        call_pos = grepl(pattern, .data$call)
      ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(keys, "sample", "tube", "gate_type")))) |>
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

# Cells each lowering rule adds over the ordinary gates, per tube, sample
# and marker. Every rule only adds positives (its positive cells include the
# ordinary ones), so the differences are the added cells, split into cells
# truly positive and truly negative for that marker.
.simLowSepAdded <- function(recovery) {
  keys <- c(.simLowSepKeys(recovery), "sample", "tube", "marker")
  ord <- recovery |>
    dplyr::filter(.data$gate_type == "ordinary") |>
    dplyr::select(dplyr::all_of(keys),
      ord_called = "n_called_pos", ord_recovered = "n_recovered", ord_contaminating = "n_contaminating")
  recovery |>
    dplyr::filter(.data$gate_type != "ordinary") |>
    dplyr::left_join(ord, by = keys) |>
    dplyr::mutate(
      n_added = .data$n_called_pos - .data$ord_called,
      n_added_true = .data$n_recovered - .data$ord_recovered,
      n_added_negative = .data$n_contaminating - .data$ord_contaminating
    ) |>
    dplyr::select(dplyr::all_of(c(keys, "gate_type", "n_cell_tube", "n_added", "n_added_true", "n_added_negative")))
}

# Whether each rule lowered each sample's gate for each marker, and by how
# much (zero when not lowered). The antimode gate is one per marker; a
# co-expression gate is the one lowered for cells positive for the other
# marker.
.simLowSepLowered <- function(gates, coex_gates) {
  keys <- c(.simLowSepKeys(gates), "sample")
  anti <- gates |>
    dplyr::transmute(dplyr::across(dplyr::all_of(keys)),
      marker = .data$marker, gate_type = "antimode",
      lowered = (.data$gate_cyt < .data$gate) %in% TRUE,
      shift = dplyr::if_else(.data$lowered, .data$gate - .data$gate_cyt, 0))
  coex <- coex_gates |>
    dplyr::transmute(dplyr::across(dplyr::all_of(keys)),
      marker = .data$b, gate_type = .data$gate_type,
      lowered = .data$lowered %in% TRUE,
      shift = dplyr::if_else(.data$lowered, .data$gateB - .data$cut, 0))
  dplyr::bind_rows(anti, coex)
}

# What each rule gates, by scenario, rule and marker: the share of samples
# whose gate was lowered, the median shift where lowered, and the cells
# added per 100,000 cells of each tube, split into cells truly positive and
# truly negative for the marker (means over samples, pooled over datasets).
.simLowSepGatingSummary <- function(gates, coex_gates, counts) {
  low <- .simLowSepLowered(gates, coex_gates)
  added <- .simLowSepAdded(.simLowSepRecovery(counts))
  scen <- .simLowSepScen(low)
  low_sum <- low |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen, "gate_type", "marker")))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_dataset = dplyr::n_distinct(.data$iter),
      n_lowered = sum(.data$lowered),
      share_lowered = mean(.data$lowered),
      median_shift_lowered = if (any(.data$lowered)) stats::median(.data$shift[.data$lowered]) else NA_real_,
      .groups = "drop"
    )
  add_sum <- added |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen, "gate_type", "marker", "tube")))) |>
    dplyr::summarise(
      added_per_1e5 = mean(1e5 * .data$n_added / .data$n_cell_tube),
      added_true_per_1e5 = mean(1e5 * .data$n_added_true / .data$n_cell_tube),
      added_negative_per_1e5 = mean(1e5 * .data$n_added_negative / .data$n_cell_tube),
      share_added_true = if (sum(.data$n_added) > 0) sum(.data$n_added_true) / sum(.data$n_added) else NA_real_,
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "tube",
      values_from = c("added_per_1e5", "added_true_per_1e5", "added_negative_per_1e5", "share_added_true")
    )
  dplyr::left_join(low_sum, add_sum, by = c(scen, "gate_type", "marker"))
}

# Background-subtracted marginal and exclusive frequencies, estimated and
# true, per stimulated sample and rule, with the tube counts behind them.
.simLowSepFrequencies <- function(counts) {
  keys <- .simLowSepKeys(counts)
  group_cols <- c(keys, "sample", "gate_type", "tube")
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
      id_cols = dplyr::all_of(c(keys, "sample", "gate_type", "quantity")),
      names_from = "tube",
      values_from = c("n_tube", "n_est", "n_truth", "prop_est", "prop_truth")
    ) |>
    dplyr::mutate(
      prop_bs_est = .data$prop_est_stim - .data$prop_est_uns,
      prop_bs_truth = .data$prop_truth_stim - .data$prop_truth_uns,
      error_bs = .data$prop_bs_est - .data$prop_bs_truth
    ) |>
    dplyr::left_join(.simLowSepQuantities, by = "quantity")
}

# Means over each scenario's samples (pooled over datasets), by rule, with
# how many samples each lowering rule moved closer to or further from the
# truth than the ordinary gates (NA for the ordinary gates themselves).
.simLowSepFrequencySummary <- function(freq) {
  keys <- .simLowSepKeys(freq)
  scen <- .simLowSepScen(freq)
  ord <- freq |>
    dplyr::filter(.data$gate_type == "ordinary") |>
    dplyr::select(dplyr::all_of(c(keys, "sample", "quantity")), error_ordinary = "error_bs")
  is_ordinary <- function() identical(as.character(dplyr::cur_group()$gate_type), "ordinary")
  freq |>
    dplyr::left_join(ord, by = c(keys, "sample", "quantity")) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen, "quantity", "quantity_type", "quantity_lab", "gate_type")))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_dataset = dplyr::n_distinct(.data$iter),
      mean_prop_bs_truth = mean(.data$prop_bs_truth),
      mean_prop_bs_est = mean(.data$prop_bs_est),
      mean_error_bs = mean(.data$error_bs),
      mean_abs_error_bs = mean(abs(.data$error_bs)),
      n_closer = if (is_ordinary()) NA_integer_ else sum(abs(.data$error_bs) < abs(.data$error_ordinary)),
      n_further = if (is_ordinary()) NA_integer_ else sum(abs(.data$error_bs) > abs(.data$error_ordinary)),
      n_unchanged = if (is_ordinary()) NA_integer_ else sum(abs(.data$error_bs) == abs(.data$error_ordinary)),
      .groups = "drop"
    )
}

# Cell-level precision, recall (sensitivity) and F1 of each quantity's calls
# in every stimulated tube: true positives are cells whose simulated label and
# call both belong to the quantity. Undefined proportions are NA.
.simLowSepF1 <- function(counts) {
  keys <- .simLowSepKeys(counts)
  stim <- dplyr::filter(counts, .data$tube == "stim")
  purrr::map_dfr(.simLowSepQuantities$quantity, function(q) {
    hit <- function(x) {
      if (q == "TNF+") grepl("TNF\\+", x) else if (q == "IFNg+") grepl("IFNg\\+", x) else x == q
    }
    stim |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(keys, "sample", "gate_type")))) |>
      dplyr::summarise(
        tp = sum(.data$n[hit(.data$truth) & hit(.data$call)]),
        n_called = sum(.data$n[hit(.data$call)]),
        n_true = sum(.data$n[hit(.data$truth)]),
        .groups = "drop"
      ) |>
      dplyr::mutate(quantity = q)
  }) |>
    dplyr::mutate(
      precision = dplyr::if_else(.data$n_called > 0, .data$tp / .data$n_called, NA_real_),
      recall = dplyr::if_else(.data$n_true > 0, .data$tp / .data$n_true, NA_real_),
      f1 = dplyr::if_else(.data$n_called + .data$n_true > 0,
        2 * .data$tp / (.data$n_called + .data$n_true), NA_real_)
    ) |>
    dplyr::left_join(.simLowSepQuantities, by = "quantity")
}

# 10th, 50th and 90th percentiles of a per-sample value `col` over every
# sample of every dataset in a scenario, by quantity and rule, with the
# median truth when `truth` names a column. Undefined values are dropped and
# counted.
.simLowSepPercentiles <- function(tbl, col, truth = NULL) {
  scen <- .simLowSepScen(tbl)
  quant_finite <- function(x, p) {
    if (any(is.finite(x))) stats::quantile(x, p, na.rm = TRUE, names = FALSE) else NA_real_
  }
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen, "quantity", "quantity_lab", "gate_type")))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      n_dataset = if ("iter" %in% names(tbl)) dplyr::n_distinct(.data$iter) else NA_integer_,
      n_defined = sum(is.finite(.data[[col]])),
      p10 = quant_finite(.data[[col]], 0.1),
      p50 = quant_finite(.data[[col]], 0.5),
      p90 = quant_finite(.data[[col]], 0.9),
      truth = if (is.null(truth)) NA_real_ else stats::median(.data[[truth]]),
      .groups = "drop"
    )
}

# Response tests --------------------------------------------------------------

# One-sided Fisher's exact test p-value that the stimulated tube has a larger
# share of positive cells than the unstimulated tube: `a` of `n_a`
# stimulated and `b` of `n_b` unstimulated cells are positive. Equal to
# stats::fisher.test(matrix(c(a, b, n_a - a, n_b - b), 2),
# alternative = "greater")$p.value, vectorised.
.simLowSepFisherGreater <- function(a, n_a, b, n_b) {
  stats::phyper(a - 1, a + b, (n_a - a) + (n_b - b), n_a, lower.tail = FALSE)
}

# One-sided p-value that paired differences `d` (stimulated minus
# unstimulated frequency, one per sample) are positive: Wilcoxon signed-rank
# (zero differences dropped, as wilcox.test() does) or t. With no non-zero
# difference (Wilcoxon) or no spread (t) there is no evidence of a response
# and the p-value is 1, or 0 for t when every difference is the same
# positive value. NA with fewer than two samples.
.simLowSepPairedP <- function(d, test = c("wilcoxon", "t")) {
  test <- match.arg(test)
  d <- d[is.finite(d)]
  if (length(d) < 2L) {
    return(NA_real_)
  }
  if (identical(test, "wilcoxon")) {
    if (all(d == 0)) {
      return(1)
    }
    return(suppressWarnings(stats::wilcox.test(d, mu = 0, alternative = "greater")$p.value))
  }
  if (stats::sd(d) == 0) {
    return(if (mean(d) > 0) 0 else 1)
  }
  stats::t.test(d, mu = 0, alternative = "greater")$p.value
}

# Per stimulated sample, rule and quantity: Fisher's exact test of the
# stimulated against the unstimulated tube's counts.
.simLowSepSampleTests <- function(freq, alpha = 0.05) {
  keys <- .simLowSepKeys(freq)
  freq |>
    dplyr::mutate(
      p_value = .simLowSepFisherGreater(.data$n_est_stim, .data$n_tube_stim, .data$n_est_uns, .data$n_tube_uns),
      reject = .data$p_value < alpha
    ) |>
    dplyr::select(dplyr::all_of(c(keys, "sample", "gate_type", "quantity", "quantity_type", "quantity_lab",
      "n_est_stim", "n_tube_stim", "n_est_uns", "n_tube_uns", "prop_bs_est", "prop_bs_truth",
      "p_value", "reject")))
}

# Per dataset, rule and quantity: paired Wilcoxon signed-rank and t tests
# across the dataset's samples of stimulated minus unstimulated frequency.
.simLowSepDatasetTests <- function(freq, alpha = 0.05) {
  keys <- .simLowSepKeys(freq)
  freq |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(keys, "gate_type", "quantity", "quantity_type", "quantity_lab")))) |>
    dplyr::summarise(
      n_sample = dplyr::n(),
      mean_prop_bs_est = mean(.data$prop_bs_est),
      p_wilcoxon = .simLowSepPairedP(.data$prop_bs_est, "wilcoxon"),
      p_t = .simLowSepPairedP(.data$prop_bs_est, "t"),
      .groups = "drop"
    ) |>
    dplyr::mutate(reject_wilcoxon = .data$p_wilcoxon < alpha, reject_t = .data$p_t < alpha)
}

# Clopper-Pearson interval for `x` rejections of `n` independent datasets.
.simLowSepBinomInterval <- function(x, n) {
  if (!is.finite(n) || n < 1L) {
    return(c(NA_real_, NA_real_))
  }
  as.numeric(stats::binom.test(x, n)$conf.int)
}

# Rejection rates at `alpha` by scenario, rule and quantity, for each test.
# Fisher's test is per sample, pooled over datasets, with a 95% dataset
# bootstrap interval (datasets are the independent unit; the same resampled
# datasets are used for every rule, quantity and mismatch setting of a
# baseline scenario). The paired tests are per dataset, with a 95%
# Clopper-Pearson interval. `response_truth` says whether the quantity has
# a net response (power) or not (false-positive rate). With `mcse = FALSE`
# the rates are unchanged and no intervals are computed.
.simLowSepRejectionRates <- function(sample_tests, dataset_tests, expected_net, mcse = TRUE) {
  scen <- .simLowSepScen(sample_tests)
  grp <- c(scen, "gate_type", "quantity", "quantity_lab")
  fisher <- sample_tests |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::group_modify(function(d, key) {
      cols <- .analysis_mcse_pooled_cols(
        as.numeric(d$reject), d$iter, mean, "rate",
        family = c("low-separation", key$base_scenario_id[[1]]), mcse = mcse
      )
      tibble::tibble(
        n_unit = sum(!is.na(d$reject)), n_reject = sum(d$reject, na.rm = TRUE),
        n_dataset = dplyr::n_distinct(d$iter),
        rate = cols$rate, rate_lower = cols$rate_lower, rate_upper = cols$rate_upper,
        interval_available = cols$interval_available_rate
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(test = "fisher", unit = "sample")
  paired <- purrr::map_dfr(c("wilcoxon", "t"), function(tst) {
    dataset_tests |>
      dplyr::mutate(.reject = .data[[paste0("reject_", tst)]]) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
      dplyr::summarise(
        n_unit = sum(!is.na(.data$.reject)), n_reject = sum(.data$.reject, na.rm = TRUE),
        n_dataset = dplyr::n_distinct(.data$iter),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        rate = dplyr::if_else(.data$n_unit > 0, .data$n_reject / .data$n_unit, NA_real_),
        rate_lower = if (isTRUE(mcse)) purrr::map2_dbl(.data$n_reject, .data$n_unit, ~ .simLowSepBinomInterval(.x, .y)[1]) else NA_real_,
        rate_upper = if (isTRUE(mcse)) purrr::map2_dbl(.data$n_reject, .data$n_unit, ~ .simLowSepBinomInterval(.x, .y)[2]) else NA_real_,
        interval_available = isTRUE(mcse) & .data$n_unit > 0,
        test = tst, unit = "dataset"
      )
  })
  dplyr::bind_rows(fisher, paired) |>
    dplyr::left_join(dplyr::select(expected_net, "scenario_id", "quantity", "net_expected", "response_truth"),
      by = c("scenario_id", "quantity")) |>
    dplyr::relocate("test", "unit", .after = "quantity_lab")
}

# Plots (return ggplot objects; no files are written here) -------------------

.simLowSepGateFactor <- function(x) {
  factor(x, levels = names(.simLowSepGateTypes))
}

.simLowSepScaleRule <- function(aesthetics = "colour") {
  scale_colour_manual(values = .simLowSepGateColours, labels = .simLowSepGateTypes, name = NULL,
    aesthetics = aesthetics, drop = TRUE)
}

# First sample of the first dataset of each scenario (chosen in advance).
.simLowSepFirst <- function(tbl) {
  out <- tbl
  if ("iter" %in% names(out)) out <- dplyr::filter(out, .data$iter == 1L)
  dplyr::filter(out, .data$sample == 1L)
}

# Hex plot of TNF against IFNg in the first sample's stimulated tube of the
# first dataset. Ordinary gates span the plot; each antimode gate is drawn
# only where the other cytokine is ordinarily positive, the only cells it
# applies to; co-expression gates are drawn from their raised cut on the
# other cytokine, where lowered. `cells` carries the scenario columns.
.simLowSepPlotHex <- function(cells, gates, coex_gates, bins = 70) {
  cells <- .simLowSepFacetLabels(cells)
  g <- .simLowSepFirst(gates) |>
    dplyr::select(dplyr::all_of(c("scenario_id", "marker", "gate", "gate_cyt"))) |>
    tidyr::pivot_wider(names_from = "marker", values_from = c("gate", "gate_cyt"))
  lim <- cells |>
    dplyr::group_by(.data$scenario_id) |>
    dplyr::summarise(
      x_max = max(.data$tnf), y_max = max(.data$ifng),
      x_min = min(.data$tnf), y_min = min(.data$ifng), .groups = "drop"
    )
  g <- dplyr::inner_join(g, lim, by = "scenario_id") |>
    dplyr::left_join(dplyr::distinct(cells, .data$scenario_id, .data$scenario_lab), by = "scenario_id")
  # Regions recovered only by the antimode gates (empty if no move).
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
      y = .data$y_min, yend = .data$y_max),
    dplyr::transmute(g, .data$scenario_lab, x = .data$x_min, xend = .data$x_max,
      y = .data$gate_IFNg, yend = .data$gate_IFNg)
  ) |>
    dplyr::mutate(gate_type = "ordinary")
  antimode <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_cyt_TNF, xend = .data$gate_cyt_TNF,
      y = .data$gate_IFNg, yend = .data$y_max),
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_TNF, xend = .data$x_max,
      y = .data$gate_cyt_IFNg, yend = .data$gate_cyt_IFNg)
  ) |>
    dplyr::mutate(gate_type = "antimode")
  # Co-expression gates: lowered TNF (vertical) among cells above the raised
  # IFNg cut, lowered IFNg (horizontal) among cells above the raised TNF cut.
  coex <- .simLowSepFirst(coex_gates) |>
    dplyr::filter(.data$lowered) |>
    dplyr::inner_join(dplyr::select(g, "scenario_id", "scenario_lab", "x_max", "y_max"), by = "scenario_id")
  coexpr <- dplyr::bind_rows(
    dplyr::transmute(dplyr::filter(coex, .data$b == "TNF"), .data$scenario_lab, .data$gate_type,
      x = .data$cut, xend = .data$cut, y = .data$condCut, yend = .data$y_max),
    dplyr::transmute(dplyr::filter(coex, .data$b == "IFNg"), .data$scenario_lab, .data$gate_type,
      x = .data$condCut, xend = .data$x_max, y = .data$cut, yend = .data$cut)
  )
  lines <- dplyr::bind_rows(ordinary, antimode, coexpr)
  # Distinct groups keep coincident lines (unmoved gates) from being merged.
  lines$group <- paste(lines$gate_type, seq_len(nrow(lines)))
  lines$gate_type <- .simLowSepGateFactor(lines$gate_type)
  ggplot(cells, aes(x = .data$tnf, y = .data$ifng)) +
    geom_hex(bins = bins) +
    scale_fill_viridis_c(trans = "log10", name = "Cells") +
    geom_rect(
      data = rect, inherit.aes = FALSE,
      aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = .simLowSepGateColours[["antimode"]], alpha = 0.18
    ) +
    geom_segment(
      data = lines, inherit.aes = FALSE,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend,
        group = .data$group, colour = .data$gate_type, linetype = .data$gate_type),
      linewidth = 0.7
    ) +
    .simLowSepScaleRule() +
    scale_linetype_manual(
      values = c(ordinary = "22", antimode = "solid", coex_screened = "solid", coex_unscreened = "42"),
      labels = .simLowSepGateTypes, name = NULL, drop = TRUE
    ) +
    facet_wrap(~ scenario_lab, ncol = 2, labeller = ggplot2::label_wrap_gen(width = 36)) +
    coord_equal() +
    labs(x = "TNF", y = "IFNγ") +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# IFNg density among all stimulated cells and among TNF+ cells (ordinary TNF
# gate), with the ordinary, antimode and (where lowered) co-expression IFNg
# gates of the first sample of the first dataset.
.simLowSepPlotConditionalDensity <- function(cells, gates, coex_gates, n_grid = 512L) {
  g <- .simLowSepFirst(gates) |>
    dplyr::select(dplyr::all_of(c("scenario_id", "marker", "gate", "gate_cyt"))) |>
    tidyr::pivot_wider(names_from = "marker", values_from = c("gate", "gate_cyt"))
  scen <- .simLowSepScen(cells)
  dens <- cells |>
    dplyr::inner_join(g, by = "scenario_id") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scen))) |>
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
    dplyr::inner_join(dplyr::distinct(dens, .data$scenario_id, .data$scenario_lab), by = "scenario_id") |>
    dplyr::mutate(
      shift = .data$gate_IFNg - .data$gate_cyt_IFNg,
      shift_lab = dplyr::if_else(
        .data$shift > 0,
        paste0("Antimode: lowered by ", formatC(.data$shift, format = "f", digits = 2), "\namong TNF+ cells"),
        "Antimode gate not lowered"
      )
    )
  coex <- .simLowSepFirst(coex_gates) |>
    dplyr::filter(.data$a == "TNF", .data$b == "IFNg", .data$lowered) |>
    dplyr::inner_join(dplyr::distinct(g, .data$scenario_id, .data$scenario_lab), by = "scenario_id")
  lines <- dplyr::bind_rows(
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_IFNg, gate_type = "ordinary"),
    dplyr::transmute(g, .data$scenario_lab, x = .data$gate_cyt_IFNg, gate_type = "antimode"),
    dplyr::transmute(coex, .data$scenario_lab, x = .data$cut, gate_type = .data$gate_type)
  ) |>
    dplyr::mutate(group = seq_len(dplyr::n()), gate_type = .simLowSepGateFactor(.data$gate_type))
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
      colour = .simLowSepGateColours[["antimode"]], linewidth = 0.8,
      arrow = grid::arrow(length = grid::unit(0.2, "cm"))
    ) +
    geom_text(
      data = g, inherit.aes = FALSE,
      aes(x = Inf, y = Inf, label = .data$shift_lab),
      hjust = 1.02, vjust = 1.3, size = 2.7
    ) +
    scale_colour_manual(values = .simLowSepGateColours, labels = c(
      ordinary = "Ordinary IFNγ gate", antimode = "Antimode IFNγ gate",
      coex_screened = "Co-expression IFNγ gate, screened (where lowered)",
      coex_unscreened = "Co-expression IFNγ gate, unscreened (where lowered)"
    ), name = NULL, drop = TRUE) +
    scale_y_sqrt() +
    facet_wrap(~ scenario_lab, ncol = 2, scales = "free_y", labeller = ggplot2::label_wrap_gen(width = 36)) +
    labs(x = "IFNγ", y = "Density (square-root scale)") +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# How far each sample's gates moved below its ordinary gate under each
# lowering rule (zero where not lowered), every sample of every dataset.
.simLowSepPlotGateShift <- function(gates, coex_gates) {
  .simLowSepLowered(gates, coex_gates) |>
    .simLowSepFacetLabels() |>
    dplyr::mutate(
      marker = factor(.data$marker, levels = unname(.simLowSepMarkers), labels = c("TNF", "IFNγ")),
      cell_lab = factor(.simLowSepCellLab(.data$n_cell), levels = .simLowSepCellLab(sort(unique(.data$n_cell)))),
      gate_type = .simLowSepGateFactor(.data$gate_type)
    ) |>
    ggplot(aes(x = .data$cell_lab, y = .data$shift, colour = .data$gate_type)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.7, seed = 1L),
      alpha = 0.5, size = 0.8
    ) +
    .simLowSepScaleRule() +
    facet_grid(marker ~ scenario_lab, labeller = ggplot2::labeller(scenario_lab = ggplot2::label_wrap_gen(width = 18))) +
    labs(x = NULL, y = "Ordinary gate minus lowered gate") +
    .analysis_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
}

# Paired per-sample recovery of true positives and contamination by negative
# cells in the stimulated tube, under the four rules.
.simLowSepPlotRecovery <- function(recovery) {
  recovery |>
    dplyr::filter(.data$tube == "stim") |>
    .simLowSepFacetLabels() |>
    tidyr::pivot_longer(
      c("sensitivity", "fdp"),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      panel = factor(
        paste0(
          ifelse(.data$marker == "TNF", "TNF", "IFNγ"), ": ",
          ifelse(.data$metric == "sensitivity", "true positives recovered", "called positives that are negative")
        ),
        levels = c(
          "TNF: true positives recovered", "IFNγ: true positives recovered",
          "TNF: called positives that are negative", "IFNγ: called positives that are negative"
        )
      ),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      unit = paste(.data$iter, .data$sample)
    ) |>
    ggplot(aes(x = .data$gate_type, y = .data$value)) +
    geom_line(aes(group = .data$unit), colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = .data$gate_type), size = 0.8) +
    .simLowSepScaleRule() +
    scale_x_discrete(labels = .simLowSepGateShort) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    .simLowSepYFloor(c(0, 0.1)) +
    facet_wrap(~ panel + scenario_lab, ncol = 4, scales = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 24, multi_line = FALSE)) +
    labs(x = NULL, y = NULL) +
    .analysis_theme() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 7))
}

# Signed error of background-subtracted frequencies (estimate minus truth),
# paired by sample, under the four rules.
.simLowSepPlotFrequencyError <- function(freq) {
  freq |>
    .simLowSepFacetLabels() |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      error_pp = 100 * .data$error_bs,
      unit = paste(.data$iter, .data$sample)
    ) |>
    ggplot(aes(x = .data$gate_type, y = .data$error_pp)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(aes(group = .data$unit), colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = .data$gate_type), size = 0.8) +
    .simLowSepScaleRule() +
    scale_x_discrete(labels = .simLowSepGateShort) +
    facet_grid(scenario_lab ~ quantity_lab, labeller = ggplot2::labeller(
      scenario_lab = ggplot2::label_wrap_gen(width = 18),
      quantity_lab = ggplot2::label_wrap_gen(width = 12)
    )) +
    labs(x = NULL, y = "Estimated minus true frequency (percentage points)") +
    .analysis_theme() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6))
}

# Percentiles over every sample of a scenario (10th to 90th as a bar, median
# as a point) by scenario, quantity and rule. With `truth = TRUE` the
# median true value is marked with a black cross (frequencies); otherwise
# the values are proportions in [0, 1] (F1). `scale` multiplies the values
# for display (100 for percentages).
.simLowSepPlotPercentiles <- function(pct, x_lab, truth = FALSE, scale = 1) {
  pct <- pct |>
    .simLowSepFacetLabels() |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      scenario_lab = factor(.data$scenario_lab, levels = rev(levels(.data$scenario_lab))),
      dplyr::across(c("p10", "p50", "p90", "truth"), ~ scale * .x)
    )
  dodge <- position_dodge(width = 0.75)
  p <- ggplot(pct, aes(y = .data$scenario_lab, colour = .data$gate_type)) +
    geom_linerange(aes(xmin = .data$p10, xmax = .data$p90), position = dodge, linewidth = 0.5) +
    geom_point(aes(x = .data$p50), position = dodge, size = 1.2)
  if (isTRUE(truth)) {
    p <- p + geom_point(
      data = dplyr::distinct(pct, .data$scenario_lab, .data$quantity_lab, .data$truth),
      aes(x = .data$truth, y = .data$scenario_lab), inherit.aes = FALSE,
      shape = 4, size = 2.2, colour = "black"
    )
  } else {
    p <- p + .simLowSepXFloor(c(0, 1))
  }
  p +
    .simLowSepScaleRule() +
    facet_wrap(~ quantity_lab, nrow = 1, scales = if (isTRUE(truth)) "free_x" else "fixed") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 40)) +
    labs(x = x_lab, y = NULL) +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# Rejection rates of one family's scenarios (rows, with cell counts) by test
# (panel rows) and quantity (panel columns). Filled points are power (the
# quantity has a net response), open points false-positive rates; the
# dashed line is the 5% level. Monte Carlo intervals are marked as such.
.simLowSepPlotRates <- function(rates) {
  rates <- rates |>
    .simLowSepFacetLabels(cells = TRUE) |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      scenario_lab = factor(.data$scenario_lab, levels = rev(levels(.data$scenario_lab))),
      test = factor(.data$test, levels = names(.simLowSepTests), labels = .simLowSepTests),
      response_truth = factor(.data$response_truth, levels = c("response", "none"))
    )
  dodge <- position_dodge(width = 0.75)
  intervals <- dplyr::filter(rates, is.finite(.data$rate_lower), is.finite(.data$rate_upper))
  ggplot(rates, aes(y = .data$scenario_lab, colour = .data$gate_type)) +
    geom_vline(xintercept = 0.05, linetype = "22", colour = "grey50") +
    .analysis_mcse_layer(geom_linerange(
      data = intervals, aes(xmin = .data$rate_lower, xmax = .data$rate_upper),
      position = dodge, linewidth = 0.4, alpha = 0.7
    )) +
    geom_point(aes(x = .data$rate, shape = .data$response_truth), position = dodge, size = 1.4) +
    scale_shape_manual(values = c(response = 16, none = 1), labels = c(
      response = "Net response (power)", none = "No net response (false-positive rate)"
    ), name = NULL, drop = FALSE) +
    .simLowSepScaleRule() +
    .simLowSepXFloor(c(0, 0.1)) +
    scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
    facet_grid(test ~ quantity_lab, labeller = ggplot2::labeller(
      test = ggplot2::label_wrap_gen(width = 20), quantity_lab = ggplot2::label_wrap_gen(width = 12)
    )) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 36)) +
    labs(x = "Tests rejecting at the 5% level", y = NULL) +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical", axis.text.x = element_text(size = 7))
}

# Cells each lowering rule adds over the ordinary gates in the stimulated
# tube (per 100,000 cells, mean over samples), split into cells truly
# positive and truly negative for the marker, with the share of samples
# whose gate was lowered printed above each bar.
.simLowSepPlotAdded <- function(gating) {
  long <- gating |>
    .simLowSepFacetLabels(cells = TRUE) |>
    dplyr::mutate(
      gate_type = .simLowSepGateFactor(.data$gate_type),
      marker = factor(.data$marker, levels = unname(.simLowSepMarkers), labels = c("TNF", "IFNγ"))
    )
  bars <- long |>
    tidyr::pivot_longer(c("added_true_per_1e5_stim", "added_negative_per_1e5_stim"),
      names_to = "class", values_to = "added") |>
    dplyr::mutate(class = factor(.data$class,
      levels = c("added_negative_per_1e5_stim", "added_true_per_1e5_stim"),
      labels = c("Truly negative for the marker", "Truly positive for the marker")))
  ggplot(bars, aes(x = .data$gate_type, y = .data$added)) +
    geom_col(aes(fill = .data$class), width = 0.7) +
    geom_text(data = long, aes(y = .data$added_per_1e5_stim,
      label = paste0(round(100 * .data$share_lowered), "%")), vjust = -0.4, size = 2.4) +
    scale_fill_manual(values = c("#A6611A", "#018571"), name = NULL) +
    scale_x_discrete(labels = .simLowSepGateShort) +
    .simLowSepYFloor(c(0, 10)) +
    facet_grid(marker ~ scenario_lab, scales = "free_y",
      labeller = ggplot2::labeller(scenario_lab = ggplot2::label_wrap_gen(width = 18))) +
    labs(x = NULL, y = "Cells added per 100,000 stimulated cells") +
    .analysis_theme() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6))
}

# Analysis 11b: figures against the mismatch -----------------------------

# Rows of one mismatch section with the x-axis factor `mismatch_x`. The SD
# inflation section includes the no-mismatch rows (mean shift 0) as its
# reference.
.simLowSepMismatchSection <- function(tbl, type) {
  if (identical(type, "mean_shift")) {
    out <- dplyr::filter(tbl, .data$mismatch_type == "mean_shift")
    lev <- sort(unique(out$stim_mean_shift))
    out$mismatch_x <- factor(ifelse(out$stim_mean_shift == 0, "None", as.character(out$stim_mean_shift)),
      levels = ifelse(lev == 0, "None", as.character(lev)))
  } else if (identical(type, "sd_inflation")) {
    out <- dplyr::filter(tbl, .data$mismatch_type == "sd_inflation" |
      (.data$mismatch_type == "mean_shift" & .data$stim_mean_shift == 0))
    out$mismatch_x <- factor(ifelse(out$mismatch_type == "sd_inflation",
      paste0("SD +", round(100 * (out$stim_sd_multiplier - 1)), "%"), "None"))
    if ("None" %in% levels(out$mismatch_x)) out$mismatch_x <- stats::relevel(out$mismatch_x, "None")
  } else {
    stop("type must be \"mean_shift\" or \"sd_inflation\".")
  }
  out
}

.simLowSepMismatchXLab <- c(
  mean_shift = "Mean shift of the stimulated TNF- IFNγ- cells (both markers)",
  sd_inflation = "SD inflation of the stimulated TNF- IFNγ- cells (both markers)"
)

# Rejection rates of one test against the mismatch, by baseline scenario
# (rows) and quantity (columns); lines join the mismatch settings of a rule.
.simLowSepPlotMismatchRates <- function(rates, type) {
  rates <- rates |>
    .simLowSepFacetLabels(mismatch = FALSE, cells = TRUE) |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      response_truth = factor(.data$response_truth, levels = c("response", "none"))
    )
  dodge <- position_dodge(width = 0.4)
  intervals <- dplyr::filter(rates, is.finite(.data$rate_lower), is.finite(.data$rate_upper))
  ggplot(rates, aes(x = .data$mismatch_x, y = .data$rate, colour = .data$gate_type, group = .data$gate_type)) +
    geom_hline(yintercept = 0.05, linetype = "22", colour = "grey50") +
    .analysis_mcse_layer(geom_linerange(
      data = intervals, aes(ymin = .data$rate_lower, ymax = .data$rate_upper),
      position = dodge, linewidth = 0.4, alpha = 0.7
    )) +
    geom_line(position = dodge, linewidth = 0.4) +
    geom_point(aes(shape = .data$response_truth), position = dodge, size = 1.4) +
    scale_shape_manual(values = c(response = 16, none = 1), labels = c(
      response = "Net response (power)", none = "No net response (false-positive rate)"
    ), name = NULL, drop = FALSE) +
    .simLowSepScaleRule() +
    .simLowSepYFloor(c(0, 0.1)) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    facet_grid(scenario_lab ~ quantity_lab, labeller = ggplot2::labeller(
      scenario_lab = ggplot2::label_wrap_gen(width = 20), quantity_lab = ggplot2::label_wrap_gen(width = 12)
    )) +
    labs(x = .simLowSepMismatchXLab[[type]], y = "Tests rejecting at the 5% level") +
    .analysis_theme() +
    theme(legend.position = "bottom", legend.box = "vertical", axis.text.x = element_text(size = 7, angle = 30, hjust = 1))
}

# Percentiles of the signed frequency error (percentage points) against the
# mismatch: median as a point, 10th to 90th percentiles as a bar, every
# sample of every dataset pooled.
.simLowSepPlotMismatchError <- function(pct, type) {
  pct |>
    .simLowSepFacetLabels(mismatch = FALSE, cells = TRUE) |>
    dplyr::mutate(
      quantity_lab = factor(.data$quantity_lab, levels = .simLowSepQuantities$quantity_lab),
      gate_type = .simLowSepGateFactor(.data$gate_type),
      dplyr::across(c("p10", "p50", "p90"), ~ 100 * .x)
    ) |>
    ggplot(aes(x = .data$mismatch_x, colour = .data$gate_type, group = .data$gate_type)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_linerange(aes(ymin = .data$p10, ymax = .data$p90), position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(aes(y = .data$p50), position = position_dodge(width = 0.6), size = 1.2) +
    .simLowSepScaleRule() +
    facet_wrap(~ scenario_lab + quantity_lab, ncol = 5, scales = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 22, multi_line = TRUE)) +
    labs(x = .simLowSepMismatchXLab[[type]], y = "Estimated minus true frequency (percentage points)") +
    .analysis_theme() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
      strip.text = element_text(size = 7))
}

# What the lowering rules gate against the mismatch: share of sample-marker
# gates lowered, and truly negative and truly positive cells added per
# 100,000 stimulated cells (summed over the two markers, so a cell newly
# positive for both counts twice).
.simLowSepPlotMismatchGating <- function(gating, type) {
  scen <- c(.simLowSepScen(gating), "mismatch_x")
  metric_labs <- c(
    share_lowered = "Sample gates lowered",
    added_negative = "Truly negative cells added per 100,000",
    added_true = "Truly positive cells added per 100,000"
  )
  gating |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(scen, "gate_type")))) |>
    dplyr::summarise(
      share_lowered = sum(.data$n_lowered) / sum(.data$n_sample),
      added_negative = sum(.data$added_negative_per_1e5_stim),
      added_true = sum(.data$added_true_per_1e5_stim),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(names(metric_labs), names_to = "metric", values_to = "value") |>
    .simLowSepFacetLabels(mismatch = FALSE, cells = TRUE) |>
    dplyr::mutate(
      metric = factor(.data$metric, levels = names(metric_labs), labels = metric_labs),
      gate_type = .simLowSepGateFactor(.data$gate_type)
    ) |>
    ggplot(aes(x = .data$mismatch_x, y = .data$value, colour = .data$gate_type, group = .data$gate_type)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.2) +
    .simLowSepScaleRule() +
    .simLowSepYFloor(c(0, 0.1)) +
    facet_wrap(~ scenario_lab + metric, ncol = 3, scales = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 26, multi_line = TRUE)) +
    labs(x = .simLowSepMismatchXLab[[type]], y = NULL) +
    .analysis_theme() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
      strip.text = element_text(size = 7))
}

# A blank layer that trains the x range without censoring data.
.simLowSepXFloor <- function(range) {
  geom_blank(data = data.frame(.x_floor = range), aes(x = .data$.x_floor), inherit.aes = FALSE)
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

# Save a table as CSV; list columns (bootstrap draws) are dropped.
.simLowSepWriteTable <- function(tbl, name, fig_key, path_root = NULL) {
  tbl <- tbl[, !vapply(tbl, is.list, logical(1)), drop = FALSE]
  path <- file.path(.simLowSepTableDir(fig_key, path_root), paste0(name, ".csv"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tbl, path, row.names = FALSE)
  invisible(path)
}

# Summaries computed at render time from the collated tables, shared by
# 11a and 11b. `mcse` controls the rejection-rate intervals.
.simLowSepSummaries <- function(results, settings, mcse = TRUE) {
  freq <- .simLowSepFrequencies(results$counts)
  expected <- .simLowSepExpectedNet(results$scenarios, settings$splits)
  sample_tests <- .simLowSepSampleTests(freq)
  dataset_tests <- .simLowSepDatasetTests(freq)
  f1 <- .simLowSepF1(results$counts)
  list(
    recovery = .simLowSepRecovery(results$counts),
    freq = freq,
    freq_summary = .simLowSepFrequencySummary(freq),
    f1 = f1,
    f1_pct = .simLowSepPercentiles(f1, "f1"),
    freq_pct = .simLowSepPercentiles(freq, "prop_bs_est", truth = "prop_bs_truth"),
    error_pct = .simLowSepPercentiles(freq, "error_bs", truth = "prop_bs_truth"),
    expected_net = expected,
    sample_tests = sample_tests,
    dataset_tests = dataset_tests,
    rates = .simLowSepRejectionRates(sample_tests, dataset_tests, expected, mcse = mcse),
    gating = .simLowSepGatingSummary(results$gates, results$coex_gates, results$counts)
  )
}

# Run plumbing ----------------------------------------------------------------

# Names of the collated tables saved at promotion (`.simLowSepCollate()`).
.simLowSepCollatedNames <- c("scenarios", "gates", "coex_gates", "counts", "cells")

# Read the promoted collated tables of a run (or read-only results) context.
.simLowSepReadResults <- function(run_ctx, required_params) {
  stats::setNames(lapply(.simLowSepCollatedNames, function(nm) {
    .analysis_read_current(
      run_ctx = run_ctx,
      relative_path = c("collated", paste0(nm, ".rds")),
      required_params = required_params
    )
  }), .simLowSepCollatedNames)
}

# Rows of this chunk: the selected grid shuffled with a fixed seed and dealt
# round-robin into `n_chunks` chunks (as in Analyses 2a and 2b).
.simLowSepChunkRows <- function(grid, shuffle_seed, n_chunks, chunk_index) {
  .analysis_with_seed(shuffle_seed, dplyr::slice_sample(grid, prop = 1)) |>
    dplyr::mutate(.sim_grid_chunk = ((dplyr::row_number() - 1L) %% n_chunks) + 1L) |>
    dplyr::filter(.data$.sim_grid_chunk == chunk_index) |>
    dplyr::select(-".sim_grid_chunk")
}
