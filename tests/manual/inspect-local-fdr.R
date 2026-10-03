# Manual inspection of the baseline local-FDR gating path
#
# Run from the repository root after:
#   devtools::load_all()
#
# Edit these values to compare controlled scenarios, then source this file again.
scenario <- list(
  seed = 364L,
  n_cell = 5000L,
  response_prop = 0.08,
  separation = 4,
  background_relative_to_response = 0,
  bandwidth = NULL
)

run_local_fdr_inspection <- function(scenario) {
  if (!exists("gateStim", mode = "function")) {
    stop("Run devtools::load_all() before sourcing this script.")
  }
  if (!requireNamespace("simcyto", quietly = TRUE)) {
    stop(
      "The manual inspection script requires the development dependency ",
      "simcyto."
    )
  }

  old_intermediate <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_)
  on.exit(
    if (is.na(old_intermediate)) {
      Sys.unsetenv("STIMGATE_INTERMEDIATE")
    } else {
      Sys.setenv(STIMGATE_INTERMEDIATE = old_intermediate)
    },
    add = TRUE
  )
  Sys.setenv(STIMGATE_INTERMEDIATE = "all")

  set.seed(scenario$seed)

  prob_background <- scenario$response_prop *
    scenario$background_relative_to_response

  sim <- simcyto::simCytExperiment(
    nSample = 1L,
    nMarker = 1L,
    nCondition = 2L,
    nCluster = 2L,
    nCellByCondition = rep(scenario$n_cell, 2L),
    transformationFunc = simcyto::simCytTransformGaussian(),
    mixtureType = "gaussianOnly",
    meanExprMat = matrix(
      c(0, scenario$separation),
      byrow = TRUE,
      ncol = 1L
    ),
    clusterLabelVec = c("gn", "gp"),
    probVecUns = c(1 - prob_background, prob_background),
    probExact = TRUE,
    probResponseVecByStimCondition = list(
      c(-scenario$response_prop, scenario$response_prop)
    ),
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    covEvMin = 1,
    covEvMax = 1
  )

  fs <- as(sim$flowFrameList, "flowSet")
  gs <- flowWorkspace::GatingSet(fs)

  path_project <- file.path(tempdir(), "stimgate-local-fdr-manual")
  if (dir.exists(path_project)) {
    unlink(path_project, recursive = TRUE, force = TRUE)
  }

  invisible(gateStim(
    .data = gs,
    pathProject = path_project,
    popGate = "root",
    batchList = list(c(1L, 2L)),
    marker = "MarkerF1",
    calcCytPosGates = FALSE,
    tolClust = NULL,
    biasUns = 0,
    bw = scenario$bandwidth,
    gateCombn = "no",
    locEnforceShapeThreshold = FALSE
  ))

  path_ind <- file.path(
    path_project,
    "intermediateData",
    "init",
    "F1",
    "ind",
    "2"
  )

  read_optional <- function(name) {
    path <- file.path(path_ind, paste0(name, ".rds"))
    if (file.exists(path)) readRDS(path) else NULL
  }

  scalar <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0L) NA_real_ else x[[1L]]
  }

  ex_stim <- read_optional("exTblStimThreshold")
  ex_uns <- read_optional("exTblUnsThreshold")
  dens_raw <- read_optional("densTblRaw")
  prob_tables <- read_optional("probTblList")
  data_mod <- read_optional("dataMod")
  trim_info <- read_optional("dataModTrimInfo")
  data_trim <- read_optional("dataModTrim")
  data_threshold <- read_optional("dataThreshold")
  detail <- read_optional("locDetailCondition")
  if (is.null(detail)) {
    detail <- read_optional("locDetailSample")
  }
  cp_ind <- read_optional("cpInd")

  bandwidth <- read_optional("bwCpUnsLoc")
  if (is.null(bandwidth)) {
    bandwidth <- read_optional("bwCpUnsLocAdaptive")
  }

  shape_info <- list(
    requested = attr(data_mod, "locShapeThresholdRequested"),
    applied = attr(data_mod, "locShapeThresholdApplied"),
    threshold_x = attr(data_mod, "locShapeThresholdX"),
    tailgate_x = attr(data_mod, "locShapeTailgateX"),
    antimode_x = attr(data_mod, "locShapeAntimodeX"),
    detail = attr(data_mod, "locShapeThresholdInfo")
  )

  final_decisions <- if (is.list(trim_info)) trim_info$final else NULL
  if (!is.list(final_decisions)) {
    final_decisions <- list()
  }

  returned_threshold <- if (is.data.frame(detail) && nrow(detail) > 0L) {
    scalar(detail$threshold)
  } else {
    scalar(cp_ind)
  }
  returned_generated <- is.data.frame(detail) &&
    nrow(detail) > 0L &&
    detail$locGenerated[[1L]] %in% TRUE
  returned_label <- if (returned_generated) {
    "final local-FDR gate"
  } else {
    "returned fallback threshold"
  }

  decision_table <- data.frame(
    decision = c(
      "preliminary modelling lower bound",
      "clear-response initial reference",
      "density-dominance reference",
      "clear-response boundary",
      "quality boundary",
      "antimode boundary",
      "final filtering boundary",
      returned_label
    ),
    x = c(
      scalar(attr(data_mod, "minProbXPos")),
      scalar(final_decisions$xClearInit),
      scalar(final_decisions$xDom),
      scalar(final_decisions$xClear),
      scalar(final_decisions$xQual),
      scalar(final_decisions$xAntimode),
      scalar(final_decisions$xSum),
      returned_threshold
    )
  )

  cat("\nScenario\n")
  print(scenario)
  cat("\nDensity bandwidth\n")
  print(bandwidth)
  cat("\nShape/lower-bound information\n")
  print(shape_info)
  cat("\nFiltering/reference decisions\n")
  print(decision_table, row.names = FALSE)
  cat("\nFiltering details\n")
  if (is.list(trim_info)) {
    print(trim_info[c("reason", "clear", "marginal", "antimode", "final")])
  } else {
    cat("Post-smoothing filtering was not reached.\n")
  }
  cat("\nReturned threshold and provenance\n")
  if (is.data.frame(detail) && nrow(detail) > 0L) {
    print(detail)
  } else {
    print(list(
      threshold = returned_threshold,
      locGenerated = FALSE,
      locReason = "Detailed threshold provenance was not saved"
    ))
  }
  cat("\nIntermediate directory\n", path_ind, "\n", sep = "")

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

  if (is.data.frame(dens_raw) && nrow(dens_raw) > 0L) {
    dens_stim <- dens_raw[dens_raw$stim == "yes", , drop = FALSE]
    dens_uns <- dens_raw[dens_raw$stim == "no", , drop = FALSE]
    dens_stim <- dens_stim[order(dens_stim$xStim), , drop = FALSE]
    dens_uns <- dens_uns[order(dens_uns$xStim), , drop = FALSE]

    graphics::plot(
      dens_stim$xStim,
      dens_stim$dens,
      type = "l",
      xlab = "F1 expression",
      ylab = "Density",
      main = "Stimulated and unstimulated densities"
    )
    graphics::lines(dens_uns$xStim, dens_uns$dens, lty = 2)
    if (is.finite(returned_threshold)) {
      graphics::abline(v = returned_threshold, lty = 3, lwd = 2)
    }
    graphics::legend(
      "topright",
      legend = c("stimulated", "unstimulated", returned_label),
      lty = c(1, 2, 3),
      lwd = c(1, 1, 2),
      bty = "n"
    )
  } else {
    graphics::plot.new()
    graphics::title("Density estimation was not reached")
  }

  if (
    is.list(prob_tables) &&
      is.data.frame(prob_tables$all) &&
      is.data.frame(prob_tables$pos) &&
      is.data.frame(data_mod) &&
      all(c("F1", "pred") %in% names(data_mod))
  ) {
    raw_prob <- prob_tables$all[
      order(prob_tables$all$xStim),
      ,
      drop = FALSE
    ]
    model_prob <- prob_tables$pos[
      order(prob_tables$pos$xStim),
      ,
      drop = FALSE
    ]
    data_mod <- data_mod[order(data_mod$F1), , drop = FALSE]

    graphics::plot(
      raw_prob$xStim,
      raw_prob$probStimNorm,
      type = "l",
      ylim = c(0, 1),
      xlab = "F1 expression",
      ylab = "Response probability",
      main = "Raw, modelled and smoothed response probability"
    )
    graphics::points(
      model_prob$xStim,
      model_prob$probStimNorm,
      pch = 16,
      cex = 0.45
    )
    graphics::lines(data_mod$F1, data_mod$pred, lty = 2, lwd = 2)

    plot_decision_names <- c(
      "preliminary modelling lower bound",
      "clear-response boundary",
      "quality boundary",
      "antimode boundary",
      "final filtering boundary",
      returned_label
    )
    plot_decisions <- decision_table[
      match(plot_decision_names, decision_table$decision),
      ,
      drop = FALSE
    ]
    plot_decisions$lty <- c(3, 4, 5, 6, 2, 1)
    plot_decisions$label <- c(
      "preliminary bound",
      "x_clear",
      "x_qual",
      "x_antimode",
      "x_sum",
      returned_label
    )
    plot_decisions <- plot_decisions[
      is.finite(plot_decisions$x),
      ,
      drop = FALSE
    ]

    if (nrow(plot_decisions) > 0L) {
      for (i in seq_len(nrow(plot_decisions))) {
        graphics::abline(
          v = plot_decisions$x[[i]],
          lty = plot_decisions$lty[[i]],
          lwd = if (plot_decisions$decision[[i]] == returned_label) 2 else 1
        )
      }
    }
    graphics::legend(
      "bottomright",
      legend = c("raw probability", "preliminary model region", "smoothed curve"),
      lty = c(1, NA, 2),
      pch = c(NA, 16, NA),
      bty = "n",
      cex = 0.8
    )
    if (nrow(plot_decisions) > 0L) {
      graphics::legend(
        "topleft",
        legend = plot_decisions$label,
        lty = plot_decisions$lty,
        bty = "n",
        cex = 0.7
      )
    }
  } else {
    graphics::plot.new()
    graphics::title("Response-probability smoothing was not reached")
  }

  if (is.data.frame(data_threshold) && nrow(data_threshold) > 0L) {
    data_threshold <- data_threshold[order(data_threshold$F1), , drop = FALSE]
    graphics::plot(
      data_threshold$F1,
      data_threshold$propBsDiff,
      type = "l",
      xlab = "Candidate F1 threshold",
      ylab = "Observed - estimated response proportion",
      main = "Final threshold matching"
    )
    graphics::abline(h = 0, lty = 2)
    if (is.finite(returned_threshold)) {
      graphics::abline(v = returned_threshold, lty = 3, lwd = 2)
    }
  } else {
    graphics::plot.new()
    graphics::title("Threshold matching was not reached")
  }

  list(
    scenario = scenario,
    path_project = path_project,
    path_intermediate = path_ind,
    expression = list(stim = ex_stim, unstim = ex_uns),
    density = dens_raw,
    bandwidth = bandwidth,
    probabilities = prob_tables,
    smoothed_probability = data_mod,
    filtered_probability = data_trim,
    filtering = trim_info,
    decisions = decision_table,
    threshold_candidates = data_threshold,
    threshold_detail = detail,
    shape = shape_info
  )
}

inspection <- run_local_fdr_inspection(scenario)
