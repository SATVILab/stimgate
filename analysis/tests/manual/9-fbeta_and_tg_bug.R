# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------

# IMPORTANT: ACS batches contain four stimulated files followed by one
# unstimulated file. The current package comparator assumes that the first
# file is unstimulated. This standalone diagnostic uses the correct ordering
# and stops if the filenames do not confirm it.

library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)

root_dir <- projr::projr_path_get("project")

devtools::load_all(root_dir)

scripts_r_dir <- file.path(root_dir, "scripts", "r")

script_vec <- c(
  "analysis-runtime.R",
  "acs_cytof-helper.R",
  "acs_cytof-preprocess.R",
  "acs_cytof-gate.R",
  "sim-misc.R",
  "sim-compare-freq_bs.R",
  "acs_cytof-methods.R"
)

purrr::walk(
  file.path(scripts_r_dir, script_vec),
  source
)

path_fcs_base <- projr::projr_path_get(
  "raw-data-large",
  "comparison_data",
  "acscytof",
  "fcs"
)

path_gs_base <- projr::projr_path_get(
  "cache",
  "acs_cytof",
  "gs"
)

path_scratch_base <- projr::projr_path_get(
  "cache",
  "acs_cytof",
  "scratch"
)

path_fbeta <- file.path(
  root_dir,
  "scripts",
  "python",
  "fbeta.py"
)

stopifnot(
  dir.exists(path_fcs_base),
  dir.exists(path_gs_base),
  file.exists(path_fbeta)
)

reticulate::py_config()


# second tranche
.tryAcsMethod <- function(code) {
  tryCatch(
    list(
      success = TRUE,
      result = force(code),
      error = NA_character_
    ),
    error = function(e) {
      list(
        success = FALSE,
        result = NULL,
        error = conditionMessage(e)
      )
    }
  )
}

.firstFinite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  x[[1L]]
}

.countAbove <- function(x, threshold) {
  if (!is.finite(threshold)) {
    return(NA_integer_)
  }

  as.integer(sum(x > threshold, na.rm = TRUE))
}

.isAcsUnstimFile <- function(path) {
  grepl(
    "(^|[^[:alnum:]])uns([^[:alnum:]]|$)",
    tolower(basename(path))
  )
}

.acsBatchLayout <- function(batchList, fcsFiles) {
  layout <- purrr::imap_dfr(
    batchList,
    function(indBatch, batchIndex) {
      tibble(
        batch = as.integer(batchIndex),
        position = seq_along(indBatch),
        ind = indBatch,
        expectedRole = if_else(
          seq_along(indBatch) == length(indBatch),
          "unstimulated",
          "stimulated"
        ),
        sample = basename(fcsFiles[indBatch]),
        filenameLooksUnstimulated = .isAcsUnstimFile(fcsFiles[indBatch])
      )
    }
  )

  invalidBatch <- layout |>
    group_by(.data$batch) |>
    summarise(
      nFiles = n(),
      finalFileIsUnstimulated = last(.data$filenameLooksUnstimulated),
      earlierFileLooksUnstimulated = any(
        head(.data$filenameLooksUnstimulated, -1L)
      ),
      .groups = "drop"
    ) |>
    filter(
      .data$nFiles != 5L |
        !.data$finalFileIsUnstimulated |
        .data$earlierFileLooksUnstimulated
    )

  if (nrow(invalidBatch) > 0L) {
    print(
      layout |>
        filter(.data$batch %in% invalidBatch$batch),
      n = Inf,
      width = Inf
    )
    stop(
      "The FCS ordering does not match four stimulated files followed by ",
      "one unstimulated file in every five-file batch."
    )
  }

  layout
}

.inspectAcsComparatorCase <- function(
  pop = "tcrgd",
  batchIndex = 1L,
  stimOffset = 1L,
  chnl = "Ho165Di",
  outputGroup = NULL,
  showTailgatePackagePlot = FALSE,
  pathFcsBase = path_fcs_base,
  pathGsBase = path_gs_base,
  pathScratchBase = path_scratch_base,
  pathFbeta = path_fbeta
) {
  channelMap <- .acsCytofChannelMap()

  if (!chnl %in% names(channelMap)) {
    stop(
      "chnl must be one of: ",
      paste(names(channelMap), collapse = ", ")
    )
  }

  paths <- .acsCytofPopulationPaths(
    pop = pop,
    pathFcsBase = pathFcsBase,
    pathGsBase = pathGsBase,
    pathScratchBase = pathScratchBase,
    outputGroup = outputGroup
  )

  if (!dir.exists(paths$gs)) {
    stop("GatingSet not found at: ", paths$gs)
  }

  gs <- flowWorkspace::load_gs(paths$gs)
  batchList <- lapply(
    seq.int(1L, length(gs), by = 5L),
    \(indStart) seq.int(indStart, length.out = 5L)
  )
  fcsFiles <- .acsCytofFcsFiles(paths$fcs)
  fcsFiles <- fcsFiles[seq_len(length(gs))]
  batchLayout <- .acsBatchLayout(
    batchList = batchList,
    fcsFiles = fcsFiles
  )

  if (batchIndex < 1L || batchIndex > length(batchList)) {
    stop("batchIndex is outside the available batches.")
  }

  if (stimOffset < 1L || stimOffset > 4L) {
    stop("stimOffset must be between 1 and 4.")
  }

  indBatch <- batchList[[batchIndex]]
  if (length(indBatch) != 5L) {
    stop("Each ACS batch must contain exactly five samples.")
  }

  # ACS ordering is four stimulated samples followed by the unstimulated
  # sample. This differs from the old comparator assumption that the first
  # sample in each batch was unstimulated.
  indStimVec <- indBatch[seq_len(4L)]
  indUns <- indBatch[[5L]]
  indStim <- indStimVec[[stimOffset]]

  cat("\nSelected five-file batch layout\n")
  print(
    batchLayout |>
      filter(.data$batch == .env$batchIndex),
    n = Inf,
    width = Inf
  )

  xUns <- .acsCytofExpressionMatrix(
    gs = gs,
    ind = indUns,
    channels = chnl
  )[, 1L]

  xStim <- .acsCytofExpressionMatrix(
    gs = gs,
    ind = indStim,
    channels = chnl
  )[, 1L]

  xUnsFinite <- xUns[is.finite(xUns)]
  xStimFinite <- xStim[is.finite(xStim)]

  sampleInfo <- tibble(
    pop = pop,
    batch = batchIndex,
    stimulusOffset = stimOffset,
    cytokine = unname(channelMap[[chnl]]),
    chnl = chnl,
    indUns = indUns,
    sampleUns = basename(fcsFiles[[indUns]]),
    indStim = indStim,
    sampleStim = basename(fcsFiles[[indStim]])
  )

  print(sampleInfo, width = Inf)

  eventSummary <- bind_rows(
    tibble(
      condition = "unstimulated",
      n = length(xUnsFinite),
      minimum = min(xUnsFinite),
      q01 = quantile(xUnsFinite, 0.01, names = FALSE),
      q25 = quantile(xUnsFinite, 0.25, names = FALSE),
      median = median(xUnsFinite),
      q75 = quantile(xUnsFinite, 0.75, names = FALSE),
      q99 = quantile(xUnsFinite, 0.99, names = FALSE),
      maximum = max(xUnsFinite)
    ),
    tibble(
      condition = "stimulated",
      n = length(xStimFinite),
      minimum = min(xStimFinite),
      q01 = quantile(xStimFinite, 0.01, names = FALSE),
      q25 = quantile(xStimFinite, 0.25, names = FALSE),
      median = median(xStimFinite),
      q75 = quantile(xStimFinite, 0.75, names = FALSE),
      q99 = quantile(xStimFinite, 0.99, names = FALSE),
      maximum = max(xStimFinite)
    )
  )

  print(eventSummary, width = Inf)

  # -----------------------------------------------------------------------
  # Run exactly the calls used by the ACS comparator
  # -----------------------------------------------------------------------

  fbetaSettings <- .acsCytofComparatorSettings("fbeta")
  tailgateSettings <- .acsCytofComparatorSettings("tailgate")

  fbetaDirect <- .tryAcsMethod(
    do.call(
      .simCompareFbetaThreshold,
      c(
        list(
          xUns = xUns,
          xStim = xStim,
          pathFbeta = pathFbeta
        ),
        fbetaSettings$params
      )
    )
  )

  tailgateDirect <- .tryAcsMethod(
    .simCompareTailgateThreshold(
      x = xStim,
      adjust = tailgateSettings$params$adjust,
      bandwidth = tailgateSettings$params$bandwidth,
      numPeaks = tailgateSettings$params$numPeaks,
      refPeak = tailgateSettings$params$refPeak,
      method = tailgateSettings$params$derivativeMethod,
      tol = tailgateSettings$params$tol,
      side = tailgateSettings$params$side,
      strict = tailgateSettings$params$strict,
      autoTol = tailgateSettings$params$autoTol
    )
  )

  # These are the wrapped calls used in .acsCytofRunComparator().
  # Any error is replaced here by a high-value fallback threshold.
  fbetaWrapped <- .acsCytofThresholdOne(
    method = "fbeta",
    xUns = xUns,
    xStim = xStim,
    settings = fbetaSettings,
    pathFbeta = pathFbeta
  )

  tailgateWrapped <- .acsCytofThresholdOne(
    method = "tailgate",
    xUns = xUns,
    xStim = xStim,
    settings = tailgateSettings
  )

  thresholdSummary <- bind_rows(
    mutate(fbetaWrapped, method = "fbeta", .before = 1L),
    mutate(tailgateWrapped, method = "tailgate", .before = 1L)
  ) |>
    mutate(
      nPosStim = map_int(
        .data$threshold,
        \(threshold) .countAbove(xStim, threshold)
      ),
      nPosUns = map_int(
        .data$threshold,
        \(threshold) .countAbove(xUns, threshold)
      ),
      propStim = .data$nPosStim / length(xStim),
      propUns = .data$nPosUns / length(xUns),
      propBs = .data$propStim - .data$propUns,
      thresholdAboveStimMaximum = .data$threshold >= max(xStimFinite)
    )

  cat("\nThresholds and resulting counts\n")
  print(thresholdSummary, width = Inf)

  if (!fbetaDirect$success) {
    message("\nDIRECT F-BETA ERROR: ", fbetaDirect$error)
  }

  if (!tailgateDirect$success) {
    message("\nDIRECT TAILGATE ERROR: ", tailgateDirect$error)
  }

  # -----------------------------------------------------------------------
  # F-beta diagnostics
  # -----------------------------------------------------------------------

  # fbeta.py constructs the bin edges from the unstimulated distribution,
  # then uses those same edges for the stimulated histogram.
  fbetaRangeCheck <- tibble(
    nBins = floor(sqrt(max(length(xUnsFinite), length(xStimFinite)))),
    minimumUns = min(xUnsFinite),
    maximumUns = max(xUnsFinite),
    minimumStim = min(xStimFinite),
    maximumStim = max(xStimFinite),
    nStimBelowUnsRange = sum(xStimFinite < min(xUnsFinite)),
    nStimAboveUnsRange = sum(xStimFinite > max(xUnsFinite)),
    propStimOutsideUnsRange = mean(
      xStimFinite < min(xUnsFinite) |
        xStimFinite > max(xUnsFinite)
    )
  )

  cat("\nF-beta histogram range check\n")
  print(fbetaRangeCheck, width = Inf)

  fbetaCurve <- NULL
  pFbetaDensity <- NULL
  pFbetaScore <- NULL

  if (fbetaDirect$success) {
    fbetaObject <- fbetaDirect$result$fbeta

    fbetaCurve <- tibble(
      x = as.numeric(fbetaObject[["pdfx"]]),
      pdfUns = as.numeric(fbetaObject[["pdfneg"]]),
      pdfStim = as.numeric(fbetaObject[["pdfpos"]]),
      fscore = as.numeric(fbetaObject[["fscores"]]),
      precision = as.numeric(fbetaObject[["precision"]]),
      recall = as.numeric(fbetaObject[["recall"]])
    )

    fbetaThreshold <- as.numeric(
      fbetaDirect$result$threshold
    )[[1L]]

    pFbetaDensity <- fbetaCurve |>
      select(.data$x, .data$pdfUns, .data$pdfStim) |>
      pivot_longer(
        cols = all_of(c("pdfUns", "pdfStim")),
        names_to = "distribution",
        values_to = "density"
      ) |>
      ggplot(aes(x = .data$x, y = .data$density, colour = .data$distribution)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_vline(
        xintercept = fbetaThreshold,
        linetype = "dashed"
      ) +
      theme_bw() +
      labs(
        title = "F-beta histograms used by fbeta.py",
        subtitle = paste("Threshold:", signif(fbetaThreshold, 5)),
        x = paste0(chnl, " expression"),
        y = "Smoothed histogram density",
        colour = NULL
      )

    pFbetaScore <- fbetaCurve |>
      select(
        .data$x,
        .data$fscore,
        .data$precision,
        .data$recall
      ) |>
      pivot_longer(
        cols = all_of(c("fscore", "precision", "recall")),
        names_to = "metric",
        values_to = "value"
      ) |>
      ggplot(aes(x = .data$x, y = .data$value, colour = .data$metric)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_vline(
        xintercept = fbetaThreshold,
        linetype = "dashed"
      ) +
      theme_bw() +
      labs(
        title = "F-beta objective",
        x = paste0(chnl, " expression"),
        y = "Metric",
        colour = NULL
      )

    path_p_fbeta_density <- projr::projr_path_get(
      "cache",
      "acs_cytof",
      "scratch",
      "debug_fbeta_and_tg",
      paste0("fbeta_density_", pop, "_", chnl, ".pdf")
    )
    ggsave(
      path_p_fbeta_density,
      plot = pFbetaDensity,
      units = "cm",
      width = 20,
      height = 20
    )
    path_p_fbeta_score <- projr::projr_path_get(
      "cache",
      "acs_cytof",
      "scratch",
      "debug_fbeta_and_tg",
      paste0("fbeta_score_", pop, "_", chnl, ".pdf")
    )
    ggsave(
      path_p_fbeta_score,
      plot = pFbetaScore,
      units = "cm",
      width = 20,
      height = 20
    )
  }

  # -----------------------------------------------------------------------
  # Tailgate internals
  # -----------------------------------------------------------------------

  tailgateInternal <- .tryAcsMethod({
    pars <- tailgateSettings$params

    bandwidth <- if (is.null(pars$bandwidth)) {
      ks::hpi(xStimFinite, deriv.order = 1L)
    } else {
      pars$bandwidth
    }

    derivative <- cytoUtils:::.deriv_density(
      x = xStimFinite,
      deriv = 1L,
      bandwidth = bandwidth,
      adjust = pars$adjust
    )

    peakObject <- openCyto:::.find_peaks(
      xStimFinite,
      num_peaks = pars$numPeaks,
      adjust = pars$adjust,
      plot = FALSE
    )

    peaks <- if (
      is.matrix(peakObject) ||
        is.data.frame(peakObject)
    ) {
      as.numeric(peakObject[, "x"])
    } else {
      as.numeric(peakObject)
    }

    referencePeak <- if (length(peaks) > 0L) {
      peaks[[min(pars$refPeak, length(peaks))]]
    } else {
      NA_real_
    }

    valleys <- with(
      derivative,
      openCyto:::.find_valleys(
        x = x,
        y = y,
        adjust = pars$adjust
      )
    )

    rightValleys <- sort(
      as.numeric(valleys)[as.numeric(valleys) > referencePeak]
    )

    selectedValley <- .firstFinite(rightValleys)

    tolerance <- if (isTRUE(pars$autoTol)) {
      0.01 * max(abs(derivative$y))
    } else {
      pars$tol
    }

    cutpointCandidates <- with(
      derivative,
      x[x > selectedValley & abs(y) < tolerance]
    )

    list(
      bandwidth = as.numeric(bandwidth)[[1L]],
      derivative = derivative,
      peaks = peaks,
      referencePeak = referencePeak,
      valleys = as.numeric(valleys),
      selectedValley = selectedValley,
      tolerance = tolerance,
      cutpointCandidates = cutpointCandidates,
      reconstructedThreshold = .firstFinite(cutpointCandidates)
    )
  })

  pTailgateDerivative <- NULL

  if (tailgateInternal$success) {
    tailgateDetails <- tailgateInternal$result

    tailgateSummary <- tibble(
      bandwidth = tailgateDetails$bandwidth,
      referencePeak = tailgateDetails$referencePeak,
      selectedDerivativeValley = tailgateDetails$selectedValley,
      automaticTolerance = tailgateDetails$tolerance,
      reconstructedThreshold = tailgateDetails$reconstructedThreshold,
      wrapperThreshold = tailgateWrapped$threshold
    )

    cat("\nTailgate internal calculation\n")
    print(tailgateSummary, width = Inf)

    tailgateLines <- tibble(
      feature = c(
        "reference peak",
        "derivative valley",
        "cut-point"
      ),
      x = c(
        tailgateDetails$referencePeak,
        tailgateDetails$selectedValley,
        tailgateDetails$reconstructedThreshold
      )
    ) |>
      filter(is.finite(.data$x))

    pTailgateDerivative <- as_tibble(
      tailgateDetails$derivative
    ) |>
      ggplot(aes(x = .data$x, y = .data$y)) +
      geom_hline(
        yintercept = c(
          -tailgateDetails$tolerance,
          tailgateDetails$tolerance
        ),
        colour = "grey50",
        linetype = "dotted"
      ) +
      geom_line(linewidth = 0.7) +
      geom_vline(
        data = tailgateLines,
        aes(xintercept = .data$x, colour = .data$feature),
        linewidth = 0.7,
        show.legend = TRUE
      ) +
      theme_bw() +
      labs(
        title = "Tailgate first derivative",
        subtitle = "Dotted lines show the automatic tolerance",
        x = paste0(chnl, " expression"),
        y = "First derivative",
        colour = NULL
      )

    print(pTailgateDerivative)
    path_p_tailgate_derivative <- projr::projr_path_get(
      "cache",
      "acs_cytof",
      "scratch",
      "debug_fbeta_and_tg",
      paste0("tailgate_derivative_", pop, "_", chnl, ".pdf")
    )
    ggsave(
      path_p_tailgate_derivative,
      plot = pTailgateDerivative,
      units = "cm",
      width = 20,
      height = 20
    )
  } else {
    message("\nTAILGATE INTERNAL ERROR: ", tailgateInternal$error)
  }

  if (isTRUE(showTailgatePackagePlot)) {
    invisible(cytoUtils:::.cytokine_cutpoint(
      x = xStimFinite,
      num_peaks = tailgateSettings$params$numPeaks,
      ref_peak = tailgateSettings$params$refPeak,
      method = "first_deriv",
      tol = tailgateSettings$params$tol,
      adjust = tailgateSettings$params$adjust,
      side = tailgateSettings$params$side,
      strict = tailgateSettings$params$strict,
      plot = TRUE,
      auto_tol = tailgateSettings$params$autoTol,
      bandwidth = tailgateSettings$params$bandwidth
    ))
  }

  # -----------------------------------------------------------------------
  # Observed event distributions with final gates
  # -----------------------------------------------------------------------

  eventTbl <- bind_rows(
    tibble(expr = xUnsFinite, condition = "unstimulated"),
    tibble(expr = xStimFinite, condition = "stimulated")
  )

  gateLines <- thresholdSummary |>
    transmute(
      method = .data$method,
      threshold = .data$threshold
    ) |>
    filter(is.finite(.data$threshold))

  pDistributionFull <- ggplot(
    eventTbl,
    aes(x = .data$expr, colour = .data$condition)
  ) +
    geom_density(linewidth = 0.8) +
    geom_vline(
      data = gateLines,
      aes(
        xintercept = .data$threshold,
        colour = .data$method
      ),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste(
        pop,
        basename(fcsFiles[[indStim]]),
        unname(channelMap[[chnl]])
      ),
      subtitle = "Dashed lines are the wrapped comparator thresholds",
      x = paste0(chnl, " expression"),
      y = "Density",
      colour = NULL
    )

  centralRange <- quantile(
    c(xUnsFinite, xStimFinite),
    c(0.001, 0.999),
    names = FALSE
  )

  pDistributionZoom <- pDistributionFull +
    coord_cartesian(xlim = centralRange) +
    labs(
      subtitle = paste(
        "Central 99.8% of events.",
        "An extreme gate may be outside this view."
      )
    )

  path_p_distribution_full <- projr::projr_path_get(
    "cache",
    "acs_cytof",
    "scratch",
    "debug_fbeta_and_tg",
    paste0("distribution_full_", pop, "_", chnl, ".pdf")
  )
  ggsave(
    path_p_distribution_full,
    plot = pDistributionFull,
    units = "cm",
    width = 20,
    height = 20
  )
  path_p_distribution_zoom <- projr::projr_path_get(
    "cache",
    "acs_cytof",
    "scratch",
    "debug_fbeta_and_tg",
    paste0("distribution_zoom_", pop, "_", chnl, ".pdf")
  )
  ggsave(
    path_p_distribution_zoom,
    plot = pDistributionZoom,
    units = "cm",
    width = 20,
    height = 20
  )

  invisible(list(
    sampleInfo = sampleInfo,
    eventSummary = eventSummary,
    thresholds = thresholdSummary,
    fbetaRangeCheck = fbetaRangeCheck,
    fbetaDirect = fbetaDirect,
    tailgateDirect = tailgateDirect,
    tailgateInternal = tailgateInternal,
    batchLayout = batchLayout,
    fbetaCurve = fbetaCurve,
    plots = list(
      distributionFull = pDistributionFull,
      distributionZoom = pDistributionZoom,
      fbetaDensity = pFbetaDensity,
      fbetaScore = pFbetaScore,
      tailgateDerivative = pTailgateDerivative
    ),
    xUns = xUns,
    xStim = xStim
  ))
}


diagnostic_tcrgd <- .inspectAcsComparatorCase(
  pop = "tcrgd",
  batchIndex = 1L,
  stimOffset = 2L,
  chnl = "Ho165Di", # IFNg
  outputGroup = NULL, # use "tester" for the tester GatingSet
  showTailgatePackagePlot = TRUE
)

diagnostic_cd4 <- .inspectAcsComparatorCase(
  pop = "cd4",
  batchIndex = 2L,
  stimOffset = 4L,
  chnl = "Nd146Di" # TNF
)
