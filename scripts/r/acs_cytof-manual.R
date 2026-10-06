.acsCytofManualPopulationMap <- function() {
  c(
    tcrgd = "TCRgd T cells",
    cd4 = "CD4 T cells",
    cd8 = "CD8 T cells",
    nk = "NK cells",
    nk_pre = NA_character_,
    b = "B cells"
  )
}

.acsCytofManualRead <- function(fn) {
  pathManual <- projr::projr_path_get(
    "raw-data-small",
    "comparison_data",
    "acscytof",
    fn,
    format = "absolute"
  )
  if (!file.exists(pathManual)) {
    stop("Manual ACS cytokine file not found at: ", pathManual)
  }

  manualRaw <- suppressMessages(suppressWarnings(readr::read_csv(
    pathManual,
    name_repair = "unique",
    show_col_types = FALSE
  )))
  idCol <- names(manualRaw)[[1]]
  sourceParts <- stringr::str_split(as.character(manualRaw[[idCol]]), "_")
  if (any(lengths(sourceParts) < 5L)) {
    stop("Could not parse SampleID and stimulus from the manual ACS file.")
  }

  emptyCol <- vapply(
    manualRaw,
    function(x) {
      xCharacter <- trimws(as.character(x))
      all(is.na(x) | is.na(xCharacter) | !nzchar(xCharacter))
    },
    logical(1)
  )
  emptyCol[[idCol]] <- FALSE

  list(data = manualRaw[, !emptyCol, drop = FALSE])
}

.comp_against_manual_cyt_format_manual <- function(
  fn,
  pop = NULL,
  cyt = NULL
) {
  manualObject <- .acsCytofManualRead(fn)

  manualTblInit <- manualObject$data |>
    dplyr::rename(fn = `...1`)
  x <- manualTblInit[[1]]
  sampleidVec <- sub("^([^_]+_[^_]+).*", "\\1", x)
  stimVec <- sub(".*_", "", x)
  manualTbl <- manualTblInit |>
    dplyr::mutate(
      SampleID = sampleidVec,
      stim = dplyr::if_else(stimVec == "mtbaux", "mtb", stimVec)
    ) |>
    dplyr::select(-fn) |>
    dplyr::select(SampleID, stim, dplyr::everything())

  names(manualTbl) <- stringr::str_remove(
    names(manualTbl),
    " FreqofParent$"
  )

  manualTbl <- manualTbl |>
    tidyr::pivot_longer(
      cols = -c(SampleID, stim),
      names_to = "popCyt",
      values_to = "freq_stim_man"
    ) |>
    tidyr::separate(
      col = "popCyt",
      into = c("pop", "cyt"),
      sep = "/",
      extra = "merge",
      fill = "right"
    ) |>
    dplyr::filter(!is.na(.data$cyt), .data$cyt != "Perf") |>
    dplyr::mutate(
      pop = dplyr::if_else(
        .data$pop == "TCRgd+",
        "TCRgd T cells",
        .data$pop
      ),
      freq_stim_man = as.numeric(.data$freq_stim_man)
    )

  manualUns <- manualTbl |>
    dplyr::filter(.data$stim == "uns") |>
    dplyr::count(SampleID, pop, cyt, name = "nUns")
  if (any(manualUns$nUns != 1L)) {
    stop("The manual file has repeated unstimulated population/cytokine rows.")
  }

  manualUns <- manualTbl |>
    dplyr::filter(.data$stim == "uns") |>
    dplyr::select(
      SampleID,
      pop,
      cyt,
      freq_uns_man = freq_stim_man
    )
  missingUns <- manualTbl |>
    dplyr::distinct(SampleID, pop, cyt) |>
    dplyr::anti_join(manualUns, by = c("SampleID", "pop", "cyt"))
  if (nrow(missingUns) > 0L) {
    stop("The manual file is missing an unstimulated population/cytokine row.")
  }

  manualTbl <- manualTbl |>
    dplyr::left_join(manualUns, by = c("SampleID", "pop", "cyt")) |>
    dplyr::mutate(
      freq_bs_man = pmax(.data$freq_stim_man - .data$freq_uns_man, 0)
    ) |>
    dplyr::select(
      SampleID,
      stim,
      pop,
      cyt,
      freq_stim_man,
      freq_uns_man,
      freq_bs_man
    )
  if (!is.null(pop)) {
    manualTbl <- dplyr::filter(manualTbl, .data$pop %in% .env$pop)
  }
  if (!is.null(cyt)) {
    manualTbl <- dplyr::filter(manualTbl, .data$cyt %in% .env$cyt)
  }
  manualTbl
}

.acsCytofManualResolvePopCodes <- function(pop, popMap) {
  if (is.null(pop)) {
    return(names(popMap))
  }

  unique(purrr::map_chr(pop, function(x) {
    if (x %in% names(popMap)) {
      return(x)
    }
    matchIndex <- which(!is.na(popMap) & popMap == x)
    if (length(matchIndex) == 1L) {
      return(names(popMap)[[matchIndex]])
    }
    stop("Unknown or unsupported ACS population: ", x)
  }))
}

.acsCytofManualSampleMapFromFcs <- function(
  pathFcsBase,
  popCode,
  sampleLookup
) {
  fcsFiles <- .acsCytofFcsFiles(file.path(pathFcsBase, popCode))
  .acsCytofMapFiles(fcsFiles, sampleLookup) |>
    dplyr::mutate(popCode = .env$popCode)

}

.acsCytofManualValidateSampleMap <- function(sampleMap, popCodes) {
  requiredColumns <- c("popCode", "ind", "SampleID", "stim")
  missingColumns <- setdiff(requiredColumns, names(sampleMap))
  if (length(missingColumns) > 0L) {
    stop(
      "sampleMap is missing required column(s): ",
      paste(missingColumns, collapse = ", ")
    )
  }

  sampleMap <- sampleMap |>
    dplyr::transmute(
      popCode = as.character(.data$popCode),
      ind = as.character(.data$ind),
      SampleID = as.character(.data$SampleID),
      stim = dplyr::if_else(
        .data$stim == "mtbaux",
        "mtb",
        as.character(.data$stim)
      )
    ) |>
    dplyr::filter(.data$popCode %in% .env$popCodes)

  duplicateMap <- sampleMap |>
    dplyr::count(popCode, ind, name = "n") |>
    dplyr::filter(.data$n != 1L)
  if (nrow(duplicateMap) > 0L) {
    stop("sampleMap has duplicate popCode/ind keys.")
  }
  sampleMap
}

.acsCytofManualMethodPath <- function(
  pathScratchBase,
  popCode,
  method,
  outputGroup = NULL
) {
  pathParts <- c(
    pathScratchBase,
    if (!is.null(outputGroup)) outputGroup,
    popCode
  )
  pathPopulation <- do.call(file.path, as.list(pathParts))

  switch(
    method,
    stimgate = file.path(pathPopulation, "stimgate"),
    tailgate = file.path(pathPopulation, "tailgate", "result.rds"),
    fbeta = file.path(pathPopulation, "fbeta", "result.rds"),
    stop("Unknown ACS method: ", method)
  )
}

.acsCytofManualReadStats <- function(path, method, gateName) {
  if (identical(method, "stimgate")) {
    if (!dir.exists(path)) {
      stop("StimGate output not found at: ", path)
    }
    statsTbl <- stimgate::getStimStats(path) |>
      dplyr::filter(.data$gateName == .env$gateName) |>
      dplyr::mutate(method = "stimgate")
    if (nrow(statsTbl) == 0L) {
      stop("No '", gateName, "' rows were found at: ", path)
    }
    return(statsTbl)
  }

  .acsCytofReadComparatorCache(path = path, method = method)$stats
}

.acsCytofCombinationToStandard <- function(cytCombn, channelMap) {
  combinationStd <- cytCombn |>
    stringr::str_replace_all("~\\+~", "+") |>
    stringr::str_replace_all("~-~", "-")
  combinationCompass <- UtilsCompassSV::convert_cyt_combn_format(
    cyt_combn = combinationStd,
    to = "compass",
    silent = TRUE,
    lab = channelMap
  )
  UtilsCompassSV::convert_cyt_combn_format(
    cyt_combn = combinationCompass,
    to = "std",
    silent = TRUE
  )
}

.acsCytofStatsSingleMarkers <- function(
  statsTbl,
  method,
  popCode,
  popLabel,
  sampleMap,
  cyt = NULL
) {
  requiredColumns <- c(
    "ind",
    "cytCombn",
    "countStim",
    "nCellStim",
    "countUns",
    "nCellUns"
  )
  missingColumns <- setdiff(requiredColumns, names(statsTbl))
  if (length(missingColumns) > 0L) {
    stop(
      method,
      " stats for '",
      popCode,
      "' are missing column(s): ",
      paste(missingColumns, collapse = ", ")
    )
  }

  channelMap <- .acsCytofChannelMap()
  channels <- names(channelMap)
  channelsToKeep <- channels
  if (!is.null(cyt)) {
    unknownCytokines <- setdiff(cyt, unname(channelMap))
    if (length(unknownCytokines) > 0L) {
      stop(
        "Unknown ACS cytokine label(s): ",
        paste(unknownCytokines, collapse = ", ")
      )
    }
    channelsToKeep <- channels[channelMap %in% cyt]
  }

  statsTbl <- statsTbl |>
    dplyr::mutate(
      ind = as.character(.data$ind),
      method = .env$method,
      popCode = .env$popCode
    )
  groupColumns <- intersect(
    c(
      "gateName",
      "method",
      "popCode",
      "pop",
      "batch",
      "ind",
      "indUns",
      "nCellStim",
      "nCellUns"
    ),
    names(statsTbl)
  )

  singleTbl <- purrr::map_dfr(channelsToKeep, function(channel) {
    UtilsCytoRSV::sum_over_markers(
      .data = statsTbl,
      grp = groupColumns,
      cmbn = "cytCombn",
      levels = c("~+~", "~-~"),
      markers_to_sum = setdiff(channels, channel),
      resp = c("countStim", "countUns")
    ) |>
      dplyr::filter(.data$cytCombn == paste0(channel, "~+~"))
  }) |>
    dplyr::mutate(
      cytCombn = .acsCytofCombinationToStandard(
        .data$cytCombn,
        channelMap = channelMap
      ),
      cyt = stringr::str_remove(.data$cytCombn, "[+-]$")
    )

  mapPopulation <- sampleMap |>
    dplyr::filter(.data$popCode == .env$popCode) |>
    dplyr::select(ind, SampleID, stim)
  unmatchedIndex <- setdiff(unique(singleTbl$ind), mapPopulation$ind)
  if (length(unmatchedIndex) > 0L) {
    stop("Unmapped ", method, " sample indices for '", popCode, "': ",
         paste(unmatchedIndex, collapse = ", "))
  }

  singleTbl |>
    dplyr::inner_join(mapPopulation, by = "ind") |>
    dplyr::mutate(
      pop = .env$popLabel,
      freq_stim_auto = .data$countStim / .data$nCellStim * 100,
      freq_uns_auto = .data$countUns / .data$nCellUns * 100,
      freq_bs_auto = pmax(.data$freq_stim_auto - .data$freq_uns_auto, 0)
    ) |>
    dplyr::select(
      method,
      SampleID,
      stim,
      popCode,
      pop,
      cyt,
      cytCombn,
      dplyr::any_of(c("gateName", "batch", "ind", "indUns")),
      countStim,
      nCellStim,
      countUns,
      nCellUns,
      freq_stim_auto,
      freq_uns_auto,
      freq_bs_auto
    )
}

.acsCytofManualAutoTable <- function(
  pathScratchBase,
  pathFcsBase,
  pop = NULL,
  cyt = NULL,
  methods = c("stimgate", "tailgate", "fbeta"),
  outputGroup = NULL,
  gateName = "loc_minClust",
  sampleMap = NULL
) {
  methods <- match.arg(
    methods,
    choices = c("stimgate", "tailgate", "fbeta"),
    several.ok = TRUE
  )
  popMap <- .acsCytofManualPopulationMap()
  popCodes <- .acsCytofManualResolvePopCodes(pop, popMap)
  noManual <- popCodes[is.na(popMap[popCodes])]
  if (length(noManual) > 0L) {
    warning(
      "No matching manual population exists for: ",
      paste(noManual, collapse = ", "),
      ". These population(s) will not be included."
    )
    popCodes <- setdiff(popCodes, noManual)
  }

  if (is.null(pop)) {
    hasStimGate <- vapply(
      popCodes,
      function(popCode) {
        dir.exists(.acsCytofManualMethodPath(
          pathScratchBase,
          popCode,
          method = "stimgate",
          outputGroup = outputGroup
        ))
      },
      logical(1)
    )
    popCodes <- popCodes[hasStimGate]
  }
  if (length(popCodes) == 0L) {
    stop("No requested ACS population has an automated and manual result.")
  }

  manifests <- list()
  result <- purrr::map_dfr(popCodes, function(popCode) {
    paths <- lapply(methods, function(method) {
      .acsCytofManualMethodPath(pathScratchBase, popCode, method, outputGroup)
    })
    names(paths) <- methods
    saved <- lapply(methods, function(method) {
      if (method == "stimgate") {
        pathManifest <- file.path(paths[[method]], "acs-manifest.rds")
        if (!file.exists(pathManifest)) stop("ACS StimGate manifest missing; re-run all methods.")
        manifest <- readRDS(pathManifest)
        if (!identical(manifest$channelSettings, stimgate::stimgateMetaReadSettingsChnls(paths[[method]]))) {
          stop("Mismatched ACS StimGate settings manifest; re-run all methods.")
        }
        manifest
      } else {
        object <- .acsCytofReadComparatorCache(paths[[method]], method)
        if (!identical(object$manifest$settings, object$settings)) {
          stop("Mismatched ACS comparator manifest settings; re-run all methods.")
        }
        object$manifest
      }
    })
    names(saved) <- methods
    .acsCytofValidateManifests(saved)
    manifests[[popCode]] <<- saved
    mapped <- saved[[1]]$context$preprocessing$sampleMap |>
      dplyr::mutate(popCode = .env$popCode)
    mapped <- .acsCytofManualValidateSampleMap(mapped, popCode)
    purrr::map_dfr(methods, function(method) {
      statsTbl <- .acsCytofManualReadStats(paths[[method]], method, gateName)
      single <- .acsCytofStatsSingleMarkers(
        statsTbl, method, popCode, unname(popMap[[popCode]]), mapped, cyt
      )
      thresholds <- if (method == "stimgate") {
        stimgate::getStimGates(paths[[method]]) |>
          dplyr::filter(.data$gateName == .env$gateName) |>
          dplyr::mutate(
            threshold = .data$gate,
            thresholdOrigin = .data$locSource,
            thresholdFallbackUsed = !(.data$locGenerated %in% TRUE),
            cyt = unname(.acsCytofChannelMap()[.data$chnl])
          )
      } else {
        .acsCytofReadComparatorCache(paths[[method]], method)$thresholds
      }
      .acsCytofJoinProvenance(single, thresholds)
    })
  })
  .acsCytofValidateCohorts(result, methods)
  attr(result, "manifest") <- manifests
  result

}

.acsCytofManualComparisonTable <- function(
  fn,
  pathScratchBase,
  pathFcsBase,
  pop = NULL,
  cyt = NULL,
  methods = c("stimgate", "tailgate", "fbeta"),
  outputGroup = NULL,
  gateName = "loc_minClust",
  sampleMap = NULL
) {
  autoTbl <- .acsCytofManualAutoTable(
    pathScratchBase = pathScratchBase,
    pathFcsBase = pathFcsBase,
    pop = pop,
    cyt = cyt,
    methods = methods,
    outputGroup = outputGroup,
    gateName = gateName,
    sampleMap = sampleMap
  )
  manualTbl <- .comp_against_manual_cyt_format_manual(
    fn = fn,
    pop = unique(autoTbl$pop),
    cyt = cyt
  ) |>
    dplyr::filter(.data$stim != "uns")

  joinBy <- c("SampleID", "stim", "pop", "cyt")
  unmatchedAuto <- autoTbl |>
    dplyr::anti_join(manualTbl, by = joinBy)
  unmatchedManual <- manualTbl |>
    dplyr::anti_join(autoTbl, by = joinBy)
  exclusions <- dplyr::bind_rows(
    dplyr::mutate(unmatchedAuto, exclusionReason = "no_manual_key"),
    dplyr::mutate(unmatchedManual, exclusionReason = "no_automated_key")
  )
  if (anyDuplicated(manualTbl[joinBy])) stop("Duplicate ACS manual comparison keys.")

  methodLevels <- c("stimgate", "tailgate", "fbeta")
  popLevels <- stats::na.omit(unname(.acsCytofManualPopulationMap()))
  cytLevels <- unname(.acsCytofChannelMap())

  result <- autoTbl |>
    dplyr::inner_join(manualTbl, by = joinBy) |>
    dplyr::mutate(
      method = factor(.data$method, levels = methodLevels),
      pop = factor(.data$pop, levels = popLevels),
      cyt = factor(.data$cyt, levels = cytLevels),
      diff = .data$freq_bs_auto - .data$freq_bs_man,
      abs_diff = abs(.data$diff),
      rel_error = dplyr::if_else(
        .data$freq_bs_man > 0,
        .data$diff / .data$freq_bs_man,
        NA_real_
      ),
      abs_rel_error = abs(.data$rel_error)
    ) |>
    dplyr::arrange(
      .data$method,
      .data$pop,
      .data$cyt,
      dplyr::desc(.data$abs_diff)
    )
  attr(result, "manifest") <- list(
    methods = attr(autoTbl, "manifest"),
    comparisonSettings = list(gateName = gateName, pop = pop, cyt = cyt, methods = methods),
    manualInputHash = .acsCytofHash(manualTbl)
  )
  attr(result, "exclusions") <- dplyr::bind_rows(
    exclusions,
    tibble::tibble(popCode = setdiff(.acsCytofManualResolvePopCodes(pop, .acsCytofManualPopulationMap()), unique(autoTbl$popCode)),
                   exclusionReason = "population_without_manual_result")
  )
  result

}

.acsCytofManualSummaryTable <- function(comparisonTbl) {
  if (!"thresholdFailed" %in% names(comparisonTbl)) comparisonTbl$thresholdFailed <- FALSE
  comparisonTbl |>
    dplyr::mutate(dplyr::across(c(freq_bs_auto, abs_diff, abs_rel_error),
                              ~ dplyr::if_else(.data$thresholdFailed, NA_real_, .x))) |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_failed = sum(.data$thresholdFailed),
      n_manual_nonpositive = sum(is.finite(.data$freq_bs_man) & .data$freq_bs_man <= 0),
      n_manual_missing = sum(!is.finite(.data$freq_bs_man)),
      n_relative = sum(is.finite(.data$abs_rel_error)),
      n_relative_excluded = dplyr::n() - n_relative,
      prop_failed = mean(.data$thresholdFailed),
      n = sum(is.finite(.data$freq_bs_auto) & is.finite(.data$freq_bs_man)),
      pcc = if (n > 1L) {
        suppressWarnings(stats::cor(
          .data$freq_bs_auto,
          .data$freq_bs_man,
          use = "complete.obs"
        ))
      } else {
        NA_real_
      },
      mae = mean(.data$abs_diff, na.rm = TRUE),
      median_abs_rel_error = stats::median(
        .data$abs_rel_error,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

# Resample donors once for the whole comparison. Every stimulated tube and its
# already-subtracted shared control travel together, across stimuli and methods.
.acsCytofManualUncertainty <- function(comparisonTbl, reps = 999L, seed = 20261005L) {
  if (!"SampleID" %in% names(comparisonTbl) || anyNA(comparisonTbl$SampleID)) {
    stop("ACS donor bootstrap requires complete SampleID values.")
  }
  if (length(reps) != 1L || !is.finite(reps) || reps < 2L || reps != as.integer(reps)) {
    stop("ACS donor bootstrap needs at least two integer replicates.")
  }
  donors <- unique(as.character(comparisonTbl$SampleID))
  weights <- .analysis_with_seed(seed, replicate(reps, tabulate(
    sample.int(length(donors), length(donors), replace = TRUE), length(donors)
  )))
  weights <- matrix(weights, nrow = length(donors))
  comparisonTbl |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::group_modify(function(rows, key) {
      failed <- if ("thresholdFailed" %in% names(rows)) rows$thresholdFailed %in% TRUE else FALSE
      valid <- !failed & is.finite(rows$freq_bs_auto) & is.finite(rows$freq_bs_man)
      absError <- abs(rows$freq_bs_auto - rows$freq_bs_man)
      relative <- valid & rows$freq_bs_man > 0
      interval <- function(values, keep) {
        if (dplyr::n_distinct(rows$SampleID[keep]) < 2L) return(c(NA_real_, NA_real_, 0L))
        sums <- counts <- numeric(length(donors))
        for (i in seq_along(donors)) {
          selected <- keep & as.character(rows$SampleID) == donors[i]
          sums[i] <- sum(values[selected])
          counts[i] <- sum(selected)
        }
        denominator <- as.numeric(crossprod(counts, weights))
        draws <- as.numeric(crossprod(sums, weights)) / denominator
        c(stats::quantile(draws[is.finite(draws)], c(0.025, 0.975), names = FALSE), sum(is.finite(draws)))
      }
      absCi <- interval(absError, valid)
      relCi <- interval(absError / rows$freq_bs_man, relative)
      tibble::tibble(
        n_donors = dplyr::n_distinct(rows$SampleID[valid]),
        n_relative_donors = dplyr::n_distinct(rows$SampleID[relative]),
        mae = if (any(valid)) mean(absError[valid]) else NA_real_,
        mae_lower = absCi[1], mae_upper = absCi[2],
        mean_abs_rel_error = if (any(relative)) mean(absError[relative] / rows$freq_bs_man[relative]) else NA_real_,
        mean_abs_rel_error_lower = relCi[1], mean_abs_rel_error_upper = relCi[2],
        n_absolute_bootstrap_finite = as.integer(absCi[3]),
        n_relative_bootstrap_finite = as.integer(relCi[3]),
        bootstrap_reps = reps
      )
    }) |>
    dplyr::ungroup()
}

.acsCytofManualPlotAbsoluteError <- function(comparisonTbl) {
  comparisonTbl$method <- .acsCytofManualMethodFactor(comparisonTbl$method)
  ggplot2::ggplot(comparisonTbl, ggplot2::aes(
    x = .data$method, y = .data$abs_diff, fill = .data$method
  )) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::facet_grid(rows = ggplot2::vars(pop), cols = ggplot2::vars(cyt), scales = "free_y") +
    ggplot2::scale_x_discrete(labels = .analysis_method_labels) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number) +
    .analysis_scale_method("fill") +
    ggplot2::labs(x = NULL, y = "Absolute frequency difference (percentage points)") +
    .analysis_theme(grid = "y") +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# Method as a factor in the standard order (StimGate, Tailgate, F-beta).
.acsCytofManualMethodFactor <- function(method) {
  method <- as.character(method)
  levels <- c(
    intersect(names(.analysis_method_labels), method),
    setdiff(unique(method), names(.analysis_method_labels))
  )
  factor(method, levels = levels)
}

.acsCytofManualPlotScatter <- function(comparisonTbl) {
  comparisonTbl$method <- .acsCytofManualMethodFactor(comparisonTbl$method)
  ggplot2::ggplot(
    comparisonTbl,
    ggplot2::aes(
      x = .data$freq_bs_man,
      y = .data$freq_bs_auto,
      colour = .data$method,
      shape = .data$stim
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "grey45",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(pop),
      cols = ggplot2::vars(cyt),
      scales = "free"
    ) +
    ggplot2::scale_x_continuous(
      labels = .analysis_label_number,
      n.breaks = 3,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number, n.breaks = 3) +
    .analysis_scale_method() +
    ggplot2::labs(
      x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency (automated gating, %)",
      colour = "Method",
      shape = "Stimulus"
    ) +
    .analysis_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7))
}

.acsCytofManualPlotRelativeError <- function(comparisonTbl) {
  comparisonTbl$method <- .acsCytofManualMethodFactor(comparisonTbl$method)
  ggplot2::ggplot(
    comparisonTbl,
    ggplot2::aes(
      x = .data$method,
      y = .data$abs_rel_error,
      fill = .data$method
    )
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(pop),
      cols = ggplot2::vars(cyt),
      scales = "free_y"
    ) +
    ggplot2::scale_x_discrete(labels = .analysis_method_labels) +
    ggplot2::scale_y_continuous(labels = .analysis_label_number) +
    .analysis_scale_method("fill") +
    ggplot2::labs(x = NULL, y = "Absolute relative error") +
    .analysis_theme(grid = "y") +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# Over- and under-estimates of the relative error (estimate - manual) / manual,
# summarised over samples for each method, population and cytokine. Uses the
# same rows as `.acsCytofManualPlotRelativeError()`: rows with a non-positive
# manual frequency have no relative error. Over-estimates plot above zero and
# under-estimates below; point size is the share of samples in that direction.
# Needs `.simBandwidthSignedErrorSummary()` and `.simBandwidthSignedErrorLayers()`
# from `sim-bandwidth-analysis-plot.R`.
.acsCytofManualPlotSignedError <- function(comparisonTbl) {
  comparisonTbl$method <- .acsCytofManualMethodFactor(comparisonTbl$method)
  summaryTbl <- .simBandwidthSignedErrorSummary(
    comparisonTbl,
    c("method", "pop", "cyt")
  )
  statCols <- c(median = "Median", q95 = "95th percentile", max = "Maximum")
  plotTbl <- summaryTbl |>
    tidyr::pivot_longer(
      dplyr::all_of(names(statCols)),
      names_to = "statistic",
      values_to = "value"
    ) |>
    dplyr::mutate(
      statistic = factor(
        .data$statistic,
        levels = names(statCols),
        labels = statCols
      )
    )
  ggplot2::ggplot(
    plotTbl,
    ggplot2::aes(
      x = .data$cyt,
      y = .data$value,
      colour = .data$method,
      size = .data$prop,
      group = interaction(.data$method, .data$direction)
    )
  ) +
    .simBandwidthSignedErrorLayers("Relative error") +
    ggplot2::geom_point(
      alpha = 0.75,
      position = ggplot2::position_dodge(width = 0.6)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$statistic),
      cols = ggplot2::vars(.data$pop),
      scales = "free_y"
    ) +
    ggplot2::scale_size_continuous(
      range = c(0.5, 3),
      limits = c(0, 1),
      labels = .analysis_label_percent
    ) +
    .analysis_scale_method() +
    ggplot2::labs(
      x = "Cytokine",
      colour = "Method",
      size = "Share of samples\nin this direction"
    ) +
    .analysis_theme() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

.acsCytofManualWrite <- function(
  comparisonTbl,
  pathDirSave,
  savePlots = TRUE
) {
  dir.create(pathDirSave, recursive = TRUE, showWarnings = FALSE)
  summaryTbl <- .acsCytofManualSummaryTable(comparisonTbl)
  uncertaintyTbl <- .acsCytofManualUncertainty(comparisonTbl)
  utils::write.csv(uncertaintyTbl, file.path(pathDirSave, "manual-comparison-donor-uncertainty.csv"), row.names = FALSE)
  exclusions <- attr(comparisonTbl, "exclusions")
  if (is.null(exclusions)) exclusions <- data.frame(exclusionReason = character())
  utils::write.csv(exclusions, file.path(pathDirSave, "manual-comparison-exclusions.csv"), row.names = FALSE)
  saveRDS(attr(comparisonTbl, "manifest"), file.path(pathDirSave, "acs-manifest.rds"))

  utils::write.csv(
    comparisonTbl,
    file.path(pathDirSave, "manual-comparison.csv"),
    row.names = FALSE
  )
  .write_rds_atomic(
    comparisonTbl,
    file.path(pathDirSave, "manual-comparison.rds")
  )
  utils::write.csv(
    summaryTbl,
    file.path(pathDirSave, "manual-comparison-summary.csv"),
    row.names = FALSE
  )

  if (!isTRUE(savePlots)) {
    return(invisible(list(summary = summaryTbl)))
  }

  scatter <- .acsCytofManualPlotScatter(comparisonTbl)
  relativeError <- .acsCytofManualPlotRelativeError(comparisonTbl)
  nPop <- max(1L, dplyr::n_distinct(comparisonTbl$pop))
  ggplot2::ggsave(
    file.path(pathDirSave, "manual-comparison-scatter.png"),
    plot = scatter,
    width = 30,
    height = max(16, 5 * nPop),
    units = "cm"
  )
  ggplot2::ggsave(
    file.path(pathDirSave, "manual-comparison-relative-error.png"),
    plot = relativeError,
    width = 30,
    height = max(16, 5 * nPop),
    units = "cm"
  )

  invisible(list(
    scatter = scatter,
    relativeError = relativeError,
    summary = summaryTbl
  ))
}

.acsCytofManualSave <- function(comparisonTbl, pathDirSave, savePlots = TRUE) {
  result <- NULL
  .acsCytofReplaceDir(pathDirSave, function(pathTmp) {
    result <<- .acsCytofManualWrite(comparisonTbl, pathTmp, savePlots)
  })
  invisible(result)
}

comp_against_manual_cyt <- function(
  fn,
  path_scratch_base,
  path_fcs_base,
  pop = NULL,
  cyt = NULL,
  methods = c("stimgate", "tailgate", "fbeta"),
  output_group = NULL,
  gate_name = "loc_minClust",
  sample_map = NULL,
  path_dir_save = NULL,
  save_plots = TRUE
) {
  comparisonTbl <- .acsCytofManualComparisonTable(
    fn = fn,
    pathScratchBase = path_scratch_base,
    pathFcsBase = path_fcs_base,
    pop = pop,
    cyt = cyt,
    methods = methods,
    outputGroup = output_group,
    gateName = gate_name,
    sampleMap = sample_map
  )
  attr(comparisonTbl, "summary") <- .acsCytofManualSummaryTable(comparisonTbl)

  if (!is.null(path_dir_save)) {
    .acsCytofManualSave(
      comparisonTbl = comparisonTbl,
      pathDirSave = path_dir_save,
      savePlots = save_plots
    )
  }
  invisible(comparisonTbl)
}

.acsCytofJoinProvenance <- function(single, thresholds) {
  required <- c("ind", "cyt", "threshold", "thresholdOrigin", "thresholdFallbackUsed")
  if (!all(required %in% names(thresholds))) stop("Missing ACS threshold provenance; re-run methods.")
  provenance <- thresholds |>
    dplyr::mutate(ind = as.character(.data$ind)) |>
    dplyr::select(dplyr::all_of(required), dplyr::any_of(c(
      "thresholdRaw", "gateCyt", "locGenerated", "locGeneratedDirect", "locSource", "locReason"
    )))
  keys <- c("ind", "cyt")
  if (anyDuplicated(provenance[keys]) ||
      nrow(dplyr::anti_join(single, provenance, by = keys))) {
    stop("Missing or duplicate ACS per-marker threshold provenance.")
  }
  single |>
    dplyr::left_join(provenance, by = keys) |>
    dplyr::mutate(
      thresholdFailed = !(.data$thresholdFallbackUsed %in% FALSE) | !is.finite(.data$threshold),
      dplyr::across(c(freq_stim_auto, freq_uns_auto, freq_bs_auto),
                    ~ dplyr::if_else(.data$thresholdFailed, NA_real_, .x))
    )
}

.acsCytofValidateCohorts <- function(table, methods) {
  keys <- c("SampleID", "stim", "pop", "cyt")
  reference <- table[as.character(table$method) == methods[[1]], keys]
  for (method in methods) {
    cohort <- table[as.character(table$method) == method, keys]
    if (anyNA(cohort) || anyDuplicated(cohort) ||
        nrow(dplyr::anti_join(reference, cohort, by = keys)) ||
        nrow(dplyr::anti_join(cohort, reference, by = keys))) {
      stop("ACS methods have different or duplicate sample/stratum keys: ", method)
    }
  }
  invisible(TRUE)
}

# One donor stream across every method/population/cytokine, including donors
# absent from a stratum as empty blocks, so resampling retains their pairing.
.acsCytofManualSignedPercentiles <- function(comparisonTbl) {
  if (!"SampleID" %in% names(comparisonTbl) || anyNA(comparisonTbl$SampleID)) {
    stop("ACS donor bootstrap requires complete SampleID values.")
  }
  donors <- sort(unique(as.character(comparisonTbl$SampleID)))
  comparisonTbl |>
    dplyr::group_by(.data$method, .data$pop, .data$cyt) |>
    dplyr::group_modify(function(rows, key) {
      missing <- setdiff(donors, as.character(rows$SampleID))
      errors <- rows$rel_error
      if ("thresholdFailed" %in% names(rows)) errors[rows$thresholdFailed %in% TRUE] <- NA_real_
      .simBandwidthSignedErrorPercentiles(
        c(errors, rep(NA_real_, length(missing))), mcse = TRUE,
        unit = c(as.character(rows$SampleID), missing), bootstrap_family = "acs-donors")
    }) |>
    dplyr::ungroup()
}
