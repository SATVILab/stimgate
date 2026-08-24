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

.acsCytofManualNormaliseId <- function(x) {
  x |>
    basename() |>
    tools::file_path_sans_ext() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^[:alnum:]]", "")
}

.acsCytofManualRead <- function(fn) {
  pathManual <- projr::projr_path_get(
    "raw-data-small",
    "comparison_data",
    "acscytof",
    fn
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
  sourceId <- as.character(manualRaw[[idCol]])
  sourceParts <- stringr::str_split(sourceId, "_")
  if (any(lengths(sourceParts) < 5L)) {
    stop("Could not parse SampleID and stimulus from the manual ACS file.")
  }

  sampleLookup <- tibble::tibble(
    sourceId = sourceId,
    SampleID = purrr::map_chr(
      sourceParts,
      function(x) paste0(x[[1]], "_", x[[2]])
    ),
    stim = purrr::map_chr(sourceParts, function(x) x[[5]]),
    pid = stringr::str_extract(
      stringr::str_to_lower(sourceId),
      "pid[12]"
    )
  ) |>
    dplyr::mutate(
      stim = dplyr::if_else(.data$stim == "mtbaux", "mtb", .data$stim),
      sourceKey = .acsCytofManualNormaliseId(.data$sourceId)
    ) |>
    dplyr::distinct()

  emptyCol <- vapply(
    manualRaw,
    function(x) {
      xCharacter <- trimws(as.character(x))
      all(is.na(x) | is.na(xCharacter) | !nzchar(xCharacter))
    },
    logical(1)
  )
  emptyCol[[idCol]] <- FALSE

  list(
    data = manualRaw[, !emptyCol, drop = FALSE],
    idCol = idCol,
    sampleLookup = sampleLookup
  )
}

.comp_against_manual_cyt_format_manual <- function(
    fn,
    pop = NULL,
    cyt = NULL) {
  manualObject <- .acsCytofManualRead(fn)
  idCol <- manualObject$idCol

  manualTbl <- manualObject$data |>
    dplyr::mutate(.sourceId = as.character(.data[[idCol]])) |>
    dplyr::left_join(
      manualObject$sampleLookup |>
        dplyr::select(.sourceId = sourceId, SampleID, stim),
      by = ".sourceId"
    ) |>
    dplyr::select(-dplyr::all_of(idCol))
  names(manualTbl) <- stringr::str_remove(
    names(manualTbl),
    " FreqofParent$"
  )

  manualTbl <- manualTbl |>
    tidyr::pivot_longer(
      cols = -c(.sourceId, SampleID, stim),
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

.acsCytofManualStimPattern <- function(stim) {
  switch(
    stim,
    mtb = "(^|[^[:alnum:]])mtb(aux)?([^[:alnum:]]|$)",
    paste0("(^|[^[:alnum:]])", stim, "([^[:alnum:]]|$)")
  )
}

.acsCytofManualMatchOneFcs <- function(fcs, sampleLookup) {
  fcsKey <- .acsCytofManualNormaliseId(fcs)
  exactIndex <- which(
    nzchar(sampleLookup$sourceKey) &
      (
        stringr::str_detect(fcsKey, stringr::fixed(sampleLookup$sourceKey)) |
          stringr::str_detect(sampleLookup$sourceKey, stringr::fixed(fcsKey))
      )
  )
  if (length(exactIndex) == 1L) {
    return(exactIndex)
  }

  fcsLower <- stringr::str_to_lower(basename(fcs))
  sampleKey <- .acsCytofManualNormaliseId(sampleLookup$SampleID)
  sampleIndex <- which(stringr::str_detect(
    fcsKey,
    stringr::fixed(sampleKey)
  ))
  if (length(sampleIndex) == 0L) {
    return(NA_integer_)
  }

  stimMatch <- vapply(
    sampleLookup$stim[sampleIndex],
    function(stim) {
      stringr::str_detect(
        fcsLower,
        stringr::regex(
          .acsCytofManualStimPattern(stim),
          ignore_case = TRUE
        )
      )
    },
    logical(1)
  )
  sampleIndex <- sampleIndex[stimMatch]

  fcsPid <- stringr::str_extract(fcsLower, "pid[12]")
  if (length(sampleIndex) > 1L && !is.na(fcsPid)) {
    pidMatch <- sampleLookup$pid[sampleIndex] == fcsPid
    pidMatch[is.na(pidMatch)] <- FALSE
    sampleIndex <- sampleIndex[pidMatch]
  }
  if (length(sampleIndex) == 1L) sampleIndex else NA_integer_
}

.acsCytofManualSampleMapFromFcs <- function(
    pathFcsBase,
    popCode,
    sampleLookup) {
  fcsFiles <- .acsCytofFcsFiles(file.path(pathFcsBase, popCode))
  lookupIndex <- vapply(
    fcsFiles,
    .acsCytofManualMatchOneFcs,
    integer(1),
    sampleLookup = sampleLookup
  )
  matched <- sampleLookup[lookupIndex, , drop = FALSE]

  tibble::tibble(
    popCode = popCode,
    ind = as.character(seq_along(fcsFiles)),
    SampleID = matched$SampleID,
    stim = matched$stim
  ) |>
    dplyr::filter(!is.na(.data$SampleID), !is.na(.data$stim))
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
    outputGroup = NULL) {
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
    cyt = NULL) {
  requiredColumns <- c(
    "ind", "cytCombn", "countStim", "nCellStim",
    "countUns", "nCellUns"
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
      "gateName", "method", "popCode", "pop", "batch",
      "ind", "indUns", "nCellStim", "nCellUns"
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
    fn,
    pop = NULL,
    cyt = NULL,
    methods = c("stimgate", "tailgate", "fbeta"),
    outputGroup = NULL,
    gateName = "loc_min",
    sampleMap = NULL) {
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
    hasStimGate <- vapply(popCodes, function(popCode) {
      dir.exists(.acsCytofManualMethodPath(
        pathScratchBase,
        popCode,
        method = "stimgate",
        outputGroup = outputGroup
      ))
    }, logical(1))
    popCodes <- popCodes[hasStimGate]
  }
  if (length(popCodes) == 0L) {
    stop("No requested ACS population has an automated and manual result.")
  }

  if (is.null(sampleMap)) {
    sampleLookup <- .acsCytofManualRead(fn)$sampleLookup
    sampleMap <- purrr::map_dfr(popCodes, function(popCode) {
      .acsCytofManualSampleMapFromFcs(
        pathFcsBase = pathFcsBase,
        popCode = popCode,
        sampleLookup = sampleLookup
      )
    })
  }
  sampleMap <- .acsCytofManualValidateSampleMap(sampleMap, popCodes)

  purrr::map_dfr(popCodes, function(popCode) {
    purrr::map_dfr(methods, function(method) {
      statsTbl <- .acsCytofManualReadStats(
        path = .acsCytofManualMethodPath(
          pathScratchBase,
          popCode,
          method,
          outputGroup
        ),
        method = method,
        gateName = gateName
      )
      .acsCytofStatsSingleMarkers(
        statsTbl = statsTbl,
        method = method,
        popCode = popCode,
        popLabel = unname(popMap[[popCode]]),
        sampleMap = sampleMap,
        cyt = cyt
      )
    })
  })
}

.acsCytofManualComparisonTable <- function(
    fn,
    pathScratchBase,
    pathFcsBase,
    pop = NULL,
    cyt = NULL,
    methods = c("stimgate", "tailgate", "fbeta"),
    outputGroup = NULL,
    gateName = "loc_min",
    sampleMap = NULL) {
  autoTbl <- .acsCytofManualAutoTable(
    pathScratchBase = pathScratchBase,
    pathFcsBase = pathFcsBase,
    fn = fn,
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
  methodLevels <- c("stimgate", "tailgate", "fbeta")
  popLevels <- stats::na.omit(unname(.acsCytofManualPopulationMap()))
  cytLevels <- unname(.acsCytofChannelMap())

  autoTbl |>
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
}

.acsCytofManualSummaryTable <- function(comparisonTbl) {
  comparisonTbl |>
    dplyr::group_by(method, pop, cyt, stim) |>
    dplyr::summarise(
      n = sum(stats::complete.cases(
        .data$freq_bs_auto,
        .data$freq_bs_man
      )),
      pcc = if (dplyr::n() > 1L) {
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

.acsCytofManualPlotScatter <- function(comparisonTbl) {
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
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(pop),
      cols = ggplot2::vars(cyt),
      scales = "free"
    ) +
    ggplot2::labs(
      x = "Background-subtracted frequency (manual gating, %)",
      y = "Background-subtracted frequency (automated gating, %)",
      colour = "Method",
      shape = "Stimulus"
    ) +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy") +
    ggplot2::theme(legend.position = "bottom")
}

.acsCytofManualPlotRelativeError <- function(comparisonTbl) {
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
    ggplot2::labs(x = NULL, y = "Absolute relative error") +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "y") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

.acsCytofManualSave <- function(
    comparisonTbl,
    pathDirSave,
    savePlots = TRUE) {
  dir.create(pathDirSave, recursive = TRUE, showWarnings = FALSE)
  summaryTbl <- .acsCytofManualSummaryTable(comparisonTbl)

  utils::write.csv(
    comparisonTbl,
    file.path(pathDirSave, "manual-comparison.csv"),
    row.names = FALSE
  )
  saveRDS(comparisonTbl, file.path(pathDirSave, "manual-comparison.rds"))
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

comp_against_manual_cyt <- function(
    fn,
    path_scratch_base,
    path_fcs_base,
    pop = NULL,
    cyt = NULL,
    methods = c("stimgate", "tailgate", "fbeta"),
    output_group = NULL,
    gate_name = "loc_min",
    sample_map = NULL,
    path_dir_save = NULL,
    save_plots = TRUE) {
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
