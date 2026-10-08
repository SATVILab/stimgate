.acsCytofPopulationPaths <- function(
  pop,
  pathFcsBase,
  pathGsBase,
  pathScratchBase,
  outputGroup = NULL
) {
  if (!is.character(pop) || length(pop) != 1L || !nzchar(pop)) {
    stop("pop must be one non-empty character value.")
  }
  if (
    !is.null(outputGroup) &&
      (!is.character(outputGroup) ||
        length(outputGroup) != 1L ||
        !nzchar(outputGroup))
  ) {
    stop("outputGroup must be NULL or one non-empty character value.")
  }

  outputParts <- if (is.null(outputGroup)) pop else c(outputGroup, pop)
  pathGs <- do.call(file.path, as.list(c(pathGsBase, outputParts)))
  pathScratch <- do.call(file.path, as.list(c(pathScratchBase, outputParts)))

  list(
    fcs = file.path(pathFcsBase, pop),
    gs = pathGs,
    scratch = pathScratch,
    gsCheck = file.path(pathScratch, "gatingset"),
    stimgate = file.path(pathScratch, "stimgate"),
    tailgate = file.path(pathScratch, "tailgate", "result.rds"),
    fbeta = file.path(pathScratch, "fbeta", "result.rds"),
    stimgateCheck = file.path(pathScratch, "stimgate_check_sample_2.pdf")
  )
}

.acsCytofValidateSampleCount <- function(nSample) {
  if (length(nSample) != 1L || !is.numeric(nSample) ||
      !is.finite(nSample) || nSample < 5L || nSample %% 5L != 0L) {
    stop("ACS sample count must be a positive multiple of five (at least 5).")
  }
  invisible(nSample)
}

.acsCytofBatchList <- function(sampleMap) {
  if (!is.data.frame(sampleMap) ||
      !all(c("SampleID", "stim", "ind") %in% names(sampleMap)) ||
      anyNA(sampleMap[c("SampleID", "stim", "ind")]) ||
      nrow(sampleMap) == 0L ||
      anyDuplicated(sampleMap$ind) ||
      any(!nzchar(as.character(sampleMap$SampleID))) ||
      anyNA(suppressWarnings(as.integer(sampleMap$ind)))) {
    stop("ACS batches require a complete mapped SampleID/stim/ind table.")
  }
  stimuli <- c("uns", "p1", "mtb", "ebv", "p4")
  batches <- split(sampleMap, sampleMap$SampleID)
  lapply(batches, function(batch) {
    if (nrow(batch) != length(stimuli) ||
        anyDuplicated(batch$stim) || !setequal(batch$stim, stimuli)) {
      stop("ACS batch '", batch$SampleID[[1]],
           "' must contain exactly one tube for each of: ",
           paste(stimuli, collapse = ", "), ". Missing or duplicate tubes.")
    }
    as.integer(batch$ind[match(stimuli, batch$stim)])
  })
}

.acsCytofFcsFiles <- function(pathFcs) {
  if (!dir.exists(pathFcs)) {
    stop("ACS CyTOF FCS directory not found at: ", pathFcs)
  }

  fcsFiles <- list.files(
    pathFcs,
    pattern = "\\.fcs$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(fcsFiles) < 10L) {
    stop(
      "Expected at least 10 ACS CyTOF FCS files in ",
      pathFcs,
      ", but found only ",
      length(fcsFiles),
      "."
    )
  }

  sort(fcsFiles)
}

.acsCytofPreprocessPopulation <- function(
  paths,
  nSample = NULL,
  runPlots = TRUE
) {
  fcsFiles <- .acsCytofFcsFiles(paths$fcs)
  if (!is.null(nSample)) {
    .acsCytofValidateSampleCount(nSample)
    if (nSample > length(fcsFiles)) {
      stop(
        "Requested ",
        nSample,
        " tester samples, but only ",
        length(fcsFiles),
        " FCS files are available."
      )
    }
  }

  create_gatingset(
    path_fcs = paths$fcs,
    path_gs = paths$gs,
    n_sample = nSample
  )

  if (isTRUE(runPlots)) {
    plot_gatingset_check(
      path_gs = paths$gs,
      path_plot_dir = paths$gsCheck
    )
  }

  invisible(paths$gs)
}

.acsCytofEnsureCurrentCheckout <- function(pathRoot = NULL) {
  if (is.null(pathRoot) || !nzchar(pathRoot)) {
    if (requireNamespace("projr", quietly = TRUE)) {
      pathRoot <- tryCatch(
        projr::projr_path_get("project", format = "absolute"),
        error = function(e) normalizePath(".", winslash = "/", mustWork = FALSE)
      )
    } else {
      pathRoot <- normalizePath(".", winslash = "/", mustWork = FALSE)
    }
  }
  pathRoot <- normalizePath(pathRoot, winslash = "/", mustWork = FALSE)

  if (requireNamespace("stimgate", quietly = TRUE)) {
    ns <- tryCatch(asNamespace("stimgate"), error = function(e) NULL)
    if (!is.null(ns)) {
      nsPath <- normalizePath(
        getNamespaceInfo(ns, "path"),
        winslash = "/",
        mustWork = FALSE
      )
      if (isTRUE(identical(nsPath, pathRoot))) {
        return(invisible(TRUE))
      }
    }
  }

  if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(pathRoot, quiet = TRUE))
  } else if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(pathRoot, quiet = TRUE))
  } else {
    stop(
      "Neither pkgload nor devtools is available to load the current ",
      "StimGate checkout."
    )
  }

  invisible(TRUE)
}

.acsCytofSetDebug <- function() {
  oldDebug <- Sys.getenv("STIMGATE_DEBUG", unset = NA_character_)
  Sys.setenv(STIMGATE_DEBUG = "TRUE")

  function() {
    if (is.na(oldDebug)) {
      Sys.unsetenv("STIMGATE_DEBUG")
    } else {
      Sys.setenv(STIMGATE_DEBUG = oldDebug)
    }
  }
}

.acsCytofRunPopulationSafe <- function(...) {
  pop <- list(...)$pop

  tryCatch(
    {
      result <- .acsCytofRunPopulation(...)
      list(
        pop = pop,
        success = TRUE,
        result = result,
        error = NULL
      )
    },
    error = function(e) {
      list(
        pop = pop,
        success = FALSE,
        result = NULL,
        error = conditionMessage(e)
      )
    }
  )
}

.acsCytofRunPopulation <- function(
  pop,
  pathFcsBase,
  pathGsBase,
  pathScratchBase,
  runPreprocessing,
  runMethods,
  runPlots,
  nSample = NULL,
  biasUns = NULL,
  biasUnsFactor = 1,
  bwMtd = "nrd0",
  bwScope = "cytokine",
  locThresholdMethod = "region",
  outputGroup = NULL,
  runPreprocessingPlots = FALSE
) {
  paths <- .acsCytofPopulationPaths(
    pop = pop,
    pathFcsBase = pathFcsBase,
    pathGsBase = pathGsBase,
    pathScratchBase = pathScratchBase,
    outputGroup = outputGroup
  )

  if (isTRUE(runPreprocessing)) {
    .acsCytofPreprocessPopulation(
      paths = paths,
      nSample = nSample,
      runPlots = runPreprocessingPlots
    )
  }

  if (!isTRUE(runMethods) && !isTRUE(runPlots)) {
    return(invisible(list(pop = pop, paths = paths)))
  }
  if (!dir.exists(paths$gs)) {
    stop(
      "Cached ACS CyTOF GatingSet not found for '",
      pop,
      "' at: ",
      paths$gs,
      ". Run preprocessing for this population first."
    )
  }

  .acsCytofEnsureCurrentCheckout()
  gs <- flowWorkspace::load_gs(paths$gs)
  nSampleActual <- length(gs)
  if (!is.null(nSample) && nSampleActual != nSample) {
    stop(
      "Cached tester GatingSet for '",
      pop,
      "' contains ",
      nSampleActual,
      " samples; expected ",
      nSample,
      ". Re-run tester preprocessing."
    )
  }
  preprocessing <- .acsCytofReadPreprocessing(paths$gs, gs)
  batchList <- .acsCytofBatchList(preprocessing$sampleMap)

  if (isTRUE(runMethods)) {
    restoreDebug <- .acsCytofSetDebug()
    on.exit(restoreDebug(), add = TRUE)

    # Gate into a temporary sibling so a failed run keeps the last good output.
    .acsCytofReplaceDir(paths$stimgate, function(pathTmp) {
      stimgate::gateStim(
        pathProject = pathTmp,
        .data = gs,
        popGate = "root",
        batchList = batchList,
        chnl = c(
          "Ho165Di",
          "Gd158Di",
          "Nd146Di",
          "Dy164Di",
          "Gd156Di",
          "Nd150Di"
        ),
        # NULL or NA biasUns: StimGate sets it to biasUnsFactor times the
        # shared bandwidth, i.e. the bandwidth itself with factor 1.
        biasUns = if (length(biasUns) == 1L && is.na(biasUns)) NULL else biasUns,
        control = stimgate::stimControl(
          biasUnsFactor = biasUnsFactor,
          bwMtd = bwMtd,
          bwScope = bwScope,
          bwNcellMax = 1e4,
          bwFallback = "auto",
          bwMin = "none",
          bwMax = "none",
          gateCombn = "min",
          clusterGates = TRUE,
          calcCytPosGates = TRUE,
          minCell = 100,
          locThresholdMethod = locThresholdMethod
        )
      )
      manifest <- list(
        context = .acsCytofManifest(preprocessing),
        settings = list(biasUns = biasUns, biasUnsFactor = biasUnsFactor,
                        clusterGates = TRUE, calcCytPosGates = TRUE,
                        locThresholdMethod = locThresholdMethod),
        channelSettings = stimgate::stimgateMetaReadSettingsChnls(pathTmp)
      )
      # Fail before the swap if any channel resolved to another method.
      .acsCytofValidateStimGateMethod(manifest, locThresholdMethod)
      saveRDS(manifest, file.path(pathTmp, "acs-manifest.rds"))
    })
  }

  if (isTRUE(runPlots)) {
    if (!dir.exists(paths$stimgate)) {
      stop(
        "StimGate output not found for '",
        pop,
        "' at: ",
        paths$stimgate,
        ". Run the methods for this population first."
      )
    }

    p <- .acsCytofPlotGateCheck(gs, paths$stimgate)
    dir.create(
      dirname(paths$stimgateCheck),
      recursive = TRUE,
      showWarnings = FALSE
    )
    ggplot2::ggsave(
      filename = paths$stimgateCheck,
      plot = p,
      width = 18,
      height = 16,
      units = "cm"
    )
  }

  invisible(list(
    pop = pop,
    nSample = nSampleActual,
    batchList = batchList,
    paths = paths
  ))
}

# Remove package diagnostic panel labels before assembling the analysis figure.
.acsCytofPlotGateCheck <- function(gs, pathProject) {
  plots <- stimgate::plotStim(
    ind = 2L, .data = gs, pathProject = pathProject,
    pop = "root", chnl = c("Ho165Di", "Nd146Di"), grid = FALSE
  )
  if (length(plots) == 0L) return(NULL)
  plots <- lapply(plots, function(p) {
    p$labels$title <- NULL
    p$labels$subtitle <- NULL
    p
  })
  cowplot::plot_grid(plotlist = plots, ncol = 2, align = "hv") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white"),
      panel.background = ggplot2::element_rect(fill = "white")
    )
}
