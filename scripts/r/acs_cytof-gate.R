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
  if (!is.null(outputGroup) &&
    (!is.character(outputGroup) ||
      length(outputGroup) != 1L ||
      !nzchar(outputGroup))) {
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
    stimgateCheck = file.path(pathScratch, "stimgate_check.pdf")
  )
}

.acsCytofBatchList <- function(nSample) {
  if (!is.numeric(nSample) ||
    length(nSample) != 1L ||
    is.na(nSample) ||
    !is.finite(nSample) ||
    nSample != as.integer(nSample) ||
    nSample < 5L) {
    stop("nSample must be one integer of at least 5.")
  }

  nSample <- as.integer(nSample)
  nPerBatch <- 5L
  if (nSample %% nPerBatch != 0L) {
    stop(
      "The ACS CyTOF sample count must be a multiple of five so that ",
      "every batch contains one unstimulated and four stimulated samples."
    )
  }

  lapply(seq.int(1L, nSample, by = nPerBatch), function(indUns) {
    seq.int(indUns, length.out = nPerBatch)
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
      "Expected at least 10 ACS CyTOF FCS files in ", pathFcs,
      ", but found only ", length(fcsFiles), "."
    )
  }

  sort(fcsFiles)
}

.acsCytofPreprocessPopulation <- function(paths, nSample = NULL, runPlots = TRUE) {
  fcsFiles <- .acsCytofFcsFiles(paths$fcs)
  if (!is.null(nSample)) {
    .acsCytofBatchList(nSample)
    if (nSample > length(fcsFiles)) {
      stop(
        "Requested ", nSample, " tester samples, but only ",
        length(fcsFiles), " FCS files are available."
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
        projr::projr_path_get("project"),
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

.acsCytofRunPopulation <- function(
  pop,
  pathFcsBase,
  pathGsBase,
  pathScratchBase,
  runPreprocessing,
  runMethods,
  runPlots,
  nSample = NULL,
  biasUns = 0.15,
  outputGroup = NULL
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
      runPlots = runPlots
    )
  }

  if (!isTRUE(runMethods) && !isTRUE(runPlots)) {
    return(invisible(list(pop = pop, paths = paths)))
  }
  if (!dir.exists(paths$gs)) {
    stop(
      "Cached ACS CyTOF GatingSet not found for '", pop, "' at: ",
      paths$gs, ". Run preprocessing for this population first."
    )
  }

  .acsCytofEnsureCurrentCheckout()
  gs <- flowWorkspace::load_gs(paths$gs)
  nSampleActual <- length(gs)
  if (!is.null(nSample) && nSampleActual != nSample) {
    stop(
      "Cached tester GatingSet for '", pop, "' contains ", nSampleActual,
      " samples; expected ", nSample, ". Re-run tester preprocessing."
    )
  }
  batchList <- .acsCytofBatchList(nSampleActual)

  if (isTRUE(runMethods)) {
    if (dir.exists(paths$stimgate)) {
      unlink(paths$stimgate, recursive = TRUE)
    }
    dir.create(paths$stimgate, recursive = TRUE, showWarnings = FALSE)

    restoreDebug <- .acsCytofSetDebug()
    on.exit(restoreDebug(), add = TRUE)

    invisible(stimgate::gateStim(
      pathProject = paths$stimgate,
      .data = gs,
      popGate = "root",
      batchList = batchList,
      chnl = c(
        "Ho165Di", "Gd158Di", "Nd146Di",
        "Dy164Di", "Gd156Di", "Nd150Di"
      ),
      biasUns = biasUns,
      bwMtd = "hpi1",
      bwNcellMax = 1e4,
      bwFallback = "auto",
      bwMin = "none",
      bwMax = "none",
      minCell = 100,
      gateCombn = "min",
      tolClust = NULL,
      calcCytPosGates = TRUE
    ))
  }

  if (isTRUE(runPlots)) {
    if (!dir.exists(paths$stimgate)) {
      stop(
        "StimGate output not found for '", pop, "' at: ", paths$stimgate,
        ". Run the methods for this population first."
      )
    }

    p <- stimgate::plotStim(
      ind = 2L,
      .data = gs,
      pathProject = paths$stimgate,
      pop = "root",
      chnl = c("Ho165Di", "Nd146Di")
    )
    dir.create(dirname(paths$stimgateCheck), recursive = TRUE, showWarnings = FALSE)
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
