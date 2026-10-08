.acs_get_expr <- function(gs, ind, chnl, pop = "root") {
  gh <- gs[[as.integer(ind)]]
  ff <- flowWorkspace::gh_pop_get_data(gh, pop)
  ex <- flowCore::exprs(ff)
  if (!chnl %in% colnames(ex)) {
    stop("Channel '", chnl, "' is not present in sample index ", ind, ".")
  }
  as.numeric(ex[, chnl])
}

# Build a directory in a temporary sibling, then swap it in. A failed build
# leaves the previous directory untouched.
.acsCytofReplaceDir <- function(path, build) {
  pathTmp <- paste0(path, ".tmp-", Sys.getpid())
  pathOld <- paste0(path, ".old-", Sys.getpid())
  unlink(c(pathTmp, pathOld), recursive = TRUE)
  on.exit(unlink(c(pathTmp, pathOld), recursive = TRUE), add = TRUE)
  dir.create(pathTmp, recursive = TRUE, showWarnings = FALSE)

  build(pathTmp)

  hadOld <- dir.exists(path)
  if (hadOld && !file.rename(path, pathOld)) {
    stop("Could not move the previous output aside: ", path)
  }
  if (!file.rename(pathTmp, path)) {
    if (hadOld) file.rename(pathOld, path)
    stop("Could not move the new output into place: ", path)
  }

  invisible(path)
}

# Run `fn` over `popVec` with a multisession plan (sequential for one worker)
# and stop with every failed population's message. `fn` must return the
# list(pop, success, error) shape of the *Safe runners.
.acsCytofMapPopulations <- function(popVec, fn, nWorkers, seed, label) {
  nWorkersUse <- min(nWorkers, length(popVec))
  oldPlan <- future::plan()

  runList <- tryCatch(
    {
      if (nWorkersUse > 1L) {
        future::plan(future::multisession, workers = nWorkersUse)
      } else {
        future::plan(future::sequential)
      }
      furrr::future_map(
        popVec,
        fn,
        .options = furrr::furrr_options(seed = seed, scheduling = Inf)
      )
    },
    finally = future::plan(oldPlan)
  )
  names(runList) <- popVec

  failed <- !vapply(runList, function(x) isTRUE(x$success), logical(1))
  if (any(failed)) {
    stop(
      label,
      " failed:\n",
      paste(
        vapply(
          runList[failed],
          function(x) paste0(x$pop, ": ", x$error),
          character(1)
        ),
        collapse = "\n"
      )
    )
  }

  invisible(runList)
}

# Identity is independent of directory location, but preserves loaded tube order.
.acsCytofHash <- function(object) {
  path <- tempfile("acs-hash-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2)
  unname(tools::md5sum(path))
}

.acsCytofMapFiles <- function(files, sampleLookup = NULL) {
  if (is.null(sampleLookup)) {
    env <- new.env()
    utils::data("fn_to_sampleid_map", package = "DataTidyACSCyTOFFAUST", envir = env)
    sampleLookup <- env$fn_to_sampleid_map
  }
  clean <- DataTidyACSCyTOFFAUST::clean_fcs_for_matching(basename(files))
  if (anyDuplicated(clean) || anyDuplicated(sampleLookup$MatchFCSName)) {
    stop("Duplicate ACS FCS filenames or filename lookup keys.")
  }
  index <- match(clean, sampleLookup$MatchFCSName)
  if (anyNA(index)) {
    stop("Unmapped ACS FCS files: ", paste(basename(files)[is.na(index)], collapse = ", "))
  }
  tibble::tibble(
    ind = as.character(seq_along(files)),
    file = basename(files),
    SampleID = as.character(sampleLookup$SampleID[index]),
    stim = sub("^mtbaux$", "mtb", as.character(sampleLookup$Stim[index]))
  )
}

# Stored beside, not inside, the GatingSet folder: flowWorkspace::load_gs()
# rejects any extra file or folder in it (an ".rds" file as a legacy archive).
.acsCytofPreprocessingFile <- function(path) {
  paste0(path, ".acs-preprocessing.rds")
}

.acsCytofReadPreprocessing <- function(path, gs = NULL) {
  file <- .acsCytofPreprocessingFile(path)
  if (!file.exists(file)) stop("ACS preprocessing manifest missing; re-run preprocessing: ", path)
  manifest <- readRDS(file)
  .acsCytofBatchList(manifest$sampleMap)
  if (!is.null(gs) && !identical(basename(flowWorkspace::sampleNames(gs)), manifest$sampleMap$file)) {
    stop("ACS GatingSet files are missing, duplicated or reordered relative to preprocessing; re-run preprocessing.")
  }
  manifest
}

.acsCytofManifest <- function(preprocessing) {
  list(version = 1L, gitSha = .git_sha(), preprocessing = preprocessing)
}

.acsCytofValidateManifests <- function(manifests) {
  contexts <- lapply(manifests, function(x) {
    context <- x$context
    context$gitSha <- NULL
    context
  })
  if (!length(contexts) || any(vapply(contexts, is.null, logical(1))) ||
      !all(vapply(contexts, identical, logical(1), contexts[[1]]))) {
    stop("Mismatched ACS result manifests (data or preprocessing). Re-run all methods together.")
  }
  invisible(TRUE)
}

# A saved StimGate result must record the requested local-FDR threshold
# method, and every channel's saved settings must have resolved to it.
# Results made before the method was recorded (matched gates) are rejected.
.acsCytofValidateStimGateMethod <- function(manifest, locThresholdMethod) {
  recorded <- manifest$settings$locThresholdMethod
  resolved <- vapply(
    manifest$channelSettings,
    function(x) {
      method <- if (is.list(x)) x[["locThresholdMethod"]] else NULL
      if (length(method) == 1L) as.character(method) else NA_character_
    },
    character(1)
  )
  if (!identical(recorded, locThresholdMethod) || !length(resolved) ||
      !all(resolved %in% locThresholdMethod)) {
    stop(
      "ACS StimGate result was not made with locThresholdMethod = '",
      locThresholdMethod, "' (recorded: '",
      paste(recorded, collapse = ", "), "'); re-run all methods."
    )
  }
  invisible(TRUE)
}

# The analysis 9 run manifest records the StimGate threshold method used for
# the saved comparison; a manifest without it (or with another) is rejected.
.acsCytofCheckRunManifestMethod <- function(
    manifest, locThresholdMethod, qmdPath) {
  if (!identical(
    manifest$stimgate_loc_threshold_method, locThresholdMethod
  )) {
    .analysis_cache_error(
      "acs_cytof",
      paste0(
        "Cached stimgate_loc_threshold_method does not match '",
        locThresholdMethod, "'."
      ),
      qmdPath
    )
  }
  invisible(manifest)
}

# Read a saved StimGate project's ACS manifest, checking that it matches the
# project's saved channel settings and the requested threshold method.
.acsCytofReadStimGateManifest <- function(path, locThresholdMethod) {
  pathManifest <- file.path(path, "acs-manifest.rds")
  if (!file.exists(pathManifest)) {
    stop("ACS StimGate manifest missing; re-run all methods.")
  }
  manifest <- readRDS(pathManifest)
  if (!identical(
    manifest$channelSettings, stimgate::stimgateMetaReadSettingsChnls(path)
  )) {
    stop("Mismatched ACS StimGate settings manifest; re-run all methods.")
  }
  .acsCytofValidateStimGateMethod(manifest, locThresholdMethod)
  manifest
}

.acsCytofValidateComparisonManifest <- function(
    table, locThresholdMethod = "region") {
  manifest <- attr(table, "manifest")
  if (is.null(manifest$methods) || !length(manifest$methods) ||
      is.null(manifest$comparisonSettings) || is.null(manifest$manualInputHash) ||
      !"thresholdFailed" %in% names(table)) {
    stop("Legacy or incomplete ACS comparison manifest; re-run analysis 9 with RUN_SIMULATIONS=true RUN_PLOTS=false.")
  }
  if (!identical(
    manifest$comparisonSettings$locThresholdMethod, locThresholdMethod
  )) {
    stop(
      "ACS comparison does not record StimGate locThresholdMethod = '",
      locThresholdMethod, "'; re-run analysis 9 with ",
      "RUN_SIMULATIONS=true RUN_PLOTS=false."
    )
  }
  for (population in manifest$methods) {
    .acsCytofValidateManifests(population)
    if (!is.null(population$stimgate)) {
      .acsCytofValidateStimGateMethod(population$stimgate, locThresholdMethod)
    }
  }
  isStimGate <- as.character(table$method) %in% "stimgate"
  if (any(isStimGate) && (!"locThresholdMethod" %in% names(table) ||
    !all(table$locThresholdMethod[isStimGate] %in% locThresholdMethod))) {
    stop(
      "ACS StimGate rows do not record locThresholdMethod = '",
      locThresholdMethod, "'; re-run analysis 9 with ",
      "RUN_SIMULATIONS=true RUN_PLOTS=false."
    )
  }
  invisible(TRUE)
}
