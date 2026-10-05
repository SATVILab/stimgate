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

.acsCytofReadPreprocessing <- function(path, gs = NULL) {
  file <- file.path(path, "acs-preprocessing.rds")
  if (!file.exists(file)) stop("ACS preprocessing manifest missing; re-run preprocessing: ", path)
  manifest <- readRDS(file)
  .acsCytofBatchList(manifest$sampleMap)
  if (!is.null(gs) && !identical(basename(flowWorkspace::sampleNames(gs)), manifest$sampleMap$file)) {
    stop("ACS GatingSet files are missing, duplicated or reordered relative to preprocessing; re-run preprocessing.")
  }
  manifest
}

.acsCytofManifest <- function(preprocessing) {
  sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  if (length(sha) != 1L || !nzchar(sha)) stop("Cannot record ACS git SHA.")
  list(version = 1L, gitSha = sha, preprocessing = preprocessing)
}

.acsCytofValidateManifests <- function(manifests) {
  contexts <- lapply(manifests, function(x) x$context)
  if (!length(contexts) || any(vapply(contexts, is.null, logical(1))) ||
      !all(vapply(contexts, identical, logical(1), contexts[[1]]))) {
    stop("Mismatched ACS result manifests (data, preprocessing or git SHA). Re-run all methods together.")
  }
  invisible(TRUE)
}
