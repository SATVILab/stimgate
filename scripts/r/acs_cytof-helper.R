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
