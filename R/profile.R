# Debug profiling helpers ----------------------------------------------------

.profileState <- new.env(parent = emptyenv())
.profileState$initialized <- FALSE
.profileState$pathProject <- NULL
.profileState$rawDir <- NULL
.profileState$context <- list(
  marker = NA_character_,
  channel = NA_character_,
  batch = NA_character_,
  sample = NA_character_,
  stage = NA_character_
)
.profileState$stack <- character()
.profileState$runTimer <- NULL

#' @keywords internal
.profileEnabled <- function() {
  tolower(trimws(Sys.getenv("STIMGATE_DEBUG"))) %in%
    c("y", "true", "yes", "1")
}

#' @keywords internal
.profilePath <- function(pathProject) {
  file.path(pathProject, "profile")
}

#' @keywords internal
.profileNormalisePath <- function(pathProject) {
  normalizePath(pathProject, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.profileContextDefault <- function() {
  list(
    marker = NA_character_,
    channel = NA_character_,
    batch = NA_character_,
    sample = NA_character_,
    stage = NA_character_
  )
}

#' @keywords internal
.profileMessage <- function(msg) {
  if (.profileEnabled()) {
    message("StimGate profiling: ", msg)
  }
  invisible(FALSE)
}

#' @keywords internal
.profileStateReset <- function() {
  .profileState$initialized <- FALSE
  .profileState$pathProject <- NULL
  .profileState$rawDir <- NULL
  .profileState$context <- .profileContextDefault()
  .profileState$stack <- character()
  .profileState$runTimer <- NULL
  invisible(TRUE)
}

#' Start a fresh profile directory for a StimGate run
#' @keywords internal
.profileInit <- function(pathProject, reset = FALSE) {
  if (!.profileEnabled()) {
    return(invisible(FALSE))
  }

  tryCatch(
    {
      pathProject <- .profileNormalisePath(pathProject)
      if (
        isTRUE(.profileState$initialized) &&
          identical(.profileState$pathProject, pathProject) &&
          !isTRUE(reset)
      ) {
        return(invisible(TRUE))
      }

      pathProfile <- .profilePath(pathProject)
      if (dir.exists(pathProfile)) {
        unlink(pathProfile, recursive = TRUE, force = TRUE)
      }
      rawDir <- file.path(pathProfile, "raw")
      dir.create(rawDir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(rawDir)) {
        stop("could not create profile/raw directory")
      }

      .profileState$initialized <- TRUE
      .profileState$pathProject <- pathProject
      .profileState$rawDir <- rawDir
      .profileState$context <- .profileContextDefault()
      .profileState$stack <- character()
      .profileState$runTimer <- .profileStart(
        level = "run",
        major = "gateStim",
        operation = "gateStim",
        pathProject = pathProject
      )
      invisible(TRUE)
    },
    error = function(e) {
      .profileStateReset()
      .profileMessage(conditionMessage(e))
      invisible(FALSE)
    }
  )
}

#' Attach profiling state without deleting existing raw records
#'
#' This is mainly a safeguard for worker processes that enter a profiled helper
#' after the main process has already created the profile directory.
#' @keywords internal
.profileAttach <- function(pathProject) {
  if (!.profileEnabled()) {
    return(invisible(FALSE))
  }

  tryCatch(
    {
      pathProject <- .profileNormalisePath(pathProject)
      rawDir <- file.path(.profilePath(pathProject), "raw")
      dir.create(rawDir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(rawDir)) {
        stop("could not attach to profile/raw directory")
      }
      .profileState$initialized <- TRUE
      .profileState$pathProject <- pathProject
      .profileState$rawDir <- rawDir
      .profileState$context <- .profileContextDefault()
      .profileState$stack <- character()
      .profileState$runTimer <- NULL
      invisible(TRUE)
    },
    error = function(e) {
      .profileStateReset()
      .profileMessage(conditionMessage(e))
      invisible(FALSE)
    }
  )
}

#' @keywords internal
.profileEnsureAttached <- function(pathProject = NULL) {
  if (!.profileEnabled()) {
    return(FALSE)
  }
  if (isTRUE(.profileState$initialized)) {
    return(TRUE)
  }
  if (is.null(pathProject) || length(pathProject) == 0L) {
    return(FALSE)
  }
  isTRUE(.profileAttach(pathProject))
}

#' Temporarily add profiling context
#' @keywords internal
.profileWithContext <- function(
    expr,
    marker = NULL,
    channel = NULL,
    batch = NULL,
    sample = NULL,
    stage = NULL) {
  if (!.profileEnabled()) {
    return(force(expr))
  }

  oldContext <- .profileState$context
  newContext <- oldContext
  values <- list(
    marker = marker,
    channel = channel,
    batch = batch,
    sample = sample,
    stage = stage
  )
  for (name in names(values)) {
    value <- values[[name]]
    if (!is.null(value) && length(value) > 0L) {
      newContext[[name]] <- as.character(value[[1L]])
    }
  }
  .profileState$context <- newContext
  on.exit({
    .profileState$context <- oldContext
  }, add = TRUE)

  value <- withVisible(force(expr))
  if (isTRUE(value$visible)) {
    value$value
  } else {
    invisible(value$value)
  }
}

#' @keywords internal
.profileStart <- function(
    level,
    major,
    minor = NA_character_,
    operation = NA_character_,
    pathProject = NULL) {
  if (!.profileEnsureAttached(pathProject)) {
    return(NULL)
  }

  rawDir <- .profileState$rawDir
  pathRecord <- tempfile(
    pattern = "profile-",
    tmpdir = rawDir,
    fileext = ".rds"
  )
  recordId <- sub("\\.rds$", "", basename(pathRecord))
  stack <- .profileState$stack
  depth <- length(stack)
  parentId <- if (depth == 0L) NA_character_ else stack[[depth]]
  context <- .profileState$context

  timer <- list(
    recordId = recordId,
    parentId = parentId,
    pathRecord = pathRecord,
    depth = depth,
    level = as.character(level[[1L]]),
    major = as.character(major[[1L]]),
    minor = as.character(minor[[1L]]),
    operation = as.character(operation[[1L]]),
    marker = context$marker,
    channel = context$channel,
    batch = context$batch,
    sample = context$sample,
    stage = context$stage,
    pid = Sys.getpid(),
    startedAt = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%OS3%z"
    ),
    startedElapsed = proc.time()[["elapsed"]]
  )
  .profileState$stack <- c(stack, recordId)
  timer
}

#' @keywords internal
.profilePop <- function(recordId) {
  stack <- .profileState$stack
  if (length(stack) == 0L) {
    return(invisible(FALSE))
  }
  matchInd <- which(stack == recordId)
  if (length(matchInd) == 0L) {
    return(invisible(FALSE))
  }
  removeInd <- matchInd[[length(matchInd)]]
  .profileState$stack <- stack[-removeInd]
  invisible(TRUE)
}

#' @keywords internal
.profileWriteRecord <- function(record, pathRecord) {
  tryCatch(
    {
      saveRDS(record, pathRecord)
      invisible(TRUE)
    },
    error = function(e) {
      .profileMessage(
        paste0("could not write timing record: ", conditionMessage(e))
      )
      invisible(FALSE)
    }
  )
}

#' Finish one profiling timer and persist it immediately
#' @keywords internal
.profileStop <- function(timer) {
  if (is.null(timer)) {
    return(invisible(NULL))
  }

  elapsed <- proc.time()[["elapsed"]] - timer$startedElapsed
  finishedAt <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%OS3%z"
  )
  .profilePop(timer$recordId)

  record <- data.frame(
    record_id = timer$recordId,
    parent_id = timer$parentId,
    depth = as.integer(timer$depth),
    level = timer$level,
    major = timer$major,
    minor = timer$minor,
    operation = timer$operation,
    marker = timer$marker,
    channel = timer$channel,
    batch = timer$batch,
    sample = timer$sample,
    stage = timer$stage,
    pid = as.integer(timer$pid),
    elapsed_sec = as.numeric(elapsed),
    started_at = timer$startedAt,
    finished_at = finishedAt,
    status = "completed",
    stringsAsFactors = FALSE
  )
  .profileWriteRecord(record, timer$pathRecord)
  invisible(record)
}

#' Discard an unfinished timer while restoring the hierarchy stack
#' @keywords internal
.profileCancel <- function(timer) {
  if (!is.null(timer)) {
    .profilePop(timer$recordId)
  }
  invisible(FALSE)
}

#' Time one expression and persist the completed timing record
#' @keywords internal
.profileTime <- function(
    expr,
    level,
    major,
    minor = NA_character_,
    operation = NA_character_,
    pathProject = NULL) {
  if (!.profileEnabled()) {
    return(force(expr))
  }

  timer <- .profileStart(
    level = level,
    major = major,
    minor = minor,
    operation = operation,
    pathProject = pathProject
  )
  if (is.null(timer)) {
    return(force(expr))
  }

  completed <- FALSE
  on.exit({
    if (!completed) {
      .profileCancel(timer)
    }
  }, add = TRUE)

  value <- withVisible(force(expr))
  .profileStop(timer)
  completed <- TRUE
  if (isTRUE(value$visible)) {
    value$value
  } else {
    invisible(value$value)
  }
}

#' @keywords internal
.profileDataBatch <- function(.data) {
  batch <- attr(.data, "batch")
  if (!is.null(batch) && length(batch) > 0L && !is.na(batch[[1L]])) {
    return(as.character(batch[[1L]]))
  }
  if (
    is.data.frame(.data) &&
      "batch" %in% names(.data) &&
      nrow(.data) > 0L &&
      !is.na(.data$batch[[1L]])
  ) {
    return(as.character(.data$batch[[1L]]))
  }
  NA_character_
}

#' @keywords internal
.profileInitialSampleActive <- function() {
  context <- .profileState$context
  identical(context$stage, "init") &&
    !is.na(context$sample) &&
    nzchar(context$sample)
}

#' Merge the incrementally saved timing records
#' @keywords internal
.profileFinalise <- function(pathProject = .profileState$pathProject) {
  if (!.profileEnabled() || is.null(pathProject)) {
    return(invisible(NULL))
  }

  tryCatch(
    {
      pathProfile <- .profilePath(.profileNormalisePath(pathProject))
      rawDir <- file.path(pathProfile, "raw")
      files <- list.files(
        rawDir,
        pattern = "\\.rds$",
        full.names = TRUE
      )
      if (length(files) == 0L) {
        return(invisible(NULL))
      }

      records <- lapply(files, function(path) {
        tryCatch(readRDS(path), error = function(e) NULL)
      })
      records <- Filter(Negate(is.null), records)
      if (length(records) == 0L) {
        return(invisible(NULL))
      }

      profileTbl <- do.call(rbind, records)
      ord <- order(
        profileTbl$started_at,
        profileTbl$depth,
        profileTbl$record_id
      )
      profileTbl <- profileTbl[ord, , drop = FALSE]
      rownames(profileTbl) <- NULL

      saveRDS(profileTbl, file.path(pathProfile, "profile.rds"))
      utils::write.csv(
        profileTbl,
        file.path(pathProfile, "profile.csv"),
        row.names = FALSE,
        na = ""
      )
      invisible(profileTbl)
    },
    error = function(e) {
      .profileMessage(
        paste0("could not collate profiling records: ", conditionMessage(e))
      )
      invisible(NULL)
    }
  )
}

#' Finish the run timer and collate the final profile
#' @keywords internal
.profileFinishRun <- function(pathProject = .profileState$pathProject) {
  if (!.profileEnabled() || !isTRUE(.profileState$initialized)) {
    return(invisible(NULL))
  }

  runTimer <- .profileState$runTimer
  if (!is.null(runTimer)) {
    .profileStop(runTimer)
    .profileState$runTimer <- NULL
  }
  profileTbl <- .profileFinalise(pathProject)
  .profileStateReset()
  invisible(profileTbl)
}
