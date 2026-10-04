#' @keywords internal
.checkStimInputPop <- function(popGate) {
  if (anyNA(popGate) || any(popGate != "root")) {
    stop('Only the "root" population exists for non-GatingSet input.')
  }
  invisible(NULL)
}

#' @keywords internal
.asStimGatingSet <- function(.data, popGate = "root") {
  if (inherits(.data, "GatingSet")) {
    return(.data)
  }
  .checkStimInputPop(popGate)

  if (inherits(.data, "cytoframe")) {
    .data <- flowWorkspace::cytoset(list(sample1 = .data))
  } else if (inherits(.data, "flowFrame")) {
    .data <- flowCore::flowSet(list(sample1 = .data))
    flowCore::sampleNames(.data) <- "sample1"
  }
  if (inherits(.data, c("flowSet", "cytoset"))) {
    return(flowWorkspace::GatingSet(.data))
  }

  if (is.character(.data)) {
    files <- .data
    if (length(files) == 1L && !is.na(files) && dir.exists(files)) {
      files <- sort(list.files(
        files, pattern = "\\.fcs$", ignore.case = TRUE, full.names = TRUE,
        recursive = FALSE
      ))
      files <- files[!dir.exists(files)]
    }
    if (length(files) == 0L) {
      stop("No FCS files found in `.data`.")
    }
    missing <- is.na(files) | !file.exists(files) | dir.exists(files)
    if (any(missing)) {
      stop("Missing FCS file(s): ", paste(files[missing], collapse = ", "))
    }
    cs <- flowWorkspace::load_cytoset_from_fcs(files)
    flowWorkspace::sampleNames(cs) <- basename(files)
    return(flowWorkspace::GatingSet(cs))
  }

  if (is.data.frame(.data) && "sample" %in% names(.data)) {
    sample <- .data$sample
    if (anyNA(sample)) {
      stop("The `sample` column must not contain missing values.")
    }
    sampleNames <- if (is.factor(sample)) {
      levels(droplevels(sample))
    } else {
      unique(as.character(sample))
    }
    channels <- .data[, names(.data) != "sample", drop = FALSE]
    .data <- stats::setNames(lapply(sampleNames, function(nm) {
      channels[as.character(sample) == nm, , drop = FALSE]
    }), sampleNames)
  }

  if (is.list(.data) && !is.data.frame(.data)) {
    if (length(.data) == 0L) {
      stop("`.data` must contain at least one sample.")
    }
    frames <- lapply(seq_along(.data), function(i) {
      x <- .data[[i]]
      if (!is.matrix(x) && !is.data.frame(x)) {
        stop("Sample ", i, " must be a numeric matrix or data frame.")
      }
      if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
        stop("Sample ", i, " channel columns must be numeric.")
      }
      m <- as.matrix(x)
      if (!is.numeric(m)) {
        stop("Sample ", i, " channel columns must be numeric.")
      }
      nms <- colnames(m)
      if (
        is.null(nms) || length(nms) == 0L || anyNA(nms) ||
          any(!nzchar(nms)) || anyDuplicated(nms) > 0L
      ) {
        stop("Sample ", i, " must have unique, non-empty column names.")
      }
      firstNames <- colnames(.data[[1]])
      if (!setequal(nms, firstNames)) {
        stop("Sample ", i, " has mismatched channel column names.")
      }
      m <- m[, firstNames, drop = FALSE]
      storage.mode(m) <- "double"
      fr <- flowCore::flowFrame(m)
      params <- flowCore::parameters(fr)
      params@data$desc <- colnames(m)
      flowCore::parameters(fr) <- params
      fr
    })
    sampleNames <- names(.data)
    if (is.null(sampleNames)) {
      sampleNames <- paste0("sample", seq_along(.data))
    }
    if (
      anyNA(sampleNames) || any(!nzchar(sampleNames)) ||
        anyDuplicated(sampleNames) > 0L
    ) {
      stop("Sample names must be unique and non-empty.")
    }
    fs <- flowCore::flowSet(frames)
    flowCore::sampleNames(fs) <- sampleNames
    return(flowWorkspace::GatingSet(fs))
  }

  stop(
    "`.data` must be a GatingSet, flowSet, cytoset, flowFrame, cytoframe, ",
    "FCS directory or file paths, a list of numeric matrices/data frames, ",
    "or a data frame with a `sample` column."
  )
}

#' @keywords internal
.resolveBatchList <- function(batchList, .data) {
  sampleNames <- flowWorkspace::sampleNames(.data)
  lapply(batchList, function(batch) {
    if (!is.character(batch)) {
      return(batch)
    }
    indices <- match(batch, sampleNames)
    if (anyNA(indices)) {
      stop(
        "Unknown sample name(s) in `batchList`: ",
        paste(batch[is.na(indices)], collapse = ", ")
      )
    }
    indices
  })
}
