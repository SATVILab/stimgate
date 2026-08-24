.acsCytofChannelMap <- function() {
  c(
    Ho165Di = "IFNg",
    Gd158Di = "IL2",
    Nd146Di = "TNF",
    Dy164Di = "IL17",
    Gd156Di = "IL6",
    Nd150Di = "IL22"
  )
}

.acsCytofComparatorSettings <- function(method = c("tailgate", "fbeta")) {
  method <- match.arg(method)

  params <- switch(
    method,
    tailgate = list(
      tailgateX = "stim",
      adjust = 1,
      bandwidth = NULL,
      numPeaks = 1L,
      refPeak = 1L,
      derivativeMethod = "firstDeriv",
      tol = 1e-2,
      side = "right",
      strict = FALSE,
      autoTol = TRUE
    ),
    fbeta = list(
      beta = 0.8,
      theta = 2,
      width = 10L,
      numBins = NULL,
      patchPy2Compat = TRUE
    )
  )

  list(
    cacheVersion = if (identical(method, "fbeta")) 2L else 1L,
    method = method,
    channelMap = .acsCytofChannelMap(),
    params = params
  )
}

.acsCytofCombinationLevels <- function(channels) {
  if (length(channels) == 0L || anyNA(channels) || any(!nzchar(channels))) {
    stop("channels must contain at least one non-empty channel name.")
  }

  levelGrid <- expand.grid(
    rep(list(c(FALSE, TRUE)), length(channels)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  apply(levelGrid, 1L, function(isPositive) {
    paste0(
      paste0(channels, ifelse(isPositive, "~+~", "~-~")),
      collapse = ""
    )
  })
}

.acsCytofCombinationLabels <- function(isPositive, channels) {
  isPositive <- as.matrix(isPositive)
  if (ncol(isPositive) != length(channels)) {
    stop("The positivity matrix must have one column per channel.")
  }
  if (nrow(isPositive) == 0L) {
    return(character())
  }

  tokenList <- lapply(seq_along(channels), function(i) {
    paste0(channels[[i]], ifelse(isPositive[, i], "~+~", "~-~"))
  })
  do.call(paste0, tokenList)
}

.acsCytofCombinationCounts <- function(
  xStim,
  xUns,
  thresholds,
  channels = names(thresholds)
) {
  xStim <- as.matrix(xStim)
  xUns <- as.matrix(xUns)
  thresholds <- as.numeric(thresholds)

  if (
    ncol(xStim) != length(channels) ||
      ncol(xUns) != length(channels) ||
      length(thresholds) != length(channels)
  ) {
    stop(
      "Expression matrices, thresholds and channels do not have matching widths."
    )
  }
  if (any(!is.finite(thresholds))) {
    stop(
      "All thresholds must be finite before combination counts are calculated."
    )
  }

  combinationLevels <- .acsCytofCombinationLevels(channels)
  labelStim <- .acsCytofCombinationLabels(
    sweep(xStim, 2L, thresholds, FUN = ">"),
    channels = channels
  )
  labelUns <- .acsCytofCombinationLabels(
    sweep(xUns, 2L, thresholds, FUN = ">"),
    channels = channels
  )

  tibble::tibble(
    cytCombn = combinationLevels,
    countStim = as.integer(table(factor(
      labelStim,
      levels = combinationLevels
    ))),
    nCellStim = nrow(xStim),
    countUns = as.integer(table(factor(labelUns, levels = combinationLevels))),
    nCellUns = nrow(xUns)
  )
}

.acsCytofExpressionMatrix <- function(gs, ind, channels) {
  exprList <- lapply(channels, function(chnl) {
    .acs_get_expr(gs = gs, ind = ind, chnl = chnl, pop = "root")
  })
  nCell <- lengths(exprList)
  if (length(unique(nCell)) != 1L) {
    stop(
      "Cytokine channels returned different cell counts for sample index ",
      ind,
      "."
    )
  }

  out <- do.call(cbind, exprList)
  colnames(out) <- channels
  out
}

.acsCytofFallbackThreshold <- function(xStim, xUns, margin = 0.05) {
  if (exists(".simCompareHighValueGate", mode = "function", inherits = TRUE)) {
    return(.simCompareHighValueGate(
      xStim = xStim,
      xUns = xUns,
      margin = margin
    ))
  }

  x <- c(xStim, xUns)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(0)
  }
  rng <- range(x)
  rng[[2]] + max(1, diff(rng)) * margin
}

.acsCytofThresholdOne <- function(
  method,
  xUns,
  xStim,
  settings,
  pathFbeta = NULL,
  fbetaEnv = NULL
) {
  thresholdObj <- if (identical(method, "fbeta")) {
    if (
      !exists(
        ".simCompareFbetaThreshold",
        mode = "function",
        inherits = TRUE
      )
    ) {
      stop("Source scripts/r/sim-compare-freq_bs.R before running F-beta.")
    }
    do.call(
      .simCompareFbetaThreshold,
      c(
        list(
          xUns = xUns,
          xStim = xStim,
          pathFbeta = pathFbeta,
          fbetaEnv = fbetaEnv
        ),
        settings$params
      )
    )
  } else {
    if (
      !exists(
        ".simCompareTailgateThreshold",
        mode = "function",
        inherits = TRUE
      )
    ) {
      stop(
        "Source scripts/r/sim-compare-freq_bs.R before running Tailgate."
      )
    }
    do.call(
      .simCompareTailgateThreshold,
      list(
        x = xStim,
        adjust = settings$params$adjust,
        bandwidth = settings$params$bandwidth,
        numPeaks = settings$params$numPeaks,
        refPeak = settings$params$refPeak,
        method = settings$params$derivativeMethod,
        tol = settings$params$tol,
        side = settings$params$side,
        strict = settings$params$strict,
        autoTol = settings$params$autoTol
      )
    )
  }

  thresholdRaw <- suppressWarnings(as.numeric(thresholdObj$threshold))[[1]]
  fallbackUsed <- !is.finite(thresholdRaw)
  threshold <- if (fallbackUsed) {
    .acsCytofFallbackThreshold(xStim = xStim, xUns = xUns)
  } else {
    thresholdRaw
  }

  thresholdMetric <- thresholdObj$thresholdMetric
  if (is.null(thresholdMetric) || length(thresholdMetric) == 0L) {
    thresholdMetric <- NA_real_
  }
  thresholdOrigin <- thresholdObj$thresholdOrigin
  if (is.null(thresholdOrigin) || length(thresholdOrigin) == 0L) {
    thresholdOrigin <- "unknown"
  }

  tibble::tibble(
    thresholdRaw = thresholdRaw,
    threshold = threshold,
    thresholdMetric = suppressWarnings(
      as.numeric(thresholdMetric)
    )[[1]],
    thresholdOrigin = as.character(thresholdOrigin)[[1]],
    thresholdFallbackUsed = fallbackUsed
  )
}

.acsCytofRunComparator <- function(
  gs,
  pop,
  method = c("tailgate", "fbeta"),
  batchList = .acsCytofBatchList(length(gs)),
  pathFbeta = NULL
) {
  method <- match.arg(method)
  settings <- .acsCytofComparatorSettings(method)
  channelMap <- settings$channelMap
  channels <- names(channelMap)
  fbetaEnv <- if (identical(method, "fbeta")) {
    if (
      !exists(
        ".simCompareFbetaEnvironment",
        mode = "function",
        inherits = TRUE
      )
    ) {
      stop("Source scripts/r/sim-compare-freq_bs.R before running F-beta.")
    }
    .simCompareFbetaEnvironment(
      pathFbeta = pathFbeta,
      patchPy2Compat = settings$params$patchPy2Compat
    )
  } else {
    NULL
  }

  resultList <- purrr::imap(batchList, function(indBatch, batchIndex) {
    indUns <- indBatch[[1]]
    xUns <- .acsCytofExpressionMatrix(
      gs = gs,
      ind = indUns,
      channels = channels
    )

    purrr::map(indBatch[-1], function(indStim) {
      xStim <- .acsCytofExpressionMatrix(
        gs = gs,
        ind = indStim,
        channels = channels
      )

      thresholdTbl <- purrr::map_dfr(seq_along(channels), function(i) {
        .acsCytofThresholdOne(
          method = method,
          xUns = xUns[, i],
          xStim = xStim[, i],
          settings = settings,
          pathFbeta = pathFbeta,
          fbetaEnv = fbetaEnv
        ) |>
          dplyr::mutate(
            method = .env$method,
            pop = .env$pop,
            batch = as.character(batchIndex),
            ind = as.character(indStim),
            indUns = as.character(indUns),
            chnl = channels[[i]],
            cyt = unname(channelMap[[i]]),
            .before = 1L
          )
      })

      thresholdVec <- stats::setNames(
        thresholdTbl$threshold,
        thresholdTbl$chnl
      )
      statsTbl <- .acsCytofCombinationCounts(
        xStim = xStim,
        xUns = xUns,
        thresholds = thresholdVec[channels],
        channels = channels
      ) |>
        dplyr::mutate(
          gateName = .env$method,
          method = .env$method,
          pop = .env$pop,
          batch = as.character(batchIndex),
          ind = as.character(indStim),
          indUns = as.character(indUns),
          .before = 1L
        ) |>
        dplyr::mutate(
          propStim = .data$countStim / .data$nCellStim,
          propUns = .data$countUns / .data$nCellUns,
          propBs = .data$propStim - .data$propUns,
          freqStim = .data$propStim * 100,
          freqUns = .data$propUns * 100,
          freqBs = .data$freqStim - .data$freqUns
        )

      list(stats = statsTbl, thresholds = thresholdTbl)
    })
  }) |>
    unlist(recursive = FALSE)

  list(
    settings = settings,
    nSample = length(gs),
    createdAt = Sys.time(),
    stats = purrr::list_rbind(lapply(
      resultList,
      function(x) x[["stats"]]
    )),
    thresholds = purrr::list_rbind(lapply(
      resultList,
      function(x) x[["thresholds"]]
    ))
  )
}

.acsCytofCacheIsCurrent <- function(
  object,
  settings,
  nSample = NULL
) {
  if (
    !is.list(object) ||
      !all(c("settings", "nSample", "stats", "thresholds") %in% names(object))
  ) {
    return(FALSE)
  }
  if (
    !isTRUE(all.equal(
      object$settings,
      settings,
      check.attributes = FALSE
    ))
  ) {
    return(FALSE)
  }
  if (
    !is.null(nSample) &&
      !identical(as.integer(object$nSample), as.integer(nSample))
  ) {
    return(FALSE)
  }

  requiredStats <- c(
    "method",
    "pop",
    "ind",
    "cytCombn",
    "countStim",
    "nCellStim",
    "countUns",
    "nCellUns"
  )
  requiredThresholds <- c(
    "method",
    "pop",
    "ind",
    "chnl",
    "threshold",
    "thresholdOrigin",
    "thresholdFallbackUsed"
  )
  all(requiredStats %in% names(object$stats)) &&
    all(requiredThresholds %in% names(object$thresholds))
}

.acsCytofWriteComparatorCache <- function(object, path) {
  if (exists(".write_rds_atomic", mode = "function", inherits = TRUE)) {
    .write_rds_atomic(object, path)
    return(invisible(path))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pathTmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, pathTmp)
  if (!file.rename(pathTmp, path)) {
    file.copy(pathTmp, path, overwrite = TRUE)
    unlink(pathTmp)
  }
  invisible(path)
}

.acsCytofRemoveComparatorResults <- function(
  paths,
  methods = c("fbeta", "tailgate")
) {
  pathResultVec <- unname(unlist(paths[methods], use.names = FALSE))
  pathResultVec <- pathResultVec[file.exists(pathResultVec)]
  if (length(pathResultVec) == 0L) {
    return(invisible(character()))
  }

  unlink(pathResultVec)
  pathRemaining <- pathResultVec[file.exists(pathResultVec)]
  if (length(pathRemaining) > 0L) {
    stop(
      "Could not remove previous comparator result(s): ",
      paste(pathRemaining, collapse = ", ")
    )
  }

  invisible(pathResultVec)
}

.acsCytofReadComparatorCache <- function(
  path,
  method,
  nSample = NULL
) {
  if (!file.exists(path)) {
    stop(
      "Cached ",
      method,
      " result not found at: ",
      path,
      ". Run the ACS comparison methods first."
    )
  }

  object <- tryCatch(readRDS(path), error = function(e) NULL)
  settings <- .acsCytofComparatorSettings(method)
  if (!.acsCytofCacheIsCurrent(object, settings, nSample = nSample)) {
    stop(
      "Cached ",
      method,
      " result is missing, incomplete or was created ",
      "with different settings. Re-run that ACS comparison method."
    )
  }
  object
}

.acsCytofRunComparisonMethodsSafe <- function(...) {
  pop <- list(...)$pop

  tryCatch(
    {
      result <- .acsCytofRunComparisonMethods(...)
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

.acsCytofRunComparisonMethods <- function(
  pop,
  pathFcsBase,
  pathGsBase,
  pathScratchBase,
  runMethods,
  nSample = NULL,
  outputGroup = NULL,
  pathFbeta = NULL
) {
  paths <- .acsCytofPopulationPaths(
    pop = pop,
    pathFcsBase = pathFcsBase,
    pathGsBase = pathGsBase,
    pathScratchBase = pathScratchBase,
    outputGroup = outputGroup
  )
  methodVec <- c("tailgate", "fbeta") |> rev()

  if (!isTRUE(runMethods)) {
    resultList <- lapply(methodVec, function(method) {
      .acsCytofReadComparatorCache(
        path = paths[[method]],
        method = method,
        nSample = nSample
      )
    })
    names(resultList) <- methodVec
    return(invisible(list(pop = pop, paths = paths, results = resultList)))
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
  batchList <- .acsCytofBatchList(nSampleActual)

  .acsCytofRemoveComparatorResults(paths = paths, methods = methodVec)

  resultList <- lapply(methodVec, function(method) {
    message("Running ", method, " for ", pop, ".")
    result <- .acsCytofRunComparator(
      gs = gs,
      pop = pop,
      method = method,
      batchList = batchList,
      pathFbeta = pathFbeta
    )
    .acsCytofWriteComparatorCache(result, paths[[method]])
    result
  })
  names(resultList) <- methodVec

  invisible(list(
    pop = pop,
    nSample = nSampleActual,
    batchList = batchList,
    paths = paths,
    results = resultList
  ))
}
