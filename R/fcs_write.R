#' @title Export stimulation-positive cells as FCS files
#' @description Select positive cells using saved or supplied gates and write
#'   one FCS file per sample with retained cells, plus `manifest.csv`.
#' @param pathProject character Project directory from [gateStim()].
#' @param .data GatingSet Cytometry data used for gating.
#' @param indBatchList list Sample indices grouped by batch, control first.
#' @param pathDirSave character Output directory; existing contents are deleted.
#' @param pop character or NULL Population to export. NULL uses the single saved
#'   population, or "root" when `gateTbl` is supplied. Default: NULL.
#' @param chnl character vector or NULL Channels used to select positive cells;
#'   NULL uses all channels in the gate table. Default: NULL.
#' @param gateTbl data.frame or NULL Gates with `chnl`, `batch`, `ind`, `gate`
#'   and, for refined gates, `gateCyt`. NULL reads saved gates. Default: NULL.
#' @param gateTypeCytPos character Positivity rule: "base" uses main gates;
#'   "cyt" also uses refined gates for cells positive for another marker.
#'   Default: "cyt".
#' @param mult logical Require positivity for at least two markers. Default: FALSE.
#' @param combnExc list or NULL Channel combinations to exclude: each vector
#'   specifies positive channels, with other selected channels negative.
#'   Default: NULL.
#' @param gateUnsMethod character Summary of stimulated gates used for missing
#'   control gates: "min", "max", "mean", "tmean" (20% trimmed mean), or "med".
#'   Default: "min".
#' @param transFn function or NULL Transformation of retained expression before
#'   export. Default: NULL.
#' @param transChnl character vector or NULL Columns to transform; NULL transforms
#'   all expression columns. Default: NULL.
#' @return Invisibly, a tibble with one row per sample and columns `ind`, `batch`,
#'   `fileName`, `nCellPos`, `written`, `reason`. The `pathDirSave` attribute holds
#'   the output path. Samples with no retained cells have no FCS file.
#' @examples
#' exampleData <- getExampleData()
#' gs <- flowWorkspace::load_gs(exampleData$pathGs)
#' pathProject <- gateStim(
#'   tempfile("stimgate_"), gs, exampleData$batchList,
#'   marker = exampleData$marker
#' )
#' manifest <- writeStimFCS(
#'   pathProject, gs, indBatchList = exampleData$batchList,
#'   pathDirSave = tempfile("positive_fcs_")
#' )
#' @export
writeStimFCS <- function(
  pathProject, # project directory
  .data, # gatingset
  pop = NULL, # population that was gated on
  indBatchList, # indices by batch
  pathDirSave, # directory to save to
  chnl = NULL, # specific channels to gate on
  gateTbl = NULL, # whether gateTbl is pre-available
  transFn = NULL, # transformation to apply
  transChnl = NULL, # columns to transform
  combnExc = NULL, # combinations of chnl to exclude
  gateTypeCytPos = "cyt", # gate type to use for cyt-pos cells # nolint
  mult = FALSE, # whether cells must be multi-positive
  gateUnsMethod = "min"
) {
  # how to calculate unstim thresholds # nolint
  popUnspecified <- is.null(pop)
  pop <- pop %||% if (is.null(gateTbl)) .gateGetPop(pathProject) else "root"
  if (is.null(pop) || length(pop) == 0) {
    stop(
      "No population provided and no populations found in project directory."
    )
  }
  if (length(pop) > 1) {
    stop(paste0(
      "Multiple populations found in project directory. ",
      "Please specify 'pop' parameter."
    ))
  }
  if (popUnspecified && is.null(gateTbl)) {
    message(paste0("Using population '", pop, "' from project directory."))
  }

  # get gates
  gateTbl <- .fcsWriteGetGateTbl(
    gateTbl = gateTbl,
    chnl = chnl,
    pop = pop,
    .data = .data,
    indBatchList = indBatchList,
    gateUnsMethod = gateUnsMethod,
    pathProject = pathProject
  )

  chnl <- chnl %||% unique(gateTbl$chnl)

  # clear and create directory to save to
  if (dir.exists(pathDirSave)) {
    unlink(pathDirSave, force = TRUE, recursive = TRUE)
  }
  dir.create(pathDirSave, recursive = TRUE)

  nFn <- length(.data)

  manifestRows <- purrr::map(seq_along(.data), function(ind) {
    txt <- paste0("Writing ", ind, " of ", nFn, " files")
    message(txt)
    .fcsWriteImpl(
      .data = .data,
      ind = ind,
      pop = pop,
      gateTbl = gateTbl,
      pathDirSave = pathDirSave,
      chnl = chnl,
      mult = mult,
      gateTypeCytPos = gateTypeCytPos,
      combnExc = combnExc,
      transFn = transFn,
      transChnl = transChnl,
      indBatchList = indBatchList
    )
  })
  manifest <- if (nFn == 0L) {
    tibble::tibble(
      ind = character(0),
      batch = character(0),
      fileName = character(0),
      nCellPos = integer(0),
      written = logical(0),
      reason = character(0)
    )
  } else {
    dplyr::bind_rows(manifestRows)
  }
  utils::write.csv(
    manifest,
    file = file.path(pathDirSave, "manifest.csv"),
    row.names = FALSE
  )
  attr(manifest, "pathDirSave") <- pathDirSave
  invisible(manifest)
}


# ================
# Get Gates
# ================

#' @keywords internal
.fcsWriteGetGateTbl <- function(
  gateTbl,
  chnl,
  pop,
  .data,
  indBatchList,
  gateUnsMethod,
  pathProject
) {
  # Get gate table if not provided
  if (is.null(gateTbl)) {
    chnlUnspecified <- is.null(chnl)
    chnl <- chnl %||% .gateGetChnl(pathProject, pop)
    if (is.null(chnl) || length(chnl) == 0) {
      stop("No channels provided and no channels found in project directory.")
    }
    if (chnlUnspecified) {
      message(paste0(
        "Using channels '",
        paste(chnl, collapse = ", "),
        "' from project directory."
      ))
    }
    gateTbl <- .gateGetGateTblAll(pop, chnl, pathProject)
  }

  # Check if gateTbl already contains all required information
  # (i.e., it has both stimulated and unstimulated gates)
  hasAllSamples <- all(unlist(indBatchList) %in% gateTbl$ind)

  if (!"marker" %in% names(gateTbl)) {
    gateTbl <- .fcsWriteGetGateTblAddMarker(
      gateTbl, chnl %||% unique(gateTbl$chnl), .data
    )
  }

  # Only process unstimulated gates if needed
  if (!hasAllSamples) {
    gateTbl <- gateTbl |>
      .fcsWriteGetGateTblAddUns(
        gateUnsMethod = gateUnsMethod,
        indBatchList = indBatchList
      )
  }

  chnl <- chnl %||% unique(gateTbl$chnl)

  # Apply remaining processing
  gateTbl <- gateTbl |>
    dplyr::filter(chnl %in% .env$chnl) |>
    .fcsWriteGetGateTblAddMarker(chnl, .data) |>
    # duplicates are not a possible issue,
    # as the gates must be the same for all duplicates
    # as we separate them if they are not
    dplyr::distinct()

  gateTbl
}

#' @keywords internal
.gateGetGateTblAll <- function(pop, chnl, pathProject) {
  purrr::map_df(chnl, function(chnlCurr) {
    pathCurr <- .gatesGetPathAll(pathProject, pop, chnlCurr, FALSE)
    if (!file.exists(pathCurr)) {
      stop(paste0("Gate table not found for channel: ", chnlCurr))
    }
    readRDS(pathCurr)
  })
}

#' @keywords internal
.fcsWriteGetGateTblAddUns <- function(
  gateTbl,
  gateUnsMethod,
  indBatchList
) {
  calcUnsGate <- switch(gateUnsMethod,
    "min" = min,
    "max" = max,
    "mean" = mean,
    "tmean" = function(x) mean(x, trim = 0.2),
    "med" = stats::median,
    stop("gateUnsMethod not recognised")
  )
  gateTblUns <- .fcsWriteGetGateTblAddUnsGetUnsImpl(
    gateTbl = gateTbl,
    calc = calcUnsGate,
    indBatchList = indBatchList
  )

  if ("gateCyt" %in% colnames(gateTbl)) {
    gateTblUns <- gateTblUns |>
      dplyr::mutate(gateCyt = pmin(gate, gateCyt)) # nolint
  }

  gateTbl |>
    dplyr::bind_rows(gateTblUns)
}

#' @keywords internal
.fcsWriteGetGateTblAddUnsGetUnsImpl <- function(
  gateTbl,
  calc,
  indBatchList
) {
  gateTblDistinct <- gateTbl |>
    dplyr::distinct(chnl, marker, batch, ind, .keep_all = TRUE)
  thresholdCols <- c("chnl", "marker", "batch", "ind", "gate", "gateCyt")
  gateTblThresholds <- gateTbl |>
    dplyr::distinct(dplyr::across(dplyr::any_of(thresholdCols)))
  if (nrow(gateTblThresholds) != nrow(gateTblDistinct)) {
    stop("Gates are not the same for all duplicates in gateTbl.")
  }
  gateTblDistinct |>
    dplyr::group_by(chnl, marker, batch) |> # nolint
    dplyr::summarise(
      indStim = paste0(ind |> sort(), collapse = "_"),
      dplyr::across(c("gate", dplyr::any_of("gateCyt")), calc),
      .groups = "drop"
    ) |>
    .fcsWriteGetGateTblAddUnsGetUnsInd(indBatchList)
}

#' @keywords internal
.fcsWriteGetGateTblAddUnsGetUnsInd <- function(
  gateTbl,
  indBatchList
) {
  indBatchVec <- lapply(indBatchList, function(x) {
    (x[-1]) |>
      sort() |>
      paste0(collapse = "_")
  }) |>
    unlist()
  indUnsVec <- lapply(indBatchList, function(x) x[[1]]) |>
    unlist()
  indVec <- vapply(gateTbl$indStim, function(indStim) {
    indMatch <- which(indBatchVec == indStim)
    stopifnot(length(indMatch) == 1L)
    as.character(indUnsVec[indMatch])
  }, character(1), USE.NAMES = FALSE)
  gateTbl |>
    dplyr::mutate(ind = indVec) |>
    dplyr::select(chnl, marker, batch, ind, dplyr::everything()) |> # nolint
    dplyr::arrange(chnl, marker, batch, ind)
}


#' @keywords internal
.fcsWriteGetGateTblAddMarker <- function(gateTbl, chnl, .data) {
  chnlLabVec <- .getLabs(.data = .data[[1]], chnlCut = chnl) # nolint

  gateTbl |>
    dplyr::mutate(marker = chnlLabVec[.data$chnl]) |> # nolint
    dplyr::select(dplyr::any_of(c(
      "chnl", "marker", "batch", "ind", "gate", "gateCyt", "gateName"
    ))) |>
    dplyr::arrange(chnl, marker, batch, ind)
}

# ===============
# Implementation
# ================

#' @keywords internal
.fcsWriteImpl <- function(
  .data,
  ind,
  pop,
  gateTbl,
  pathDirSave,
  chnl,
  mult,
  gateTypeCytPos,
  combnExc,
  transFn,
  transChnl,
  indBatchList = NULL
) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]], y = pop)
  if (inherits(fr, "cytoframe")) {
    fr <- flowWorkspace::cytoframe_to_flowFrame(fr)
  }
  guid <- flowCore::keyword(fr)[["GUID"]]
  fileName <- if (!is.null(guid) && length(guid) > 0L && !is.na(guid[[1]])) {
    basename(as.character(guid[[1]]))
  } else {
    basename(as.character(flowCore::identifier(fr)))
  }
  batch <- .fcsWriteGetBatch(ind, indBatchList)

  ex <- flowCore::exprs(fr) |> tibble::as_tibble()

  if (nrow(ex) == 0L || (nrow(ex) == 1L && is.na(ex[[chnl[1]]][1]))) {
    return(tibble::tibble(
      ind = as.character(ind),
      batch = as.character(batch),
      fileName = fileName,
      nCellPos = 0L,
      written = FALSE,
      reason = "empty_sample"
    ))
  }

  gateTblInd <- gateTbl |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .dataGetExCytPosInc(
    ex = ex,
    gateTblInd = gateTblInd,
    mult = mult,
    chnl = chnl,
    gateTypeCytPos = gateTypeCytPos
  )

  if (nrow(ex) == 0L) {
    message("No stimulation-positive cells. No FCS file written.")
    return(tibble::tibble(
      ind = as.character(ind),
      batch = as.character(batch),
      fileName = fileName,
      nCellPos = 0L,
      written = FALSE,
      reason = "no_positive_cells"
    ))
  }

  ex <- .dataGetExCytPosExc(
    ex = ex,
    combnExc = combnExc,
    gateTblInd = gateTblInd,
    chnlGate = chnl,
    gateTypeCytPos = gateTypeCytPos
  )

  if (nrow(ex) == 0L) {
    message(
      "No cells after excluding particular combinations. No FCS file written."
    )
    return(tibble::tibble(
      ind = as.character(ind),
      batch = as.character(batch),
      fileName = fileName,
      nCellPos = 0L,
      written = FALSE,
      reason = "none_after_exclusion"
    ))
  }

  nCellPos <- as.integer(nrow(ex))
  ex <- .dataGetExTrans(ex, transFn, transChnl)

  .fcsWriteImplWrite(ex, fr, pathDirSave)
  tibble::tibble(
    ind = as.character(ind),
    batch = as.character(batch),
    fileName = fileName,
    nCellPos = nCellPos,
    written = TRUE,
    reason = "written"
  )
}

#' @keywords internal
.fcsWriteGetBatch <- function(ind, indBatchList) {
  bNames <- names(indBatchList)
  if (!is.list(indBatchList) || is.null(bNames)) {
    return(NA_character_)
  }
  for (i in seq_along(indBatchList)) {
    items <- indBatchList[[i]]
    if (as.character(ind) %in% as.character(items)) {
      b <- bNames[i]
      if (!is.na(b) && nzchar(b)) {
        return(as.character(b))
      }
      return(NA_character_)
    }
  }
  NA_character_
}

#' @keywords internal
.fcsWriteImplWrite <- function(ex, fr, pathDirSave) {
  flowCore::exprs(fr) <- as.matrix(ex)
  fn <- flowCore::keyword(fr)[["GUID"]] |> basename()
  fnOut <- file.path(pathDirSave, fn)
  flowCore::write.FCS(x = fr, filename = fnOut)
  txt <- paste0("Wrote ", fn)
  message(txt)
  invisible(TRUE)
}
