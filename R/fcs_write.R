#' @title Write FCS files of marker-positive FCS files
#'
#' @description
#' Uses the gates to write FCS files of marker-positive FCS files.
#'
#' @param pathProject character. Path to project directory.
#' @param .data GatingSet. GatingSet object containing the flow cytometry data.
#' @param indBatchList list. List of indices grouped by batch.
#' @param pathDirSave character. Directory path to save the FCS files to.
#' @param pop character. Population that was gated on.
#' @param chnl character vector. Specific channels to gate on.
#' @param gateTbl data.frame. Pre-computed gate table, if available.
#' @param transFn function. Transformation function to apply.
#' @param transChnl character vector. Columns to transform.
#' @param combnExc list. Combinations of channels to exclude.
#' @param gateTypeCytPos character. Gate type to use for cytokine-positive cells.
#' @param gateTypeSinglePos character. Gate type to use for single-positive cells.
#' @param mult logical. Whether cells must be multi-positive.
#' @param gateUnsMethod character. Method to calculate unstimulated thresholds.
#'
#' @details
#' This function processes flow cytometry data to identify and export cytokine-positive
#' cells to FCS files. It requires that gates have been pre-computed using
#' \code{\link{stimgate_gate}} or that a complete gate table is provided.
#'
#' The function will create the output directory and write FCS files for samples
#' that contain cytokine-positive cells. If no positive cells are found in a sample,
#' no FCS file will be written for that sample.
#'
#' @examples
#' \dontrun{
#' # Complete workflow example
#' library(stimgate)
#'
#' # Load your GatingSet (gs) and define batch structure
#' # batch_list <- list(batch1 = c(1, 2, 3), batch2 = c(4, 5, 6))
#' # where the last element in each batch is the unstimulated sample
#'
#' # First, run gating to create gates
#' pathProject <- tempfile("stimgate_project")
#' # stimgate_gate(
#' #   .data = gs,
#' #   pathProject = pathProject,
#' #   popGate = "root",
#' #   batch_list = batch_list,
#' #   marker = c("IL2", "IFNg")  # your cytokine markers
#' # )
#'
#' # Then write FCS files of cytokine-positive cells
#' path_output <- tempfile("fcs_output")
#' # writeStimFCS(
#' #   pathProject = pathProject,
#' #   .data = gs,
#' #   indBatchList = batch_list,
#' #   pathDirSave = path_output,
#' #   chnl = c("IL2", "IFNg")
#' # )
#'
#' # Alternative: provide your own gate table
#' # gateTbl <- data.frame(
#' #   chnl = c("IL2", "IFNg"),
#' #   marker = c("IL2", "IFNg"),
#' #   batch = c(1, 1),
#' #   ind = c(1, 1),
#' #   gate = c(0.5, 0.3),
#' #   gateName = c("gate", "gate")
#' # )
#' # writeStimFCS(
#' #   pathProject = pathProject,
#' #   .data = gs,
#' #   indBatchList = batch_list,
#' #   pathDirSave = path_output,
#' #   chnl = c("IL2", "IFNg"),
#' #   gateTbl = gateTbl
#' # )
#' }
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
  gateTypeSinglePos = "single", # gate type to use for single-pos cells # nolint
  mult = FALSE, # whether cells must be multi-positive
  gateUnsMethod = "min"
) {
  # how to calculate unstim thresholds # nolint
  # clear and create directory to save to
  if (dir.exists(pathDirSave)) {
    unlink(pathDirSave, force = TRUE, recursive = TRUE)
  }
  dir.create(pathDirSave, recursive = TRUE)

  # get gates
  gateTbl <- .fcsWriteGetGateTbl(
    gateTbl = gateTbl,
    chnl = chnl,
    pop = pop,
    .data = .data,
    indBatchList = indBatchList,
    gateUnsMethod = gateUnsMethod,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos,
    pathProject = pathProject
  )

  nFn <- length(.data)

  purrr::walk(seq_along(.data), function(ind) {
    txt <- paste0("Writing ", ind, " of ", nFn, " files")
    message(txt)
    .fcsWriteImpl(
      .data = .data,
      ind = ind,
      gateTbl = gateTbl,
      pathDirSave = pathDirSave,
      chnl = chnl,
      mult = mult,
      gateTypeCytPos = gateTypeCytPos,
      gateTypeSinglePos = gateTypeSinglePos,
      combnExc = combnExc,
      transFn = transFn,
      transChnl = transChnl
    )
  })
  invisible(pathDirSave)
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
  gateTypeCytPos,
  gateTypeSinglePos,
  pathProject
) {
  # Get gate table if not provided
  if (is.null(gateTbl)) {
    popUnspecified <- is.null(pop)
    pop <- pop %||% .gateGetPop(pathProject)
    if (is.null(pop) || length(pop) == 0) {
      stop(
        "No population provided and no populations found in project directory."
      )
    }
    if (length(pop) > 1) {
      stop(
        "Multiple populations found in project directory. Please specify 'pop' parameter."
      )
    }
    if (popUnspecified) {
      message(paste0("Using population '", pop, "' from project directory."))
    }
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
    gateTbl <- .gateGetGateTblAll(NULL, pop, chnl, pathProject)
  }

  # Check if gateTbl already contains all required information
  # (i.e., it has both stimulated and unstimulated gates)
  hasAllSamples <- length(unique(gateTbl$ind)) >=
    length(unlist(indBatchList))

  # Only process unstimulated gates if needed
  if (!hasAllSamples) {
    gateTbl <- gateTbl |>
      .fcsWriteGetGateTblAddUns(
        gateUnsMethod = gateUnsMethod,
        gateTypeCytPos = gateTypeCytPos,
        gateTypeSinglePos = gateTypeSinglePos,
        indBatchList = indBatchList
      )
  }

  # Apply remaining processing
  gateTbl <- gateTbl |>
    .fcsWriteGetGateTblFilterChnl(chnl) |>
    .fcsWriteGetGateTblAddMarker(chnl, .data) |>
    # duplicates are not a possible issue,
    # as the gates must be the same for all duplicates
    # as we separate them if they are not
    dplyr::distinct()

  gateTbl
}

#' @keywords internal
.gateGetGateTblAll <- function(gateTbl, pop, chnl, pathProject) {
  if (!is.null(gateTbl)) {
    return(gateTbl)
  }
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
  gateTypeCytPos,
  gateTypeSinglePos,
  indBatchList
) {
  gateTblUns <- .fcsWriteGetGateTblAddUnsGetUns(
    gateTbl = gateTbl,
    gateUnsMethod = gateUnsMethod,
    indBatchList = indBatchList
  )

  if ("gate_cyt" %in% colnames(gateTbl)) {
    gateTblUns <- gateTblUns |>
      dplyr::mutate(gate_cyt = pmin(gate, gate_cyt)) # nolint
  } else if ("gate_cyt" %in% colnames(gateTblUns)) {
    gateTblUns <- gateTblUns |> dplyr::select(-gate_cyt) # nolint
  }
  if ("gate_single" %in% colnames(gateTbl)) {
    gateTblUns <- gateTblUns |>
      dplyr::mutate(gate_single = pmax(gate, gate_single)) # nolint
  } else if ("gate_single" %in% colnames(gateTblUns)) {
    gateTblUns <- gateTblUns |> dplyr::select(-gate_single) # nolint
  }

  gateTbl |>
    dplyr::bind_rows(gateTblUns)
}

#' @keywords internal
.fcsWriteGetGateTblAddUnsGetUns <- function(
  gateTbl,
  gateUnsMethod,
  indBatchList
) {
  calcUnsGate <- .fcsWriteGetGateTblAddUnsGetUnsCalc(
    gateUnsMethod = gateUnsMethod
  )

  .fcsWriteGetGateTblAddUnsGetUnsImpl(
    gateTbl = gateTbl,
    calc = calcUnsGate,
    indBatchList = indBatchList
  )
}

#' @keywords internal
.fcsWriteGetGateTblAddUnsGetUnsCalc <- function(gateUnsMethod) {
  switch(
    gate_uns_method,
    "min" = min,
    "max" = max,
    "mean" = mean,
    "tmean" = function(x) mean(x, trim = 0.2),
    "med" = median,
    stop("gateUnsMethod not recognised")
  )
}

#' @keywords internal
.fcsWriteGetGateTblAddUnsGetUnsImpl <- function(
  gateTbl,
  calc,
  indBatchList
) {
  gateTblDistinct <- gateTbl |>
    dplyr::distinct(chnl, marker, batch, ind, .keep_all = TRUE)
  if (nrow(gateTblDistinct) != nrow(gateTbl)) {
    # check that gates are the same for all duplicates
    cnVec <- c("chnl", "marker", "batch", "ind")
    cnVecConcate <- NULL
    for (i in seq_along(cnVec)) {
      cnVecConcate <- paste(cnVecConcate, gateTbl[[cnVec[i]]], sep = "_")
    }
    gateTbl$concat <- cnVecConcate
    gateTbl <- gateTbl |>
      dplyr::group_by(concat) |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::ungroup()
    gateVec <- paste0(gateTbl$gate, gateTbl$gate_cyt, gateTbl$gate_single)
    gateTbl$gateConcat <- gateVec
    isError <- nrow(
      gateTbl |>
        dplyr::group_by(chnl, marker, batch, ind) |>
        dplyr::filter(length(unique(gateConcat)) > 1) |>
        dplyr::ungroup()
    ) >
      0
    if (isError) {
      stop("Gates are not the same for all duplicates in gateTbl.")
    }
  }
  gateTblDistinct |>
    dplyr::group_by(chnl, marker, batch) |> # nolint
    dplyr::summarise(
      indStim = paste0(ind |> sort(), collapse = "_"),
      gate = calc(gate), # nolint
      gate_cyt = calc(gate_cyt), # nolint
      gate_single = calc(gate_single), # nolint
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
    (x[-length(x)]) |>
      sort() |>
      paste0(collapse = "_")
  }) |>
    unlist()
  indUnsVec <- lapply(indBatchList, function(x) x[length(x)]) |>
    unlist()
  indVec <- lapply(seq_len(nrow(gateTbl)), function(x) {
    indMatch <- which(indBatchVec == gateTbl$indStim[[x]])
    stopifnot(length(indMatch) == 1L)
    indUnsVec[indMatch]
  }) |>
    unlist() |>
    as.character()
  gateTbl |>
    dplyr::mutate(ind = indVec) |>
    dplyr::select(chnl, marker, batch, ind, everything()) |> # nolint
    dplyr::arrange(chnl, marker, batch, ind)
}


#' @keywords internal
.fcsWriteGetGateTblFilterChnl <- function(gateTbl, chnl) {
  if (is.null(chnl)) {
    return(gateTbl)
  }
  gateTbl |> dplyr::filter(chnl %in% .env$chnl)
}

#' @keywords internal
.fcsWriteGetGateTblAddMarker <- function(gateTbl, chnl, .data) {
  chnlLabVec <- .getLabs(.data = .data[[1]], chnlCut = chnl) # nolint

  # Base columns that should always be present
  baseCols <- c("chnl", "marker", "batch", "ind", "gate")

  # Optional columns that may or may not be present
  optionalCols <- c("gate_cyt", "gate_single", "gate_name")

  # Only select columns that exist
  colsToSelect <- c(
    baseCols,
    optionalCols[optionalCols %in% colnames(gateTbl)]
  )

  gateTbl |>
    dplyr::mutate(marker = chnlLabVec[.data$chnl]) |> # nolint
    dplyr::select(dplyr::any_of(colsToSelect)) |>
    dplyr::arrange(chnl, marker, batch, ind)
}

# ===============
# Implementation
# ================

#' @keywords internal
.fcsWriteImpl <- function(
  .data,
  ind,
  gateTbl,
  pathDirSave,
  chnl,
  mult,
  gateTypeCytPos,
  gateTypeSinglePos,
  combnExc,
  transFn,
  transChnl
) {
  fr <- .fcsWriteImplLoad(.data, ind)
  ex <- flowCore::exprs(fr) |> tibble::as_tibble()

  if (is.na(ex[1, chnl[1]]) && nrow(ex) == 1) {
    return(invisible(FALSE))
  }

  gateTblInd <- gateTbl |>
    dplyr::filter(.data$ind == .env$ind) # nolint

  ex <- .dataGetExCytPosInc(
    ex = ex,
    gateTblInd = gateTblInd,
    mult = mult,
    chnl = chnl,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos
  )

  if (nrow(ex) == 0) {
    message("No stimulation-positive cells. No FCS file written.")
    return(invisible(FALSE))
  }

  ex <- .dataGetExCytPosExc(
    ex = ex,
    combnExc = combnExc,
    gateTblInd = gateTblInd,
    chnl = chnl,
    gateTypeCytPos = gateTypeCytPos,
    gateTypeSinglePos = gateTypeSinglePos
  )

  if (nrow(ex) == 0) {
    message(
      "No cells after excluding particular combinations. No FCS file written."
    )
    return(invisible(FALSE))
  }

  ex <- .dataGetExTrans(ex, transFn, transChnl)

  .fcsWriteImplWrite(ex, fr, pathDirSave)
  invisible(TRUE)
}

#' @keywords internal
.fcsWriteImplLoad <- function(.data, ind) {
  fr <- flowWorkspace::gh_pop_get_data(.data[[ind]])
  if (inherits(fr, "cytoframe")) {
    fr <- flowWorkspace::cytoframe_to_flowFrame(fr)
  }
  fr
}


#' @keywords internal
.fcsWriteImplWrite <- function(ex, fr, pathDirSave) {
  flowCore::exprs(fr) <- as.matrix(ex)
  fn <- flowCore::keyword(fr)[["GUID"]] |> basename()
  fnOut <- file.path(pathDirSave, fn)
  if (file.exists(fnOut)) {
    invisible(file.remove(fnOut))
  }
  flowCore::write.FCS(x = fr, filename = fnOut)
  txt <- paste0("Wrote ", fn)
  message(txt)
  invisible(TRUE)
}
