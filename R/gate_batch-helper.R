#' @keywords internal
.gateBatchTbl <- function(gateList, batch) {
  rows <- list()
  for (i in seq_along(gateList)) {
    .debug("gate list index", i) # nolint
    gateType <- names(gateList)[i]
    gateListElem <- gateList[[i]]
    cpList <- if ("cp" %in% names(gateListElem)) {
      gateListElem[["cp"]]
    } else {
      gateListElem
    }

    for (j in seq_along(cpList)) {
      .debug("gate list sub-index", j) # nolint
      gateCombn <- names(cpList)[[j]]
      gateVec <- cpList[[j]]
      meta <- .getCpUnsLocMetaFromCp(gateVec)
      gateUse <- if (length(gateType) > 0 && any(grepl("tgCtrl_", gateType))) {
        "ctrl"
      } else {
        "gate"
      }

      rows[[length(rows) + 1L]] <- tibble::tibble(
        gateName = paste0(gateType, "_", gateCombn),
        gateType = gateType,
        gateCombn = gateCombn,
        batch = batch,
        ind = as.character(names(gateVec)),
        gate = gateVec,
        gateUse = gateUse
      ) |>
        dplyr::bind_cols(dplyr::select(meta, -"ind"))
    }
  }

  dplyr::bind_rows(rows)
}
