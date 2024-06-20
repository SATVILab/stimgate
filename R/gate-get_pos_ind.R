#' @title Get logical indicator for whether each cell is positive
#' for any of the channels specified in chnl.
#'
#' @description Returns a logical vector specifying whether each
#' cell produces at least one of the cytokines specified, but only
#' using a single type of threshold (i.e. base, cytokine-positive or single).
#'
#' @inheritParams .get_pos_ind_simple
#' @param chnl character vector. Channels for which cell should be
#' positive for at least one. If \code{NULL}, then set to
#' \code{unique(gate_tbl$chnl)}, i.e. all channels for which we
#' have gates. Default is \code{chnl}.
#' @param gate_type "base", "cyt" or 'single'. Determines which gate
#' type is used. If \code{"base"}, then initial threshold is used.
#' If \code{"cyt"}, then cytokine-positive threshold is used. If
#' \code{'single'}, then single-positive threshold is used. No default is set.
.get_pos_ind_simple <- function(ex, gate_tbl, chnl = NULL, gate_type) {
  if (length(unique(gate_tbl$gate_name)) != 1) {
    stop(paste0(
      "gate_tbl$gate_name has ", length(unique(gate_tbl$gate_name)),
      " unique entries when it should have 1."
    ))
  }
  if (length(unique(gate_tbl$ind)) != 1) {
    stop(paste0(
      "gate_tbl$ind has ", length(unique(gate_tbl$ind)),
      " unique entries when it should have 1."
    ))
  }
  if (is.null(chnl)) chnl <- unique(gate_tbl$chnl)
  pos_vec <- rep(FALSE, nrow(ex))
  for (chnl_curr in chnl) {
    # print(chnl_curr)
    gate_tbl_chnl_ind <- which(gate_tbl$chnl == chnl_curr)
    gate_tbl_col_name <- switch(gate_type,
      "base" = "gate",
      "cyt" = "gate_cyt",
      "single" = "gate_single",
      stop(paste0(
        "gate_type ", gate_type,
        " not  recognised in .get_pos_ind_simple"
      ))
    )
    pos_vec <- pos_vec |
      ex[[chnl_curr]] > gate_tbl[[gate_tbl_col_name]][[gate_tbl_chnl_ind]]
  }
  pos_vec
}

#' @title Identify cells that are positive for at least two cytokines.
#'
#' @description Identify cells that express multiple cytokines. May
#' force each cytokine combination to contain at least one of a set of
#' cytokines (using \code{chnl} parameter) and specify which cytokines a cell
#' must express at least one of (using \code{chnl_alt} parameter). The
#' default is to consider a cell multi-functional if it expresses any cytokine
#' in conjunction with any other cytokine.
#'
#' @inheritParams .get_pos_ind_but_single_pos_for_one_cyt
#'
#' @param chnl character vector. Cytokines that a cell must be positive for
#' in any combination (exclusively, i.e. a cell need not be positive for both
#' simultaneously if both are placed here). If \code{NULL},
#' then set to \code{unique(gate_tbl$chnl)}. Default is \code{NULL}.
#' @param chnl_alt character vector. Cytokine(s) that a cell must be positive
#' for in conjunction with at least one cytokine in \code{chnl}. If \code{NULL},
#' then set to whatever \code{chnl} is, which is either specified
#' by the parameter or set to  \code{setdiff(unique(gate_tbl$chnl), <chnl_crr>)}
#' if \code{chnl} parameter is NULL, where \code{<chnl_curr>} is one of the
#' channels in \code{chnl} parameter. Default is \code{NULL}.
#' @param gate_type_cyt_pos "base" or "cyt". If \code{"base"}, then the initial thresholds are
#' used even for cells that are cytokine positive for at least one other cytokine.
#' If \code{"cyt"}, then the cytokine+ thresholds are used for cells that are cytokine
#' positive for at least one other cytokine.
#'
#' @return A logical vector with value \code{TRUE} if a cell is positive for
#' at least one cytokine in \code{chnl} and any cytokine in \code{chnl_alt}, and
#' \code{FALSE} otherwise.
#'
#' @example
#'
.get_pos_ind_mult <- function(ex, gate_tbl, chnl = NULL, chnl_alt = NULL,
                              gate_type_cyt_pos) {
  # must specify types of gates to use for cyt+ cells
  if (!gate_type_cyt_pos %in% c("base", "cyt")) {
    stop(paste0(
      "gate_type_cyt_pos value of ", ifelse(missing(gate_type_cyt_pos),
        "blank", gate_type_cyt_pos
      ),
      ' not either "cyt" or "base" in function .get_pos_ind_mult.'
    ))
  }

  # set defaults
  if (is.null(chnl)) chnl <- unique(gate_tbl$chnl)
  if (is.null(chnl_alt)) chnl_alt <- chnl

  # if not using cyt+ gates, simply count any cell as multifunctional
  # if it is above base threshold for at least two cyts
  if (gate_type_cyt_pos == "base") {
    pos_vec_cyt_pos_mult <- rep(0, nrow(ex))
    for (chnl_curr in chnl) {
      pos_vec_cyt_pos_curr <- .get_pos_ind_simple(
        ex = ex,
        gate_tbl = gate_tbl,
        chnl = chnl_curr,
        gate_type = "base"
      )
      pos_vec_cyt_pos_mult <- pos_vec_cyt_pos_mult + pos_vec_cyt_pos_curr
    }
    return(pos_vec_cyt_pos_mult >= 2)
  }

  # if using cyt+ gates, then for each cytokine in chnl, determine
  # which cells are positive for all other cytokines in chnl_alt and
  # then determine which of these are also positive for the current chnl.
  # Each of these cells is then positive for this chnl and at least
  # one other chnl. A cell is multifunctional if it is so positive
  # for at least one of the cytokines in chnl.
  if (gate_type_cyt_pos == "cyt") {
    pos_vec_cyt_pos_mult <- rep(FALSE, nrow(ex))
    for (chnl_curr in chnl) {
      # find which cells are positive for any cell but current cell,
      # based on "base" thresholds
      # based on current chnl being positive, using "base"
      # find which cells are positive for the current chnl,
      # based on cyt+ threshold
      pos_vec_curr <- .get_pos_ind_simple(
        ex = ex,
        gate_tbl = gate_tbl,
        chnl = chnl_curr,
        gate_type = "base"
      )
      pos_vec_alt <- .get_pos_ind_simple(
        ex = ex,
        gate_tbl = gate_tbl,
        chnl = setdiff(c(chnl_alt, chnl), chnl_curr),
        gate_type = "cyt"
      )

      # positive for chnl_curr using base threshold and any other using
      # cyt+ thresholds
      pos_vec_cyt_pos_mult_curr <- pos_vec_alt & pos_vec_curr
      pos_