# Get logical indicator for cytokine-positive cells
# Returns a logical vector indicating whether each cell is positive
# for any of the specified channels using a single threshold type
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

# Identify cells that are positive for at least two cytokines.
#
# Identify cells that express multiple cytokines. May
# Get logical indicator for multi-functional cells
# Identifies cells positive for multiple cytokines with customizable requirements
#
# Example usage would require data setup

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
      pos_vec_cyt_pos_mult <- pos_vec_cyt_pos_mult | pos_vec_cyt_pos_mult_curr

      for (chnl_alt_curr in setdiff(c(chnl, chnl_alt), chnl_curr)) {
        pos_vec_curr <- .get_pos_ind_simple(
          ex = ex,
          gate_tbl = gate_tbl,
          chnl = chnl_curr,
          gate_type = "cyt"
        )

        pos_vec_alt_base <- .get_pos_ind_simple(
          ex = ex,
          gate_tbl = gate_tbl,
          chnl = chnl_alt_curr,
          gate_type = "base"
        )

        # pos_vec_alt_cyt <- .get_pos_ind_simple(ex = ex,
        #                                        gate_tbl = gate_tbl,
        #                                        chnl = setdiff(c(chnl, chnl_alt), c(chnl_curr, chnl_alt_curr),
        #                                        gate_type = "cyt")
        # positive for chnl_curr using base threshold and any other using cyt+ thresholds
        pos_vec_cyt_pos_mult_curr <- pos_vec_alt_base & pos_vec_curr # & pos_vec_alt_cyt
        pos_vec_cyt_pos_mult <- pos_vec_cyt_pos_mult | pos_vec_cyt_pos_mult_curr
      }

      # find which cells are positive for any cell but current cell,
      # based on "base" thresholds
      # pos_vec_alt <- .get_pos_ind_simple(ex = ex,
      #                                   gate_tbl = gate_tbl,
      #                                   chnl = setdiff(c(chnl_alt, chnl), chnl_curr),
      #                                   gate_type = gate_type_cyt_pos)

      # find which cells are positive for the current chnl,
      # based on cyt+ threshold
      # pos_vec_curr <- .get_pos_ind_simple(ex = ex,
      #                                    gate_tbl = gate_tbl,
      #                                    chnl = chnl_curr,
      #                                    gate_type = "base")
      # find
      # pos_vec_cyt_pos_mult_curr <- pos_vec_alt & pos_vec_curr
      # pos_vec_cyt_pos_mult <- pos_vec_cyt_pos_mult | pos_vec_cyt_pos_mult_curr
    }
  }

  pos_vec_cyt_pos_mult
}

# Identify cells positive for all cytokines except one
# Finds cells positive for every cytokine except one specified channel
.get_pos_ind_but_single_pos_for_one_cyt <- function(ex,
                                                    gate_tbl,
                                                    chnl_single_exc,
                                                    chnl = NULL,
                                                    gate_type_cyt_pos,
                                                    gate_type_single_pos) {
  .get_pos_ind_but_single_pos_for_one_cyt_check(
    gate_type_single_pos = gate_type_single_pos
  )
  if (is.null(chnl)) chnl <- unique(gate_tbl$chnl)
  chnl <- c(chnl_single_exc, chnl)
  chnl <- unique(chnl)

  # cells positive for any cytokine except
  # current cytokine using base or single threshold
  pos_vec_single_ind_any_cyt_but_curr <- .get_pos_ind_simple(
    ex = ex, gate_tbl = gate_tbl,
    chnl = setdiff(chnl, chnl_single_exc),
    gate_type = gate_type_single_pos
  )

  # above specifies all cells that are polyfunctional,
  # if no adjusted thresholds are
  # used
  if (gate_type_cyt_pos == "base" && gate_type_single_pos == "base") {
    return(pos_vec_single_ind_any_cyt_but_curr)
  }

  # cells positive for any two cytokines
  pos_vec_multi_cyt <- .get_pos_ind_mult(
    ex = ex, gate_tbl = gate_tbl,
    chnl = chnl,
    chnl_alt = chnl, gate_type_cyt_pos = gate_type_cyt_pos
  )

  pos_vec_single_ind_any_cyt_but_curr | pos_vec_multi_cyt
}

.get_pos_ind_but_single_pos_for_one_cyt_check <- function(gate_type_single_pos) {
  # must specify types of gates to use for single+ cells
  if (missing(gate_type_single_pos)) {
    stop(paste0(
      "gate_type_single_pos missing"
    ))
  }
  if (!gate_type_single_pos %in% c("base", "single")) {
    stop(paste0(
      "gate_type_single_pos value of ", gate_type_single_pos,
      ' not either "cyt" or "base" in function .get_pos_ind_mult.'
    ))
  }
  invisible(TRUE)
}

# Identify cells that express at least one cytokine
# Returns a logical vector indicating cytokine-positive cells using flexible thresholds
.get_pos_ind <- function(ex, gate_tbl, chnl, chnl_alt = NULL, gate_type_cyt_pos,
                         gate_type_single_pos) {
  # must specify types of gates to use for single+ cells
  if (!gate_type_single_pos %in% c("base", "single")) {
    stop(paste0(
      "gate_type_cyt_pos value of ", ifelse(missing(gate_type_single_pos),
        "blank", gate_type_single_pos
      ),
      " not either 'single' or 'base' in function .get_pos_ind"
    ))
  }

  if (is.null(chnl)) chnl <- unique(gate_tbl$chnl)
  if (is.null(chnl_alt)) chnl_alt <- setdiff(unique(gate_tbl$chnl), chnl)
  chnl_alt <- setdiff(chnl_alt, chnl)

  # if only base thresholds are used, then it's sufficient to
  # look only for cells that are positive for the required
  # channels using base thresholds
  if (gate_type_cyt_pos == "base" && gate_type_single_pos == "base") {
    pos_ind_vec <- .get_pos_ind_simple(
      ex = ex, gate_tbl = gate_tbl,
      chnl = chnl, gate_type = "base"
    )
    return(pos_ind_vec)
  }

  # if an adjusted threshold is used (either single or cyt_pos or both) and every channel in chnl_alt
  # is required, then this chunk calculates positivity


  # cells positive for any cytokine that is required
  pos_vec_single <- .get_pos_ind_simple(
    ex = ex, gate_tbl = gate_tbl,
    chnl = chnl, gate_type = gate_type_single_pos
  )

  # cells positive for at least two cytokines, for some cytokines that are required
  pos_vec_multi <- .get_pos_ind_mult(
    ex = ex, gate_tbl = gate_tbl, chnl = chnl,
    chnl_alt = chnl_alt, gate_type_cyt_pos = gate_type_cyt_pos
  )

  # cells positive for either one of the require cytokines (possibly on its own)
  # or any other cytokine together with the required cytokine
  pos_ind_vec <- pos_vec_single | pos_vec_multi

  pos_ind_vec
}


.get_pos_ind_cyt_combn <- function(ex, gate_tbl, chnl_pos, chnl_neg, chnl_alt, gate_type_cyt_pos, gate_type_single_pos) {
  # must specify types of gates to use for single+ cells
  if (!gate_type_single_pos %in% c("base", "single")) {
    stop(paste0(
      "gate_type_cyt_pos value of ", ifelse(missing(gate_type_single_pos),
        "blank", gate_type_single_pos
      ),
      " not either 'single' or 'base' in function .get_pos_ind"
    ))
  }
  chnl <- unique(c(chnl_pos, chnl_neg, chnl_alt))
  chnl_pos_ind_vec_pos <- rep(TRUE, nrow(ex))

  # get all cells that are positive for any marker
  for (chnl_curr in chnl_pos) {
    chnl_pos_ind_vec_pos_curr <- .get_pos_ind(
      ex = ex, gate_tbl = gate_tbl, chnl = chnl_curr,
      chnl_alt = setdiff(chnl, chnl_curr),
      gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )
    chnl_pos_ind_vec_pos <- chnl_pos_ind_vec_pos & chnl_pos_ind_vec_pos_curr
  }

  if (length(chnl_neg) > 0) {
    # get all cells that are negative for any marker that they should be negative
    chnl_neg_ind_vec_pos <- .get_pos_ind(
      ex = ex, gate_tbl = gate_tbl, chnl = chnl_neg,
      chnl_alt = chnl,
      gate_type_cyt_pos = gate_type_cyt_pos,
      gate_type_single_pos = gate_type_single_pos
    )

    # all cells positive only for positive markers are those that
    # are positive for those markers and negative for all other markers
    chnl_pos_ind_vec_pos <- chnl_pos_ind_vec_pos & !chnl_neg_ind_vec_pos
  }

  chnl_pos_ind_vec_pos
}
