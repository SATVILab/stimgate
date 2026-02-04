# Simulate cytokine expression data for multiple channels sequentially.
# This is the main orchestrator that processes a list of channel specifications,
# simulating expression data for each channel while carrying forward the updated
# flowSet to the next channel.
#
# Arguments:
#   args_list - Named list where each element contains simulation parameters for one channel
#   fs - Initial flowSet to use as template
#
# Returns:
#   A named list with one element per channel, each containing:
#     - fs: updated flowSet with simulated expression
#     - ind_list: indices of positive cells per sample
#     - resp_tbl: metadata about proportions and cell counts
#     - batch_list: batch groupings
sample_fs <- function(args_list,
                      fs) {
  # Initialize output list with channel names
  chnl_obj_list <- lapply(seq_along(args_list), function(x) NULL) |>
    stats::setNames(names(args_list))
  ind_list_pos_already <- NULL
  # Process each channel sequentially
  for (i in seq_along(args_list)) {
    args <- args_list[[i]]
    # Simulate this channel
    chnl_obj_list[[i]] <- sample_chnl(
      fs = fs,
      batch_list = args$batch_list,
      chnl = args$chnl,
      ind_list_pos_already = ind_list_pos_already,
      prop_mean_pos = args$prop_mean_pos,
      prop_sd_pos = args$prop_sd_pos,
      prop_mean_neg = args$prop_mean_neg,
      prop_sd_neg = args$prop_sd_neg,
      expr_mean_neg = args$expr_mean_neg,
      expr_mean_pos = args$expr_mean_pos,
      expr_sd_pos = args$expr_sd_pos,
      expr_sd_neg = args$expr_sd_neg,
      prop_mean_pos_ctrb_from_pos = args$prop_mean_pos_ctrb_from_pos,
      prop_sd_pos_ctrb_from_pos = args$prop_sd_pos_ctrb_from_pos
    )
    # Carry forward updated flowSet to next channel
    fs <- chnl_obj_list[[i]]$fs
    ind_list_pos_already <- chnl_obj_list[[i]]$ind_list
  }
  chnl_obj_list
}

# Simulate cytokine expression data for a single channel across all samples.
# Generates bimodal expression data (positive and negative populations) for each
# sample in the flowSet, with different parameters for stimulated vs unstimulated.
#
# Arguments:
#   fs - flowSet to modify with simulated data
#   batch_list - List defining which sample indices belong to each batch
#   chnl - Name of the channel/marker to simulate
#   prop_mean_pos - Mean proportion of positive cells in stimulated samples
#   prop_sd_pos - SD of positive proportion in stimulated samples
#   prop_mean_neg - Mean proportion of positive cells in unstimulated samples
#   prop_sd_neg - SD of positive proportion in unstimulated (default: same as prop_sd_pos)
#   expr_mean_neg - Mean expression value for negative cells
#   expr_mean_pos - Mean expression value for positive cells
#   expr_sd_pos - SD of expression for positive cells
#   expr_sd_neg - SD of expression for negative cells (default: same as expr_sd_pos)
#
# Returns:
#   A list containing:
#     - fs: flowSet with updated expression values for this channel
#     - ind_list: list of indices identifying positive cells in each sample
#     - resp_tbl: tibble with metadata (batch, n_cell, prop_pos per sample)
#     - batch_list: the input batch_list (passed through)
sample_chnl <- function(fs,
                        batch_list,
                        chnl,
                        ind_list_pos_already = NULL,
                        prop_mean_pos,
                        prop_sd_pos,
                        prop_mean_neg,
                        prop_sd_neg = NULL,
                        expr_mean_neg,
                        expr_mean_pos,
                        expr_sd_pos,
                        expr_sd_neg = NULL,
                        prop_mean_pos_ctrb_from_pos = NULL,
                        prop_sd_pos_ctrb_from_pos = NULL) {
  # Initialize storage for positive cell indices
  ind_list <- lapply(seq_along(fs), function(x) NULL)
  ind_list_pos_already <- if (!is.null(ind_list_pos_already)) {
    ind_list_pos_already
  } else {
    lapply(seq_along(fs), function(x) NULL)
  }
  # Initialize response tibble to track simulation outcomes
  resp_tbl <- tibble::tibble(
    chnl = chnl,
    sample_ind = seq_along(fs),
    batch = rep(NA_integer_, length(fs)),
    n_cell = rep(NA_integer_, length(fs)),
    prop_pos = rep(NA_real_, length(fs)),
  )
  # Process each batch
  for (i in seq_along(batch_list)) {
    # Last sample in each batch is unstimulated
    ind_uns <- batch_list[[i]][length(batch_list[[i]])]
    # Other samples are stimulated
    ind_stim_vec <- setdiff(batch_list[[i]], ind_uns)
    # Process stimulated samples
    for (ind_stim in ind_stim_vec) {
      fr_stim <- fs[[ind_stim]]
      # indices of cells already positive from previous channels
      ind_pos_already <- ind_list_pos_already[[ind_stim]]
      ex_mat_stim <- flowCore::exprs(fr_stim)
      n_cell_stim <- nrow(ex_mat_stim)
      # Generate simulated response for stimulated sample
      response_stim <- sample_response(
        n = n_cell_stim,
        prop_mean = prop_mean_pos,
        prop_sd = prop_sd_pos,
        expr_mean_neg = expr_mean_neg,
        expr_mean_pos = expr_mean_pos,
        expr_sd_pos = expr_sd_pos,
        expr_sd_neg = expr_sd_neg,
        ind_pos_already = ind_pos_already,
        prop_mean_ctrb_from_pos = prop_mean_pos_ctrb_from_pos,
        prop_sd_ctrb_from_pos = prop_sd_pos_ctrb_from_pos
      )
      # Update expression values in flowFrame
      ex_mat_stim[, chnl] <- response_stim$expr_vec
      flowCore::exprs(fr_stim) <- ex_mat_stim
      fs[[ind_stim]] <- fr_stim
      # Record metadata
      resp_tbl$batch[ind_stim] <- i
      resp_tbl$n_cell[ind_stim] <- nrow(ex_mat_stim)
      resp_tbl$prop_pos[ind_stim] <- response_stim$n_pos / n_cell_stim
      ind_list[[ind_stim]] <- response_stim$ind_pos
    }
    # Process unstimulated sample
    fr_uns <- fs[[ind_uns]]
    ind_pos_already <- ind_list_pos_already[[ind_uns]]
    ex_mat_uns <- flowCore::exprs(fr_uns)
    n_cell_uns <- nrow(ex_mat_uns)
    # Use same SD as stimulated if not specified
    prop_sd_neg <- if (is.null(prop_sd_neg)) prop_sd_pos else prop_sd_neg
    # Generate simulated response for unstimulated sample
    response_uns <- sample_response(
      n = n_cell_uns,
      prop_mean = prop_mean_neg,
      prop_sd = prop_sd_neg,
      expr_mean_neg = expr_mean_neg,
      expr_mean_pos = expr_mean_pos,
      expr_sd_pos = expr_sd_pos,
      expr_sd_neg = expr_sd_neg,
      ind_pos_already = ind_pos_already,
      prop_mean_ctrb_from_pos = prop_mean_pos_ctrb_from_pos,
      prop_sd_ctrb_from_pos = prop_sd_pos_ctrb_from_pos
    )
    # Update expression values
    ex_mat_uns[, chnl] <- response_uns$expr_vec
    flowCore::exprs(fr_uns) <- ex_mat_uns
    fs[[ind_uns]] <- fr_uns
    # Record metadata
    resp_tbl$batch[ind_uns] <- i
    resp_tbl$n_cell[ind_uns] <- nrow(ex_mat_uns)
    resp_tbl$prop_pos[ind_uns] <- response_uns$n_pos / n_cell_uns
    ind_list[[ind_uns]] <- response_uns$ind_pos
  }
  # Return all components
  list(
    fs = fs,
    ind_list = ind_list,
    resp_tbl = resp_tbl,
    batch_list = batch_list
  )
}


# Generate a simulated response (positive cells and expression values) for one sample.
# This is a helper function that orchestrates the sampling of positive cell count,
# indices, and expression values.
#
# Arguments:
#   n - Total number of cells in the sample
#   prop_mean - Mean proportion of positive cells
#   prop_sd - SD of positive proportion
#   expr_mean_neg - Mean expression for negative cells
#   expr_mean_pos - Mean expression for positive cells
#   expr_sd_pos - SD of expression for positive cells
#   expr_sd_neg - SD of expression for negative cells (default: same as expr_sd_pos)
#
# Returns:
#   A list containing:
#     - n_pos: number of positive cells
#     - ind_pos: indices of positive cells
#     - expr_vec: vector of expression values for all cells
sample_response <- function(n,
                            prop_mean,
                            prop_sd,
                            expr_mean_neg,
                            expr_mean_pos,
                            expr_sd_pos,
                            expr_sd_neg = NULL,
                            ind_pos_already = NULL,
                            prop_mean_ctrb_from_pos = NULL,
                            prop_sd_ctrb_from_pos = NULL) {
  ind_pos_already <- ind_pos_already %||% integer(0L)
  n_pos_already <- length(ind_pos_already)
  prop_mean_ctrb_from_pos <- prop_mean_ctrb_from_pos %||% 0L
  prop_sd_ctrb_from_pos <- prop_sd_ctrb_from_pos %||% 0L
  # get desired proportion that are positive overall
  prop_pos <- sample_prop_mean(prop_mean, prop_sd)
  n_pos <- round(n * prop_pos)
  # get how many cells this implies must be positive
  # amongst those already positive and those
  # not yet positive
  # must first sample what the contribution is,
  # so the split isn't the same every time:
  prop_ctrb_from_pos <- sample_prop_mean(
    prop_mean_ctrb_from_pos, prop_sd_ctrb_from_pos
  )
  # now we have those counts, and that's actually
  # all we need since we're returning indices]
  # as we need indices of already positive cells
  # and that has all the information needed
  # to calculate proportions later on
  n_pos_among_already <- round(n_pos * prop_ctrb_from_pos)
  n_pos_among_not_yet <- n_pos - n_pos_among_already
  # Randomly select which already-positive cells
  # are positive for this marker as well
  ind_pos_among_already_ind <- sample_ind_pos(
    n = n_pos_already,
    n_cell_pos = n_pos_among_already
  )
  # this is actually a subset of ind_pos_already,
  # as at most all those can be positive
  # and we're choosing from there
  ind_pos_among_already <- sample_ind_pos_map(
    ind_pos_among_already_ind, ind_pos_already
  )
  # Randomly select from not-yet-positive cells
  # which will be positive for this marker
  ind_pos_not_yet <- setdiff(seq_len(n), ind_pos_already)
  ind_pos_new_ind <- sample_ind_pos(
    n = length(ind_pos_not_yet),
    n_cell_pos = n_pos_among_not_yet
  )
  ind_pos_new <- sample_ind_pos_map(ind_pos_new_ind, ind_pos_not_yet)

  # Combine new positive indices with already positive indices
  if (length(intersect(ind_pos_new, ind_pos_among_already)) > 0L) {
    stop("Overlap between new positive indices and already positive indices")
  }
  ind_pos <- sort(c(ind_pos_new, ind_pos_among_already))
  if (length(ind_pos) > n) {
    stop("Number of positive indices cannot exceed total number of cells")
  }
  # re-calculate this as the numbers could be off due to rounding
  n_pos <- length(ind_pos)
  # Generate expression values for positive and negative cells
  # could make multi-positive cells have different expression
  # distribution here, but for now we keep it simple
  expr_vec <- sample_expr(
    n = n,
    n_pos = n_pos,
    ind_pos = ind_pos,
    mean_neg = expr_mean_neg,
    mean_pos = expr_mean_pos,
    sd_pos = expr_sd_pos,
    sd_neg = expr_sd_neg
  )
  list(
    n_pos = n_pos,
    ind_pos = ind_pos,
    expr_vec = expr_vec
  )
}

# Sample the number of positive cells using a beta distribution.
# The beta distribution is parameterized to achieve the desired mean proportion
# and variance, while constraining the SD to be valid for a proportion.
#
# Arguments:
#   n - Total number of cells
#   prop_mean - Desired mean proportion of positive cells
#   prop_sd - Desired SD of positive proportion
#   eps - Small value to ensure prop_sd < theoretical maximum (default: 1e-8)
#
# Returns:
#   Integer count of positive cells (rounded from beta-sampled proportion * n)
sample_n_pos <- function(n, prop_mean, prop_sd, eps = 1e-8) {
  prop_mean <- sample_prop_mean(prop_mean, prop_sd)
  if (prop_mean == 0L) {
    return(0L)
  }
  # Sample proportion and convert to cell count
  round(n * prop_mean)
}

sample_prop_mean <- function(prop_mean, prop_sd) {
  if (prop_mean == 0L) {
    return(0L)
  }
  if (prop_mean >= 1L) {
    return(1L)
  }
  # Compute the maximum SD and force prop_sd < max_sd
  max_sd <- sqrt(prop_mean * (1 - prop_mean))
  prop_sd <- min(prop_sd, max_sd * (1 - 1e-8))
  if (prop_sd <= 0L) {
    return(prop_mean)
  }
  # Convert mean and SD to beta distribution parameters
  nu    <- prop_mean * (1 - prop_mean) / prop_sd^2 - 1
  alpha <- prop_mean * nu
  beta  <- (1 - prop_mean) * nu
  # Sample proportion
  rbeta(n = 1, shape1 = alpha, shape2 = beta)
}

# Randomly select indices of positive cells from the total cell population.
#
# Arguments:
#   n - Total number of cells
#   n_cell_pos - Number of positive cells to select
#
# Returns:
#   Integer vector of indices for positive cells
sample_ind_pos <- function(n, n_cell_pos) {
  if (n_cell_pos == 0L) {
    return(integer(0))
  }
  if (n_cell_pos > n) {
    stop("n_cell_pos cannot be greater than n")
  }
  sample.int(n, n_cell_pos)
}

# Generate expression values for positive and negative cell populations.
# Draws values from normal distributions with specified means and SDs.
#
# Arguments:
#   n - Total number of cells
#   n_pos - Number of positive cells
#   ind_pos - Indices of positive cells
#   mean_neg - Mean expression for negative cells
#   mean_pos - Mean expression for positive cells
#   sd_pos - SD of expression for positive cells
#   sd_neg - SD of expression for negative cells (default: same as sd_pos)
#
# Returns:
#   Numeric vector of expression values for all n cells
sample_expr <- function(n,
                        n_pos,
                        ind_pos,
                        mean_neg,
                        mean_pos,
                        sd_pos,
                        sd_neg = NULL) {
  n_neg <- n - n_pos
  expr_vec <- rep(NA_real_, n)
  # Generate expression values for positive cells if any exist
  if (n_pos > 0) {
    expr_vec[ind_pos] <- rnorm(n_pos, mean = mean_pos, sd = sd_pos)
  }
  # Generate expression values for negative cells if any exist
  if (n_neg > 0) {
    if (is.null(sd_neg)) {
      sd_neg <- sd_pos
    }
    expr_vec[setdiff(seq_len(n), ind_pos)] <-
      rnorm(n_neg, mean = mean_neg, sd = sd_neg)
  }
  expr_vec
}

# Calculate actual proportions for all marker combinations from simulated data.
# This function computes the proportion of cells positive for each possible combination
# of marker positivity/negativity patterns (e.g., BC1+BC2+, BC1+BC2-, etc.).
#
# Arguments:
#   chnl_list - List of channel objects containing ind_list (positive cell indices)
#
# Returns:
#   A tibble with columns:
#     - n_pos_chnl: number of positive markers in this combination
#     - chnl: marker combination string (e.g., "BC1(La139)Dd~+~BC2(Pr141)Dd~-~")
#     - batch: batch number
#     - sample_ind: sample index
#     - n_cell: total cells in sample
#     - ind: list of indices for cells matching this combination
#     - prop_pos: proportion of cells matching this combination
#     - prop_bs: background-subtracted proportion (stim - unstim)
calc_resp_combn <- function(chnl_list) {
  # Get all possible positive/negative combinations
  combn_pos_list <- calc_resp_combn_get_combn_pos(chnl_list)
  # Calculate proportions for each combination
  resp_tbl <- lapply(combn_pos_list, function(chnl_pos) {
    calc_resp_combn_ind(chnl_pos, chnl_list)
  }) |>
    dplyr::bind_rows()
  # Add background subtraction and metadata
  resp_tbl |>
    dplyr::group_by(chnl, batch) |>
    dplyr::mutate(
      prop_bs = prop_pos - prop_pos[length(prop_pos)]
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_pos_chnl = stringr::str_count(chnl, stringr::fixed("~+~")),
    ) |>
    dplyr::select(
      n_pos_chnl, chnl, batch, sample_ind, n_cell, ind, prop_pos, prop_bs,
      everything()
    )
 }

# Generate all possible marker combinations (positive/negative patterns).
# For n markers, generates 2^n combinations (including all-negative).
#
# Arguments:
#   chnl_list - List of channel objects
#
# Returns:
#   A list where each element is a character vector of marker names that are
#   positive in that combination
calc_resp_combn_get_combn_pos <- function(chnl_list) {
  chnl_vec <- names(chnl_list)
  # Get combinations for each number of positive markers (0 to n)
  combn_mat_list_pos_init <- lapply(seq_along(chnl_vec), function(n_pos) {
    utils::combn(chnl_vec, m = n_pos, simplify = FALSE)
  })
  # Flatten nested list
  combn_mat_list_pos <- list()
  for (i in seq_along(combn_mat_list_pos_init)) {
    combn_mat_list_pos_init_curr <- combn_mat_list_pos_init[[i]]
    for (j in seq_along(combn_mat_list_pos_init_curr)) {
      combn_mat_list_pos <- combn_mat_list_pos |>
        append(combn_mat_list_pos_init_curr[j])
    }
  }
  combn_mat_list_pos
}

# Calculate proportion for a specific marker combination pattern.
# Finds cells that are positive for specified markers AND negative for all others.
#
# Arguments:
#   chnl_pos - Character vector of markers that should be positive
#   chnl_list - List of channel objects containing ind_list
#
# Returns:
#   A tibble with proportion of cells matching this exact combination pattern
calc_resp_combn_ind <- function(chnl_pos, chnl_list) {
  # Identify negative markers
  chnl_neg <- setdiff(names(chnl_list), chnl_pos)
  chnl_pos_ind <- which(names(chnl_list) %in% chnl_pos)
  # Build combination name string
  sign_ind <- rep("~-~", length(chnl_list))
  sign_ind[chnl_pos_ind] <- "~+~"
  combn_nm <- paste0(names(chnl_list), sign_ind, collapse = "")
  n_sample <- length(chnl_list[[chnl_pos[1]]]$ind_list)
  # Use resp_tbl from first channel as scaffold
  resp_tbl <- chnl_list[[chnl_pos[1]]]$resp_tbl
  resp_tbl$chnl <- combn_nm
  ind_list_pos <- lapply(seq_len(n_sample), function(x) NULL)
  # For each sample, find cells matching this combination
  for (i in seq_len(n_sample)) {
    # Get positive cell indices for each channel in this sample
    ind_list_sample <- lapply(seq_along(chnl_list), function(chnl) {
      chnl_list[[chnl]]$ind_list[[i]]
    }) |>
      stats::setNames(names(chnl_list))
    # Find cells positive for ALL required markers (intersection)
    pos_vec_ind <- Reduce(intersect, lapply(chnl_pos, function(chnl) {
      ind_list_sample[[chnl]]
    }))
    # Remove cells that are positive for any negative markers
    if (length(chnl_neg) > 0L && length(pos_vec_ind) > 0L) {
      # Union of negative markers (remove if positive for ANY)
      neg_vec_ind <- Reduce(union, lapply(chnl_neg, function(chnl) {
        ind_list_sample[[chnl]]
      }))
      # Subtract negative indices from positive
      pos_vec_ind <- setdiff(pos_vec_ind, neg_vec_ind)
    }
    ind_list_pos[[i]] <- pos_vec_ind
    resp_tbl$prop_pos[i] <- length(pos_vec_ind) / resp_tbl$n_cell[i]
  }
  resp_tbl$ind <- ind_list_pos
  resp_tbl
}

  # Map back to original indices
  sample_ind_pos_map <- function(ind_subset, ind_full) {
    if (length(ind_subset) == 0L) {
      return(integer(0))
    }
    if (length(ind_subset) > length(ind_full)) {
      stop("ind_subset cannot be longer than ind_full")
    }
    if (length(ind_full) == 0L) {
      stop("ind_full cannot be empty")
    }
    ind_full[ind_subset]
  }