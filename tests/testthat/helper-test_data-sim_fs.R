# Create a channel list with simulated cytokine expression data for BC1 and BC2 markers.
# This function sets up simulation parameters for two channels (BC1 and BC2) with 
# different positive proportions and expression levels, then calls sample_chnls to 
# generate the simulated data.
#
# Arguments:
#   fs - A flowSet object to use as the template for simulation
#
# Returns:
#   A list with one element per channel, each containing:
#     - fs: updated flowSet with simulated expression values
#     - ind_list: list of indices for positive cells in each sample
#     - resp_tbl: tibble with response metadata (n_cell, prop_pos, etc.)
#     - batch_list: list defining which samples belong to each batch
sim_fs_default <- function(fs, n_ind = 8) {
  # Create batch list: 8 batches, each with 2 samples (indices)
  batch_list <- lapply(seq_len(n_ind), function(i) seq((i - 1) * 2 + 1, i * 2))
  # Define simulation parameters for BC1 channel
  # BC1 has higher positive proportion (5%) with strong separation
  args_list_bc1  <- list(
    batch_list = batch_list,
    chnl = "BC1(La139)Dd",
    prop_mean_pos = 0.05,     # Mean proportion of positive cells
    prop_sd_pos = 0.01,       # SD of positive proportion
    prop_mean_neg = 0.005,    # Mean proportion positive in unstimulated
    prop_sd_neg = 0.0075,     # SD of unstimulated positive proportion
    expr_mean_neg = 0,        # Mean expression for negative cells
    expr_mean_pos = 2,        # Mean expression for positive cells(good separation)
    expr_sd_pos = 0.1         # SD of positive expression
  )
  # Define simulation parameters for BC2 channel
  # BC2 has lower positive proportion (1%) with strong separation
  args_list_bc2 <- list(
    batch_list = batch_list,
    chnl = "BC2(Pr141)Dd",
    prop_mean_pos = 0.01,     # Lower mean proportion of positive cells
    prop_sd_pos = 0.01, 
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 2,        # Same strong separation as BC1
    expr_sd_pos = 0.1
  )
  # Combine into named list
  args_list <- list(
    "BC1(La139)Dd" = args_list_bc1,
    "BC2(Pr141)Dd" = args_list_bc2
  )
  # Generate simulated data for all channels
  sample_fs(args_list = args_list, fs = fs)
}

# Create "easy" scenario channel list with well-separated positive and negative populations.
# This is an alias for get_chnl_list, used in simulation comparisons to represent
# an easier gating scenario.
#
# Arguments:
#   fs - A flowSet object to use as the template for simulation
#
# Returns:
#   Same as get_chnl_list: a list with simulated channel data
sim_fs_easy <- function(fs, n_ind = 8) {
  batch_list <- lapply(seq_len(n_ind), function(i) seq((i - 1) * 2 + 1, i * 2))
  args_list_bc1  <- list(
    batch_list = batch_list,
    chnl = "BC1(La139)Dd",
    prop_mean_pos = 0.05,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 2,
    expr_sd_pos = 0.1
  )
  args_list_bc2 <- list(
    batch_list = batch_list,
    chnl = "BC2(Pr141)Dd",
    prop_mean_pos = 0.01,
    prop_sd_pos = 0.01, 
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 2,
    expr_sd_pos = 0.1
  )
  args_list <- list(
    "BC1(La139)Dd" = args_list_bc1, "BC2(Pr141)Dd" = args_list_bc2
  )

  sample_fs(args_list = args_list, fs = fs)
}

# Create "poor separation" scenario channel list with less distinct populations.
# Unlike the easy scenario, this uses lower mean expression for positive cells (0.5)
# and higher variance, making it harder to distinguish positive from negative cells.
#
# Arguments:
#   fs - A flowSet object to use as the template for simulation
#
# Returns:
#   Same as get_chnl_list: a list with simulated channel data
sim_fs_poor_separation <- function(fs, n_ind = 8) {
  batch_list <- lapply(seq_len(n_ind), function(i) seq((i - 1) * 2 + 1, i * 2))
  args_list_bc1  <- list(
    batch_list = batch_list,
    chnl = "BC1(La139)Dd",
    prop_mean_pos = 0.05,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 0.5,      # Lower expression = poor separation
    expr_sd_pos = 0.1
  )
  args_list_bc2 <- list(
    batch_list = batch_list,
    chnl = "BC2(Pr141)Dd",
    prop_mean_pos = 0.01,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 0.5,      # Lower expression = poor separation
    expr_sd_pos = 0.2         # Higher variance makes it even harder
  )
  args_list <- list(
    "BC1(La139)Dd" = args_list_bc1, "BC2(Pr141)Dd" = args_list_bc2
  )

  sample_fs(args_list = args_list, fs = fs)
}

# Create scenario where cells are more likely to be
# positive if they're already positive in previous c
# Arguments:
#   fs - A flowSet object to use as the template for simulation
#
# Returns:
#   Same as get_chnl_list: a list with simulated channel data
sim_fs_cyt_pos <- function(fs, n_ind = 8) {
  batch_list <- lapply(seq_len(n_ind), function(i) seq((i - 1) * 2 + 1, i * 2))
  args_list_bc1  <- list(
    batch_list = batch_list,
    chnl = "BC1(La139)Dd",
    prop_mean_pos = 0.05,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 0.5,      # Lower expression = poor separation
    expr_sd_pos = 0.1
  )
  args_list_bc2 <- list(
    batch_list = batch_list,
    chnl = "BC2(Pr141)Dd",
    prop_mean_pos = 0.01,
    prop_sd_pos = 0.01,
    prop_mean_neg = 0.005,
    prop_sd_neg = 0.0075,
    expr_mean_neg = 0,
    expr_mean_pos = 0.5,      # Lower expression = poor separation
    expr_sd_pos = 0.2,         # Higher variance makes it even harder
    # 50% of proportion that are positive for
    # this subset, are already positive in previous channel(s)
    prop_mean_pos_ctrb_from_pos = 0.5,
    prop_sd_pos_ctrb_from_pos = 0.1
  )
  args_list <- list(
    "BC1(La139)Dd" = args_list_bc1,
    "BC2(Pr141)Dd" = args_list_bc2
  )

  sample_fs(args_list = args_list, fs = fs)
}