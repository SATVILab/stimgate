#' @importFrom rlang .data .env
#' @importFrom stats setNames
#' @importFrom  utils head

# Global variable bindings to avoid R CMD check notes
# These are primarily used in dplyr and ggplot2 contexts
globalVariables(c(
  # Variables used in dplyr operations
  "marker", "batch", "ind", "gate", "gate_cyt", "gate_single", "gate_name",
  "chnl", "gate_use", "gate_type", "gate_combn", "pop_gate", "gate_tbl",
  "chnl_cut", "tol", "ind_in_batch_gate", "tol_clust_single", "ind_batch_gate",
  "cp", "grp", "cp_join_lse_orig_mean_tg", "cp_orig_quant_min", "cp_join",
  "cp_join_lse", "cp_join_lse_orig", "cp_join_lse_orig_mean", "cp_join_tg_orig",
  "cp_join_tg_orig_mean", "prop_bs_orig", 
  "prop_bs_cp_diff", "prop_bs_cp_diff_sd",
  "prop_bs_cp", "prop_l1se", "pred", "der", "cp_orig", "max_expr", "gate_05",
  "prop_bs_cp_diff_sd_max", "grp_level", "ind_vec", "x1", "x", "y", "x_ind",
  "count_stim", "n_cell_stim", "count_uns", "n_cell_uns", "prop_stim_pos",
  "prop_uns_pos", "prop_stim_sd", "prop_uns_sd", "x512", "x_vec", "freq_bs",
  "freq_stim", "pop_gate_curr", "cp_join_tg", "cp_orig_quant_min", "lse_orig",
  "cp_tg_ctrl", "cp_join_lse_orig_mean", "cp_join_tg_orig", "chnl_pos",
  "dir_save", "is.null_gate_tbl", "path_project", "exc_min", "prop_bs_diff",
  "prop_stim", "prop_uns", "prop_bs", "prob_smooth", "n_row", "y_stim", 
  "y_uns", "stim", "x_stim", "prob", "x_uns", "prop_lab", "type", "line_id",
  "dens", "no", "yes", "prob_stim", "prob_stim_norm", "prop_pos",
  # Variables used in plots and ggplot2 context
  ".debug", "prop_l1se"
))

.debug_msg <- function(.debug, msg, val = NULL) {
  if (!.debug) {
    return(invisible(FALSE))
  }
  if (!is.null(val)) {
    msg <- paste0(msg, ": ", val)
  }
  message(msg)
  invisible(TRUE)
}
