# QMD 7 presentation and stable extension of its original response grid.
.simCompareQmd7SeedGrid <- function(grid, simulation_seed) {
  # Existing response levels retain their row order and sorted biological IDs.
  # New levels are appended after the complete legacy grid, before profile filters.
  grid |>
    dplyr::mutate(.new_response = .data$prob_response == 0.01) |>
    dplyr::arrange(.data$.new_response) |>
    dplyr::group_by(
      .data$.new_response, .data$transformation, .data$prob_response,
      .data$n_cell, .data$mean_pos_setting, .data$mean_pos,
      .data$sample_perturbation_sd, .data$condition_perturbation_sd,
      .data$cluster_perturbation_sd, .data$background_relative_to_response,
      .data$n_cell_uns_relative_to_stim
    ) |>
    dplyr::mutate(
      base_scenario_id = dplyr::cur_group_id(),
      sim_seed = as.integer(simulation_seed + .data$base_scenario_id - 1L)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(sim_id = dplyr::row_number()) |>
    dplyr::select(-".new_response")
}

# Keep the plot's existing statistical subtitle and add its loop settings.
.simCompareQmd7Subtitle <- function(plot, data, method_col = "method") {
  methods <- unique(as.character(data[[method_col]]))
  method_set <- if ("tailgate" %in% methods) "All methods" else "Without Tailgate"
  settings <- c(method_set,
    paste0("Mean position: ", paste(unique(data$mean_pos_setting), collapse = ", ")))
  for (col in c("prob_response", "condition_perturbation_sd", "n_cell")) {
    values <- unique(data[[col]])
    if (length(values) == 1L) {
      label <- switch(col, prob_response = "Response probability",
        condition_perturbation_sd = "Condition perturbation SD", n_cell = "Stimulated cells")
      value <- if (col == "prob_response") .analysis_label_percent(values) else values
      settings <- c(settings, paste0(label, ": ", value))
    }
  }
  plot + ggplot2::labs(
    subtitle = paste(c(plot$labels$subtitle, settings), collapse = " | ")
  )
}

.simCompareQmd7FallbackTable <- function(summary) {
  counts <- c("n_threshold_fallback", "n_run_error", "n_no_cutpoint")
  summary |>
    dplyr::group_by(.data$method, .data$transformation,
      .data$condition_perturbation_sd, .data$mean_pos_setting) |>
    dplyr::summarise(
      n_samples = sum(.data$n),
      dplyr::across(dplyr::all_of(counts), sum), .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(counts),
      ~ paste0(.x, " / ", .data$n_samples, " (",
        .analysis_label_percent(.x / .data$n_samples), ")")))
}
