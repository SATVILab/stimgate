#' @title Create an R Markdown HTML doc based on output from get_cp and plot_cp
render_gate <- function(params_knit, pop_sub, pathProject) {
  # --------------------------
  # Preparation
  # --------------------------

  # Directories
  # ------------------

  # get directory for RMD render doc
  dir_rmd_in <- here::here("data-raw/gating/gate-render_gate.Rmd")

  # get base directory
  dir_base <- stimgate_dir_base_create(
    dir_base_init = pathProject,
    params = params_knit,
    empty_dir = FALSE
  )

  # Output file name
  # ------------------

  # get node
  popGate <- params_knit$popGate
  popGate_last_slash_mat <- stringr::str_locate_all(popGate, "/")[[1]]
  popGate_last_slash_loc <- popGate_last_slash_mat[, "end"][nrow(popGate_last_slash_mat)]
  popGate_final <- stringr::str_sub(popGate, popGate_last_slash_loc + 1)

  # file save name
  dir_html_out <- file.path(dir_base, paste0(popGate_final, " - ", params_knit$chnl_lab[params_knit$cut], ".html"))

  # --------------------------
  # Create RMD
  # --------------------------

  rmarkdown::render(
    input = system.file("extdata/gate-render_gate.Rmd", package = ""),
    output_file = dir_html_out,
    params = params_knit |> append(list(pop_sub = pop_sub)),
    quiet = TRUE
  )
}
