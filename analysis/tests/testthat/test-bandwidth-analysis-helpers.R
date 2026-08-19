root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

script_bw_io <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-io.R")
script_bw_plot <- file.path(root_dir, "scripts", "r", "sim-bandwidth-analysis-plot.R")

.load_bw_analysis_env <- function() {
  env <- new.env(parent = getNamespace("stimgate"))

  for (f in c(script_bw_io, script_bw_plot)) {
    if (!file.exists(f)) {
      stop("Expected analysis helper not found: ", f)
    }
    source(f, local = env)
  }

  env
}

test_that("output discovery finds per-simulation files in a temporary output directory", {
  env <- .load_bw_analysis_env()
  tmp_dir <- tempfile("bw-output")
  output_dir <- file.path(tmp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  file_1 <- file.path(output_dir, "bw_list_raw-chunk_001-of_002-sim_id_000001.rds")
  file_2 <- file.path(output_dir, "bw_list_raw-chunk_001-of_002-sim_id_000002.rds")
  saveRDS(1, file_1)
  saveRDS(2, file_2)

  found <- env$.find_bw_list_output_files(output_dir = output_dir)
  expect_setequal(basename(found), c(
    "bw_list_raw-chunk_001-of_002-sim_id_000001.rds",
    "bw_list_raw-chunk_001-of_002-sim_id_000002.rds"
  ))
})

test_that("bandwidth plot helpers create stable labels and scales", {
  env <- .load_bw_analysis_env()

  bw_vec <- c(0.1, 0.25, NA_real_)
  expect_equal(env$format_bw_lab(bw_vec), c("0.1", "0.25", "NA"))
  expect_equal(names(env$make_bw_colour_values(bw_vec)), env$format_bw_lab(bw_vec))

  lty <- env$make_bw_linetype_scale(bw_vec)
  expect_equal(length(lty), length(bw_vec))
  expect_true(all(names(lty) %in% c("0.1", "0.25", "NA")))

  out <- env$add_bw_labs(data.frame(bw = c(0.1, 0.25)))
  expect_true("bw_lab" %in% names(out))
  expect_equal(out$bw_lab, c("0.1", "0.25"))
})
