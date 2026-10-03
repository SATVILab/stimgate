root_dir <- normalizePath(
  file.path(testthat::test_path(), "../../.."),
  mustWork = TRUE
)

test_that("analysis QMDs do not overwrite sourced helper functions", {
  helper_env <- new.env(parent = getNamespace("stimgate"))
  for (file in c(
    "analysis-runtime.R",
    "sim-misc.R",
    "sim-bandwidth.R",
    "sim-bandwidth-analysis-io.R",
    "sim-bandwidth-analysis-plot.R",
    "sim-compare-freq_bs.R",
    "sim-trans.R"
  )) {
    source(file.path(root_dir, "scripts", "r", file), local = helper_env)
  }

  helper_names <- ls(helper_env, all.names = TRUE)
  helper_names <- helper_names[vapply(
    helper_names,
    function(name) is.function(get(name, helper_env, inherits = FALSE)),
    logical(1)
  )]

  violations <- character()
  qmd_paths <- list.files(
    file.path(root_dir, "analysis"),
    pattern = "[.]qmd$",
    full.names = TRUE
  )

  for (qmd_path in qmd_paths) {
    lines <- readLines(qmd_path, warn = FALSE)
    for (i in seq_along(lines)) {
      match <- regmatches(
        lines[[i]],
        regexec(
          "^\\s*([.A-Za-z][A-Za-z0-9._]*)\\s*<-\\s*(.*)$",
          lines[[i]],
          perl = TRUE
        )
      )[[1]]
      if (length(match) == 0L) {
        next
      }

      lhs <- match[[2]]
      rhs <- trimws(match[[3]])
      if (
        lhs %in% helper_names &&
          !grepl("^function\\s*\\(", rhs, perl = TRUE)
      ) {
        violations <- c(
          violations,
          paste0(basename(qmd_path), ":", i, ": ", trimws(lines[[i]]))
        )
      }
    }
  }

  expect_identical(violations, character())
})

test_that("analysis 8 uses deterministic scenario seeds and full-grid promotion", {
  qmd_path <- file.path(
    root_dir,
    "analysis",
    "8-sim-compare-freq_bs-batch.qmd"
  )
  content <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("simulation_seed:\\s*1", content))
  expect_true(grepl("comparison_semantics_version", content, fixed = TRUE))
  expect_true(grepl(
    "sim_seed = as.integer(simulation_seed + sim_id - 1L)",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "set.seed(as.integer(row$sim_seed[[1]]))",
    content,
    fixed = TRUE
  ))
  expect_true(grepl(
    "path_progress_file <- run_ctx$progress_file",
    content,
    fixed = TRUE
  ))
  expect_true(grepl("recursive\\s*=\\s*TRUE", content))
  expect_true(grepl("expected_sim_ids", content, fixed = TRUE))
  expect_true(grepl("Refusing to promote analysis 8", content, fixed = TRUE))
  expect_true(grepl(".analysis_current_file", content, fixed = TRUE))
})
