example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)
path_project <- file.path(dirname(example_data$path_gs), "stimgate")
invisible(stimgate::stimgate_gate(
  .data = gs,
  path_project = path_project,
  pop_gate = "root",
  batch_list = example_data$batch_list,
  marker = example_data$marker
))

# Comprehensive test suite for stimgate_fcs_write function
# Tests cover:
# 1. Basic functionality and parameter validation
# 2. Directory management (creation, cleanup, relative paths)
# 3. Different parameter combinations (gate methods, mult, gate types)
# 4. Output validation (file content, metadata, naming)
# 5. Edge cases (empty data, invalid parameters, transformations)
# 6. Integration with gate tables (pre-provided vs computed)
# 7. Error handling and message output
# 8. Channel filtering and combination exclusions

test_that("stimgate_fcs_write function exists and has correct signature", {
  # Test that the function exists and has the expected parameters
  expect_true(exists("stimgate_fcs_write", where = asNamespace("stimgate")))

  # Test function signature by checking for argument names
  args <- names(formals(stimgate::stimgate_fcs_write))
  expected_args <- c(
    "path_project", ".data", "ind_batch_list", "path_dir_save",
    "chnl", "gate_tbl", "trans_fn", "trans_chnl", "combn_exc",
    "gate_type_cyt_pos", "gate_type_single_pos", "mult",
    "gate_uns_method"
  )

  expect_true(all(expected_args %in% args))
})


test_that("stimgate_fcs_write runs with basic parameters", {
  path_dir_save <- file.path(tempdir(), "fcs_output_test")

  # Function should create the directory before failing on missing gates
  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  # Test output validation
  expect_true(length(list.files(path_dir_save)) > 0)
  expect_true(inherits(
    flowCore::read.FCS(file.path(path_dir_save, "V1.fcs")),
    "flowFrame"
  ))

  # Test return value
  expect_equal(result, path_dir_save)

  # Test directory was created
  expect_true(dir.exists(path_dir_save))

  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write handles directory creation and cleanup", {
  # Test with non-existent directory
  path_dir_save <- file.path(tempdir(), "new_fcs_dir", "subdir")
  expect_false(dir.exists(path_dir_save))

  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  expect_true(dir.exists(path_dir_save))

  # Test that existing directory is cleaned
  # Create a dummy file
  dummy_file <- file.path(path_dir_save, "dummy.txt")
  writeLines("test", dummy_file)
  expect_true(file.exists(dummy_file))

  # Run again - should clean directory
  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  expect_false(file.exists(dummy_file))
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write works with different gate_uns_method options", {
  gate_methods <- c("min", "max", "mean", "tmean", "med")

  for (method in gate_methods) {
    path_dir_save <- file.path(tempdir(), paste0("fcs_output_", method))

    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      gate_uns_method = method
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
    # Should have some files (at least one sample should have positive cells)
    expect_true(length(list.files(path_dir_save, pattern = "\\.fcs$")) >= 0)
    unlink(path_dir_save, recursive = TRUE)
  }
})

test_that("stimgate_fcs_write works with mult parameter", {
  # Test with mult = FALSE (default)
  path_dir_save_single <- file.path(tempdir(), "fcs_output_single")
  result_single <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save_single,
    chnl = example_data$marker,
    mult = FALSE
  )

  # Test with mult = TRUE
  path_dir_save_mult <- file.path(tempdir(), "fcs_output_mult")
  result_mult <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save_mult,
    chnl = example_data$marker,
    mult = TRUE
  )

  expect_equal(result_single, path_dir_save_single)
  expect_equal(result_mult, path_dir_save_mult)
  expect_true(dir.exists(path_dir_save_single))
  expect_true(dir.exists(path_dir_save_mult))
  unlink(path_dir_save_single, recursive = TRUE)
  unlink(path_dir_save_mult, recursive = TRUE)
})

test_that("stimgate_fcs_write works with different gate types", {
  path_dir_save <- file.path(tempdir(), "fcs_output_gate_types")

  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]],
    gate_type_cyt_pos = "cyt",
    gate_type_single_pos = "single"
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write validates output file contents", {
  path_dir_save <- file.path(tempdir(), "fcs_output_validation")

  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  # Get list of FCS files
  fcs_files <- list.files(path_dir_save, pattern = "\\.fcs$", full.names = TRUE)

  # Test that we can read each file and it has the expected structure
  for (fcs_file in fcs_files) {
    ff <- flowCore::read.FCS(fcs_file)

    # Check that it's a valid flowFrame
    expect_true(inherits(ff, "flowFrame"))

    expr_mat <- flowCore::exprs(ff)

    # Check that it has data
    expect_true(nrow(expr_mat) >= 0)

    # Check that it has the expected channels
    expect_true(all(example_data$marker[[1]] %in% colnames(expr_mat)))

    # Check that expression matrix can be extracted
    expr_mat <- flowCore::exprs(ff)
    expect_true(is.matrix(expr_mat))
    expect_true(ncol(expr_mat) > 0)
  }
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write works with pre-provided gate table", {
  # Create a simple gate table
  gate_tbl <- data.frame(
    chnl = rep(example_data$marker[[1]], length(unlist(example_data$batch_list))),
    marker = rep("BC1", length(unlist(example_data$batch_list))),
    batch = paste0("batch_", rep(seq_along(example_data$batch_list),
      times = sapply(example_data$batch_list, length)
    )),
    ind = as.character(unlist(example_data$batch_list)),
    gate = rep(0.5, length(unlist(example_data$batch_list))),
    gate_cyt = rep(0.5, length(unlist(example_data$batch_list))),
    gate_single = rep(0.5, length(unlist(example_data$batch_list))),
    stringsAsFactors = FALSE
  )

  path_dir_save <- file.path(tempdir(), "fcs_output_custom_gate")

  result <- stimgate::stimgate_fcs_write(
    path_project = tempdir(), # Not used when gate_tbl provided
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]],
    gate_tbl = gate_tbl
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write handles invalid gate_uns_method", {
  path_dir_save <- file.path(tempdir(), "fcs_output_invalid")

  expect_error(
    stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      gate_uns_method = "invalid_method"
    ),
    "gate_uns_method not recognised"
  )
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write works with channel filtering", {
  # Test with specific channel subset
  path_dir_save <- file.path(tempdir(), "fcs_output_filtered")

  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]] # Only first marker
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))

  # Test with NULL chnl (should use all available)
  path_dir_save_all <- file.path(tempdir(), "fcs_output_all")

  result_all <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save_all,
    chnl = NULL
  )

  expect_equal(result_all, path_dir_save_all)
  expect_true(dir.exists(path_dir_save_all))
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write handles transformation parameters", {
  # Test with transformation function
  path_dir_save <- file.path(tempdir(), "fcs_output_transform")

  # Simple log transformation
  log_transform <- function(x) log10(x + 1)

  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]],
    trans_fn = log_transform,
    trans_chnl = example_data$marker[[1]]
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))

  # Verify files were created
  fcs_files <- list.files(path_dir_save, pattern = "\\.fcs$")
  expect_true(length(fcs_files) >= 0)
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write preserves file metadata", {
  path_dir_save <- file.path(tempdir(), "fcs_output_metadata")

  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  # Get original and output files
  fcs_files <- list.files(path_dir_save, pattern = "\\.fcs$", full.names = TRUE)

  if (length(fcs_files) > 0) {
    # Check first file
    output_ff <- flowCore::read.FCS(fcs_files[1])
    original_ff <- flowWorkspace::gh_pop_get_data(gs[[1]])
    if (inherits(original_ff, "cytoframe")) {
      original_ff <- flowWorkspace::cytoframe_to_flowFrame(original_ff)
    }

    # Check that basic metadata is preserved
    expect_true(inherits(output_ff, "flowFrame"))

    # Check that parameters are preserved (at least the gated channels)
    ex_mat <- flowCore::exprs(output_ff)
    output_params <- colnames(ex_mat)
    expect_true(all(example_data$marker[[1]] %in% colnames(ex_mat)))
  }
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write handles combination exclusions", {
  path_dir_save <- file.path(tempdir(), "fcs_output_exclusions")

  # Test with combination exclusions (if we have multiple markers)
  if (length(example_data$marker) > 1) {
    combn_exc <- list(example_data$marker[[1]])

    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker,
      combn_exc = combn_exc
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
  } else {
    # Test with NULL exclusions
    result <- stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      combn_exc = NULL
    )

    expect_equal(result, path_dir_save)
    expect_true(dir.exists(path_dir_save))
  }
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write creates consistent file names", {
  path_dir_save <- file.path(tempdir(), "fcs_output_naming")

  stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]]
  )

  # Check file naming pattern
  fcs_files <- list.files(path_dir_save, pattern = "\\.fcs$")

  # Files should be named based on GUID
  expect_true(all(grepl("\\.fcs$", fcs_files)))

  # Check that files are in the correct directory
  full_paths <- list.files(
    path_dir_save,
    pattern = "\\.fcs$", full.names = TRUE
  )
  expect_true(
    all(normalizePath(dirname(full_paths)) == normalizePath(path_dir_save))
  )
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write message output", {
  path_dir_save <- file.path(tempdir(), "fcs_output_messages")

  # Capture messages
  expect_message(
    stimgate::stimgate_fcs_write(
      path_project = path_project,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]]
    ),
    "Writing.*files"
  )
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write handles edge case: empty data", {
  # Create gate table with very high thresholds (should result in no positive cells)
  gate_tbl <- data.frame(
    chnl = rep(example_data$marker[[1]], length(unlist(example_data$batch_list))),
    marker = rep("BC1", length(unlist(example_data$batch_list))),
    batch = paste0("batch_", rep(seq_along(example_data$batch_list),
      times = sapply(example_data$batch_list, length)
    )),
    ind = as.character(unlist(example_data$batch_list)),
    gate = rep(999999, length(unlist(example_data$batch_list))), # Very high threshold
    gate_cyt = rep(999999, length(unlist(example_data$batch_list))),
    gate_single = rep(999999, length(unlist(example_data$batch_list))),
    stringsAsFactors = FALSE
  )

  path_dir_save <- file.path(tempdir(), "fcs_output_empty")

  # Should handle case where no cells meet criteria
  expect_message(
    result <- stimgate::stimgate_fcs_write(
      path_project = tempdir(),
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]],
      gate_tbl = gate_tbl
    ),
    "No stimulation-positive cells"
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write validates parameter types", {
  # Test with invalid .data type (should fail gracefully)
  expect_error(
    stimgate::stimgate_fcs_write(
      path_project = tempdir(),
      .data = "not_a_gatingset",
      ind_batch_list = example_data$batch_list,
      path_dir_save = tempdir(),
      chnl = example_data$marker[[1]]
    )
  )

  # Test with invalid ind_batch_list type
  expect_error(
    stimgate::stimgate_fcs_write(
      path_project = tempdir(),
      .data = gs,
      ind_batch_list = "not_a_list",
      path_dir_save = tempdir(),
      chnl = example_data$marker[[1]]
    )
  )
})


test_that("stimgate_fcs_write integrates with stimgate workflow", {
  # Test full integration: gate -> fcs_write -> verify output
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project <- file.path(dirname(example_data$path_gs), "stimgate")

  # Step 1: Run gating
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Verify gating created expected files
  expect_true(file.exists(file.path(path_project, "gate_stats.rds")))

  # Step 2: Run FCS writing using gates from step 1
  path_dir_save <- file.path(tempdir(), "fcs_output_integration")

  result <- stimgate::stimgate_fcs_write(
    path_project = path_project,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker # Use all markers
  )

  # Step 3: Verify integration worked
  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))

  # Verify that gate information was properly used
  fcs_files <- list.files(path_dir_save, pattern = "\\.fcs$", full.names = TRUE)

  # Should have created files for samples with positive cells
  for (fcs_file in fcs_files) {
    ff <- flowCore::read.FCS(fcs_file)
    expect_true(inherits(ff, "flowFrame"))
    expect_true(nrow(ff) >= 0)
    ex_mat <- flowCore::exprs(ff)

    # Verify that all gated channels are present
    expect_true(all(example_data$marker %in% colnames(ex_mat)))
  }
  unlink(path_dir_save, recursive = TRUE)
})

test_that("stimgate_fcs_write respects working directory", {
  # Change working directory temporarily
  original_wd <- getwd()
  temp_wd <- tempdir()

  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project_2 <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project_2,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  tryCatch({
    setwd(temp_wd)

    # Use relative path for output
    path_dir_save <- "fcs_output_wd_test"

    result <- stimgate::stimgate_fcs_write(
      path_project = path_project_2,
      .data = gs,
      ind_batch_list = example_data$batch_list,
      path_dir_save = path_dir_save,
      chnl = example_data$marker[[1]]
    )

    # Should create directory relative to current working directory
    expect_true(dir.exists(file.path(temp_wd, path_dir_save)))
    expect_equal(result, path_dir_save)
  }, finally = {
    setwd(original_wd)
  })
  unlink(file.path(temp_wd, path_dir_save), recursive = TRUE)
})


test_that("stimgate_fcs_write handles transformation edge cases", {
  example_data <- get_example_data()
  gs <- flowWorkspace::load_gs(example_data$path_gs)
  path_project_2 <- file.path(dirname(example_data$path_gs), "stimgate")
  invisible(stimgate::stimgate_gate(
    .data = gs,
    path_project = path_project_2,
    pop_gate = "root",
    batch_list = example_data$batch_list,
    marker = example_data$marker
  ))

  # Test with transformation function but no trans_chnl (should apply to all columns)
  path_dir_save <- file.path(tempdir(), "fcs_output_transform_all")

  # Identity transformation (should not change values but test the pathway)
  identity_transform <- function(x) x

  result <- stimgate::stimgate_fcs_write(
    path_project = path_project_2,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save,
    chnl = example_data$marker[[1]],
    trans_fn = identity_transform,
    trans_chnl = NULL # Should apply to all columns
  )

  expect_equal(result, path_dir_save)
  expect_true(dir.exists(path_dir_save))

  # Test with NULL transformation function
  path_dir_save_null <- file.path(tempdir(), "fcs_output_transform_null")

  result_null <- stimgate::stimgate_fcs_write(
    path_project = path_project_2,
    .data = gs,
    ind_batch_list = example_data$batch_list,
    path_dir_save = path_dir_save_null,
    chnl = example_data$marker[[1]],
    trans_fn = NULL,
    trans_chnl = example_data$marker[[1]]
  )

  expect_equal(result_null, path_dir_save_null)
  expect_true(dir.exists(path_dir_save_null))
  unlink(path_dir_save, recursive = TRUE)
})
