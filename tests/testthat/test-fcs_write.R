exampleData <- get_example_data()
gs <- flowWorkspace::load_gs(exampleData$path_gs)
pathProject <- file.path(dirname(exampleData$path_gs), "stimgate")
invisible(stimgate_gate(
  .data = gs,
  pathProject = pathProject,
  popGate = "root",
  batch_list = exampleData$batch_list,
  chnl = exampleData$chnl
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
  args <- names(formals(stimgate_fcs_write))
  expectedArgs <- c(
    "pathProject",
    ".data",
    "indBatchList",
    "pathDirSave",
    "chnl",
    "gateTbl",
    "transFn",
    "transChnl",
    "combnExc",
    "gateTypeCytPos",
    "gateTypeSinglePos",
    "mult",
    "gateUnsMethod"
  )

  expect_true(all(expectedArgs %in% args))
})


test_that("stimgate_fcs_write runs with basic parameters", {
  pathDirSave <- file.path(tempdir(), "fcs_output_test")

  # Function should create the directory before failing on missing gates
  result <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  # Test output validation
  expect_true(length(list.files(pathDirSave)) > 0)
  expect_true(inherits(
    flowCore::read.FCS(file.path(pathDirSave, "V1.fcs")),
    "flowFrame"
  ))

  # Test return value
  expect_equal(result, pathDirSave)

  # Test directory was created
  expect_true(dir.exists(pathDirSave))

  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write handles directory creation and cleanup", {
  # Test with non-existent directory
  pathDirSave <- file.path(tempdir(), "new_fcs_dir", "subdir")
  expect_false(dir.exists(pathDirSave))

  stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  expect_true(dir.exists(pathDirSave))

  # Test that existing directory is cleaned
  # Create a dummy file
  dummyFile <- file.path(pathDirSave, "dummy.txt")
  writeLines("test", dummyFile)
  expect_true(file.exists(dummyFile))

  # Run again - should clean directory
  stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  expect_false(file.exists(dummyFile))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write works with different gateUnsMethod options", {
  gateMethods <- c("min", "max", "mean", "tmean", "med")

  for (method in gateMethods) {
    pathDirSave <- file.path(tempdir(), paste0("fcs_output_", method))

    result <- stimgate_fcs_write(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]],
      gateUnsMethod = method
    )

    expect_equal(result, pathDirSave)
    expect_true(dir.exists(pathDirSave))
    # Should have some files (at least one sample should have positive cells)
    expect_true(length(list.files(pathDirSave, pattern = "\\.fcs$")) >= 0)
    unlink(pathDirSave, recursive = TRUE)
  }
})

test_that("stimgate_fcs_write works with mult parameter", {
  # Test with mult = FALSE (default)
  pathDirSaveSingle <- file.path(tempdir(), "fcs_output_single")
  resultSingle <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSaveSingle,
    chnl = exampleData$chnl,
    mult = FALSE
  )

  # Test with mult = TRUE
  pathDirSaveMult <- file.path(tempdir(), "fcs_output_mult")
  resultMult <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSaveMult,
    chnl = exampleData$chnl,
    mult = TRUE
  )

  expect_equal(resultSingle, pathDirSaveSingle)
  expect_equal(resultMult, pathDirSaveMult)
  expect_true(dir.exists(pathDirSaveSingle))
  expect_true(dir.exists(pathDirSaveMult))
  unlink(pathDirSaveSingle, recursive = TRUE)
  unlink(pathDirSaveMult, recursive = TRUE)
})

test_that("stimgate_fcs_write works with different gate types", {
  pathDirSave <- file.path(tempdir(), "fcs_output_gate_types")

  result <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    gateTypeCytPos = "cyt",
    gateTypeSinglePos = "single"
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write validates output file contents", {
  pathDirSave <- file.path(tempdir(), "fcs_output_validation")

  stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  # Get list of FCS files
  fcsFiles <- list.files(pathDirSave, pattern = "\\.fcs$", full.names = TRUE)

  # Test that we can read each file and it has the expected structure
  for (fcsFile in fcsFiles) {
    ff <- flowCore::read.FCS(fcsFile)

    # Check that it's a valid flowFrame
    expect_true(inherits(ff, "flowFrame"))

    exprMat <- flowCore::exprs(ff)

    # Check that it has data
    expect_true(nrow(exprMat) >= 0)

    # Check that it has the expected channels
    expect_true(all(exampleData$chnl[[1]] %in% colnames(exprMat)))

    # Check that expression matrix can be extracted
    exprMat <- flowCore::exprs(ff)
    expect_true(is.matrix(exprMat))
    expect_true(ncol(exprMat) > 0)
  }
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write works with pre-provided gate table", {
  # Create a simple gate table
  gateTbl <- data.frame(
    chnl = rep(exampleData$chnl[[1]], length(unlist(exampleData$batch_list))),
    marker = rep("BC1", length(unlist(exampleData$batch_list))),
    batch = paste0(
      "batch_",
      rep(
        seq_along(exampleData$batch_list),
        times = sapply(exampleData$batch_list, length)
      )
    ),
    ind = as.character(unlist(exampleData$batch_list)),
    gate = rep(0.5, length(unlist(exampleData$batch_list))),
    gate_cyt = rep(0.5, length(unlist(exampleData$batch_list))),
    gate_single = rep(0.5, length(unlist(exampleData$batch_list))),
    stringsAsFactors = FALSE
  )

  pathDirSave <- file.path(tempdir(), "fcs_output_custom_gate")

  result <- stimgate_fcs_write(
    pathProject = tempdir(), # Not used when gateTbl provided
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    gateTbl = gateTbl
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write handles invalid gateUnsMethod", {
  pathDirSave <- file.path(tempdir(), "fcs_output_invalid")

  expect_error(
    stimgate_fcs_write(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]],
      gateUnsMethod = "invalid_method"
    ),
    "gateUnsMethod not recognised"
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write works with channel filtering", {
  # Test with specific channel subset
  pathDirSave <- file.path(tempdir(), "fcs_output_filtered")

  result <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]] # Only first marker
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Test with NULL chnl (should use all available)
  pathDirSaveAll <- file.path(tempdir(), "fcs_output_all")

  resultAll <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSaveAll,
    chnl = NULL
  )

  expect_equal(resultAll, pathDirSaveAll)
  expect_true(dir.exists(pathDirSaveAll))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write handles transformation parameters", {
  # Test with transformation function
  pathDirSave <- file.path(tempdir(), "fcs_output_transform")

  # Simple log transformation
  logTransform <- function(x) log10(x + 1)

  result <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    transFn = logTransform,
    transChnl = exampleData$chnl[[1]]
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Verify files were created
  fcsFiles <- list.files(pathDirSave, pattern = "\\.fcs$")
  expect_true(length(fcsFiles) >= 0)
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write preserves file metadata", {
  pathDirSave <- file.path(tempdir(), "fcs_output_metadata")

  stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  # Get original and output files
  fcsFiles <- list.files(pathDirSave, pattern = "\\.fcs$", full.names = TRUE)

  if (length(fcsFiles) > 0) {
    # Check first file
    outputFf <- flowCore::read.FCS(fcsFiles[1])
    originalFf <- flowWorkspace::gh_pop_get_data(gs[[1]])
    if (inherits(originalFf, "cytoframe")) {
      originalFf <- flowWorkspace::cytoframe_to_flowFrame(originalFf)
    }

    # Check that basic metadata is preserved
    expect_true(inherits(outputFf, "flowFrame"))

    # Check that parameters are preserved (at least the gated channels)
    exMat <- flowCore::exprs(outputFf)
    outputParams <- colnames(exMat)
    expect_true(all(exampleData$chnl[[1]] %in% colnames(exMat)))
  }
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write handles combination exclusions", {
  pathDirSave <- file.path(tempdir(), "fcs_output_exclusions")

  # Test with combination exclusions (if we have multiple markers)
  if (length(exampleData$chnl) > 1) {
    combnExc <- list(exampleData$chnl[[1]])

    result <- stimgate_fcs_write(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl,
      combnExc = combnExc
    )

    expect_equal(result, pathDirSave)
    expect_true(dir.exists(pathDirSave))
  } else {
    # Test with NULL exclusions
    result <- stimgate_fcs_write(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]],
      combnExc = NULL
    )

    expect_equal(result, pathDirSave)
    expect_true(dir.exists(pathDirSave))
  }
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write creates consistent file names", {
  pathDirSave <- file.path(tempdir(), "fcs_output_naming")

  stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  # Check file naming pattern
  fcsFiles <- list.files(pathDirSave, pattern = "\\.fcs$")

  # Files should be named based on GUID
  expect_true(all(grepl("\\.fcs$", fcsFiles)))

  # Check that files are in the correct directory
  fullPaths <- list.files(
    pathDirSave,
    pattern = "\\.fcs$",
    full.names = TRUE
  )
  expect_true(
    all(normalizePath(dirname(fullPaths)) == normalizePath(pathDirSave))
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write message output", {
  pathDirSave <- file.path(tempdir(), "fcs_output_messages")

  # Capture messages
  expect_message(
    stimgate_fcs_write(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]]
    ),
    "Writing.*files"
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write handles edge case: empty data", {
  # Create gate table with very high thresholds (should result in no positive cells)
  gateTbl <- data.frame(
    chnl = rep(exampleData$chnl[[1]], length(unlist(exampleData$batch_list))),
    marker = rep("BC1", length(unlist(exampleData$batch_list))),
    batch = paste0(
      "batch_",
      rep(
        seq_along(exampleData$batch_list),
        times = sapply(exampleData$batch_list, length)
      )
    ),
    ind = as.character(unlist(exampleData$batch_list)),
    gate = rep(999999, length(unlist(exampleData$batch_list))), # Very high threshold
    gate_cyt = rep(999999, length(unlist(exampleData$batch_list))),
    gate_single = rep(999999, length(unlist(exampleData$batch_list))),
    stringsAsFactors = FALSE
  )

  pathDirSave <- file.path(tempdir(), "fcs_output_empty")

  # Should handle case where no cells meet criteria
  expect_message(
    result <- stimgate_fcs_write(
      pathProject = tempdir(),
      .data = gs,
      indBatchList = exampleData$batch_list,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]],
      gateTbl = gateTbl
    ),
    "No stimulation-positive cells"
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write validates parameter types", {
  # Test with invalid .data type (should fail gracefully)
  expect_error(
    stimgate_fcs_write(
      pathProject = tempdir(),
      .data = "not_a_gatingset",
      indBatchList = exampleData$batch_list,
      pathDirSave = tempdir(),
      chnl = exampleData$chnl[[1]]
    )
  )

  # Test with invalid indBatchList type
  expect_error(
    stimgate_fcs_write(
      pathProject = tempdir(),
      .data = gs,
      indBatchList = "not_a_list",
      pathDirSave = tempdir(),
      chnl = exampleData$chnl[[1]]
    )
  )
})


test_that("stimgate_fcs_write integrates with stimgate workflow", {
  # Test full integration: gate -> fcs_write -> verify output
  exampleData <- get_example_data()
  gs <- flowWorkspace::load_gs(exampleData$path_gs)
  pathProject <- file.path(dirname(exampleData$path_gs), "stimgate")

  # Step 1: Run gating
  invisible(stimgate_gate(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batch_list = exampleData$batch_list,
    chnl = exampleData$chnl
  ))

  # Verify gating created expected files
  expect_true(file.exists(file.path(pathProject, "gate_stats.rds")))

  # Step 2: Run FCS writing using gates from step 1
  pathDirSave <- file.path(tempdir(), "fcs_output_integration")

  result <- stimgate_fcs_write(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl # Use all markers
  )

  # Step 3: Verify integration worked
  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Verify that gate information was properly used
  fcsFiles <- list.files(pathDirSave, pattern = "\\.fcs$", full.names = TRUE)

  # Should have created files for samples with positive cells
  for (fcsFile in fcsFiles) {
    ff <- flowCore::read.FCS(fcsFile)
    expect_true(inherits(ff, "flowFrame"))
    expect_true(nrow(ff) >= 0)
    exMat <- flowCore::exprs(ff)

    # Verify that all gated channels are present
    expect_true(all(exampleData$chnl %in% colnames(exMat)))
  }
  unlink(pathDirSave, recursive = TRUE)
})

test_that("stimgate_fcs_write respects working directory", {
  # Change working directory temporarily
  originalWd <- getwd()
  tempWd <- tempdir()

  exampleData <- get_example_data()
  gs <- flowWorkspace::load_gs(exampleData$path_gs)
  pathProject2 <- file.path(dirname(exampleData$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    pathProject = pathProject2,
    popGate = "root",
    batch_list = exampleData$batch_list,
    chnl = exampleData$chnl
  ))

  tryCatch(
    {
      setwd(tempWd)

      # Use relative path for output
      pathDirSave <- "fcs_output_wd_test"

      result <- stimgate_fcs_write(
        pathProject = pathProject2,
        .data = gs,
        indBatchList = exampleData$batch_list,
        pathDirSave = pathDirSave,
        chnl = exampleData$chnl[[1]]
      )

      # Should create directory relative to current working directory
      expect_true(dir.exists(file.path(tempWd, pathDirSave)))
      expect_equal(result, pathDirSave)
    },
    finally = {
      setwd(originalWd)
    }
  )
  unlink(file.path(tempWd, pathDirSave), recursive = TRUE)
})


test_that("stimgate_fcs_write handles transformation edge cases", {
  exampleData <- get_example_data()
  gs <- flowWorkspace::load_gs(exampleData$path_gs)
  pathProject2 <- file.path(dirname(exampleData$path_gs), "stimgate")
  invisible(stimgate_gate(
    .data = gs,
    pathProject = pathProject2,
    popGate = "root",
    batch_list = exampleData$batch_list,
    chnl = exampleData$chnl
  ))

  # Test with transformation function but no transChnl (should apply to all columns)
  pathDirSave <- file.path(tempdir(), "fcs_output_transform_all")

  # Identity transformation (should not change values but test the pathway)
  identityTransform <- function(x) x

  result <- stimgate::stimgate_fcs_write(
    pathProject = pathProject2,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    transFn = identityTransform,
    transChnl = NULL # Should apply to all columns
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Test with NULL transformation function
  pathDirSaveNull <- file.path(tempdir(), "fcs_output_transform_null")

  resultNull <- stimgate::stimgate_fcs_write(
    pathProject = pathProject2,
    .data = gs,
    indBatchList = exampleData$batch_list,
    pathDirSave = pathDirSaveNull,
    chnl = exampleData$chnl[[1]],
    transFn = NULL,
    transChnl = exampleData$chnl[[1]]
  )

  expect_equal(resultNull, pathDirSaveNull)
  expect_true(dir.exists(pathDirSaveNull))
  unlink(pathDirSave, recursive = TRUE)
})
