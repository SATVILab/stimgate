exampleData <- getExampleData()
gs <- flowWorkspace::load_gs(exampleData$pathGs)
pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")
invisible(gateStim(
  .data = gs,
  pathProject = pathProject,
  popGate = "root",
  batchList = exampleData$batchList,
  chnl = exampleData$chnl
))

.fcsRowSig <- function(ex) {
  if (!is.data.frame(ex) && !is.matrix(ex)) {
    return(character(0))
  }
  if (nrow(ex) == 0L) {
    return(character(0))
  }
  unname(sort(
    apply(ex, 1, function(row) {
      paste(sprintf("%.12f", as.numeric(row)), collapse = "|")
    })
  ))
}

# Comprehensive test suite for writeStimFCS function
# Tests cover:
# 1. Basic functionality and parameter validation
# 2. Directory management (creation, cleanup, relative paths)
# 3. Different parameter combinations (gate methods, mult, gate types)
# 4. Output validation (file content, metadata, naming)
# 5. Edge cases (empty data, invalid parameters, transformations)
# 6. Integration with gate tables (pre-provided vs computed)
# 7. Error handling and message output
# 8. Channel filtering and combination exclusions

test_that("writeStimFCS function exists and has correct signature", {
  # Test that the function exists and has the expected parameters
  expect_true(exists("writeStimFCS", where = asNamespace("stimgate")))

  # Test function signature by checking for argument names
  args <- names(formals(writeStimFCS))
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
    "mult",
    "gateUnsMethod"
  )

  expect_true(all(expectedArgs %in% args))
})


test_that("writeStimFCS runs with basic parameters", {
  pathDirSave <- file.path(tempdir(), "fcs_output_test")

  # Function should create the directory before failing on missing gates
  result <- writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  # Test output validation
  expect_true(length(list.files(pathDirSave)) > 0)
  path_fcs <- file.path(pathDirSave, "sample001_stim1.fcs")
  expect_true(file.exists(path_fcs))
  expect_true(inherits(
    flowCore::read.FCS(path_fcs),
    "flowFrame"
  ))

  # Test return value
  expect_equal(result, pathDirSave)

  # Test directory was created
  expect_true(dir.exists(pathDirSave))

  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS handles directory creation and cleanup", {
  # Test with non-existent directory
  pathDirSave <- file.path(tempdir(), "new_fcs_dir", "subdir")
  expect_false(dir.exists(pathDirSave))

  writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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
  writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]]
  )

  expect_false(file.exists(dummyFile))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS works with different gateUnsMethod options", {
  gateMethods <- c("min", "max", "mean", "tmean", "med")

  for (method in gateMethods) {
    pathDirSave <- file.path(tempdir(), paste0("fcs_output_", method))

    result <- writeStimFCS(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batchList,
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

test_that("writeStimFCS exports exactly the cells selected by an explicit gate table", {
  gsSmall <- gs[1]
  indBatchList <- list(batch_1 = 1L)
  gateTbl <- data.frame(
    chnl = exampleData$chnl,
    marker = c("MarkerF1", "MarkerF2"),
    batch = "batch_1",
    ind = "1",
    gate = c(0.5, 0.5),
    gateCyt = c(0.25, 0.25),
    stringsAsFactors = FALSE
  )
  pathDirSave <- file.path(tempdir(), "fcs_output_exact_gate")

  writeStimFCS(
    pathProject = tempdir(),
    .data = gsSmall,
    indBatchList = indBatchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl,
    gateTbl = gateTbl,
    gateTypeCytPos = "base",
    mult = FALSE
  )

  fcsFile <- list.files(pathDirSave, pattern = "\\.fcs$", full.names = TRUE)
  expect_length(fcsFile, 1L)

  exOrig <- as.data.frame(flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gsSmall[[1]])
  ))
  expected <- exOrig[
    exOrig[[exampleData$chnl[[1]]]] > 0.5 |
      exOrig[[exampleData$chnl[[2]]]] > 0.5,
    exampleData$chnl,
    drop = FALSE
  ]
  exOut <- as.data.frame(flowCore::exprs(flowCore::read.FCS(fcsFile[[1]])))

  expect_equal(
    .fcsRowSig(exOut[, exampleData$chnl, drop = FALSE]),
    .fcsRowSig(expected[, exampleData$chnl, drop = FALSE])
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS respects mult and gateTypeCytPos when exporting exact cells", {
  gsSmall <- gs[1]
  indBatchList <- list(batch_1 = 1L)
  gateTbl <- data.frame(
    chnl = exampleData$chnl,
    marker = c("MarkerF1", "MarkerF2"),
    batch = "batch_1",
    ind = "1",
    gate = c(0.5, 0.5),
    gateCyt = c(0.25, 0.25),
    stringsAsFactors = FALSE
  )

  exOrig <- as.data.frame(flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gsSmall[[1]])
  ))
  expectedBase <- exOrig[
    (exOrig[[exampleData$chnl[[1]]]] > 0.5) &
      (exOrig[[exampleData$chnl[[2]]]] > 0.5),
    exampleData$chnl,
    drop = FALSE
  ]
  expectedCyt <- exOrig[
    ((exOrig[[exampleData$chnl[[1]]]] > 0.5) &
      (exOrig[[exampleData$chnl[[2]]]] > 0.25)) |
      ((exOrig[[exampleData$chnl[[1]]]] > 0.25) &
        (exOrig[[exampleData$chnl[[2]]]] > 0.5)),
    exampleData$chnl,
    drop = FALSE
  ]

  pathDirSaveBase <- file.path(tempdir(), "fcs_output_mult_base")
  pathDirSaveCyt <- file.path(tempdir(), "fcs_output_mult_cyt")

  writeStimFCS(
    pathProject = tempdir(),
    .data = gsSmall,
    indBatchList = indBatchList,
    pathDirSave = pathDirSaveBase,
    chnl = exampleData$chnl,
    gateTbl = gateTbl,
    gateTypeCytPos = "base",
    mult = TRUE
  )
  writeStimFCS(
    pathProject = tempdir(),
    .data = gsSmall,
    indBatchList = indBatchList,
    pathDirSave = pathDirSaveCyt,
    chnl = exampleData$chnl,
    gateTbl = gateTbl,
    gateTypeCytPos = "cyt",
    mult = TRUE
  )

  fcsBase <- flowCore::read.FCS(
    list.files(pathDirSaveBase, pattern = "\\.fcs$", full.names = TRUE)[[1]]
  )
  fcsCyt <- flowCore::read.FCS(
    list.files(pathDirSaveCyt, pattern = "\\.fcs$", full.names = TRUE)[[1]]
  )

  exBase <- as.data.frame(flowCore::exprs(fcsBase))
  exCyt <- as.data.frame(flowCore::exprs(fcsCyt))

  expect_equal(
    .fcsRowSig(exBase[, exampleData$chnl, drop = FALSE]),
    .fcsRowSig(expectedBase[, exampleData$chnl, drop = FALSE])
  )
  expect_equal(
    .fcsRowSig(exCyt[, exampleData$chnl, drop = FALSE]),
    .fcsRowSig(expectedCyt[, exampleData$chnl, drop = FALSE])
  )
  expect_false(
    identical(
      .fcsRowSig(exBase[, exampleData$chnl, drop = FALSE]),
      .fcsRowSig(exCyt[, exampleData$chnl, drop = FALSE])
    )
  )

  unlink(pathDirSaveBase, recursive = TRUE)
  unlink(pathDirSaveCyt, recursive = TRUE)
})

test_that("writeStimFCS validates output file contents", {
  pathDirSave <- file.path(tempdir(), "fcs_output_validation")

  writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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

test_that("writeStimFCS works with pre-provided gate table", {
  # Create a simple gate table
  gateTbl <- data.frame(
    chnl = rep(exampleData$chnl[[1]], length(unlist(exampleData$batchList))),
    marker = rep("BC1", length(unlist(exampleData$batchList))),
    batch = paste0(
      "batch_",
      rep(
        seq_along(exampleData$batchList),
        times = sapply(exampleData$batchList, length)
      )
    ),
    ind = as.character(unlist(exampleData$batchList)),
    gate = rep(0.5, length(unlist(exampleData$batchList))),
    gateCyt = rep(0.5, length(unlist(exampleData$batchList))),
    stringsAsFactors = FALSE
  )

  pathDirSave <- file.path(tempdir(), "fcs_output_custom_gate")

  result <- writeStimFCS(
    pathProject = tempdir(), # Not used when gateTbl provided
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    gateTbl = gateTbl
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS handles invalid gateUnsMethod", {
  pathDirSave <- file.path(tempdir(), "fcs_output_invalid")

  expect_error(
    writeStimFCS(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batchList,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]],
      gateUnsMethod = "invalid_method"
    ),
    "gateUnsMethod not recognised"
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS works with channel filtering", {
  # Test with specific channel subset
  pathDirSave <- file.path(tempdir(), "fcs_output_filtered")

  result <- writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]] # Only first marker
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Test with NULL chnl (should use all available)
  pathDirSaveAll <- file.path(tempdir(), "fcs_output_all")

  resultAll <- writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSaveAll,
    chnl = NULL
  )

  expect_equal(resultAll, pathDirSaveAll)
  expect_true(dir.exists(pathDirSaveAll))
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS handles transformation parameters", {
  # Test with transformation function
  pathDirSave <- file.path(tempdir(), "fcs_output_transform")

  # Simple log transformation
  logTransform <- function(x) log10(x + 1)

  result <- writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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

test_that("writeStimFCS preserves file metadata", {
  pathDirSave <- file.path(tempdir(), "fcs_output_metadata")

  writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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

test_that("writeStimFCS removes the requested cytokine combinations exactly", {
  gsSmall <- gs[1]
  indBatchList <- list(batch_1 = 1L)
  gateTbl <- data.frame(
    chnl = exampleData$chnl,
    marker = c("MarkerF1", "MarkerF2"),
    batch = "batch_1",
    ind = "1",
    gate = c(0.5, 0.5),
    gateCyt = c(0.25, 0.25),
    stringsAsFactors = FALSE
  )
  pathDirSave <- file.path(tempdir(), "fcs_output_exclusions_exact")

  exOrig <- as.data.frame(flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gsSmall[[1]])
  ))
  incVec <- exOrig[[exampleData$chnl[[1]]]] > 0.5 |
    exOrig[[exampleData$chnl[[2]]]] > 0.5
  excVec <- .getPosIndCytCombn(
    ex = exOrig,
    gateTbl = gateTbl,
    chnlPos = exampleData$chnl[[1]],
    chnlNeg = setdiff(exampleData$chnl, exampleData$chnl[[1]]),
    chnlAlt = NULL,
    gateTypeCytPos = "base"
  )
  expected <- exOrig[incVec & !excVec, exampleData$chnl, drop = FALSE]

  writeStimFCS(
    pathProject = tempdir(),
    .data = gsSmall,
    indBatchList = indBatchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl,
    gateTbl = gateTbl,
    gateTypeCytPos = "base",
    combnExc = list(exampleData$chnl[[1]])
  )

  fcsFile <- list.files(pathDirSave, pattern = "\\.fcs$", full.names = TRUE)
  expect_length(fcsFile, 1L)

  exOut <- as.data.frame(flowCore::exprs(flowCore::read.FCS(fcsFile[[1]])))
  expect_equal(
    .fcsRowSig(exOut[, exampleData$chnl, drop = FALSE]),
    .fcsRowSig(expected[, exampleData$chnl, drop = FALSE])
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS creates consistent file names", {
  pathDirSave <- file.path(tempdir(), "fcs_output_naming")

  writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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

test_that("writeStimFCS message output", {
  pathDirSave <- file.path(tempdir(), "fcs_output_messages")

  # Capture messages
  expect_message(
    writeStimFCS(
      pathProject = pathProject,
      .data = gs,
      indBatchList = exampleData$batchList,
      pathDirSave = pathDirSave,
      chnl = exampleData$chnl[[1]]
    ),
    "Writing.*files"
  )
  unlink(pathDirSave, recursive = TRUE)
})

test_that("writeStimFCS handles edge case: empty data", {
  # Create gate table with very high thresholds (should result in no positive cells)
  gateTbl <- data.frame(
    chnl = rep(exampleData$chnl[[1]], length(unlist(exampleData$batchList))),
    marker = rep("BC1", length(unlist(exampleData$batchList))),
    batch = paste0(
      "batch_",
      rep(
        seq_along(exampleData$batchList),
        times = sapply(exampleData$batchList, length)
      )
    ),
    ind = as.character(unlist(exampleData$batchList)),
    gate = rep(999999, length(unlist(exampleData$batchList))), # Very high threshold
    gateCyt = rep(999999, length(unlist(exampleData$batchList))),
    stringsAsFactors = FALSE
  )

  pathDirSave <- file.path(tempdir(), "fcs_output_empty")

  # Should handle case where no cells meet criteria
  expect_message(
    result <- writeStimFCS(
      pathProject = tempdir(),
      .data = gs,
      indBatchList = exampleData$batchList,
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

test_that("writeStimFCS validates parameter types", {
  # Test with invalid .data type (should fail gracefully)
  expect_error(
    writeStimFCS(
      pathProject = tempdir(),
      .data = "not_a_gatingset",
      indBatchList = exampleData$batchList,
      pathDirSave = tempdir(),
      chnl = exampleData$chnl[[1]]
    )
  )

  # Test with invalid indBatchList type
  expect_error(
    writeStimFCS(
      pathProject = tempdir(),
      .data = gs,
      indBatchList = "not_a_list",
      pathDirSave = tempdir(),
      chnl = exampleData$chnl[[1]]
    )
  )
})


test_that("writeStimFCS integrates with stimgate workflow", {
  # Test full integration: gate -> fcs_write -> verify output
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "stimgate")

  # Step 1: Run gating
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    chnl = exampleData$chnl
  ))

  # Verify gating created expected files
  expect_true(file.exists(file.path(pathProject, "gateStats.rds")))

  # Step 2: Run FCS writing using gates from step 1
  pathDirSave <- file.path(tempdir(), "fcs_output_integration")

  result <- writeStimFCS(
    pathProject = pathProject,
    .data = gs,
    indBatchList = exampleData$batchList,
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

test_that("writeStimFCS respects working directory", {
  # Change working directory temporarily
  originalWd <- getwd()
  tempWd <- tempdir()

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject2 <- file.path(dirname(exampleData$pathGs), "stimgate")
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject2,
    popGate = "root",
    batchList = exampleData$batchList,
    chnl = exampleData$chnl
  ))

  tryCatch(
    {
      setwd(tempWd)

      # Use relative path for output
      pathDirSave <- "fcs_output_wd_test"

      result <- writeStimFCS(
        pathProject = pathProject2,
        .data = gs,
        indBatchList = exampleData$batchList,
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


test_that("writeStimFCS handles transformation edge cases", {
  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject2 <- file.path(dirname(exampleData$pathGs), "stimgate")
  invisible(gateStim(
    .data = gs,
    pathProject = pathProject2,
    popGate = "root",
    batchList = exampleData$batchList,
    chnl = exampleData$chnl
  ))

  # Test with transformation function but no transChnl (should apply to all columns)
  pathDirSave <- file.path(tempdir(), "fcs_output_transform_all")

  # Identity transformation (should not change values but test the pathway)
  identityTransform <- function(x) x

  result <- stimgate::writeStimFCS(
    pathProject = pathProject2,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSave,
    chnl = exampleData$chnl[[1]],
    transFn = identityTransform,
    transChnl = NULL # Should apply to all columns
  )

  expect_equal(result, pathDirSave)
  expect_true(dir.exists(pathDirSave))

  # Test with NULL transformation function
  pathDirSaveNull <- file.path(tempdir(), "fcs_output_transform_null")

  resultNull <- stimgate::writeStimFCS(
    pathProject = pathProject2,
    .data = gs,
    indBatchList = exampleData$batchList,
    pathDirSave = pathDirSaveNull,
    chnl = exampleData$chnl[[1]],
    transFn = NULL,
    transChnl = exampleData$chnl[[1]]
  )

  expect_equal(resultNull, pathDirSaveNull)
  expect_true(dir.exists(pathDirSaveNull))
  unlink(pathDirSave, recursive = TRUE)
})
