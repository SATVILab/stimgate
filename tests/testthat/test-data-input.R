testthat::test_that("cytometry sets and individual frames preserve expression", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  cs <- flowWorkspace::gs_pop_get_data(gs, y = "root")
  matrices <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]], y = "root"))
  })
  frames <- lapply(seq_along(gs), function(i) {
    fr <- flowWorkspace::gh_pop_get_data(gs[[i]], y = "root")
    flowCore::flowFrame(matrices[[i]], parameters = flowCore::parameters(fr))
  })
  fs <- flowCore::flowSet(frames)
  flowCore::sampleNames(fs) <- flowWorkspace::sampleNames(gs)

  testthat::expect_identical(.asStimGatingSet(gs), gs)
  for (input in list(fs, cs)) {
    converted <- .asStimGatingSet(input)
    testthat::expect_s4_class(converted, "GatingSet")
    testthat::expect_length(converted, length(gs))
    testthat::expect_identical(
      flowWorkspace::sampleNames(converted), flowWorkspace::sampleNames(gs)
    )
    for (i in seq_along(gs)) {
      testthat::expect_identical(
        flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]])),
        matrices[[i]]
      )
    }
  }
  for (input in list(frames[[1]], cs[[1, returnType = "cytoframe"]])) {
    converted <- .asStimGatingSet(input)
    testthat::expect_s4_class(converted, "GatingSet")
    testthat::expect_length(converted, 1L)
    testthat::expect_identical(flowWorkspace::sampleNames(converted), "sample1")
    testthat::expect_identical(
      flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[1]])),
      matrices[[1]]
    )
  }
})

testthat::test_that("matrix and data-frame lists align channels and name samples", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  matrices <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
  })
  reordered <- matrices
  reordered[[2]] <- reordered[[2]][, rev(colnames(reordered[[2]])), drop = FALSE]
  for (input in list(matrices, lapply(matrices, as.data.frame), reordered)) {
    converted <- .asStimGatingSet(input)
    testthat::expect_s4_class(converted, "GatingSet")
    testthat::expect_length(converted, length(gs))
    testthat::expect_identical(
      flowWorkspace::sampleNames(converted), paste0("sample", seq_along(gs))
    )
    # flowCore stores `desc` as AsIs; compare the labels themselves.
    testthat::expect_identical(
      unclass(chnlLab(converted)),
      stats::setNames(colnames(matrices[[1]]), colnames(matrices[[1]]))
    )
    for (i in seq_along(gs)) {
      testthat::expect_identical(
        flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]])),
        matrices[[i]]
      )
    }
  }
  named <- stats::setNames(matrices, flowWorkspace::sampleNames(gs))
  converted <- .asStimGatingSet(named)
  testthat::expect_identical(flowWorkspace::sampleNames(converted), names(named))
  for (i in seq_along(gs)) {
    testthat::expect_identical(
      flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]])), matrices[[i]]
    )
  }

  # Missing descriptions already fall back to channel names in chnlLab().
  fr <- flowCore::flowFrame(matrices[[1]])
  params <- flowCore::parameters(fr)
  params@data$desc <- NA_character_
  flowCore::parameters(fr) <- params
  testthat::expect_identical(
    chnlLab(fr), stats::setNames(colnames(matrices[[1]]), colnames(matrices[[1]]))
  )
})

testthat::test_that("long data frames use first appearance or observed factor levels", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  matrices <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
  })
  sampleNames <- paste0("tube", seq_along(gs))
  appearance <- rev(seq_along(gs))
  long <- do.call(rbind, lapply(appearance, function(i) {
    data.frame(sample = sampleNames[[i]], matrices[[i]], check.names = FALSE)
  }))
  for (factorSample in c(FALSE, TRUE)) {
    input <- long
    sampleOrder <- appearance
    if (factorSample) {
      input$sample <- factor(input$sample, levels = c(sampleNames, "unused"))
      sampleOrder <- seq_along(gs)
    }
    converted <- .asStimGatingSet(input)
    testthat::expect_s4_class(converted, "GatingSet")
    testthat::expect_length(converted, length(gs))
    testthat::expect_identical(
      flowWorkspace::sampleNames(converted), sampleNames[sampleOrder]
    )
    for (i in seq_along(sampleOrder)) {
      # Long-format row names are presentation metadata, not expression values.
      actual <- flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]]))
      rownames(actual) <- NULL
      expected <- matrices[[sampleOrder[[i]]]]
      rownames(expected) <- NULL
      testthat::expect_identical(actual, expected)
    }
  }
})

testthat::test_that("FCS directories sort filenames and explicit paths preserve order", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathFcs <- tempfile("stimgate_input_fcs_")
  dir.create(pathFcs)
  withr::defer(unlink(pathFcs, recursive = TRUE))
  pathFcs <- normalizePath(pathFcs, winslash = "/")
  filenames <- paste0("tube", rev(seq_along(gs)), c(".fcs", ".FCS"))
  files <- file.path(pathFcs, filenames)
  matrices <- lapply(seq_along(gs), function(i) {
    fr <- flowWorkspace::gh_pop_get_data(gs[[i]])
    flowCore::write.FCS(fr, files[[i]], what = "double")
    flowCore::exprs(fr)
  })
  # Ignore non-FCS files and FCS files in subdirectories.
  file.create(file.path(pathFcs, "ignore.txt"))
  dir.create(file.path(pathFcs, "nested"))
  file.copy(files[[1]], file.path(pathFcs, "nested", "hidden.fcs"))
  for (directory in c(FALSE, TRUE)) {
    converted <- .asStimGatingSet(if (directory) pathFcs else files)
    sampleOrder <- if (directory) order(files) else seq_along(files)
    testthat::expect_s4_class(converted, "GatingSet")
    testthat::expect_length(converted, length(gs))
    testthat::expect_identical(
      flowWorkspace::sampleNames(converted), filenames[sampleOrder]
    )
    for (i in seq_along(sampleOrder)) {
      testthat::expect_identical(
        flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]])),
        matrices[[sampleOrder[[i]]]]
      )
    }
  }
  converted <- withr::with_dir(pathFcs, .asStimGatingSet(c(filenames[[1]], files[-1])))
  testthat::expect_identical(flowWorkspace::sampleNames(converted), filenames)
  for (i in seq_along(gs)) {
    testthat::expect_identical(
      flowCore::exprs(flowWorkspace::gh_pop_get_data(converted[[i]])), matrices[[i]]
    )
  }
})

testthat::test_that("invalid inputs and non-root populations fail clearly", {
  m <- matrix(seq_len(8), ncol = 2, dimnames = list(NULL, c("A", "B")))
  missingFile <- tempfile("missing_input_")
  testthat::expect_error(
    .asStimGatingSet(missingFile), basename(missingFile), fixed = TRUE
  )
  emptyDir <- tempfile("empty_input_")
  dir.create(emptyDir)
  withr::defer(unlink(emptyDir, recursive = TRUE))
  testthat::expect_error(.asStimGatingSet(emptyDir), "No FCS files")
  testthat::expect_error(.asStimGatingSet(list(m, m[, 1, drop = FALSE])), "mismatched")
  testthat::expect_error(.asStimGatingSet(list(data.frame(A = 1, B = "x"))), "numeric")
  testthat::expect_error(
    .asStimGatingSet(list(matrix("x", dimnames = list(NULL, "A")))), "numeric"
  )
  testthat::expect_error(.asStimGatingSet(42), "GatingSet, flowSet, cytoset")
  testthat::expect_error(.asStimGatingSet(data.frame(A = 1)), "sample.*column")
  testthat::expect_error(.asStimGatingSet(list()), "at least one sample")
  testthat::expect_error(.asStimGatingSet(list(m), "CD3"), 'Only the "root"')
  testthat::expect_error(.asStimGatingSet(list(m), c("root", "CD3")), 'Only the "root"')
  testthat::expect_error(.asStimGatingSet(list(unname(m))), "column names")
  testthat::expect_error(
    .asStimGatingSet(stats::setNames(list(m, m), c("same", "same"))), "Sample names"
  )
  testthat::expect_error(
    .asStimGatingSet(data.frame(sample = NA_character_, A = 1)), "missing values"
  )

  pathProject <- tempfile("invalid_input_project_")
  withr::defer(unlink(pathProject, recursive = TRUE))
  testthat::expect_error(
    gateStim(pathProject, list(m, m), list(c(1L, 2L)), chnl = "A", popGate = "CD3"),
    'Only the "root"'
  )
  testthat::expect_error(
    gateStim(
      pathProject, list(m, m), list(c(1L, 2L)), chnl = "A",
      markerControl = list(A = list(popGate = "CD3"))
    ),
    'Only the "root"'
  )
  testthat::expect_error(
    gateStim(pathProject, list(m, m), list(c("sample1", "absent")), chnl = "A"),
    "Unknown sample name.*absent"
  )
  testthat::expect_false(dir.exists(file.path(pathProject, "metaData")))
  testthat::expect_error(
    plotStim(1, list(m), pathProject, chnl = "A", pop = "CD3"), 'Only the "root"'
  )
  testthat::expect_error(
    getStimExpr(pathProject, .data = list(m), pop = "CD3", ind = "1", chnl = "A"),
    'Only the "root"'
  )
  testthat::expect_error(
    writeStimFCS(
      pathProject, list(m), pop = "CD3", indBatchList = list(1L),
      pathDirSave = file.path(pathProject, "fcs")
    ),
    'Only the "root"'
  )
})

testthat::test_that("matrix gating and named batches preserve GatingSet thresholds", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  matrices <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
  })
  matrices <- stats::setNames(matrices, flowWorkspace::sampleNames(gs))
  paths <- file.path(
    dirname(exampleData$pathGs), c("gs-results", "matrices", "names")
  )
  namedBatches <- lapply(exampleData$batchList, function(x) names(matrices)[x])
  inputs <- list(gs, matrices, matrices)
  for (i in seq_along(paths)) {
    withr::with_seed(42, gateStim(
      pathProject = paths[[i]], .data = inputs[[i]],
      batchList = if (i == 3L) namedBatches else exampleData$batchList,
      chnl = exampleData$chnl
    ))
  }
  reference <- getStimGates(paths[[1]])
  testthat::expect_gt(nrow(reference), 0L)
  keys <- c("chnl", "batch", "ind", "gateName")
  reference <- dplyr::arrange(reference, !!!rlang::syms(keys))
  for (path in paths[-1]) {
    actual <- dplyr::arrange(getStimGates(path), !!!rlang::syms(keys))
    testthat::expect_equal(actual[, c(keys, "gate")], reference[, c(keys, "gate")])
  }
  savedBatches <- stimgateMetaReadBatchList(paths[[3]])
  expectedBatches <- lapply(namedBatches, function(x) match(x, names(matrices)))
  if (is.null(names(expectedBatches))) {
    names(expectedBatches) <- paste0("batch", seq_along(expectedBatches))
  }
  testthat::expect_identical(savedBatches, expectedBatches)
})

testthat::test_that("getStimExpr reads matrix input on a cache miss", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  matrices <- lapply(seq_along(gs), function(i) {
    flowCore::exprs(flowWorkspace::gh_pop_get_data(gs[[i]]))
  })
  pathProject <- file.path(dirname(exampleData$pathGs), "expression")
  actual <- getStimExpr(
    pathProject, .data = matrices, pop = "root", ind = "2", chnl = exampleData$chnl
  )
  testthat::expect_equal(
    as.matrix(actual[, exampleData$chnl]), matrices[[2]][, exampleData$chnl, drop = FALSE]
  )
  cached <- getStimExpr(pathProject, pop = "root", ind = "2", chnl = exampleData$chnl)
  testthat::expect_identical(actual, cached)
})
