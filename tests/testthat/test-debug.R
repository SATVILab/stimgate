test_that("debug mode is inactive when STIMGATE_DEBUG is disabled", {
  withr::local_envvar(c(STIMGATE_DEBUG = "false"))
  .debugStateReset()
  pathProject <- tempfile("stimgate-debug-off-")
  dir.create(pathProject, recursive = TRUE)
  withr::defer({
    .debugStateReset()
    unlink(pathProject, recursive = TRUE, force = TRUE)
  })

  expect_false(.debugInit(pathProject))
  expect_false(dir.exists(file.path(pathProject, "debug")))
  expect_false(.debug("test message"))
})

test_that(
  "debugInit creates pathProject/debug/debug.txt and resets old debug dir",
  {
    withr::local_envvar(c(STIMGATE_DEBUG = "true"))
    .debugStateReset()
    pathProject <- tempfile("stimgate-debug-init-")
    withr::defer({
      .debugStateReset()
      unlink(pathProject, recursive = TRUE, force = TRUE)
    })

    dirDebugOld <- file.path(pathProject, "debug")
    dir.create(dirDebugOld, recursive = TRUE)
    pathOldFile <- file.path(dirDebugOld, "old.txt")
    writeLines("stale content", pathOldFile)

    expect_true(.debugInit(pathProject, reset = TRUE))
    expect_false(file.exists(pathOldFile))

    pathDebugFile <- file.path(pathProject, "debug", "debug.txt")
    expect_true(file.exists(pathDebugFile))
    expect_identical(.debugState$file, pathDebugFile)
    expect_true(.debugState$initialized)

    expect_true(.debug("First line", 42))
    expect_true(.debug("Second line"))
    lines <- readLines(pathDebugFile)
    expect_equal(lines, c("First line: 42", "Second line"))

    .debugStateReset()
    expect_null(.debugState$file)
    expect_false(.debugState$initialized)
  }
)

test_that(
  "gateStim in debug mode creates debug.txt, clears old run, tempdir clean",
  {
    skip_if_not_installed("flowWorkspace")
    skip_if_not_installed("flowCore")

    withr::local_envvar(c(STIMGATE_DEBUG = "true"))
    .debugStateReset()
    .profileStateReset()

    exampleData <- getExampleData()
    gs <- flowWorkspace::load_gs(exampleData$pathGs)
    pathProject <- file.path(tempdir(), "test_debug_gatestim_run")

    withr::defer({
      .debugStateReset()
      .profileStateReset()
      if (dir.exists(pathProject)) {
        unlink(pathProject, recursive = TRUE, force = TRUE)
      }
      if (dir.exists(exampleData$pathGs)) {
        unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
      }
    })

    # Record tempdir contents before run
    tempFilesBefore <- list.files(tempdir(), full.names = TRUE)

    # 1. Run gateStim with STIMGATE_DEBUG enabled
    resPath <- gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      calcCytPosGates = FALSE,
      tolClust = NULL
    )
    expect_equal(resPath, pathProject)

    pathDebugFile <- file.path(pathProject, "debug", "debug.txt")
    expect_true(file.exists(pathDebugFile))
    debugLines1 <- readLines(pathDebugFile)
    expect_true(length(debugLines1) > 0L)

    # Verify no new stimgate directory or files created in tempdir
    expect_false(dir.exists(file.path(tempdir(), "stimgate")))

    # Debug state is cleaned up after gateStim exit
    expect_null(.debugState$file)
    expect_false(.debugState$initialized)

    # 2. Run gateStim a second time and confirm debug.txt was wiped & recreated
    resPath2 <- gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      calcCytPosGates = FALSE,
      tolClust = NULL
    )
    expect_equal(resPath2, pathProject)
    expect_true(file.exists(pathDebugFile))
    debugLines2 <- readLines(pathDebugFile)
    expect_true(length(debugLines2) > 0L)
  }
)

test_that("non-debug gateStim does not create or reset debug directory", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  withr::local_envvar(c(STIMGATE_DEBUG = "false"))
  .debugStateReset()

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "test_non_debug_gatestim_run")

  withr::defer({
    .debugStateReset()
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE, force = TRUE)
    }
    if (dir.exists(exampleData$pathGs)) {
      unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
    }
  })

  gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    calcCytPosGates = FALSE,
    tolClust = NULL
  )

  expect_false(dir.exists(file.path(pathProject, "debug")))
})

test_that("debug write errors never fail gateStim", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .debugStateReset()

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "test_debug_error_resilience")

  withr::defer({
    .debugStateReset()
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE, force = TRUE)
    }
    if (dir.exists(exampleData$pathGs)) {
      unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
    }
  })

  # Set an invalid unwritable path in .debugState to simulate write failures
  expect_false(.debug("test error message"))

  testthat::with_mocked_bindings(
    .debugInit = function(...) FALSE,
    {
      expect_no_error(
        resPath <- gateStim(
          .data = gs,
          pathProject = pathProject,
          popGate = "root",
          batchList = exampleData$batchList,
          marker = exampleData$marker,
          calcCytPosGates = FALSE,
          tolClust = NULL
        )
      )
      expect_equal(resPath, pathProject)
      expect_true(is.data.frame(getStimGates(pathProject)))
    }
  )
})

test_that("stimgate_debug_copy and debug_print are removed from exports", {
  exports <- getNamespaceExports("stimgate")
  expect_false("stimgate_debug_copy" %in% exports)
  expect_false("stimgate_debug_print" %in% exports)
})
