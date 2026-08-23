test_that("profiling is disabled outside debug mode", {
  withr::local_envvar(c(STIMGATE_DEBUG = "false"))
  .profileStateReset()
  pathProject <- tempfile("stimgate-profile-off-")
  dir.create(pathProject, recursive = TRUE)
  withr::defer({
    .profileStateReset()
    unlink(pathProject, recursive = TRUE, force = TRUE)
  })

  expect_false(.profileInit(pathProject))
  expect_false(dir.exists(file.path(pathProject, "profile")))
})

test_that("a debug run resets and incrementally saves profiling records", {
  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()
  pathProject <- tempfile("stimgate-profile-")
  withr::defer({
    .profileStateReset()
    unlink(pathProject, recursive = TRUE, force = TRUE)
  })
  pathOldRaw <- file.path(pathProject, "profile", "raw")
  dir.create(pathOldRaw, recursive = TRUE)
  pathOld <- file.path(pathOldRaw, "old.rds")
  saveRDS(data.frame(old = TRUE), pathOld)

  expect_true(.profileInit(pathProject))
  expect_false(file.exists(pathOld))

  timer <- .profileStart(
    level = "major",
    major = "initial_gating",
    operation = "test_major",
    pathProject = pathProject
  )
  expect_false(is.null(timer))
  .profileStop(timer)

  rawFiles <- list.files(
    file.path(pathProject, "profile", "raw"),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  expect_length(rawFiles, 1L)
  raw <- readRDS(rawFiles[[1L]])
  expect_equal(raw$operation, "test_major")
  expect_equal(raw$status, "completed")
  expect_true(is.numeric(raw$elapsed_sec))

  .profileFinishRun(pathProject)
  expect_true(file.exists(file.path(pathProject, "profile", "profile.rds")))
  expect_true(file.exists(file.path(pathProject, "profile", "profile.csv")))
  expect_true(dir.exists(file.path(pathProject, "profile", "raw")))
})

test_that("profiling records explicit hierarchy and sample context", {
  withr::local_envvar(c(STIMGATE_DEBUG = "yes"))
  .profileStateReset()
  pathProject <- tempfile("stimgate-profile-hierarchy-")
  withr::defer({
    .profileStateReset()
    unlink(pathProject, recursive = TRUE, force = TRUE)
  })
  expect_true(.profileInit(pathProject))

  .profileWithContext(
    .profileTime(
      .profileTime(
        1 + 1,
        level = "sample_detail",
        major = "initial_gating",
        minor = "local_fdr",
        operation = "marginal_filtering",
        pathProject = pathProject
      ),
      level = "sample",
      major = "initial_gating",
      minor = "local_fdr",
      operation = "sample_initial_gating",
      pathProject = pathProject
    ),
    marker = "IL2",
    channel = "IL2-A",
    batch = "participant_17",
    sample = "104",
    stage = "init"
  )

  profileTbl <- .profileFinishRun(pathProject)
  sampleRow <- profileTbl[
    profileTbl$operation == "sample_initial_gating", ,
    drop = FALSE
  ]
  detailRow <- profileTbl[
    profileTbl$operation == "marginal_filtering", ,
    drop = FALSE
  ]
  runRow <- profileTbl[profileTbl$operation == "gateStim", , drop = FALSE]

  expect_equal(nrow(sampleRow), 1L)
  expect_equal(nrow(detailRow), 1L)
  expect_equal(nrow(runRow), 1L)
  expect_equal(sampleRow$marker, "IL2")
  expect_equal(sampleRow$channel, "IL2-A")
  expect_equal(sampleRow$batch, "participant_17")
  expect_equal(sampleRow$sample, "104")
  expect_equal(sampleRow$stage, "init")
  expect_equal(detailRow$parent_id, sampleRow$record_id)
  expect_equal(sampleRow$parent_id, runRow$record_id)
  expect_gt(detailRow$depth, sampleRow$depth)
})

test_that("profiling instrumentation preserves wrapped function arguments", {
  expect_identical(
    names(formals(.gateInit)),
    names(formals(.profileOriginalGateInit))
  )
  expect_identical(
    names(formals(.gateCytPos)),
    names(formals(.profileOriginalGateCytPos))
  )
  expect_identical(
    names(formals(.gateStats)),
    names(formals(.profileOriginalGateStats))
  )
  expect_identical(
    names(formals(.gateChnl)),
    names(formals(.profileOriginalGateChnl))
  )
  expect_identical(
    names(formals(.gateBatch)),
    names(formals(.profileOriginalGateBatch))
  )
  expect_identical(
    names(formals(.getCpUnsLocCondition)),
    names(formals(.profileOriginalGetCpUnsLocCondition))
  )
  expect_identical(
    names(formals(.getCpUnsLocGetProb)),
    names(formals(.profileOriginalGetCpUnsLocGetProb))
  )
  expect_identical(
    names(formals(.getCpUnsLocGetCp)),
    names(formals(.profileOriginalGetCpUnsLocGetCp))
  )
})

test_that("gateStim with STIMGATE_DEBUG produces full profiling records and cleans up", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "test_profile_gatestim_full")

  withr::defer({
    .profileStateReset()
    if (dir.exists(pathProject)) unlink(pathProject, recursive = TRUE, force = TRUE)
    if (dir.exists(exampleData$pathGs)) unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
  })

  # 1. Run gateStim in debug mode
  resPath <- gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    calcCytPosGates = TRUE,
    tolClust = NULL
  )

  expect_equal(resPath, pathProject)

  # Verify standard gating outputs exist
  gates <- getStimGates(pathProject)
  expect_true(is.data.frame(gates) && nrow(gates) > 0L)
  stats <- getStimStats(pathProject)
  expect_true(is.data.frame(stats) && nrow(stats) > 0L)

  # Verify profile outputs exist
  pathProfileRds <- file.path(pathProject, "profile", "profile.rds")
  pathProfileCsv <- file.path(pathProject, "profile", "profile.csv")
  expect_true(file.exists(pathProfileRds))
  expect_true(file.exists(pathProfileCsv))

  profileTbl <- readRDS(pathProfileRds)
  expect_true(is.data.frame(profileTbl) && nrow(profileTbl) > 0L)

  # Check root record
  rootRow <- profileTbl[profileTbl$level == "run", , drop = FALSE]
  expect_equal(nrow(rootRow), 1L)
  expect_equal(rootRow$operation, "gateStim")
  expect_equal(rootRow$depth, 0L)
  expect_true(is.na(rootRow$parent_id))
  expect_equal(rootRow$status, "completed")
  expect_true(rootRow$elapsed_sec > 0)

  # Check expected major stages
  majorOps <- unique(profileTbl$operation[profileTbl$level == "major"])
  expect_true("initial_gating" %in% majorOps)
  expect_true("cytokine_positive_gating" %in% majorOps)
  expect_true("final_statistics" %in% majorOps)

  # Check marker, batch and sample context on sample/detail records
  sampleRows <- profileTbl[profileTbl$level == "sample", , drop = FALSE]
  expect_true(nrow(sampleRows) > 0L)
  expect_true(all(!is.na(sampleRows$marker)))
  expect_true(all(!is.na(sampleRows$batch)))
  expect_true(all(!is.na(sampleRows$sample)))
  expect_true(all(!is.na(sampleRows$stage)))

  # Check that each non-root record parent chain reaches root gateStim record
  nonRootRows <- profileTbl[profileTbl$level != "run", , drop = FALSE]
  expect_true(nrow(nonRootRows) > 0L)
  expect_true(all(!is.na(nonRootRows$parent_id)))
  expect_true(all(nonRootRows$parent_id %in% profileTbl$record_id))

  # Trace each record to root
  rootId <- rootRow$record_id[[1L]]
  parentMap <- stats::setNames(profileTbl$parent_id, profileTbl$record_id)
  for (recId in nonRootRows$record_id) {
    curr <- recId
    visited <- character()
    while (!is.na(curr) && curr != rootId) {
      expect_false(curr %in% visited) # no cycles
      visited <- c(visited, curr)
      curr <- parentMap[[curr]]
    }
    expect_equal(curr, rootId)
  }

  # In-memory state is clean
  expect_false(.profileState$initialized)
  expect_null(.profileState$runTimer)

  # 2. Run a second time and confirm records from previous run were wiped
  firstRunRecordIds <- profileTbl$record_id
  resPath2 <- gateStim(
    .data = gs,
    pathProject = pathProject,
    popGate = "root",
    batchList = exampleData$batchList,
    marker = exampleData$marker,
    calcCytPosGates = FALSE,
    tolClust = NULL
  )
  profileTbl2 <- readRDS(pathProfileRds)
  rootRow2 <- profileTbl2[profileTbl2$level == "run", , drop = FALSE]
  expect_equal(nrow(rootRow2), 1L)
  # No intersection with previous run IDs (fresh run)
  expect_equal(length(intersect(firstRunRecordIds, profileTbl2$record_id)), 0L)
})

test_that(
  "gateStim with STIMGATE_DEBUG handles failure and persists failed status",
  {
    skip_if_not_installed("flowWorkspace")
    skip_if_not_installed("flowCore")

    withr::local_envvar(c(STIMGATE_DEBUG = "true"))
    .profileStateReset()

    exampleData <- getExampleData()
    gs <- flowWorkspace::load_gs(exampleData$pathGs)
    pathProject <- file.path(tempdir(), "test_profile_failed_run")

    withr::defer({
      .profileStateReset()
      if (dir.exists(pathProject)) {
        unlink(pathProject, recursive = TRUE, force = TRUE)
      }
      if (dir.exists(exampleData$pathGs)) {
        unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
      }
    })

    # Validation failure (invalid batchList)
    expect_error(
      gateStim(
        .data = gs,
        pathProject = pathProject,
        popGate = "root",
        batchList = "not_a_list",
        marker = exampleData$marker
      )
    )

    # Profile directory and collated records exist despite failure
    pathProfileRds <- file.path(pathProject, "profile", "profile.rds")
    pathProfileCsv <- file.path(pathProject, "profile", "profile.csv")
    expect_true(file.exists(pathProfileRds))
    expect_true(file.exists(pathProfileCsv))

    profileTbl <- readRDS(pathProfileRds)
    rootRow <- profileTbl[profileTbl$level == "run", , drop = FALSE]
    expect_equal(nrow(rootRow), 1L)
    expect_equal(rootRow$operation, "gateStim")
    expect_equal(rootRow$status, "failed")
    expect_true(rootRow$elapsed_sec >= 0)

    # In-memory profiling state is clean
    expect_false(.profileState$initialized)
    expect_null(.profileState$runTimer)
  }
)

test_that("profiling errors never cause gateStim to fail", {
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")

  withr::local_envvar(c(STIMGATE_DEBUG = "true"))
  .profileStateReset()

  exampleData <- getExampleData()
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(tempdir(), "test_profile_io_error")

  withr::defer({
    .profileStateReset()
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE, force = TRUE)
    }
    if (dir.exists(exampleData$pathGs)) {
      unlink(exampleData$pathGs, recursive = TRUE, force = TRUE)
    }
  })

  # Mock .profileWriteRecord or .profileFinalise to throw an error
  testthat::with_mocked_bindings(
    .profileWriteRecord = function(record, pathRecord) {
      stop("Simulated profiling disk write error")
    },
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
