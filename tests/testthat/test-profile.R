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
    profileTbl$operation == "sample_initial_gating",
    ,
    drop = FALSE
  ]
  detailRow <- profileTbl[
    profileTbl$operation == "marginal_filtering",
    ,
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
