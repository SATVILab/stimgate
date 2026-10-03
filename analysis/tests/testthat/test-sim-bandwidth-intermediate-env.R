test_that("bandwidth frequency simulations restore intermediate saving after success and errors", {
  withr::local_seed(42)
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  source(file.path(root, "scripts", "r", "sim-bandwidth.R"), local = env)
  env$.simMiscGetTrans <- function(...) identity

  frame <- flowCore::flowFrame(matrix(as.numeric(1:4), ncol = 1, dimnames = list(NULL, "F1")))
  experiment <- list(
    flowFrameList = list(frame, frame),
    labelsList = list(rep("gn", 4), c("gn", "gn", "gp", "gp"))
  )
  failure <- "none"
  gating_calls <- 0L
  env$gateStim <- function(pathProject, control, ...) {
    expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE"), "TRUE")
    expect_s3_class(control, "stimControl")
    expect_identical(control$clusterGates, FALSE)
    expect_identical(control$bw, 0.1)
    expect_true(all(names(list(...)) %in% names(formals(stimgate::gateStim))))
    gating_calls <<- gating_calls + 1L
    if (failure == "gate") stop("mock gating failure")
    saveRDS(tibble::tibble(), file.path(pathProject, "gateStats.rds"))
    dir.create(file.path(pathProject, "intermediateData", "init", "F1"), recursive = TRUE)
  }
  env$.simBandwidthReadLocDetails <- function(...) {
    expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE"), "TRUE")
    if (failure == "details") stop("mock detail failure")
    tibble::tibble(
      sample = "1", ind = "2", chnl = "F1", method = "loc_sample",
      propRespEst = 0.25, detailLevel = "sample"
    )
  }
  run <- function() {
    env$.simBandwidthBsFreq(
      nSample = 1L, nMarker = 1L, nCondition = 2L, nCluster = 2L,
      nIter = 2L, biasUns = 0, bw = 0.1, nCellStim = 4L,
      probResponse = 0.5, meanPos = 8, transformation = "gaussian",
      samplePerturbationSd = 0, conditionPerturbationSd = 0,
      clusterPerturbationSd = 0, backgroundRelativeToResponse = 0,
      ncellUnsRelativeToStim = 1
    )
  }

  testthat::with_mocked_bindings(
    simCytExperiment = function(...) experiment,
    .package = "simcyto",
    testthat::with_mocked_bindings(
      GatingSet = function(fs) fs,
      .package = "flowWorkspace",
      {
        for (initial in c(NA_character_, "FALSE", "TRUE", "", "caller-setting")) {
          withr::with_envvar(c(STIMGATE_INTERMEDIATE = initial), {
            # Windows treats an empty environment value as unset.
            before <- Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_)
            failure <- "none"
            gating_calls <- 0L
            result <- run()
            expect_identical(gating_calls, 2L)
            expect_equal(result$propRespTruth, rep(0.5, 6))
            expect_equal(result$propRespEst[result$method == "loc_sample"], c(0.25, 0.25))
            expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_), before)

            for (failure in c("gate", "details")) {
              expect_error(run(), paste("mock", if (failure == "gate") "gating" else "detail", "failure"))
              expect_identical(Sys.getenv("STIMGATE_INTERMEDIATE", unset = NA_character_), before)
            }
          })
        }
      }
    )
  )
})
