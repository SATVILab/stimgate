test_that("gateStim validates parallel as a logical scalar", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "validation")
  for (parallel in list("true", 1, NA, logical(), c(TRUE, FALSE))) {
    expect_error(gateStim(
      .data = gs,
      pathProject = pathProject,
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      parallel = parallel
    ), "`parallel` must be TRUE or FALSE", fixed = TRUE)
  }
})

test_that("parallel FALSE preserves lapply RNG and never creates futures", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  withr::local_preserve_seed()
  oldPlan <- future::plan(future::multisession, workers = 2)
  withr::defer(future::plan(oldPlan))
  testthat::local_mocked_bindings(
    future_lapply = function(...) stop("unexpected future"),
    .package = "future.apply"
  )
  testthat::local_mocked_bindings(
    .gateCacheChnl = function(...) stop("unexpected cache warmup"),
    .gateInitChnl = function(...) sample.int(100000L, 1L)
  )
  settings <- list(first = list(), second = list())
  set.seed(42)
  expected <- lapply(settings, function(...) sample.int(100000L, 1L))
  expectedSeed <- .Random.seed
  set.seed(42)
  actual <- .gateMapChnl(settings, NULL, list(), tempdir(), parallel = FALSE)
  expect_identical(actual, expected)
  expect_identical(.Random.seed, expectedSeed)
})

test_that("one channel uses lapply even when parallel is requested", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  oldPlan <- future::plan(future::multisession, workers = 2)
  withr::defer(future::plan(oldPlan))
  testthat::local_mocked_bindings(
    future_lapply = function(...) stop("unexpected future"),
    .package = "future.apply"
  )
  testthat::local_mocked_bindings(
    .gateCacheChnl = function(...) invisible(NULL),
    .gateInitChnl = function(...) "sequential"
  )
  expect_identical(
    .gateMapChnl(
      list(marker = list()), NULL, list(), tempdir(),
      parallel = TRUE
    ),
    list(marker = "sequential")
  )
})

test_that("cache warmup reads each sample once and supports NULL GatingSets", {
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  pathProject <- file.path(dirname(exampleData$pathGs), "cache")
  settings <- lapply(exampleData$chnl, function(x) {
    list(chnlCut = x, popGate = "root")
  })
  nReads <- 0L
  getData <- flowWorkspace::gh_pop_get_data
  testthat::local_mocked_bindings(
    gh_pop_get_data = function(...) {
      nReads <<- nReads + 1L
      getData(...)
    },
    .package = "flowWorkspace"
  )
  .gateCacheChnl(gs, exampleData$batchList, settings, pathProject)
  expect_equal(nReads, length(unique(unlist(exampleData$batchList))))
  .gateCacheChnl(NULL, exampleData$batchList, settings, pathProject)
  expect_equal(nReads, length(unique(unlist(exampleData$batchList))))
  indBatch <- exampleData$batchList[[1L]]
  for (chnl in exampleData$chnl) {
    cached <- .getExList(
      NULL, indBatch, "batch1", "root", chnl,
      pathProject = pathProject
    )
    withData <- .getExList(
      gs, indBatch, "batch1", "root", chnl,
      pathProject = pathProject
    )
    expect_identical(cached, withData)
  }
  unlink(.getExChnlPath(
    exampleData$chnl[[1L]], indBatch[[1L]], "root", pathProject
  ))
  expect_error(
    .getExList(
      NULL, indBatch, "batch1", "root", exampleData$chnl[[1L]],
      pathProject = pathProject
    ),
    "Incomplete expression cache"
  )
})

test_that("multisession gates are close to sequential and reproducible", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  # Multisession workers load the installed stimgate, not a load_all() copy,
  # so run this only against an installed package (R CMD check).
  skip_if(exists(
    ".__DEVTOOLS__",
    envir = asNamespace("stimgate"), inherits = FALSE
  ))
  withr::local_preserve_seed()
  withr::local_envvar(
    c(STIMGATE_DEBUG = "true", STIMGATE_INTERMEDIATE = "all")
  )
  oldPlan <- future::plan(future::multisession, workers = 2)
  withr::defer(future::plan(oldPlan))
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  paths <- file.path(
    dirname(exampleData$pathGs), c("sequential", "parallel", "repeat")
  )

  # Reject attempts to export an external pointer, including the GatingSet.
  withr::local_options(future.globals.onReference = "error")
  for (i in seq_along(paths)) {
    set.seed(42)
    gateStim(
      .data = gs,
      pathProject = paths[[i]],
      popGate = "root",
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      parallel = i > 1L
    )
  }
  sequential <- getStimGates(paths[[1L]])
  parallel <- getStimGates(paths[[2L]])
  repeated <- getStimGates(paths[[3L]])
  expect_gt(nrow(parallel), 0L)
  expect_identical(names(parallel), names(sequential))
  keys <- c("pop", "marker", "chnl", "batch", "ind", "gateName")
  expect_equal(parallel[, keys], sequential[, keys])
  expect_true(all(is.finite(sequential$gate)))
  expect_true(all(is.finite(parallel$gate)))
  # Random per-bin thinning can slightly move the fitted threshold.
  expect_lt(max(abs(parallel$gate - sequential$gate)), 0.5)
  expect_equal(parallel, repeated)

  profile <- readRDS(file.path(paths[[2L]], "profile", "profile.rds"))
  workerRows <- profile[profile$operation == "marker_total", , drop = FALSE]
  expect_equal(nrow(workerRows), length(exampleData$chnl))
  expect_true(all(workerRows$pid != Sys.getpid()))
  expect_true(all(workerRows$status == "completed"))
  expect_true(file.exists(file.path(paths[[2L]], "debug", "debug.txt")))
  for (chnl in exampleData$chnl) {
    expect_true(length(list.files(
      file.path(paths[[2L]], "intermediateData", "init", chnl),
      recursive = TRUE, pattern = "\\.rds$"
    )) > 0L)
  }
})
