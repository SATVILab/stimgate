test_that("bwScope sets shared local-FDR bandwidths in the channel settings", {
  skip_if_not_installed("flowWorkspace")
  exampleData <- getExampleData()
  withr::defer(unlink(dirname(exampleData$pathGs), recursive = TRUE))
  gs <- flowWorkspace::load_gs(exampleData$pathGs)
  indAll <- as.character(unlist(exampleData$batchList))

  runSettings <- function(..., minCell = 1e2) {
    pathProject <- withr::local_tempdir(.local_envir = parent.frame())
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      batchList = exampleData$batchList,
      marker = exampleData$marker,
      control = stimControl(
        calcCytPosGates = FALSE, clusterGates = FALSE,
        minCell = minCell, ...
      )
    ))
    stimgateMetaReadSettingsChnls(pathProject)
  }

  # Default: one trimmed-mean bandwidth per channel
  for (settings in runSettings()) {
    expect_identical(settings$bwScope, "cytokine")
    expect_false(isTRUE(settings$clusterGates))
    expect_length(settings$bwShared, 1L)
    expect_true(is.finite(settings$bwShared) && settings$bwShared > 0)
    expect_null(settings$bwSharedTbl)
  }

  for (settings in runSettings(bwScope = "cluster")) {
    tbl <- settings$bwSharedTbl
    expect_setequal(tbl$ind, indAll)
    # Fewer than 100 tubes, so every tube's bandwidth is estimated
    expect_true(all(is.finite(tbl$bwEst)))
    bwGrp <- c(tapply(tbl$bwEst, tbl$grp, stats::median))
    expect_equal(tbl$bw, unname(bwGrp[tbl$grp]))
    expect_true(is.finite(settings$bwShared))
  }

  # Tubes below minCell are excluded, leaving only the fallback
  for (settings in runSettings(bwScope = "cluster", minCell = 1e7)) {
    expect_null(settings$bwSharedTbl)
    expect_identical(settings$bwShared, settings$bwFallback)
  }
  for (settings in runSettings(minCell = 1e7)) {
    expect_identical(settings$bwShared, settings$bwFallback)
  }

  for (settings in runSettings(bwScope = "sample")) {
    expect_null(settings$bwShared)
  }
  for (settings in runSettings(bwAdaptive = TRUE)) {
    expect_null(settings$bwShared)
  }
})

test_that("local-FDR uses the smaller shared bandwidth of the sample's tubes", {
  set.seed(1)
  exStim <- data.frame(x = stats::rnorm(200))
  exUns <- data.frame(x = stats::rnorm(200))
  attr(exStim, "chnlCut") <- attr(exUns, "chnlCut") <- "x"
  attr(exStim, "ind") <- "2"
  attr(exUns, "ind") <- "1"
  getBw <- function(chnlSettings) {
    .getCpUnsLocGetDensRawDensitiesBw(exStim, exUns, chnlSettings)
  }

  expect_identical(getBw(list(bwShared = 0.3)), 0.3)

  tbl <- tibble::tibble(ind = c("1", "2"), grp = c("1", "2"), bw = c(0.2, 0.4))
  expect_identical(getBw(list(bwShared = 0.3, bwSharedTbl = tbl)), 0.2)

  # A tube missing from the cluster table (e.g. a prejoined sample) uses the
  # channel-level shared bandwidth.
  tbl$ind[1] <- "9"
  expect_identical(getBw(list(bwShared = 0.3, bwSharedTbl = tbl)), 0.3)

  # A fixed bandwidth takes precedence
  expect_identical(getBw(list(bw = 0.1, bwShared = 0.3)), 0.1)

  # Without a shared bandwidth (bwScope = "sample"), each sample's bandwidth is
  # still estimated from its own stim and unstim tubes.
  exUns$x <- exUns$x * 3
  chnlSettings <- list(
    bwMtd = "hpi1", bwAdj = 1, bwMin = 0, bwMax = Inf, bwFallback = 99
  )
  bwStim <- .getCpUnsLocGetDensRawDensitiesBwInit(exStim$x, chnlSettings)
  bwUns <- .getCpUnsLocGetDensRawDensitiesBwInit(exUns$x, chnlSettings)
  expect_lt(bwStim, bwUns)
  expect_identical(getBw(chnlSettings), bwStim)
})

test_that("tube clustering separates differently shaped backgrounds", {
  set.seed(1)
  xList <- c(
    lapply(1:4, function(i) stats::rnorm(2000, sd = 0.3)),
    lapply(1:4, function(i) stats::rnorm(2000, sd = 1.5))
  ) |>
    stats::setNames(as.character(1:8))

  clusterTbl <- .bwSharedCluster(xList, bw = 0.2)

  expect_identical(clusterTbl$ind, names(xList))
  expect_length(unique(clusterTbl$grp[1:4]), 1L)
  expect_length(unique(clusterTbl$grp[5:8]), 1L)
  expect_false(clusterTbl$grp[1] == clusterTbl$grp[5])

  expect_identical(nrow(.bwSharedCluster(list(a = 1, b = 1:2), bw = 0.2)), 0L)
})

test_that("threshold sharing uses bwCluster, then the shared bandwidth", {
  getBw <- function(chnlSettings) {
    .getCpClusterLocCommonBw(character(), list(), chnlSettings)
  }
  expect_identical(getBw(list(bwCluster = 0.2, bwShared = 0.3)), 0.2)
  expect_identical(getBw(list(bwShared = 0.3, bw = 0.5)), 0.3)
  expect_identical(getBw(list(bw = 0.5)), 0.5)
})

test_that("shared bandwidths prefer large tubes, then the highest size bands", {
  withr::local_seed(1)
  nCell <- c(a = 12000, b = 9500, c = 9100, d = 8200, e = 5500, f = 4000)
  # Enough tubes with at least bwNcellMax cells: only those are used
  expect_identical(.bwSharedSelect(c(nCell, g = 2e4), 2L, 1e4), c("a", "g"))
  # Otherwise the shortfall comes from the 9-10k band before 8-9k
  sel <- .bwSharedSelect(nCell, 3L, 1e4)
  expect_setequal(sel, c("a", "b", "c"))
  expect_setequal(.bwSharedSelect(nCell, 4L, 1e4), c("a", "b", "c", "d"))
  # Bands stop at half of bwNcellMax once some tube has been found
  expect_setequal(.bwSharedSelect(nCell, 10L, 1e4), c("a", "b", "c", "d", "e"))
  # With no tube from half of bwNcellMax upwards, smaller tubes are used
  expect_identical(.bwSharedSelect(c(f = 4000, h = 300), 1L, 1e4), "f")
  expect_setequal(.bwSharedSelect(c(f = 4000, h = 300), 5L, 1e4), c("f", "h"))
  # Within a band the choice is random
  band <- stats::setNames(rep(9500, 20), letters[1:20])
  draws <- replicate(20, .bwSharedSelect(band, 1L, 1e4))
  expect_gt(length(unique(draws)), 1L)
  # Without bwNcellMax, tubes are spread across the batches
  expect_identical(.bwSharedSelect(nCell, 2L, NULL), c("a", "f"))
})

test_that("shared bandwidths are estimated on bwNcellMax cells", {
  calls <- list()
  testthat::local_mocked_bindings(
    .getCpUnsLocGetDensRawDensitiesBwInit = function(.data, chnlSettings) {
      calls[[length(calls) + 1L]] <<- chnlSettings[c("bwNcellMin", "bwNcellMax")]
      0.3
    }
  )
  exList <- list(t1 = stats::rnorm(6000), t2 = stats::rnorm(12000))
  testthat::local_mocked_bindings(
    .bwSharedReadBatch = function(i, ...) exList[i]
  )
  settings <- .completeChnlSettingsBwShared(
    list(bwScope = "cytokine", bw = NULL, minCell = 100, bwNcellMin = 100,
      bwNcellMax = 1e4, bwFallback = 0.5, marker = "m"),
    indBatchList = list(b1 = "t1", b2 = "t2"), .data = NULL, pathProject = NULL
  )
  expect_identical(settings$bwShared, 0.3)
  # Both tubes are selected (6,000 cells is in the 6-7k band) and each is
  # upsampled or downsampled to bwNcellMax cells
  expect_length(calls, 2L)
  for (call in calls) {
    expect_identical(call$bwNcellMin, 1e4)
    expect_identical(call$bwNcellMax, 1e4)
  }
})
