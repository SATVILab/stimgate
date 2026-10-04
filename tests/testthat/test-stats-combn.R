local({
  # Cached vectors exercise the production streaming path without a GatingSet.
  .statsCombnFixture <- function(exList, gates, gateType = "base") {
    project <- tempfile("stats-combn-")
    withr::defer(unlink(project, recursive = TRUE), envir = parent.frame())
    chnl <- names(exList[[1]])
    for (ind in seq_along(exList)) {
      for (channel in chnl) {
        path <- .getExChnlPath(channel, ind, "root", project)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(exList[[ind]][[channel]], path)
      }
    }
    .getStats(
      gateTbl = gates, chnl = chnl, chnlLab = stats::setNames(chnl, chnl),
      gateName = "g", gateTypeCytPosCalc = gateType, popGate = "root",
      .data = NULL, indBatchList = list(batch = seq_along(exList)),
      pathProject = project
    )
  }

  .statsCombnGates <- function(ind = 2L, chnl = c("A", "B", "C")) {
    tibble::tibble(
      ind = as.character(ind), gateName = "g", chnl = chnl,
      gate = 0.5, gateCyt = 0.25
    )
  }

  .statsCombnExpression <- function(code) {
    tibble::tibble(
      A = as.double(bitwAnd(code, 1L) > 0L),
      B = as.double(bitwAnd(code, 2L) > 0L),
      C = as.double(bitwAnd(code, 4L) > 0L)
    )
  }

  .statsCombnReference <- function(ex, gates, gateType) {
    chnl <- names(ex)
    matrices <- .getStatsCombnMatListGet(length(chnl))
    counts <- unlist(lapply(matrices, function(mat) {
      vapply(seq_len(nrow(mat)), function(i) {
        positive <- chnl[mat[i, , drop = TRUE]]
        as.integer(sum(.getPosIndCytCombn(
          ex = ex, gateTbl = gates, chnlPos = positive,
          chnlNeg = setdiff(chnl, positive), gateTypeCytPos = gateType
        )))
      }, integer(1))
    }), use.names = FALSE)
    c(counts, as.integer(nrow(ex) - sum(counts)))
  }

  test_that("streamed exact combinations retain counts, order, types and frequencies", {
    stim <- .statsCombnExpression(c(0:7, 7L))
    uns <- .statsCombnExpression(c(0L, 0L, 1L, 3L, 4L))
    actual <- .statsCombnFixture(list(uns, stim), .statsCombnGates())
    expect_identical(actual$countStim, c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L))
    expect_identical(actual$countUns, c(1L, 0L, 1L, 1L, 0L, 0L, 0L, 2L))
    expect_identical(actual$nCellStim, rep(9L, 8L))
    expect_identical(actual$nCellUns, rep(5L, 8L))
    expect_identical(sum(actual$countStim), 9L)
    expect_identical(sum(actual$countUns), 5L)
    expect_identical(actual$cytCombn, c(
      "A~+~B~-~C~-~", "A~-~B~+~C~-~", "A~-~B~-~C~+~",
      "A~+~B~+~C~-~", "A~+~B~-~C~+~", "A~-~B~+~C~+~",
      "A~+~B~+~C~+~", "A~-~B~-~C~-~"
    ))
    expect_identical(names(actual), c(
      "gateName", "ind", "cytCombn", "countStim", "nCellStim",
      "countUns", "nCellUns", "propStim", "propUns", "propBs",
      "freqStim", "freqUns", "freqBs"
    ))
    expect_equal(actual$propStim, actual$countStim / 9)
    expect_equal(actual$propUns, actual$countUns / 5)
    expect_equal(actual$freqBs, 100 * (actual$countStim / 9 - actual$countUns / 5))
    expect_lt(actual$freqBs[[1]], 0)
    expect_gt(actual$freqBs[[7]], 0)
  })

  test_that("shared raw unstim expression is classified with each stim sample's gates", {
    ex <- .statsCombnExpression(0:7)
    gates <- dplyr::bind_rows(
      .statsCombnGates(),
      dplyr::mutate(.statsCombnGates(3L), gate = 2, gateCyt = 2)
    )
    actual <- .statsCombnFixture(list(ex, ex, ex), gates)
    expect_identical(actual$countUns[actual$ind == "2"], rep(1L, 8L))
    expect_identical(actual$countUns[actual$ind == "3"], c(rep(0L, 7L), 8L))
    expect_identical(actual$countStim, actual$countUns)
  })

  test_that("missing channel gates are FALSE and missing sample gates yield integer NAs", {
    ex <- .statsCombnExpression(0:7)
    ex$C <- rep(NA_real_, nrow(ex))
    gates <- .statsCombnGates(chnl = c("A", "B"))
    actual <- .statsCombnFixture(list(ex, ex, ex), gates)
    expect_identical(actual$countStim[actual$ind == "2"], c(2L, 2L, 0L, 2L, 0L, 0L, 0L, 2L))
    expect_identical(actual$countUns[actual$ind == "2"], actual$countStim[actual$ind == "2"])
    expect_identical(actual$countStim[actual$ind == "3"], rep(NA_integer_, 8L))
    expect_identical(actual$countUns[actual$ind == "3"], rep(NA_integer_, 8L))
    expect_identical(actual$nCellStim, rep(8L, 16L))
    expect_identical(actual$nCellUns, rep(8L, 16L))
    expect_true(all(is.na(actual$freqBs[actual$ind == "3"])))
  })

  test_that("streaming preserves context-dependent cytPos gates", {
    ex <- tibble::tibble(A = c(11, 9, 9, 0), B = c(9, 11, 9, 0), C = c(0, 0, 11, 0))
    gates <- dplyr::mutate(.statsCombnGates(), gate = 10, gateCyt = 8)
    base <- .statsCombnFixture(list(ex, ex), gates)
    cyt <- .statsCombnFixture(list(ex, ex), gates, "cyt")
    expect_identical(base$countStim, c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L))
    expect_identical(cyt$countStim, c(0L, 0L, 0L, 2L, 0L, 0L, 1L, 1L))
    expect_identical(cyt$countUns, cyt$countStim)
  })

  test_that("NA fallback agrees with the original per-combination Reduce semantics", {
    stim <- tibble::tibble(A = c(NA_real_, 1, 0), B = c(0, NA_real_, 1), C = c(0, 1, 0))
    uns <- tibble::tibble(A = c(0, 1), B = c(1, 0), C = c(0, 1))
    gates <- .statsCombnGates()
    for (gateType in c("base", "cyt")) {
      actual <- .statsCombnFixture(list(uns, stim), gates, gateType)
      expect_identical(actual$countStim, .statsCombnReference(stim, gates, gateType))
      expect_identical(actual$countUns, .statsCombnReference(uns, gates, gateType))
      expect_true(anyNA(actual$countStim))
      expect_identical(actual$countStim[[8]], NA_integer_)
      expect_identical(actual$nCellStim, rep(3L, 8L))
      expect_identical(sum(actual$countUns), 2L)
      # Other known positives rule out C-only despite NA values in A/B.
      expect_identical(actual$countStim[[3]], 0L)
      reversed <- .statsCombnFixture(list(stim, uns), gates, gateType)
      expect_identical(reversed$countUns, actual$countStim)
      expect_identical(reversed$countStim, actual$countUns)
    }
  })

  test_that("fast counts match the original algorithm on seeded random tubes", {
    exList <- withr::with_seed(981, lapply(c(53L, 71L), function(n) {
      tibble::as_tibble(stats::setNames(
        lapply(seq_len(3L), function(k) stats::rnorm(n)), c("A", "B", "C")
      ))
    }))
    gates <- dplyr::mutate(.statsCombnGates(), gate = c(0, 0.4, -0.3), gateCyt = gate - 0.5)
    for (gateType in c("base", "cyt")) {
      actual <- .statsCombnFixture(exList, gates, gateType)
      expect_identical(actual$countStim, .statsCombnReference(exList[[2]], gates, gateType))
      expect_identical(actual$countUns, .statsCombnReference(exList[[1]], gates, gateType))
      expect_identical(sum(actual$countStim), 71L)
      expect_identical(sum(actual$countUns), 53L)
    }
  })

  test_that("combination codes handle empty and single-channel tubes and reject over 30 channels", {
    ex <- tibble::tibble(A = c(0, 1, 1))
    actual <- .statsCombnFixture(list(ex, ex), .statsCombnGates(chnl = "A"))
    expect_identical(actual$countStim, c(2L, 1L))
    empty <- .statsCombnFixture(list(ex[0, ], ex[0, ]), .statsCombnGates(chnl = "A"))
    expect_identical(empty$countStim, c(0L, 0L))
    expect_identical(empty$countUns, c(0L, 0L))
    expect_error(.getStatsCombnMatListGet(31L), "at most 30 channels")
  })

  test_that("streaming retains the incomplete-cache error without a GatingSet", {
    project <- tempfile("stats-combn-missing-")
    withr::defer(unlink(project, recursive = TRUE))
    expect_error(.getStats(
      gateTbl = .statsCombnGates(), chnl = c("A", "B", "C"),
      chnlLab = c(A = "A", B = "B", C = "C"), gateName = "g",
      gateTypeCytPosCalc = "base", popGate = "root", .data = NULL,
      indBatchList = list(batch = 1:2), pathProject = project
    ), "Incomplete expression cache for sample 2")
  })
})
