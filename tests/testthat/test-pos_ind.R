pkg_ns <- asNamespace("stimgate")

.getPosIndCache <- get(".getPosIndCache", envir = pkg_ns, mode = "function")
.getPosIndCacheGet <- get(".getPosIndCacheGet", envir = pkg_ns, mode = "function")
.getPosIndCacheCount <- get(".getPosIndCacheCount", envir = pkg_ns, mode = "function")
.getPosIndCacheAnyFromCount <- get(".getPosIndCacheAnyFromCount", envir = pkg_ns, mode = "function")
.getPosIndMult <- get(".getPosIndMult", envir = pkg_ns, mode = "function")
.getPosIndByChnl <- get(".getPosIndByChnl", envir = pkg_ns, mode = "function")
.getPosInd <- get(".getPosInd", envir = pkg_ns, mode = "function")
.getPosIndButSinglePosForOneCyt <- get(".getPosIndButSinglePosForOneCyt", envir = pkg_ns, mode = "function")
.getPosIndCytCombn <- get(".getPosIndCytCombn", envir = pkg_ns, mode = "function")
.getStatsBatchGnFilterMasks <- get(".getStatsBatchGnFilterMasks", envir = pkg_ns, mode = "function")
.getStatsBatchGnFilterOrNonCombn <- get(".getStatsBatchGnFilterOrNonCombn", envir = pkg_ns, mode = "function")

test_that(".getPosIndCache records and validates per-channel threshold comparisons", {
  ex <- data.frame(
    A = c(11, 8, 9, NA_real_),
    B = c(21, 19, 2, 30),
    stringsAsFactors = FALSE
  )

  gateTbl <- data.frame(
    chnl = c("A", "B"),
    gate = c(10, 20),
    gateCyt = c(5, 10),
    stringsAsFactors = FALSE
  )

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = c("A", "B")
  )

  expect_equal(posCache$n, nrow(ex))
  expect_equal(posCache$base[["A"]], c(TRUE, FALSE, FALSE, NA))
  expect_equal(posCache$base[["B"]], c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(posCache$cyt[["A"]], c(TRUE, TRUE, TRUE, NA))
  expect_equal(posCache$cyt[["B"]], c(TRUE, TRUE, FALSE, TRUE))

  expect_error(
    .getPosIndCacheGet(
      posCache = posCache,
      chnl = "C",
      gateType = "base"
    ),
    "No base positivity vector available for channel C"
  )
})

test_that(".getPosIndMult and .getPosIndByChnl preserve NA semantics and combination logic", {
  ex <- data.frame(
    A = c(11, 8, 9, NA_real_),
    B = c(21, 19, 2, 30),
    stringsAsFactors = FALSE
  )

  gateTbl <- data.frame(
    chnl = c("A", "B"),
    gate = c(10, 20),
    gateCyt = c(5, 10),
    stringsAsFactors = FALSE
  )

  posCache <- .getPosIndCache(
    ex = ex,
    gateTbl = gateTbl,
    chnl = c("A", "B")
  )

  countBase <- .getPosIndCacheCount(
    posCache = posCache,
    chnl = c("A", "B"),
    gateType = "base"
  )

  expect_equal(
    .getPosIndCacheAnyFromCount(countBase),
    c(TRUE, FALSE, FALSE, TRUE)
  )
  expect_equal(
    .getPosIndMult(
      ex = ex,
      gateTbl = gateTbl,
      chnl = c("A", "B"),
      gateTypeCytPos = "base"
    ),
    c(TRUE, FALSE, FALSE, NA)
  )

  posByChnl <- .getPosIndByChnl(
    ex = ex,
    gateTbl = gateTbl,
    chnl = c("A", "B"),
    gateTypeCytPos = "cyt",
    posCache = posCache
  )

  expect_equal(posByChnl[["A"]], c(TRUE, FALSE, FALSE, NA))
  expect_equal(posByChnl[["B"]], c(TRUE, FALSE, FALSE, TRUE))
})

test_that(".getPosInd, .getPosIndButSinglePosForOneCyt and .getPosIndCytCombn match the cytokine rules", {
  ex <- data.frame(
    A = c(11, 8, 9, NA_real_),
    B = c(21, 19, 2, 30),
    stringsAsFactors = FALSE
  )

  gateTbl <- data.frame(
    chnl = c("A", "B"),
    gate = c(10, 20),
    gateCyt = c(5, 10),
    stringsAsFactors = FALSE
  )

  expect_equal(
    .getPosInd(
      ex = ex,
      gateTbl = gateTbl,
      chnl = c("A", "B"),
      gateTypeCytPos = "base"
    ),
    c(TRUE, FALSE, FALSE, TRUE)
  )

  expect_equal(
    .getPosInd(
      ex = ex,
      gateTbl = gateTbl,
      chnl = "A",
      chnlAlt = "B",
      gateTypeCytPos = "cyt"
    ),
    c(TRUE, FALSE, FALSE, NA)
  )

  expect_equal(
    .getPosIndButSinglePosForOneCyt(
      ex = ex,
      gateTbl = gateTbl,
      chnlSingleExc = "B",
      chnl = "A",
      gateTypeCytPos = "base"
    ),
    c(TRUE, FALSE, FALSE, NA)
  )

  expect_equal(
    .getPosIndCytCombn(
      ex = ex,
      gateTbl = gateTbl,
      chnlPos = "A",
      chnlNeg = "B",
      chnlAlt = character(0),
      gateTypeCytPos = "cyt"
    ),
    c(FALSE, FALSE, FALSE, FALSE)
  )
})

test_that(".getStatsBatchGnFilterMasks keeps only cells free of other cytokines and respects cyt gating", {
  exUns <- data.frame(
    A = c(3, 15, 3, 17, 11),
    B = c(4, 3, 15, 20, 9),
    C = c(5, 4, 4, 4, 4),
    stringsAsFactors = FALSE
  )
  exStim <- data.frame(
    A = c(3, 15, 3, 17, 11, 25),
    B = c(4, 3, 15, 20, 9, 7),
    C = c(5, 4, 4, 4, 4, 8),
    stringsAsFactors = FALSE
  )

  attr(exUns, "ind") <- "uns"
  attr(exStim, "ind") <- "stim1"

  gateTbl <- tibble::tibble(
    ind = rep("stim1", 3L),
    gateName = "g1",
    chnl = c("A", "B", "C"),
    gate = c(10, 10, 10),
    gateCyt = c(8, 8, 8)
  )

  masksBase <- .getStatsBatchGnFilterMasks(
    exList = list(exUns, exStim),
    gateTblGn = gateTbl,
    chnl = c("A", "B", "C"),
    gateTypeCytPosFilter = "base"
  )
  expect_equal(masksBase[[1]]$stim$A, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(masksBase[[1]]$stim$B, c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(masksBase[[1]]$uns$A, c(FALSE, FALSE, TRUE, TRUE, FALSE))

  masksCyt <- .getStatsBatchGnFilterMasks(
    exList = list(exUns, exStim),
    gateTblGn = gateTbl,
    chnl = c("A", "B", "C"),
    gateTypeCytPosFilter = "cyt"
  )
  expect_equal(masksCyt[[1]]$stim$A, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(masksCyt[[1]]$uns$A, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_true(!identical(masksBase[[1]]$stim$A, masksCyt[[1]]$stim$A))
})

test_that(".getStatsBatchGnFilterOrNonCombn reports the filtered numerator and denominator for the current cytokine", {
  exUns <- data.frame(
    A = c(3, 15, 3, 17, 11),
    B = c(4, 3, 15, 20, 9),
    C = c(5, 4, 4, 4, 4),
    stringsAsFactors = FALSE
  )
  exStim <- data.frame(
    A = c(3, 15, 3, 17, 11, 25),
    B = c(4, 3, 15, 20, 9, 7),
    C = c(5, 4, 4, 4, 4, 8),
    stringsAsFactors = FALSE
  )

  attr(exUns, "ind") <- "uns"
  attr(exStim, "ind") <- "stim1"

  gateTbl <- tibble::tibble(
    ind = rep("stim1", 3L),
    gateName = "g1",
    chnl = c("A", "B", "C"),
    gate = c(10, 10, 10),
    gateCyt = c(8, 8, 8)
  )

  statBase <- .getStatsBatchGnFilterOrNonCombn(
    exList = list(exUns, exStim),
    indBatch = c("uns", "stim1"),
    gateTblGn = gateTbl,
    gn = "g1",
    chnl = c("A", "B", "C"),
    filterOtherCytPos = TRUE,
    gateTypeCytPosFilter = "base"
  )
  statBaseA <- statBase[statBase$chnl == "A", ]
  expect_equal(statBaseA$countStim, 3L)
  expect_equal(statBaseA$nCellStim, 4L)
  expect_equal(statBaseA$countUns, 2L)
  expect_equal(statBaseA$nCellUns, 3L)

  statCyt <- .getStatsBatchGnFilterOrNonCombn(
    exList = list(exUns, exStim),
    indBatch = c("uns", "stim1"),
    gateTblGn = gateTbl,
    gn = "g1",
    chnl = c("A", "B", "C"),
    filterOtherCytPos = TRUE,
    gateTypeCytPosFilter = "cyt"
  )
  statCytA <- statCyt[statCyt$chnl == "A", ]
  expect_equal(statCytA$countStim, 2L)
  expect_equal(statCytA$nCellStim, 3L)
  expect_equal(statCytA$countUns, 1L)
  expect_equal(statCytA$nCellUns, 2L)
})
