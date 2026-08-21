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
