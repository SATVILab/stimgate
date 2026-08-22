# Deterministic regressions for the paired-density sample clustering route.

make_density_ex <- function(values) {
  ex <- data.frame(expr = as.numeric(values))
  attr(ex, "chnlCut") <- "expr"
  ex
}

make_pair_fixture <- function(ind, batch, stim_shift = 0, uns_shift = 0) {
  x1 <- c(seq(0.03, 0.20, length.out = 200), seq(0.45, 0.75, length.out = 200))
  x2 <- x1 + stim_shift
  x3 <- x1 + uns_shift
  list(
    ind = ind,
    batch = batch,
    stim = make_density_ex(x2),
    uns = make_density_ex(x3)
  )
}

rng_snapshot <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

test_that("paired density clustering uses both stimulated and unstimulated blocks on one shared grid", {
  densityGrid <- seq(0, 1, length.out = 32)
  exLookup <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    C = make_pair_fixture("C", "batch2", stim_shift = 0.58, uns_shift = 0.04),
    D = make_pair_fixture("D", "batch2", stim_shift = 0.10, uns_shift = 0.62)
  )

  featureTbl <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    densityGrid = densityGrid,
    bw = 0.08
  )

  unsCols <- grep("^uns_x\\d+$", names(featureTbl), value = TRUE)
  stimCols <- grep("^stim_x\\d+$", names(featureTbl), value = TRUE)

  expect_equal(length(unsCols), length(densityGrid))
  expect_equal(length(stimCols), length(densityGrid))
  expect_equal(unsCols, sprintf("uns_x%03d", seq_along(unsCols)))
  expect_equal(stimCols, sprintf("stim_x%03d", seq_along(stimCols)))
  expect_identical(
    names(featureTbl)[-(1:2)],
    c(unsCols, stimCols)
  )
})

test_that("near-identical paired densities co-cluster while clearly separated ones do not", {
  exLookup <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    C = make_pair_fixture("C", "batch2", stim_shift = 0.58, uns_shift = 0.04),
    D = make_pair_fixture("D", "batch2", stim_shift = 0.10, uns_shift = 0.62)
  )

  featureTbl <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    densityGrid = seq(0, 1, length.out = 32),
    bw = 0.08
  )

  clusterObj <- .getCpClusterLocClusters(
    featureTbl = featureTbl,
    control = list(
      maxClusters = 6L,
      gapBootstraps = 20L,
      kmeansNstart = 5L,
      seed = 1L
    )
  )

  grp <- clusterObj$clusterTbl
  expect_equal(grp$grp[match(c("A", "B"), grp$ind)], c("1", "1"))
  expect_equal(grp$grp[match(c("C", "D"), grp$ind)], c("2", "3"))
  expect_true(grp$grp[match("A", grp$ind)] != grp$grp[match("C", grp$ind)])
  expect_true(grp$grp[match("B", grp$ind)] != grp$grp[match("D", grp$ind)])
})

test_that("changing only the stimulated distribution changes the stimulated feature block while leaving the unstimulated block fixed", {
  densityGrid <- seq(0, 1, length.out = 32)
  exLookupStim <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch1", stim_shift = 0.38, uns_shift = 0.03)
  )
  exLookupUns <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch2", stim_shift = 0.12, uns_shift = 0.39)
  )

  featureStim <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookupStim,
    densityGrid = densityGrid,
    bw = 0.08
  )
  featureUns <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookupUns,
    densityGrid = densityGrid,
    bw = 0.08
  )

  unsCols <- grep("^uns_x\\d+$", names(featureStim), value = TRUE)
  stimCols <- grep("^stim_x\\d+$", names(featureStim), value = TRUE)

  expect_equal(featureStim[1, unsCols], featureStim[2, unsCols])
  expect_false(isTRUE(all.equal(featureStim[1, stimCols], featureStim[2, stimCols])))

  expect_equal(featureUns[1, stimCols], featureUns[2, stimCols])
  expect_false(isTRUE(all.equal(featureUns[1, unsCols], featureUns[2, unsCols])))
})

test_that("stimulated samples that share an unstimulated control in the same batch reuse the same unstimulated feature block", {
  densityGrid <- seq(0, 1, length.out = 32)
  exLookup <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch1", stim_shift = 0.50, uns_shift = 0.03),
    C = make_pair_fixture("C", "batch2", stim_shift = 0.35, uns_shift = 0.60)
  )

  featureTbl <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    densityGrid = densityGrid,
    bw = 0.08
  )

  unsCols <- grep("^uns_x\\d+$", names(featureTbl), value = TRUE)
  stimCols <- grep("^stim_x\\d+$", names(featureTbl), value = TRUE)

  expect_equal(featureTbl[match("A", featureTbl$ind), unsCols],
               featureTbl[match("B", featureTbl$ind), unsCols])
  expect_false(isTRUE(all.equal(featureTbl[match("A", featureTbl$ind), stimCols],
                                featureTbl[match("B", featureTbl$ind), stimCols])))
})

test_that("a degenerate one-structure fixture yields a single cluster without changing the caller RNG state", {
  densityGrid <- seq(0, 1, length.out = 32)
  exLookup <- list(
    A = make_pair_fixture("A", "batch1", stim_shift = 0.12, uns_shift = 0.03),
    B = make_pair_fixture("B", "batch1", stim_shift = 0.12, uns_shift = 0.03)
  )

  featureTbl <- .getCpClusterLocJointFeatureTbl(
    exLookup = exLookup,
    densityGrid = densityGrid,
    bw = 0.08
  )

  rngBefore <- rng_snapshot()
  clusterObj1 <- .getCpClusterLocClusters(
    featureTbl = featureTbl,
    control = list(
      maxClusters = 6L,
      gapBootstraps = 20L,
      kmeansNstart = 5L,
      seed = 1L
    )
  )
  clusterObj2 <- .getCpClusterLocClusters(
    featureTbl = featureTbl,
    control = list(
      maxClusters = 6L,
      gapBootstraps = 20L,
      kmeansNstart = 5L,
      seed = 1L
    )
  )
  rngAfter <- rng_snapshot()

  expect_equal(clusterObj1$nInitialClusters, 1L)
  expect_equal(clusterObj1$nFinalClusters, 1L)
  expect_equal(clusterObj1$clusterTbl$grp, c("1", "1"))
  expect_identical(clusterObj1$clusterTbl, clusterObj2$clusterTbl)

  if (is.null(rngBefore)) {
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  } else {
    expect_identical(rngAfter, rngBefore)
  }
})

test_that(".getCpClusterDensTblGetBatchPrepExListFilterInd filters other cytokine positive cells", {
  exTbl <- tibble::tibble(
    IFNg = c(10, 2, 12, 1),
    IL2 = c(1, 15, 12, 2)
  )
  attr(exTbl, "ind") <- "2"
  attr(exTbl, "batch") <- "batch1"

  gateTbl <- tibble::tibble(
    ind = c("2", "2"),
    chnl = c("IFNg", "IL2"),
    gate = c(5, 5),
    gateCyt = c(8, 8)
  )

  resBase <- .getCpClusterDensTblGetBatchPrepExListFilterInd(
    exTbl = exTbl,
    gateTbl = gateTbl,
    chnlCut = "IFNg",
    calcCytPosGates = FALSE
  )
  expect_equal(resBase$IFNg, c(10, 1))
  expect_equal(resBase$IL2, c(1, 2))
  expect_equal(attr(resBase, "ind"), "2")

  resIl2 <- .getCpClusterDensTblGetBatchPrepExListFilterInd(
    exTbl = exTbl,
    gateTbl = gateTbl,
    chnlCut = "IL2",
    calcCytPosGates = FALSE
  )
  expect_equal(resIl2$IFNg, c(2, 1))
  expect_equal(resIl2$IL2, c(15, 2))

  gateTblCyt <- tibble::tibble(
    ind = c("2", "2"),
    chnl = c("IFNg", "IL2"),
    gate = c(5, 5),
    gateCyt = c(10, 10)
  )
  resCyt <- .getCpClusterDensTblGetBatchPrepExListFilterInd(
    exTbl = exTbl,
    gateTbl = gateTblCyt,
    chnlCut = "IFNg",
    calcCytPosGates = TRUE
  )
  expect_equal(resCyt$IFNg, c(10, 1))
})

test_that(".getCpClusterDensTblGetBatchPrepExListFilter filters a list of samples", {
  exUnstim <- tibble::tibble(IFNg = c(1, 2), IL2 = c(1, 1))
  attr(exUnstim, "ind") <- "1"

  exStim1 <- tibble::tibble(IFNg = c(10, 2), IL2 = c(1, 15))
  attr(exStim1, "ind") <- "2"

  exStim2 <- tibble::tibble(IFNg = c(12, 1), IL2 = c(12, 2))
  attr(exStim2, "ind") <- "3"

  exList <- list("1" = exUnstim, "2" = exStim1, "3" = exStim2)

  gateTbl <- tibble::tibble(
    ind = c("2", "2", "3", "3"),
    chnl = c("IFNg", "IL2", "IFNg", "IL2"),
    gate = c(5, 5, 5, 5)
  )

  resList <- .getCpClusterDensTblGetBatchPrepExListFilter(
    exList = exList,
    chnlCut = "IFNg",
    gateTbl = gateTbl,
    calcCytPosGates = FALSE
  )

  expect_named(resList, c("2", "3"))
  expect_equal(nrow(resList[["2"]]), 1)
  expect_equal(nrow(resList[["3"]]), 1)
})
