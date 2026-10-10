.omip111ParityEnv <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  env <- new.env(parent = getNamespace("stimgate"))
  source(file.path(root, "scripts", "r", "omip111-debug.R"), local = env)
  env
}

.omip111ParityFixture <- function() {
  gates <- data.frame(
    gateName = c("loc_minClust", "loc_minClust", "loc"),
    ind = c(2, 2, 2), marker = c("IFNg", "TNF", "IFNg"),
    gate = c(4.3070608490964881, 1.5, 9), gateCyt = c(4.3070608490964881, NA, 9)
  )
  stats <- data.frame(
    gateName = "loc_minClust", ind = 2,
    cytCombn = c("IFNg~+~TNF~+~", "IFNg~-~TNF~+~"),
    countStim = c(10, 3), countUns = c(1, 0)
  )
  list(gates = gates, stats = stats)
}

test_that("OMIP-111 diagnostic parity tolerates last-bit gate differences only", {
  env <- .omip111ParityEnv()
  x <- .omip111ParityFixture()
  rerun <- x$gates
  # As observed between OpenBLAS thread counts.
  rerun$gate[[1]] <- 4.3070608490964908
  rerun$gateCyt[[1]] <- 4.3070608490964908
  expect_true(env$.omip111DebugCheckParity(
    rerun[3:1, ], x$gates, x$stats[2:1, ], x$stats, "C57 CD8"
  ))

  moved <- x$gates
  moved$gate[[1]] <- moved$gate[[1]] + 1e-6
  expect_error(
    env$.omip111DebugCheckParity(moved, x$gates, x$stats, x$stats, "C57 CD8"),
    "differs from saved Analysis 14 gates: C57 CD8"
  )
  noCyt <- x$gates
  noCyt$gateCyt[[1]] <- NA
  expect_error(
    env$.omip111DebugCheckParity(noCyt, x$gates, x$stats, x$stats, "C57 CD8"),
    "differs from saved"
  )
  expect_error(
    env$.omip111DebugCheckParity(x$gates[-2, ], x$gates, x$stats, x$stats, "C57 CD8"),
    "differs from saved"
  )

  counts <- x$stats
  counts$countStim[[1]] <- 11
  expect_error(
    env$.omip111DebugCheckParity(x$gates, x$gates, counts, x$stats, "C57 CD8"),
    "changes saved Analysis 14 combination counts: C57 CD8"
  )
})
