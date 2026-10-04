test_that("statistics combinations retain their order and integer matrices", {
  for (n_chnl in seq_len(6L)) {
    expected <- stats::setNames(lapply(seq_len(n_chnl), function(n_pos) {
      t(utils::combn(n_chnl, n_pos))
    }), seq_len(n_chnl))
    expect_identical(.getStatsCombnMatListGet(n_chnl), expected)
  }
  combinations <- .getStatsCombnMatListGet(2L)
  expect_identical(
    .getStatsCytCombnVecListGet(combinations, c("A", "B")),
    list(`1` = c("A~+~B~-~", "A~-~B~+~"), `2` = "A~+~B~+~")
  )
})

test_that("statistics disk fallback retains gate filters and channel labels", {
  project <- tempfile("stats-fallback-")
  withr::defer(unlink(project, recursive = TRUE))
  path <- .gatesGetPathAll(project, "root", "BC1", FALSE)
  dir.create(dirname(path), recursive = TRUE)
  gates <- tibble::tibble(
    gateName = c("g", "gClust"), batch = "batch_1", ind = "2",
    gate = c(1, 2), gateCyt = c(0.5, 1)
  )
  saveRDS(gates, path)
  actual <- .getStatsGateTblGet(
    gateTbl = NULL, chnlLab = c(BC1 = "IL2"), pathProject = project,
    popGate = "root", gateName = "gClust", tolClust = TRUE
  )
  expect_identical(actual$gateName, "gClust")
  expect_identical(actual$chnl, "BC1")
  expect_identical(unname(actual$marker), "IL2")
  expect_identical(actual$gate, 2)
  expect_identical(actual$gateCyt, 1)
})


test_that("statistics gate filtering accepts multiple gate names without recycling", {
  project <- tempfile("stats-gate-names-")
  withr::defer(unlink(project, recursive = TRUE))
  gates <- tibble::tibble(
    chnl = "A", marker = "MarkerA", gateName = c("g", "other", "gClust"),
    batch = "batch", ind = "2", gate = c(1, 2, 3), gateCyt = c(0, 1, 2)
  )
  .getStatsGateTblSave(
    gateTbl = gates, pathProject = project, popGate = "root",
    chnlLab = c(A = "MarkerA"), chnl = "A", save = TRUE
  )
  actual <- .getStatsGateTblGet(
    gateTbl = NULL, chnlLab = c(A = "MarkerA"), pathProject = project,
    popGate = "root", gateName = c("gClust", "g")
  )
  expect_identical(actual$gateName, c("g", "gClust"))
  expect_identical(actual$gate, c(1, 3))
})
