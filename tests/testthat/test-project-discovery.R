test_that(".gateGetPop and .gateGetChnl handle single and multiple populations correctly", {
  tmpProj <- withr::local_tempdir()
  dir.create(file.path(tmpProj, "gates", "poproot", "chnlBC1"), recursive = TRUE)
  dir.create(file.path(tmpProj, "gates", "poproot", "chnlBC2"), recursive = TRUE)

  expect_equal(.gateGetPop(tmpProj), "root")
  expect_equal(sort(.gateGetChnl(tmpProj, "root")), c("BC1", "BC2"))

  dir.create(file.path(tmpProj, "gates", "popCD4"), recursive = TRUE)
  expect_equal(sort(.gateGetPop(tmpProj)), c("CD4", "root"))
})

test_that("expression discovery preserves filtering, ordering and empty-path behaviour", {
  project <- withr::local_tempdir()
  expect_identical(.getExProjectPop(project), character())

  sample_dir <- file.path(project, "sampleData")
  dir.create(sample_dir)
  expect_error(.getExProjectPop(project), "Expected a non-empty character vector")

  root_dir <- file.path(sample_dir, "pop_root")
  dir.create(root_dir)
  dir.create(file.path(sample_dir, "unrelated"))
  saveRDS(1, file.path(sample_dir, "pop_file"))
  expect_identical(.getExProjectPop(project), "root")
  expect_error(
    .getExProjectInd(project, "root"), "Expected a non-empty character vector"
  )
  expect_identical(.getExProjectInd(project, "missing"), character())

  for (ind in c("2", "10")) {
    dir.create(file.path(root_dir, paste0("ind_", ind)))
  }
  dir.create(file.path(root_dir, "unrelated"))
  saveRDS(1, file.path(root_dir, "ind_file"))
  expect_identical(.getExProjectInd(project), c("10", "2"))

  channel_dir <- .getExChnlPathDir("10", "root", project)
  expect_error(.getExProjectChnl(project), "Expected a non-empty character vector")
  for (chnl in c("Z", "A")) {
    saveRDS(1:3, file.path(channel_dir, paste0("chnl_", chnl, ".rds")))
  }
  saveRDS(1, file.path(channel_dir, "unrelated.rds"))
  saveRDS(1, file.path(channel_dir, "chnl_ignore.txt"))
  expect_identical(.getExProjectChnl(project), c("A", "Z"))
  expect_identical(.getExProjectChnl(project, "root", "missing"), character())
})

test_that("plotStim error handling for multiple populations and empty inputs", {
  tmpProj <- withr::local_tempdir()
  dir.create(file.path(tmpProj, "gates", "poproot"), recursive = TRUE)
  dir.create(file.path(tmpProj, "gates", "popCD4"), recursive = TRUE)

  expect_error(
    plotStim(ind = c(1, 2), .data = NULL, pathProject = tmpProj, marker = "IL2"),
    "Cannot plot gates for multiple populations"
  )

  tmpProjEmpty <- file.path(withr::local_tempdir(), "missing")
  expect_error(
    plotStim(ind = c(1, 2), .data = NULL, pathProject = tmpProjEmpty, marker = "IL2"),
    "No population found for plotting gates"
  )

  expect_null(.plotGateUvMarker(
    marker = "IL2", chnl = "BC1", pop = "root", ind = numeric(0),
    excMin = FALSE, indLab = NULL, axisLab = NULL,
    showGate = FALSE, pathProject = tmpProjEmpty, minCell = 10,
    exArgs = list()
  ))
})
