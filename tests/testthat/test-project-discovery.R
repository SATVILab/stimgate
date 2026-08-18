test_that(".gateGetPop and .gateGetChnl handle single and multiple populations correctly", {
  tmpProj <- file.path(tempdir(), paste0("proj_disc_test_", as.numeric(Sys.time())))
  dir.create(file.path(tmpProj, "gates", "poproot", "chnlBC1"), recursive = TRUE)
  dir.create(file.path(tmpProj, "gates", "poproot", "chnlBC2"), recursive = TRUE)

  expect_equal(.gateGetPop(tmpProj), "root")
  expect_equal(sort(.gateGetChnl(tmpProj, "root")), c("BC1", "BC2"))

  dir.create(file.path(tmpProj, "gates", "popCD4"), recursive = TRUE)
  expect_equal(sort(.gateGetPop(tmpProj)), c("CD4", "root"))

  unlink(tmpProj, recursive = TRUE)
})

test_that(".getExProjectPop, .getExProjectInd and .getExProjectChnl discover sampleData correctly", {
  tmpProj <- file.path(tempdir(), paste0("sample_disc_test_", as.numeric(Sys.time())))
  dir.create(file.path(tmpProj, "sampleData", "pop_root", "ind_1"), recursive = TRUE)
  dir.create(file.path(tmpProj, "sampleData", "pop_root", "ind_2"), recursive = TRUE)
  saveRDS(1:5, file.path(tmpProj, "sampleData", "pop_root", "ind_1", "chnl_BC1.rds"))

  expect_equal(.getExProjectPop(tmpProj), "root")
  expect_equal(sort(.getExProjectInd(tmpProj, "root")), c("1", "2"))
  expect_equal(.getExProjectChnl(tmpProj, "root", "1"), "BC1")

  unlink(tmpProj, recursive = TRUE)
})

test_that("plotStim error handling for multiple populations and empty inputs", {
  tmpProj <- file.path(tempdir(), paste0("plot_disc_test_", as.numeric(Sys.time())))
  dir.create(file.path(tmpProj, "gates", "poproot"), recursive = TRUE)
  dir.create(file.path(tmpProj, "gates", "popCD4"), recursive = TRUE)

  expect_error(plotStim(ind = c(1, 2), .data = NULL, pathProject = tmpProj, marker = "IL2"),
               "Cannot plot gates for multiple populations")

  unlink(tmpProj, recursive = TRUE)

  tmpProjEmpty <- file.path(tempdir(), paste0("plot_disc_empty_", as.numeric(Sys.time())))
  expect_error(plotStim(ind = c(1, 2), .data = NULL, pathProject = tmpProjEmpty, marker = "IL2"),
               "No population found for plotting gates")

  expect_null(.plotGateUvMarker(marker = "IL2", chnl = "BC1", pop = "root", ind = numeric(0),
                                 .data = NULL, excMin = FALSE, indLab = NULL, axisLab = NULL,
                                 showGate = FALSE, pathProject = tmpProjEmpty, minCell = 10,
                                 bias = FALSE, combnExc = NULL, chnlGate = NULL, markerGate = NULL,
                                 gateTypeCytPos = "cyt", mult = FALSE,
                                 gateUnsMethod = "min"))

  unlink(tmpProjEmpty, recursive = TRUE)
})
