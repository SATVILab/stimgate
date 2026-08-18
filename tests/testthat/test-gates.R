test_that("project discovery helpers return clean population and channel vectors", {
  tmpProj <- file.path(tempdir(), paste0("stimgate_disc_test_", sample(10000, 1)))
  on.exit(unlink(tmpProj, recursive = TRUE), add = TRUE)

  # 1. Single population in gates/
  dir.create(file.path(tmpProj, "gates", "poproot", "chnlIL2"), recursive = TRUE)

  pops <- stimgate:::.gateGetPop(tmpProj)
  expect_equal(pops, "root")
  expect_false("" %in% pops)

  chnls <- stimgate:::.gateGetChnl(tmpProj, "root")
  expect_equal(chnls, "IL2")
  expect_false("" %in% chnls)

  # 2. Multiple populations in gates/
  dir.create(file.path(tmpProj, "gates", "popCD4"), recursive = TRUE)
  dir.create(file.path(tmpProj, "gates", "popCD8"), recursive = TRUE)

  popsMulti <- stimgate:::.gateGetPop(tmpProj)
  expect_true(all(c("root", "CD4", "CD8") %in% popsMulti))
  expect_false("" %in% popsMulti)

  # 3. Sample data discovery in sampleData/
  dir.create(file.path(tmpProj, "sampleData", "pop_root", "ind_1"), recursive = TRUE)
  saveRDS(c(1, 2, 3), file.path(tmpProj, "sampleData", "pop_root", "ind_1", "chnl_IL2.rds"))

  samplePops <- stimgate:::.getExProjectPop(tmpProj)
  expect_equal(samplePops, "root")
  expect_false("" %in% samplePops)

  sampleInds <- stimgate:::.getExProjectInd(tmpProj, "root")
  expect_equal(sampleInds, "1")
  expect_false("" %in% sampleInds)

  sampleChnls <- stimgate:::.getExProjectChnl(tmpProj, "root", "1")
  expect_equal(sampleChnls, "IL2")
  expect_false("" %in% sampleChnls)
})
