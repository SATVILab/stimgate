test_that("stimgate_data_get_ex reads saved channel data and filters correctly", {
  tmp <- tempfile("stimgate_ex_")
  dir.create(file.path(tmp, "sample_data", "POP1", "ind_1"), recursive = TRUE)
  dir.create(file.path(tmp, "sample_data", "POP1", "ind_2"), recursive = TRUE)

  saveRDS(c(1, 2, 3), file = file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC1.rds"))
  saveRDS(c(4, 5, 6), file = file.path(tmp, "sample_data", "POP1", "ind_1", "chnl_BC2.rds"))
  saveRDS(c(7, 8), file = file.path(tmp, "sample_data", "POP1", "ind_2", "chnl_BC1.rds"))
  saveRDS(c(9, 10), file = file.path(tmp, "sample_data", "POP1", "ind_2", "chnl_BC2.rds"))
  res <- stimgate_data_get_ex(tmp)
  expect_equal(nrow(res), 5)
  expect_true(all(c("pop", "ind", "BC1", "BC2") %in% names(res)))
  expect_equal(unique(res$pop), "POP1")
  expect_equal(sum(res$ind == "1"), 3)
  expect_equal(sum(res$ind == "2"), 2)

  res_bc1 <- stimgate_data_get_ex(tmp, chnl = "BC1")
  expect_true("BC1" %in% names(res_bc1))
  expect_false("BC2" %in% names(res_bc1))

  res_ind1 <- stimgate_data_get_ex(tmp, ind = "1")
  expect_equal(nrow(res_ind1), 3)

  expect_error(stimgate_data_get_ex(""))
})