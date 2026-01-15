test_that("stimgate_data_get_ex reads saved channel data and filters correctly", {
  tmp <- tempfile("stimgate_ex_")
  dir.create(
    file.path(tmp, "sample_data", "pop_POP1", "ind_1"),
    recursive = TRUE
  )
  dir.create(
    file.path(tmp, "sample_data", "pop_POP1", "ind_2"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )
  saveRDS(
    c(7, 8),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_2", "chnl_BC1.rds")
  )
  saveRDS(
    c(9, 10),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_2", "chnl_BC2.rds")
  )
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

test_that("stimgate_data_get_ex applies bias only to unstim sample", {
  tmp <- tempfile("stimgate_ex_bias_")
  dir.create(
    file.path(tmp, "sample_data", "pop_POP1", "ind_1"),
    recursive = TRUE
  )
  dir.create(
    file.path(tmp, "sample_data", "pop_POP1", "ind_2"),
    recursive = TRUE
  )

  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 5, 6),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )
  saveRDS(
    c(7, 8),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_2", "chnl_BC1.rds")
  )
  saveRDS(
    c(9, 10),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_2", "chnl_BC2.rds")
  )

  # Create meta_data with matching names so chnl lookup works
  chnl_list <- list(BC1 = list(bias_uns = 10), BC2 = list(bias_uns = -2))
  chnl_lab <- c(BC1 = "BC1", BC2 = "BC2")
  batch_list <- list(batch1 = c("1", "2"))

  dir.create(file.path(tmp, "meta_data"), showWarnings = FALSE)
  saveRDS(chnl_list, file.path(tmp, "meta_data", "chnl_list.rds"))
  saveRDS(chnl_lab, file.path(tmp, "meta_data", "chnl_lab.rds"))
  saveRDS(batch_list, file.path(tmp, "meta_data", "batch_list.rds"))

  # unstim is ind 2 -> expect bias added
  res_uns <- stimgate_data_get_ex(tmp, ind = "2", bias = TRUE)
  expect_equal(res_uns$BC1, c(7 + 10, 8 + 10))
  expect_equal(res_uns$BC2, c(9 - 2, 10 - 2))

  # stim is ind 1 -> bias should not be applied
  res_stim <- stimgate_data_get_ex(tmp, ind = "1", bias = TRUE)
  expect_equal(res_stim$BC1, c(1, 2, 3))
  expect_equal(res_stim$BC2, c(4, 5, 6))
})

test_that("stimgate_data_get_ex excludes minimum observed values when exc_min = TRUE", {
  tmp <- tempfile("stimgate_ex_excmin_")
  dir.create(
    file.path(tmp, "sample_data", "pop_POP1", "ind_1"),
    recursive = TRUE
  )

  # BC1 min is 1, BC2 min is 4; after exclusion only the row (3,6) should remain
  saveRDS(
    c(1, 2, 3),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC1.rds")
  )
  saveRDS(
    c(4, 4, 6),
    file = file.path(tmp, "sample_data", "pop_POP1", "ind_1", "chnl_BC2.rds")
  )

  res_noexc <- stimgate_data_get_ex(tmp, ind = "1", exc_min = FALSE)
  expect_equal(nrow(res_noexc), 3)

  res_exc <- stimgate_data_get_ex(tmp, ind = "1", exc_min = TRUE)
  expect_equal(nrow(res_exc), 1)
  expect_equal(res_exc$BC1, 3)
  expect_equal(res_exc$BC2, 6)
})
