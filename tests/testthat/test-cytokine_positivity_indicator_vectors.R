library(testthat)


test_that(".get_pos_ind_simple works correctly", {
  ex <- tibble::tibble(c1 = 1:2, c2 = 2:2)
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2"), gate = c(3, 1),
    gate_cyt = c(1.5, 3), gate_single = c(1.5, 2),
    ind = 1, gate_name = "a"
  )
  expect_identical(
    .get_pos_ind_simple(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, gate_type = "base"
    ),
    c(TRUE, TRUE)
  )
  expect_identical(
    .get_pos_ind_simple(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, gate_type = "single"
    ),
    c(FALSE, TRUE)
  )
  expect_identical(
    .get_pos_ind_simple(
      ex = ex |> dplyr::mutate(c1 = 2:1), gate_tbl = gate_tbl,
      chnl = NULL, gate_type = "cyt"
    ),
    c(TRUE, FALSE)
  )
  gate_tbl_2 <- gate_tbl |> dplyr::mutate(ind = 1:2)
  expect_error(.get_pos_ind_simple(
    ex = ex, gate_tbl = gate_tbl_2,
    chnl = NULL, gate_type = "single"
  ))
  gate_tbl_3 <- gate_tbl |> dplyr::mutate(gate_name = c("a", "b"))
  expect_error(.get_pos_ind_simple(
    ex = ex, gate_tbl = gate_tbl_3,
    chnl = NULL, gate_type = "single"
  ))
})

test_that(".get_pos_ind_mult works correctly", {
  ex <- tibble::tibble(c1 = 1:2, c2 = c(2, 2.1))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2"), gate = c(3, 2),
    gate_cyt = c(1.5, 1), gate_single = c(4, 2.5),
    ind = 1, gate_name = "a"
  )
  expect_identical(
    .get_pos_ind_mult(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, chnl_alt = NULL,
      gate_type_cyt_pos = "cyt"
    ),
    c(FALSE, TRUE)
  )
  expect_identical(
    .get_pos_ind_mult(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, chnl_alt = NULL,
      gate_type_cyt_pos = "base"
    ),
    c(FALSE, FALSE)
  )

  ex <- tibble::tibble(c1 = 1:2, c2 = c(2, 2.1), c3 = c(1.5, 1.5))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind_mult(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, chnl_alt = NULL,
      gate_type_cyt_pos = "base"
    ),
    c(FALSE, TRUE)
  )
  expect_identical(
    .get_pos_ind_mult(
      ex = ex, gate_tbl = gate_tbl,
      chnl = NULL, chnl_alt = NULL,
      gate_type_cyt_pos = "cyt"
    ),
    c(TRUE, TRUE)
  )
  expect_identical(
    .get_pos_ind_mult(
      ex = ex, gate_tbl = gate_tbl,
      chnl = c("c1", "c2"), chnl_alt = NULL,
      gate_type_cyt_pos = "cyt"
    ),
    c(FALSE, TRUE)
  )
})

test_that(".get_pos_ind_but_single_pos_for_one_cyt works correctly", {
  ex <- tibble::tibble(c1 = 1:2, c2 = c(2, 2.1), c3 = c(1.5, 1.5))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    c(TRUE, TRUE)
  )

  ex <- tibble::tibble(c1 = 1:2, c2 = c(2, 1.5), c3 = c(1.5, 1))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    c(TRUE, FALSE)
  )

  ex <- tibble::tibble(c1 = c(1, 10), c2 = c(2, 1.5), c3 = c(1.5, 1))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    c(TRUE, FALSE)
  )

  # single-positive only for cell to be excluded
  ex <- tibble::tibble(c1 = c(10), c2 = c(1), c3 = c(0.25))
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )
  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "cyt",
      gate_type_single_pos = "base"
    ),
    FALSE
  )

  # multi-positive, but only using cyt-pos because cell to be excluded
  # is base-positive
  ex <- tibble::tibble(c1 = c(3.1), c2 = c(1.2), c3 = 0.75)
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "cyt",
      gate_type_single_pos = "base"
    ),
    TRUE
  )

  expect_identical(
    .get_pos_ind_but_single_pos_for_one_cyt(
      ex = ex, gate_tbl = gate_tbl,
      chnl_single_exc = "c1",
      chnl = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    FALSE
  )
})


test_that(".get_pos_ind works correctly", {
  ex <- tibble::tibble(c1 = 1, c2 = 2, c3 = 1.5)
  gate_tbl <- tibble::tibble(
    chnl = c("c1", "c2", "c3"), gate = c(3, 2, 1),
    gate_cyt = c(1.5, 1, 0.5), gate_single = c(4, 2.5, 2),
    ind = 1, gate_name = "a"
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = NULL,
      chnl_alt = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    TRUE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c1", "c2"),
      chnl_alt = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    FALSE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c1", "c2"),
      chnl_alt = c("c3"),
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    FALSE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c3"),
      chnl_alt = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    TRUE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c3"),
      chnl_alt = NULL,
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "single"
    ),
    FALSE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c3"),
      chnl_alt = "",
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "single"
    ),
    FALSE
  )

  expect_identical(
    .get_pos_ind(
      ex = ex,
      gate_tbl = gate_tbl,
      chnl = c("c3"),
      chnl_alt = "",
      gate_type_cyt_pos = "base",
      gate_type_single_pos = "base"
    ),
    TRUE
  )
})
