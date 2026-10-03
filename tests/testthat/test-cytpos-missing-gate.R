test_that("cytokine refinement skips channels without a sample gate", {
  ex <- data.frame(A = c(0, 2), B = c(0, 2))
  gates <- tibble::tibble(chnl = "B", gate = 1)
  expect_identical(
    .getCpPosGatesChnl(
      chnlCurr = "A", ex = ex, gateTblInd = gates,
      basePos = .getCytPosBasePos(ex, gates), bwMin = 0.1,
      ind = 2L, stage = "cytPos", pathProject = tempdir()
    ),
    NA_real_
  )
})
