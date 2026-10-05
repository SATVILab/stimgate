test_that("sample frequencies match applied gates when expression ties the gate", {
  ex_uns <- data.frame(CD4 = c(0, 2, 2, 3))
  ex_stim <- data.frame(CD4 = c(1, 2, 2, 3, 4))
  attr(ex_uns, "chnlCut") <- "CD4"
  attr(ex_stim, "chnlCut") <- "CD4"
  gate_tbl <- data.frame(chnl = "CD4", gate = 2)

  detail <- stimgate:::.getCpUnsLocSampleDetailTbl(
    cpVec = c(stim = 2, uns = 2),
    exListOrig = list(uns = ex_uns, stim = ex_stim),
    indUns = "uns",
    indStim = "stim",
    stage = "init",
    chnl = "CD4"
  )
  pos_stim <- stimgate:::.getPosInd(
    ex = ex_stim, gateTbl = gate_tbl, chnl = "CD4", gateTypeCytPos = "base"
  )
  pos_uns <- stimgate:::.getPosInd(
    ex = ex_uns, gateTbl = gate_tbl, chnl = "CD4", gateTypeCytPos = "base"
  )

  expect_identical(pos_stim, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(pos_uns, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(detail$propStim, 2 / 5)
  expect_equal(detail$propUns, 1 / 4)
  expect_equal(detail$propStim, sum(pos_stim) / nrow(ex_stim))
  expect_equal(detail$propUns, sum(pos_uns) / nrow(ex_uns))
  expect_equal(detail$propBs, mean(pos_stim) - mean(pos_uns))
})
