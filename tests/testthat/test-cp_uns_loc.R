test_that(".getCpUnsLocSample re-checks unstim cells after removal", {
  withr::local_preserve_seed()
  mk <- function(x, ind) {
    tbl <- tibble::tibble(marker = sort(x), ind = ind)
    attr(tbl, "chnlCut") <- "marker"
    attr(tbl, "ind") <- ind
    attr(tbl, "indUns") <- "uns"
    tbl
  }
  set.seed(1)
  ex_list <- list(
    uns = mk(stats::rnorm(100), "uns"),
    stim = mk(stats::rnorm(100, 1), "stim")
  )
  condition_called <- FALSE
  local_mocked_bindings(
    # cytokine-positive removal leaves only five unstim cells
    .getCpUnsLocSampleUnsRmCytPos = function(...) {
      list(...)[["exTblUnsBias"]][seq_len(5), , drop = FALSE]
    },
    .getCpUnsLocCondition = function(...) {
      condition_called <<- TRUE
      stop("condition-level gating should not run with too few unstim cells")
    },
    .getCpUnsLocOutput = function(...) list(...)[["cpUnsLocObjList"]],
    .package = "stimgate"
  )
  out <- stimgate:::.getCpUnsLocSample(
    exListOrig = ex_list,
    exListNoMinStim = ex_list[-1],
    exTblUnsBias = ex_list[[1]],
    chnlSettings = list(chnlCut = "marker", minCell = 20, cpMin = 0),
    bias = 0,
    pathProject = tempdir(),
    stage = "cytPos"
  )
  expect_false(condition_called)
  expect_identical(out$stim$locReason, "Too few cells")
  expect_false(out$stim$locGenerated)
  expect_true(is.finite(out$stim$cp))
})
