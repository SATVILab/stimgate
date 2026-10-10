local({
  root <- normalizePath(file.path(testthat::test_path(), "../../.."))
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "r", "analysis-plot-style.R"), env)
  sys.source(file.path(root, "scripts", "r", "acs_cytof-explore-cytpos.R"), env)

  dat <- list(
    stim = data.frame(IFNg = c(0, 1, 3, 3, 0.5, 2), TNF = c(0, 3, 0.2, 3, 0, 1),
      IL2 = c(0, 0, 0, 0, 2, 0)),
    uns = data.frame(IFNg = c(0, 0.2, 0.4, 0, 2.5), TNF = c(0, 0.1, 0, 0, 0),
      IL2 = c(0, 0, 0, 0, 0)),
    gate = c(IFNg = 1.5, TNF = 2, IL2 = NA),
    gateCyt = c(IFNg = 0.8, TNF = 2, IL2 = NA)
  )

  testthat::test_that("other-cytokine positivity is strict and ignores missing gates", {
    testthat::expect_equal(
      env$.acsCytposOtherPos(dat$stim, dat$gate, "IFNg"),
      c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
    )
    # TNF exactly at a gate is not positive.
    edge <- data.frame(IFNg = 0, TNF = 2, IL2 = 0)
    testthat::expect_false(env$.acsCytposOtherPos(edge, dat$gate, "IFNg"))
  })

  testthat::test_that("density table keeps non-zero values and both cell sets", {
    tbl <- env$.acsCytposDensityTbl(dat, "d1", "IFNg")
    testthat::expect_true(all(tbl$x > 0))
    stimOther <- tbl[tbl$tube == "stimulated" & tbl$cells != "all cells", ]
    testthat::expect_equal(sort(stimOther$x), c(1, 3))
    testthat::expect_equal(unique(tbl$zeroShare[tbl$tube == "stimulated" & tbl$cells == "all cells"]), 1 / 6)
    testthat::expect_equal(unname(env$.acsCytposDonorLabels(tbl)), "d1 (n other+ = 2)")
  })

  testthat::test_that("plots build with gate lines for both gate types", {
    gates <- list(d1 = dat[c("gate", "gateCyt")])
    hex <- env$.acsCytposHexPlot(env$.acsCytposHexTbl(dat, "d1"), gates)
    testthat::expect_s3_class(ggplot2::ggplot_build(hex), "ggplot_built")
    dens <- env$.acsCytposDensityPlot(env$.acsCytposDensityTbl(dat, "d1", "IFNg"), gates, "IFNg")
    built <- ggplot2::ggplot_build(dens)
    xs <- unlist(lapply(built$data[2:3], `[[`, "xintercept"))
    testthat::expect_setequal(xs, c(1.5, 0.8))
  })
})
