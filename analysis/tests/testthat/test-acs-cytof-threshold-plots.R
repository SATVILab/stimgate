.acs_threshold_env <- function() {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), winslash = "/")
  env <- new.env(parent = getNamespace("stimgate"))
  for (file in c("analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-methods.R", "acs_cytof-threshold-plots.R")) {
    source(file.path(root, "scripts", "r", file), local = env)
  }
  env
}

test_that("ACS threshold tables keep all donors, errors and final base-gate provenance", {
  env <- .acs_threshold_env()
  channels <- names(env$.acsCytofChannelMap())
  map <- data.frame(ind = c("1", "2", "3"), SampleID = c("a", "a", "b"), stim = c("uns", "p1", "p1"))
  gates <- expand.grid(ind = c("2", "3"), chnl = channels, stringsAsFactors = FALSE)
  gates$gateName <- "loc_minClust"
  gates$gate <- 2
  gates$locSource <- "direct"
  gates$locGenerated <- TRUE
  gates$locGenerated[1] <- FALSE
  comparator <- transform(gates, threshold = 3, thresholdOrigin = "estimated", thresholdFallbackUsed = FALSE)
  comparator$threshold[1] <- NA_real_
  out <- env$.acsCytofThresholdTable(gates, list(fbeta = comparator, tailgate = comparator), map)
  expect_equal(nrow(out), 36L)
  expect_setequal(out$SampleID, c("a", "b"))
  expect_true(all(out$stim == "p1"))
  expect_equal(sum(is.na(out$threshold)), 2L)
  expect_equal(sum(out$thresholdFallbackUsed), 1L)
  expect_equal(out$threshold[out$method == "stimgate"], rep(2, 12))
  expect_error(env$.acsCytofThresholdTable(gates[-1, ], list(fbeta = comparator), map), "Missing or duplicate")
  expect_error(env$.acsCytofThresholdTable(rbind(gates, gates[1, ]), list(), map), "Missing or duplicate")
})

test_that("ACS threshold plots use common bins and a peak-scaled real density", {
  env <- .acs_threshold_env()
  gates <- data.frame(
    method = rep(c("stimgate", "fbeta", "tailgate"), each = 3),
    threshold = c(1, 2, 3, 2, 3, NA, 4, 4, 4), thresholdFallbackUsed = c(TRUE, rep(FALSE, 8))
  )
  cells <- data.frame(ind = rep(c("1", "2"), each = 100), expression = seq(0, 6, length.out = 200))
  plot <- env$.acsCytofThresholdPlot(gates, cells, "cd4 / p1 / IFNg")
  expect_no_error(ggplot2::ggplotGrob(plot))
  expect_match(plot$labels$subtitle, "F-beta (zeros removed): 2/3 finite", fixed = TRUE)
  expect_match(plot$labels$subtitle, "StimGate: 3/3 finite, 1 fallback")
  expect_equal(max(plot$layers[[1]]$data$y), max(plot$layers[[2]]$data$y))
  steps <- plot$layers[[2]]$data
  expect_equal(length(unique(steps$x)), 21L)
  expect_match(plot$labels$caption, "200 cells, 2 tubes")
  gates$threshold <- NA_real_
  expect_no_error(ggplot2::ggplotGrob(env$.acsCytofThresholdPlot(gates, cells, "all failed")))
  cells$expression <- 0
  gates$threshold <- 0
  expect_no_error(ggplot2::ggplotGrob(env$.acsCytofThresholdPlot(gates, cells, "constant")))
})

test_that("ACS tube sampling is reproducible, bounded and preserves zeros and RNG", {
  env <- .acs_threshold_env()
  testthat::skip_if_not_installed("flowWorkspace")
  channels <- names(env$.acsCytofChannelMap())
  frames <- lapply(1:3, function(i) {
    flowCore::flowFrame(matrix(
      rep(c(0, seq_len(19)), length(channels)),
      ncol = length(channels),
      dimnames = list(NULL, channels)
    ))
  })
  fs <- flowCore::flowSet(frames)
  gs <- flowWorkspace::GatingSet(fs)
  map <- data.frame(ind = c("1", "2", "3"), SampleID = c("a", "b", "c"), stim = c("uns", "p1", "p1"))
  withr::local_preserve_seed()
  set.seed(5)
  before <- .Random.seed
  sampled <- env$.acsCytofThresholdCells(gs, map, nCell = 10L)
  expect_identical(.Random.seed, before)
  expect_identical(sampled, env$.acsCytofThresholdCells(gs, map, nCell = 10L))
  expect_equal(nrow(sampled), 120L)
  expect_setequal(sampled$ind, c("2", "3"))
  all <- env$.acsCytofThresholdCells(gs, map, nCell = 100L)
  expect_equal(sum(all$expression == 0), 12L)
  expect_equal(nrow(all), 240L)
})
