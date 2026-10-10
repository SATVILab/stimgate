root_dir <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

.omipShiftedPeakTestEnv <- function() {
  env <- new.env(parent = getNamespace("stimgate"))
  for (f in c(
    "analysis-runtime.R", "analysis-plot-style.R", "acs_cytof-helper.R",
    "acs_cytof-gate.R", "sim-misc.R", "sim-compare-freq_bs.R",
    "acs_cytof-methods.R", "sim-debug-loc.R", "acs_cytof-debug.R",
    "omip016-prepare.R", "omip016-methods.R", "omip-shifted-peak.R"
  )) {
    source(file.path(root_dir, "scripts", "r", f), local = env)
  }
  env
}

test_that("shifted-peak settings change only the rule and the semantics", {
  env <- .omipShiftedPeakTestEnv()
  default <- list(semantics = "omip111-v4", seed = 1L, control = list(bwMtd = "nrd0"))
  shifted <- env$.omipShiftedPeakSettings(default)
  expect_identical(shifted$semantics, "omip111-v4-shifted-peak")
  expect_true(shifted$control$locShiftedPeakRef)
  expect_identical(shifted$control$bwMtd, "nrd0")
  expect_identical(shifted$seed, 1L)
  # The settings must build a valid stimControl().
  expect_true(do.call(stimgate::stimControl, shifted$control)$locShiftedPeakRef)
})

test_that("OMIP-016 settings keep Analysis 15's saved form when the rule is off", {
  env <- .omipShiftedPeakTestEnv()
  default <- env$.omip016StimGateSettings()
  expect_false("locShiftedPeakRef" %in% names(default))
  shifted <- env$.omip016StimGateSettings(locShiftedPeakRef = TRUE)
  expect_true(shifted$locShiftedPeakRef)
  expect_identical(shifted[names(default)], default)
})

test_that("OMIP-111 default and shifted-peak outcomes are paired by tube", {
  env <- .omipShiftedPeakTestEnv()
  samples <- data.frame(
    sample = c("u1", "s1", "u2", "s2"), strain = "C57",
    mouse = c("M1", "M1", "M2", "M2"), condition = c("uns", "stim", "uns", "stim")
  )
  base <- expand.grid(
    mouse = c("M1", "M2"), method = c("StimGate", "F-beta"),
    stringsAsFactors = FALSE
  )
  base$strain <- "C57"
  base$population <- "CD4"
  base$marker <- "TNF"
  base$sampleStim <- ifelse(base$mouse == "M1", "s1", "s2")
  default <- base
  default$threshold <- c(5, 5, 2, 2)
  default$propBs <- c(0.01, 0.02, 0.8, 0.8)
  default$errorPp <- c(-80, -79, 0, 1)
  shifted <- default
  shifted$threshold[1:2] <- c(2, 5)
  shifted$propBs[1:2] <- c(0.8, 0.02)
  shifted$errorPp[1:2] <- c(0, -79)
  fired <- data.frame(
    strain = "C57", population = "CD4", ind = c("2", "4"), chnl = "TNF",
    marker = "TNF", shiftedPeakRef = c(TRUE, FALSE)
  )
  paired <- env$.omip111ShiftedPeakPaired(default, shifted, fired, samples)
  sg <- paired[paired$method == "StimGate", ]
  expect_identical(sg$shiftedPeakRef[order(sg$mouse)], c(TRUE, FALSE))
  expect_false(any(paired$shiftedPeakRef[paired$method != "StimGate"]))
  expect_true(env$.omipShiftedPeakComparatorsUnchanged(paired))
  summary <- env$.omipShiftedPeakErrorSummary(paired, "marker")
  expect_identical(summary$rule_applied, 1L)
  expect_equal(summary$mean_error_default_pp, -79.5)
  expect_equal(summary$mean_error_shifted_pp, -39.5)
  expect_identical(summary$changed, 1L)
  expect_s3_class(env$.omipShiftedPeakErrorPlot(paired, "strain"), "ggplot")

  paired$propBs[paired$method == "F-beta"][1] <- 0.5
  expect_false(env$.omipShiftedPeakComparatorsUnchanged(paired))
  expect_error(
    env$.omip111ShiftedPeakPaired(default[-1, ], shifted, fired, samples),
    "cohorts differ"
  )
})
