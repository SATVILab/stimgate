pkg_ns <- asNamespace("stimgate")

.getCpUnsLocGetProbSmooth <- get(
  ".getCpUnsLocGetProbSmooth",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothCheckNCell <- get(
  ".getCpUnsLocGetProbSmoothCheckNCell",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothCheckNCellOut <- get(
  ".getCpUnsLocGetProbSmoothCheckNCellOut",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothObjPred <- get(
  ".getCpUnsLocGetProbSmoothObjPred",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothAttachDeriv <- get(
  ".getCpUnsLocGetProbSmoothAttachDeriv",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothNewData <- get(
  ".getCpUnsLocGetProbSmoothNewData",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothFitEval <- get(
  ".getCpUnsLocGetProbSmoothFitEval",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothFitEvalCheck <- get(
  ".getCpUnsLocGetProbSmoothFitEvalCheck",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothActualCheck <- get(
  ".getCpUnsLocGetProbSmoothActualCheck",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothActualThird <- get(
  ".getCpUnsLocGetProbSmoothActualThird",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetProbSmoothDerivativeTbl <- get(
  ".getCpUnsLocGetProbSmoothDerivativeTbl",
  envir = pkg_ns,
  mode = "function"
)

test_that(".getCpUnsLocGetProbSmooth produces monotone fit", {
  tmp_dir <- file.path(tempdir(), "test_prob_smooth_monotone")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  x_vals <- seq(0, 5, length.out = 40)
  prob_vals <- 1 / (1 + exp(-(x_vals - 2.5) * 2))

  data_mod <- data.frame(
    IFNg = x_vals,
    probSmooth = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "IFNg"
  attr(data_mod, "idxMod") <- seq_along(x_vals)
  attr(data_mod, "ind") <- 1L
  attr(data_mod, "locDensityBw") <- 0.35
  attr(data_mod, "locStimDensity") <- list(x = seq(0, 5, length.out = 100))
  attr(data_mod, "locDensityComparison") <- list(diff = rep(0, 100))
  attr(data_mod, "locPeakX") <- 1.2
  attr(data_mod, "locWindowWidth") <- 1.8

  smooth_out <- .getCpUnsLocGetProbSmooth(
    dataMod = data_mod,
    stage = "init",
    pathProject = tmp_dir,
    chnl = "IFNg",
    chnlSettings = list()
  )

  # Output structure
  expect_s3_class(smooth_out, "data.frame")
  expect_equal(nrow(smooth_out), length(x_vals))
  expect_true("pred" %in% names(smooth_out))

  # Predictions are finite and bounded in [0, 1]
  expect_true(all(is.finite(smooth_out$pred)))
  expect_true(all(smooth_out$pred >= 0 & smooth_out$pred <= 1))

  # Fitted response is monotone increasing to numerical tolerance
  diff_pred <- diff(smooth_out$pred[order(smooth_out$IFNg)])
  expect_true(all(diff_pred >= -1e-6))

  # Smoothing method metadata
  expect_equal(attr(smooth_out, "locProbSmoothMethod"), "scam_mpi")

  # Derivative table contract
  deriv_tbl <- attr(smooth_out, "locProbDerivTbl")
  expect_s3_class(deriv_tbl, "tbl_df")
  expect_named(deriv_tbl, c("x", "pred", "deriv"))
  expect_equal(nrow(deriv_tbl), 512L)
  expect_true(all(is.finite(deriv_tbl$x)))
  expect_true(all(is.finite(deriv_tbl$pred)))
  expect_true(all(is.finite(deriv_tbl$deriv)))
  expect_true(all(deriv_tbl$deriv >= 0))
  expect_equal(min(deriv_tbl$x), min(x_vals))
  expect_equal(max(deriv_tbl$x), max(x_vals))

  # Retained attributes survive
  expect_equal(attr(smooth_out, "locDensityBw"), 0.35)
  expect_equal(attr(smooth_out, "locPeakX"), 1.2)
  expect_equal(attr(smooth_out, "locWindowWidth"), 1.8)
  expect_equal(
    attr(smooth_out, "locStimDensity"),
    attr(data_mod, "locStimDensity")
  )
  expect_equal(
    attr(smooth_out, "locDensityComparison"),
    attr(data_mod, "locDensityComparison")
  )
})

test_that(".getCpUnsLocGetProbSmooth handles small cell counts", {
  tmp_dir <- file.path(tempdir(), "test_prob_smooth_small_cells")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  x_vals <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  prob_vals <- c(0.1, 0.2, 0.4, 0.7, 0.9)

  data_mod <- data.frame(
    CD4 = x_vals,
    probSmooth = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "CD4"
  attr(data_mod, "idxMod") <- seq_along(x_vals)
  attr(data_mod, "ind") <- 2L
  attr(data_mod, "locDensityBw") <- 0.4
  attr(data_mod, "locPeakX") <- 2.0
  attr(data_mod, "locWindowWidth") <- 1.5

  # Check cell helper directly
  expect_false(.getCpUnsLocGetProbSmoothCheckNCell(data_mod))
  expect_false(.getCpUnsLocGetProbSmoothCheckNCell(NULL))
  expect_false(.getCpUnsLocGetProbSmoothCheckNCell(data.frame()))
  expect_true(
    .getCpUnsLocGetProbSmoothCheckNCell(
      data.frame(x = seq_len(15), probSmooth = seq_len(15))
    )
  )

  smooth_out <- .getCpUnsLocGetProbSmooth(
    dataMod = data_mod,
    stage = "init",
    pathProject = tmp_dir,
    chnl = "CD4"
  )

  # Check prediction assigned directly from probSmooth - 1e-4
  expect_equal(smooth_out$pred, prob_vals - 1e-4)
  expect_null(attr(smooth_out, "locProbDerivTbl"))
  expect_null(attr(smooth_out, "locProbSmoothMethod"))

  # Retained attributes survive
  expect_equal(attr(smooth_out, "locDensityBw"), 0.4)
  expect_equal(attr(smooth_out, "locPeakX"), 2.0)
  expect_equal(attr(smooth_out, "locWindowWidth"), 1.5)
})

test_that(".getCpUnsLocGetProbSmooth falls back when fit is rejected", {
  tmp_dir <- file.path(tempdir(), "test_prob_smooth_fallback")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  x_vals <- seq(0, 5, length.out = 30)
  # All probabilities close to 1 (> 0.99), triggering fit rejection
  prob_vals <- rep(0.999, length(x_vals))

  data_mod <- data.frame(
    IL2 = x_vals,
    probSmooth = prob_vals,
    stringsAsFactors = FALSE
  )
  attr(data_mod, "chnlCut") <- "IL2"
  attr(data_mod, "idxMod") <- seq_along(x_vals)
  attr(data_mod, "ind") <- 3L
  attr(data_mod, "locDensityBw") <- 0.25

  smooth_out <- .getCpUnsLocGetProbSmooth(
    dataMod = data_mod,
    stage = "init",
    pathProject = tmp_dir,
    chnl = "IL2"
  )

  # Method metadata reflects fallback route
  expect_equal(
    attr(smooth_out, "locProbSmoothMethod"),
    "probSmooth_fallback"
  )
  expect_equal(smooth_out$pred, prob_vals - 0.0001)
  expect_null(attr(smooth_out, "locProbDerivTbl"))
  expect_equal(attr(smooth_out, "locDensityBw"), 0.25)
})

test_that("pure smoothing helpers validate inputs and prediction structures", {
  # .getCpUnsLocGetProbSmoothNewData
  new_data <- .getCpUnsLocGetProbSmoothNewData("TNF", c(0.5, 1.5, 2.5))
  expect_s3_class(new_data, "data.frame")
  expect_named(new_data, "TNF")
  expect_equal(new_data$TNF, c(0.5, 1.5, 2.5))

  # .getCpUnsLocGetProbSmoothCheckNCellOut
  df_in <- data.frame(probSmooth = c(0.2, 0.6))
  df_out <- .getCpUnsLocGetProbSmoothCheckNCellOut(df_in)
  expect_equal(df_out$pred, c(0.2 - 1e-4, 0.6 - 1e-4))
  expect_null(.getCpUnsLocGetProbSmoothCheckNCellOut(NULL))

  # .getCpUnsLocGetProbSmoothObjPred
  expect_equal(
    .getCpUnsLocGetProbSmoothObjPred(list(pred = c(0.1, 0.2))),
    c(0.1, 0.2)
  )
  expect_equal(
    .getCpUnsLocGetProbSmoothObjPred(c(0.3, 0.4)),
    c(0.3, 0.4)
  )

  # .getCpUnsLocGetProbSmoothAttachDeriv
  df_base <- data.frame(x = 1:3)
  deriv_mock <- tibble::tibble(
    x = 1:3,
    pred = c(0.1, 0.2, 0.3),
    deriv = c(0, 0, 0)
  )
  df_attached <- .getCpUnsLocGetProbSmoothAttachDeriv(
    dataMod = df_base,
    smoothObj = list(derivTbl = deriv_mock, method = "scam_mpi")
  )
  expect_equal(attr(df_attached, "locProbDerivTbl"), deriv_mock)
  expect_equal(attr(df_attached, "locProbSmoothMethod"), "scam_mpi")

  # Non-list smoothObj returns unmodified dataMod
  expect_identical(
    .getCpUnsLocGetProbSmoothAttachDeriv(df_base, NULL),
    df_base
  )

  # .getCpUnsLocGetProbSmoothActualThird fallback constructor
  df_third <- data.frame(probSmooth = c(0.1, 0.5, 0.9))
  third_out <- .getCpUnsLocGetProbSmoothActualThird(df_third, stage = "init")
  expect_type(third_out, "list")
  expect_equal(third_out$pred, df_third$probSmooth - 0.0001)
  expect_true(is.na(third_out$meanAbsError))
  expect_null(third_out$derivTbl)
  expect_equal(third_out$method, "probSmooth_fallback")
})

test_that("fit evaluation and acceptance checks enforce quality thresholds", {
  # .getCpUnsLocGetProbSmoothFitEval returns NULL on invalid fit
  expect_null(.getCpUnsLocGetProbSmoothFitEval(NULL, data.frame(x = 1:5)))
  expect_null(
    .getCpUnsLocGetProbSmoothFitEval(
      try(stop("failed fit"), silent = TRUE),
      data.frame(x = 1:5)
    )
  )

  # .getCpUnsLocGetProbSmoothFitEvalCheck handles NULL
  expect_false(.getCpUnsLocGetProbSmoothFitEvalCheck(NULL))

  # Rejects when all predictions are > 0.99
  fit_eval_high <- list(
    pred = rep(0.995, 10),
    meanAbsError = 0.02
  )
  expect_false(.getCpUnsLocGetProbSmoothFitEvalCheck(fit_eval_high))

  # Rejects when meanAbsError exceeds 0.3
  fit_eval_poor <- list(
    pred = seq(0.1, 0.9, length.out = 10),
    meanAbsError = 0.35
  )
  expect_false(.getCpUnsLocGetProbSmoothFitEvalCheck(fit_eval_poor))

  # Accepts valid fit with reasonable predictions and error <= 0.3
  fit_eval_good <- list(
    pred = seq(0.05, 0.85, length.out = 10),
    meanAbsError = 0.08
  )
  expect_true(.getCpUnsLocGetProbSmoothFitEvalCheck(fit_eval_good))
})

test_that(".getCpUnsLocGetProbSmoothDerivativeTbl handles edge cases", {
  # Null/short expression vectors return NULL
  df_short <- data.frame(TNF = c(1.0, 2.0))
  attr(df_short, "chnlCut") <- "TNF"
  expect_null(
    .getCpUnsLocGetProbSmoothDerivativeTbl(
      fit = list(),
      dataMod = df_short,
      chnlSettings = list()
    )
  )

  # Identical expression values (zero range) return NULL
  df_flat_x <- data.frame(TNF = rep(2.0, 10))
  attr(df_flat_x, "chnlCut") <- "TNF"
  expect_null(
    .getCpUnsLocGetProbSmoothDerivativeTbl(
      fit = list(),
      dataMod = df_flat_x,
      chnlSettings = list()
    )
  )

  # Custom locFlatDerivGridN setting with fitted smoother
  x_vals <- seq(0, 4, length.out = 30)
  prob_vals <- 1 / (1 + exp(-(x_vals - 2) * 2))
  data_mod <- data.frame(TNF = x_vals, probSmooth = prob_vals)
  attr(data_mod, "chnlCut") <- "TNF"
  attr(data_mod, "idxMod") <- seq_along(x_vals)

  fit <- pkg_ns$.getCpUnsLocGetProbSmoothActualFirst(data_mod, stage = "init")
  expect_false(inherits(fit, "try-error"))

  deriv_custom <- .getCpUnsLocGetProbSmoothDerivativeTbl(
    fit = fit,
    dataMod = data_mod,
    chnlSettings = list(locFlatDerivGridN = 100L)
  )

  expect_s3_class(deriv_custom, "tbl_df")
  expect_equal(nrow(deriv_custom), 100L)
  expect_true(all(deriv_custom$deriv >= 0))
  expect_true(all(is.finite(deriv_custom$deriv)))
})
