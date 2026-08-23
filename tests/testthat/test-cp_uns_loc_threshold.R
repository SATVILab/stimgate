pkg_ns <- asNamespace("stimgate")

.getCpUnsLocTailPropAtThresholds <- get(
  ".getCpUnsLocTailPropAtThresholds",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetCpDataThresholdActual <- get(
  ".getCpUnsLocGetCpDataThresholdActual",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocGetCpActual <- get(
  ".getCpUnsLocGetCpActual",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocThresholdOrigin <- get(
  ".getCpUnsLocThresholdOrigin",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocPropBsAtCp <- get(
  ".getCpUnsLocPropBsAtCp",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocSelectedThresholdRow <- get(
  ".getCpUnsLocSelectedThresholdRow",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocConditionDetailRow <- get(
  ".getCpUnsLocConditionDetailRow",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocSampleDetailTbl <- get(
  ".getCpUnsLocSampleDetailTbl",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocSampleCpRep <- get(
  ".getCpUnsLocSampleCpRep",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocMetaFromCp <- get(
  ".getCpUnsLocMetaFromCp",
  envir = pkg_ns,
  mode = "function"
)
.getCpUnsLocCombineCpWithMeta <- get(
  ".getCpUnsLocCombineCpWithMeta",
  envir = pkg_ns,
  mode = "function"
)

test_that(".getCpUnsLocTailPropAtThresholds handles ties and boundaries", {
  x <- c(1.0, 2.0, 2.0, 3.0, 4.0)
  thresholds <- c(0.5, 1.0, 2.0, 2.5, 4.0, 4.5, NA_real_)

  props <- .getCpUnsLocTailPropAtThresholds(
    x = x,
    thresholds = thresholds,
    denominator = length(x)
  )

  # Expected counts:
  # >= 0.5: 5/5 = 1.0
  # >= 1.0: 5/5 = 1.0 (exact match included)
  # >= 2.0: 4/5 = 0.8 (both tied 2.0 values included)
  # >= 2.5: 2/5 = 0.4 (3.0, 4.0)
  # >= 4.0: 1/5 = 0.2 (exact match included)
  # >= 4.5: 0/5 = 0.0
  # NA: NA
  expect_equal(props[1:6], c(1.0, 1.0, 0.8, 0.4, 0.2, 0.0))
  expect_true(is.na(props[[7]]))

  # Vector with NA returns all NAs
  x_na <- c(1.0, 2.0, NA_real_)
  props_na <- .getCpUnsLocTailPropAtThresholds(
    x = x_na,
    thresholds = c(1.0, 2.0),
    denominator = 3
  )
  expect_true(all(is.na(props_na)))
})

test_that("empirical threshold selection matches response estimate", {
  # Candidate expressions and pred probabilities
  data_count <- data.frame(
    IFNg = c(2.0, 3.0, 4.0),
    probSmooth = c(0.4, 0.7, 0.9),
    pred = c(0.4, 0.7, 0.9)
  )
  attr(data_count, "chnlCut") <- "IFNg"

  # Expression matrices
  df_stim <- data.frame(IFNg = c(1.0, 2.0, 3.0, 3.5, 4.0))
  attr(df_stim, "chnlCut") <- "IFNg"
  df_uns <- data.frame(IFNg = c(0.5, 1.0, 2.0, 2.5, 3.0))
  attr(df_uns, "chnlCut") <- "IFNg"

  # Background-subtracted response estimate = 0.40
  prop_bs_est <- 0.40
  data_thresh <- .getCpUnsLocGetCpDataThresholdActual(
    dataCount = data_count,
    propBsEst = prop_bs_est,
    exTblStimOrig = df_stim,
    exTblUnsBias = df_uns,
    exTblUnsOrig = df_uns,
    bias = 0
  )

  expect_s3_class(data_thresh, "data.frame")
  expect_named(
    data_thresh,
    c(
      "IFNg", "probSmooth", "pred", "propStim", "propUns",
      "propBs", "propBsDiff"
    )
  )

  # For t = 2.0: propStim = 4/5 = 0.8, propUns = 3/5 = 0.6 -> propBs = 0.2
  # For t = 3.0: propStim = 3/5 = 0.6, propUns = 1/5 = 0.2 -> propBs = 0.4
  # For t = 4.0: propStim = 1/5 = 0.2, propUns = 0/5 = 0.0 -> propBs = 0.2
  expect_equal(data_thresh$propStim, c(0.8, 0.6, 0.2))
  expect_equal(data_thresh$propUns, c(0.6, 0.2, 0.0))
  expect_equal(data_thresh$propBs, c(0.2, 0.4, 0.2))
  expect_equal(data_thresh$propBsDiff, c(-0.2, 0.0, -0.2))

  # Selection selects t = 3.0 where propBsDiff is exactly 0.0
  cp_res <- .getCpUnsLocGetCpActual(
    dataThreshold = data_thresh,
    exTblStimNoMin = df_stim,
    exTblUnsBias = df_uns,
    cpMin = NULL,
    stage = "init"
  )

  expect_equal(cp_res$cp, 3.0)
  expect_true(cp_res$locGenerated)
  expect_true(cp_res$locGeneratedDirect)
  expect_equal(cp_res$locSource, "direct")
  expect_equal(cp_res$locReason, "local_fdr_threshold_selected")

  # Empty dataThreshold returns non-local fallback
  empty_thresh <- data.frame(
    IFNg = numeric(0),
    propBsDiff = numeric(0)
  )
  attr(empty_thresh, "chnlCut") <- "IFNg"
  attr(df_stim, "ind") <- 1L

  cp_fallback <- .getCpUnsLocGetCpActual(
    dataThreshold = empty_thresh,
    exTblStimNoMin = df_stim,
    exTblUnsBias = df_uns,
    cpMin = 1.5,
    stage = "init"
  )
  expect_false(cp_fallback$locGenerated)
  expect_false(cp_fallback$locGeneratedDirect)
  expect_equal(cp_fallback$locSource, "not_calculated")
})

test_that(".getCpUnsLocThresholdOrigin labels all provenance pathways", {
  expect_equal(
    .getCpUnsLocThresholdOrigin(TRUE, TRUE, "direct"),
    "condition_detected_response"
  )
  expect_equal(
    .getCpUnsLocThresholdOrigin(TRUE, FALSE, "combined"),
    "sample_imputed_from_other_stim_conditions"
  )
  expect_equal(
    .getCpUnsLocThresholdOrigin(TRUE, FALSE, "prejoin"),
    "prejoin_generated_from_joined_stim_conditions"
  )
  expect_equal(
    .getCpUnsLocThresholdOrigin(TRUE, FALSE, "unstim_summary"),
    "unstim_summary_from_generated_thresholds"
  )
  expect_equal(
    .getCpUnsLocThresholdOrigin(FALSE, FALSE, "not_calculated"),
    "not_generated_fallback"
  )
  expect_equal(
    .getCpUnsLocThresholdOrigin(TRUE, FALSE, "custom", "tailgate_boundary"),
    "generated_tailgate_boundary"
  )
})

test_that(".getCpUnsLocSampleCpRep averages only generated thresholds", {
  tmp_dir <- file.path(tempdir(), "test_sample_cp_rep")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  # Stim 1: generated threshold = 4.0
  # Stim 2: fallback threshold = 1.0 (locGenerated = FALSE)
  obj1 <- list(
    cp = 4.0,
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "selected"
  )
  obj2 <- list(
    cp = 1.0,
    locGenerated = FALSE,
    locGeneratedDirect = FALSE,
    locSource = "not_calculated",
    locReason = "fallback"
  )
  obj_list <- list("stim1" = obj1, "stim2" = obj2)

  rep_out <- .getCpUnsLocSampleCpRep(
    stage = "init",
    cpUnsLocObjList = obj_list,
    indUns = "uns",
    indStim = c("stim1", "stim2"),
    pathProject = tmp_dir,
    chnl = "IFNg"
  )

  # Threshold vector values: stim1 = 4.0, stim2 = 1.0, uns = 4.0
  # Note: uns is the mean of generated thresholds only (NOT (4+1)/2 = 2.5)
  expect_named(rep_out, c("stim1", "stim2", "uns"))
  expect_equal(as.numeric(rep_out["stim1"]), 4.0)
  expect_equal(as.numeric(rep_out["stim2"]), 1.0)
  expect_equal(as.numeric(rep_out["uns"]), 4.0)

  meta <- .getCpUnsLocMetaFromCp(rep_out)
  expect_equal(nrow(meta), 3L)

  # Stim 1 metadata
  expect_true(meta$locGenerated[meta$ind == "stim1"])
  expect_true(meta$locGeneratedDirect[meta$ind == "stim1"])
  expect_equal(meta$locSource[meta$ind == "stim1"], "direct")

  # Stim 2 metadata
  expect_false(meta$locGenerated[meta$ind == "stim2"])
  expect_false(meta$locGeneratedDirect[meta$ind == "stim2"])
  expect_equal(meta$locSource[meta$ind == "stim2"], "not_calculated")

  # Unstim summary metadata
  expect_true(meta$locGenerated[meta$ind == "uns"])
  expect_false(meta$locGeneratedDirect[meta$ind == "uns"])
  expect_equal(meta$locSource[meta$ind == "uns"], "unstim_summary")
  expect_equal(
    meta$locReason[meta$ind == "uns"],
    "mean_of_generated_local_fdr_thresholds"
  )

  # When NO stimulated condition is generated, unstim threshold is NA
  obj_none <- list("stim1" = obj2, "stim2" = obj2)
  rep_none <- .getCpUnsLocSampleCpRep(
    stage = "init",
    cpUnsLocObjList = obj_none,
    indUns = "uns",
    indStim = c("stim1", "stim2"),
    pathProject = tmp_dir,
    chnl = "IFNg"
  )

  expect_true(is.na(rep_none["uns"]))
  meta_none <- .getCpUnsLocMetaFromCp(rep_none)
  expect_false(meta_none$locGenerated[meta_none$ind == "uns"])
  expect_equal(
    meta_none$locReason[meta_none$ind == "uns"],
    "no_generated_local_fdr_thresholds"
  )
})

test_that(".getCpUnsLocCombineCpWithMeta propagates combined metadata", {
  cp_vec <- c("stim1" = 4.0, "stim2" = 1.0, "uns" = 4.0)
  attr(cp_vec, "locGenerated") <- c(TRUE, FALSE, TRUE)
  attr(cp_vec, "locGeneratedDirect") <- c(TRUE, FALSE, FALSE)
  attr(cp_vec, "locSource") <- c("direct", "not_calculated", "unstim_summary")
  attr(cp_vec, "locReason") <- c("selected", "fallback", "mean")

  combined_out <- .getCpUnsLocCombineCpWithMeta(cp_vec, gateCombn = "mean")
  cp_mean <- combined_out[["mean"]]

  # Both stim conditions now have threshold 4.0
  expect_equal(as.numeric(cp_mean["stim1"]), 4.0)
  expect_equal(as.numeric(cp_mean["stim2"]), 4.0)
  expect_equal(as.numeric(cp_mean["uns"]), 4.0)

  meta_comb <- .getCpUnsLocMetaFromCp(cp_mean)
  # stim1 remains direct
  expect_true(meta_comb$locGenerated[meta_comb$ind == "stim1"])
  expect_true(meta_comb$locGeneratedDirect[meta_comb$ind == "stim1"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "stim1"], "direct")

  # stim2 is now imputed from stim1 (combined)
  expect_true(meta_comb$locGenerated[meta_comb$ind == "stim2"])
  expect_false(meta_comb$locGeneratedDirect[meta_comb$ind == "stim2"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "stim2"], "combined")
  expect_equal(
    meta_comb$locReason[meta_comb$ind == "stim2"],
    "combined_from_generated_local_fdr_thresholds"
  )
})

test_that(".getCpUnsLocCombineCpWithMeta uses min of generated thresholds and ignores fallback cutpoints", {
  # Multi-condition batch:
  # stim1 generated 3.5, stim2 generated 2.8, stim3 is fallback at 0.5 (lower than 2.8), uns is 3.15
  cp_vec <- c("stim1" = 3.5, "stim2" = 2.8, "stim3" = 0.5, "uns" = 3.15)
  attr(cp_vec, "locGenerated") <- c(TRUE, TRUE, FALSE, TRUE)
  attr(cp_vec, "locGeneratedDirect") <- c(TRUE, TRUE, FALSE, FALSE)
  attr(cp_vec, "locSource") <- c("direct", "direct", "not_calculated", "unstim_summary")
  attr(cp_vec, "locReason") <- c("detected", "detected", "fallback_above_range", "mean")

  combined_out <- .getCpUnsLocCombineCpWithMeta(cp_vec, gateCombn = "min")
  cp_min <- combined_out[["min"]]

  # All conditions receive min of generated thresholds (2.8), ignoring fallback (0.5)
  expect_equal(as.numeric(cp_min["stim1"]), 2.8)
  expect_equal(as.numeric(cp_min["stim2"]), 2.8)
  expect_equal(as.numeric(cp_min["stim3"]), 2.8)
  expect_equal(as.numeric(cp_min["uns"]), 2.8)

  meta_comb <- .getCpUnsLocMetaFromCp(cp_min)

  # stim2 matches the minimum and remains direct
  expect_true(meta_comb$locGenerated[meta_comb$ind == "stim2"])
  expect_true(meta_comb$locGeneratedDirect[meta_comb$ind == "stim2"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "stim2"], "direct")

  # stim1 was lowered from 3.5 to 2.8 -> combined
  expect_true(meta_comb$locGenerated[meta_comb$ind == "stim1"])
  expect_false(meta_comb$locGeneratedDirect[meta_comb$ind == "stim1"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "stim1"], "combined")
  expect_equal(
    meta_comb$locReason[meta_comb$ind == "stim1"],
    "combined_from_generated_local_fdr_thresholds"
  )

  # stim3 was non-generated fallback -> combined
  expect_true(meta_comb$locGenerated[meta_comb$ind == "stim3"])
  expect_false(meta_comb$locGeneratedDirect[meta_comb$ind == "stim3"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "stim3"], "combined")
  expect_equal(
    meta_comb$locReason[meta_comb$ind == "stim3"],
    "combined_from_generated_local_fdr_thresholds"
  )

  # uns receives unstim_summary
  expect_true(meta_comb$locGenerated[meta_comb$ind == "uns"])
  expect_false(meta_comb$locGeneratedDirect[meta_comb$ind == "uns"])
  expect_equal(meta_comb$locSource[meta_comb$ind == "uns"], "unstim_summary")
  expect_equal(
    meta_comb$locReason[meta_comb$ind == "uns"],
    "summary_of_combined_generated_local_fdr_thresholds"
  )
})

test_that(".getCpUnsLocCombineCpWithMeta handles all-fallback conditions cleanly under min", {
  cp_vec <- c("stim1" = 5.0, "stim2" = 6.0, "uns" = NA_real_)
  attr(cp_vec, "locGenerated") <- c(FALSE, FALSE, FALSE)
  attr(cp_vec, "locGeneratedDirect") <- c(FALSE, FALSE, FALSE)
  attr(cp_vec, "locSource") <- c("not_calculated", "not_calculated", "unstim_summary")
  attr(cp_vec, "locReason") <- c("fallback", "fallback", "no_generated_local_fdr_thresholds")

  combined_out <- .getCpUnsLocCombineCpWithMeta(cp_vec, gateCombn = "min")
  cp_min <- combined_out[["min"]]

  expect_equal(as.numeric(cp_min["stim1"]), 5.0)
  expect_equal(as.numeric(cp_min["stim2"]), 5.0)
  expect_equal(as.numeric(cp_min["uns"]), 5.0)

  meta_comb <- .getCpUnsLocMetaFromCp(cp_min)
  expect_false(any(meta_comb$locGenerated))
  expect_false(any(meta_comb$locGeneratedDirect))
  expect_true(all(meta_comb$locSource == "not_calculated"))
  expect_true(all(meta_comb$locReason == "no_generated_local_fdr_threshold_to_combine"))
})

test_that("diagnostic tables compute frequencies matching threshold cutoffs", {
  df_uns <- data.frame(CD4 = c(0.5, 1.0, 1.5, 2.0, 2.5))
  attr(df_uns, "chnlCut") <- "CD4"
  df_stim <- data.frame(CD4 = c(1.0, 2.0, 3.0, 3.5, 4.0))
  attr(df_stim, "chnlCut") <- "CD4"
  attr(df_stim, "ind") <- "stim1"

  # .getCpUnsLocPropBsAtCp
  freq_tbl <- .getCpUnsLocPropBsAtCp(
    cp = 3.0,
    exTblStim = df_stim,
    exTblUns = df_uns
  )
  expect_equal(freq_tbl$nCellStim, 5L)
  expect_equal(freq_tbl$nCellUns, 5L)
  expect_equal(freq_tbl$propStim, 3 / 5)
  expect_equal(freq_tbl$propUns, 0 / 5)
  expect_equal(freq_tbl$propBs, 3 / 5)

  # .getCpUnsLocConditionDetailRow with direct threshold
  cp_direct <- list(
    cp = 3.0,
    locGenerated = TRUE,
    locGeneratedDirect = TRUE,
    locSource = "direct",
    locReason = "local_fdr_threshold_selected"
  )
  data_thresh_row <- data.frame(
    CD4 = 3.0,
    propStim = 3 / 5,
    propUns = 0 / 5,
    propBs = 3 / 5,
    propBsDiff = 0.05
  )
  attr(data_thresh_row, "chnlCut") <- "CD4"

  detail_cond <- .getCpUnsLocConditionDetailRow(
    cpObj = cp_direct,
    dataThreshold = data_thresh_row,
    exTblStimOrig = df_stim,
    exTblUnsOrig = df_uns,
    exTblStimNoMin = df_stim,
    bias = 0.1,
    stage = "init",
    chnl = "CD4"
  )

  expect_equal(detail_cond$detailLevel, "condition")
  expect_equal(detail_cond$threshold, 3.0)
  expect_equal(detail_cond$thresholdOrigin, "condition_detected_response")
  expect_equal(detail_cond$propBsEst, (3 / 5) - 0.05)
  expect_equal(detail_cond$propBsDiff, 0.05)
  expect_equal(detail_cond$propStim, 3 / 5)
  expect_equal(detail_cond$propUns, 0 / 5)

  # .getCpUnsLocSampleDetailTbl across conditions
  cp_vec <- c("stim1" = 3.0, "uns" = 3.0)
  attr(cp_vec, "locGenerated") <- c(TRUE, TRUE)
  attr(cp_vec, "locGeneratedDirect") <- c(TRUE, FALSE)
  attr(cp_vec, "locSource") <- c("direct", "unstim_summary")
  attr(cp_vec, "locReason") <- c(
    "selected",
    "mean_of_generated_local_fdr_thresholds"
  )

  sample_detail <- .getCpUnsLocSampleDetailTbl(
    cpVec = cp_vec,
    exListOrig = list("uns" = df_uns, "stim1" = df_stim),
    indUns = "uns",
    indStim = "stim1",
    stage = "init",
    chnl = "CD4"
  )

  expect_equal(nrow(sample_detail), 1L)
  expect_equal(sample_detail$detailLevel, "sample")
  expect_equal(sample_detail$ind, "stim1")
  expect_equal(sample_detail$threshold, 3.0)
  expect_equal(sample_detail$thresholdOrigin, "condition_detected_response")
  expect_equal(sample_detail$propStim, 3 / 5)
  expect_equal(sample_detail$propUns, 0 / 5)
  expect_equal(sample_detail$propBs, 3 / 5)
})
