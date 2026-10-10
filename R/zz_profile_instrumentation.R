# Debug profiling instrumentation -------------------------------------------
#
# This file is deliberately prefixed with `zz_` so the implementation
# functions it wraps have already been defined. The wrappers add only debug
# profiling and delegate unchanged to the original implementation otherwise.

.profileOriginalGateInit <- .gateInit
.profileOriginalGateCytPos <- .gateCytPos
.profileOriginalGateStats <- .gateStats
.profileOriginalGateChnl <- .gateChnl
.profileOriginalGateBatch <- .gateBatch
.profileOriginalGetCpUnsLoc <- .getCpUnsLoc
.profileOriginalGetCpCluster <- .getCpCluster
.profileOriginalGetCpUnsLocCondition <- .getCpUnsLocCondition
.profileOriginalGetCpUnsLocGetProb <- .getCpUnsLocGetProb
.profileOriginalGetCpUnsLocGetDensRawDensities <-
  .getCpUnsLocGetDensRawDensities
.profileOriginalGetCpUnsLocAntimodeDensity <- .getCpUnsLocAntimodeDensity
.profileOriginalGetCpUnsLocGetCp <- .getCpUnsLocGetCp
.profileOriginalGetCpUnsLocFilterMarginal <- .getCpUnsLocFilterMarginal
.profileOriginalCompleteChnlSettingsBwShared <- .completeChnlSettingsBwShared

.completeChnlSettingsBwShared <- function(
  chnlSettings,
  indBatchList,
  .data,
  pathProject
) {
  .profileTime(
    .profileOriginalCompleteChnlSettingsBwShared(
      chnlSettings = chnlSettings,
      indBatchList = indBatchList,
      .data = .data,
      pathProject = pathProject
    ),
    level = "major",
    major = "shared_bandwidth",
    operation = paste0("shared_bandwidth_", chnlSettings$bwScope),
    pathProject = pathProject
  )
}

.gateInit <- function(
  chnlSettings,
  .data,
  indBatchList,
  pathProject,
  parallel = FALSE
) {
  .profileInit(pathProject)
  .profileWithContext(
    .profileTime(
      .profileOriginalGateInit(
        chnlSettings = chnlSettings,
        .data = .data,
        indBatchList = indBatchList,
        pathProject = pathProject,
        parallel = parallel
      ),
      level = "major",
      major = "initial_gating",
      operation = "initial_gating",
      pathProject = pathProject
    ),
    stage = "init"
  )
}

.gateCytPos <- function(
  chnlSettings,
  indBatchList,
  .data,
  calcCytPos = TRUE,
  stage,
  pathProject
) {
  .profileWithContext(
    .profileTime(
      .profileOriginalGateCytPos(
        chnlSettings = chnlSettings,
        indBatchList = indBatchList,
        .data = .data,
        calcCytPos = calcCytPos,
        stage = stage,
        pathProject = pathProject
      ),
      level = "major",
      major = "cytokine_positive_gating",
      operation = "cytokine_positive_gating",
      pathProject = pathProject
    ),
    stage = stage
  )
}

.gateStats <- function(
  .data,
  gateTbl = NULL,
  calcCytPosGates,
  chnlSettings,
  indBatchList,
  pathProject
) {
  out <- .profileWithContext(
    .profileTime(
      .profileOriginalGateStats(
        .data = .data,
        gateTbl = gateTbl,
        calcCytPosGates = calcCytPosGates,
        chnlSettings = chnlSettings,
        indBatchList = indBatchList,
        pathProject = pathProject
      ),
      level = "major",
      major = "final_statistics",
      operation = "final_statistics",
      pathProject = pathProject
    ),
    stage = "stats"
  )
  if (.profileEnabled()) out else invisible(out)
}

.gateChnl <- function(
  .data,
  indBatchList,
  chnlSettings,
  calcCytPosGates,
  pathProject,
  stage
) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGateChnl(
      .data = .data,
      indBatchList = indBatchList,
      chnlSettings = chnlSettings,
      calcCytPosGates = calcCytPosGates,
      pathProject = pathProject,
      stage = stage
    ))
  }

  marker <- chnlSettings$marker %||% chnlSettings$chnlCut
  .profileWithContext(
    .profileTime(
      .profileOriginalGateChnl(
        .data = .data,
        indBatchList = indBatchList,
        chnlSettings = chnlSettings,
        calcCytPosGates = calcCytPosGates,
        pathProject = pathProject,
        stage = stage
      ),
      level = "minor",
      major = "initial_gating",
      minor = "marker",
      operation = "marker_total",
      pathProject = pathProject
    ),
    marker = marker,
    channel = chnlSettings$chnlCut,
    stage = stage
  )
}

.gateBatch <- function(
  .data,
  indBatch,
  chnlSettings,
  batch,
  stage,
  pathProject
) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGateBatch(
      .data = .data,
      indBatch = indBatch,
      chnlSettings = chnlSettings,
      batch = batch,
      stage = stage,
      pathProject = pathProject
    ))
  }

  .profileWithContext(
    .profileTime(
      .profileOriginalGateBatch(
        .data = .data,
        indBatch = indBatch,
        chnlSettings = chnlSettings,
        batch = batch,
        stage = stage,
        pathProject = pathProject
      ),
      level = "minor",
      major = "initial_gating",
      minor = "batch",
      operation = "batch_total",
      pathProject = pathProject
    ),
    batch = batch,
    stage = stage
  )
}

.getCpUnsLoc <- function(
  exList,
  .data,
  chnlSettings,
  stage,
  pathProject
) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGetCpUnsLoc(
      exList = exList,
      .data = .data,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLoc(
      exList = exList,
      .data = .data,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject
    ),
    level = "minor",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "local_fdr_initial_gating",
    pathProject = pathProject
  )
}

.getCpCluster <- function(
  .data,
  gateTbl,
  chnlSettings,
  stage,
  pathProject,
  control = list(),
  filterOtherCytPos,
  calcCytPosGates,
  indBatchList,
  exLookup = NULL
) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGetCpCluster(
      .data = .data,
      gateTbl = gateTbl,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject,
      control = control,
      filterOtherCytPos = filterOtherCytPos,
      calcCytPosGates = calcCytPosGates,
      indBatchList = indBatchList,
      exLookup = exLookup
    ))
  }

  .profileTime(
    .profileOriginalGetCpCluster(
      .data = .data,
      gateTbl = gateTbl,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject,
      control = control,
      filterOtherCytPos = filterOtherCytPos,
      calcCytPosGates = calcCytPosGates,
      indBatchList = indBatchList,
      exLookup = exLookup
    ),
    level = "minor",
    major = "initial_gating",
    minor = "cross_sample",
    operation = "cluster_refinement",
    pathProject = pathProject
  )
}

.getCpUnsLocCondition <- function(
  exTblUnsBias,
  exTblStimNoMin,
  chnlSettings,
  exTblStimOrig,
  exTblUnsOrig,
  bias,
  pathProject,
  stage
) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGetCpUnsLocCondition(
      exTblUnsBias = exTblUnsBias,
      exTblStimNoMin = exTblStimNoMin,
      chnlSettings = chnlSettings,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig,
      bias = bias,
      pathProject = pathProject,
      stage = stage
    ))
  }

  batch <- .profileDataBatch(exTblStimNoMin)
  if (is.na(batch)) {
    batch <- NULL
  }
  sample <- as.character(.getInd(exTblStimNoMin))
  marker <- chnlSettings$marker %||% chnlSettings$chnlCut

  .profileWithContext(
    .profileTime(
      .profileOriginalGetCpUnsLocCondition(
        exTblUnsBias = exTblUnsBias,
        exTblStimNoMin = exTblStimNoMin,
        chnlSettings = chnlSettings,
        exTblStimOrig = exTblStimOrig,
        exTblUnsOrig = exTblUnsOrig,
        bias = bias,
        pathProject = pathProject,
        stage = stage
      ),
      level = "sample",
      major = "initial_gating",
      minor = "local_fdr",
      operation = "sample_initial_gating",
      pathProject = pathProject
    ),
    marker = marker,
    channel = chnlSettings$chnlCut,
    batch = batch,
    sample = sample,
    stage = stage
  )
}

.getCpUnsLocGetProb <- function(
  exTblStimNoMin,
  exTblStimThreshold,
  exTblUnsThreshold,
  exTblUnsBias,
  bias,
  exTblUnsOrig,
  stage,
  pathProject,
  chnlSettings
) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocGetProb(
      exTblStimNoMin = exTblStimNoMin,
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      exTblUnsBias = exTblUnsBias,
      bias = bias,
      exTblUnsOrig = exTblUnsOrig,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocGetProb(
      exTblStimNoMin = exTblStimNoMin,
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      exTblUnsBias = exTblUnsBias,
      bias = bias,
      exTblUnsOrig = exTblUnsOrig,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ),
    level = "sample_minor",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "probability_model",
    pathProject = pathProject
  )
}

.getCpUnsLocGetDensRawDensities <- function(
  exTblStimThreshold,
  exTblUnsThreshold,
  stage,
  pathProject,
  chnlSettings
) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocGetDensRawDensities(
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocGetDensRawDensities(
      exTblStimThreshold = exTblStimThreshold,
      exTblUnsThreshold = exTblUnsThreshold,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ),
    level = "sample_detail",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "density_bandwidth",
    pathProject = pathProject
  )
}

.getCpUnsLocAntimodeDensity <- function(expr) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocAntimodeDensity(expr = expr))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocAntimodeDensity(expr = expr),
    level = "sample_detail",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "antimode"
  )
}

.getCpUnsLocGetCp <- function(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsOrig,
  exTblUnsBias,
  bias,
  cpMin,
  stage,
  pathProject,
  chnlSettings = list()
) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocGetCp(
      dataMod = dataMod,
      exTblStimOrig = exTblStimOrig,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsOrig = exTblUnsOrig,
      exTblUnsBias = exTblUnsBias,
      bias = bias,
      cpMin = cpMin,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocGetCp(
      dataMod = dataMod,
      exTblStimOrig = exTblStimOrig,
      exTblStimNoMin = exTblStimNoMin,
      exTblUnsOrig = exTblUnsOrig,
      exTblUnsBias = exTblUnsBias,
      bias = bias,
      cpMin = cpMin,
      stage = stage,
      pathProject = pathProject,
      chnlSettings = chnlSettings
    ),
    level = "sample_minor",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "filtering_threshold",
    pathProject = pathProject
  )
}

.getCpUnsLocFilterMarginal <- function(
  dataMod,
  chnlSettings,
  probCol,
  antimodeX = NULL,
  threshold = NULL,
  dominance = NULL,
  globalLowerBoundX = NA_real_,
  shapeLowerBoundX = NA_real_,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL
) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocFilterMarginal(
      dataMod = dataMod,
      chnlSettings = chnlSettings,
      probCol = probCol,
      antimodeX = antimodeX,
      threshold = threshold,
      dominance = dominance,
      globalLowerBoundX = globalLowerBoundX,
      shapeLowerBoundX = shapeLowerBoundX,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocFilterMarginal(
      dataMod = dataMod,
      chnlSettings = chnlSettings,
      probCol = probCol,
      antimodeX = antimodeX,
      threshold = threshold,
      dominance = dominance,
      globalLowerBoundX = globalLowerBoundX,
      shapeLowerBoundX = shapeLowerBoundX,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig
    ),
    level = "sample_detail",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "marginal_filtering"
  )
}
