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

.gateInit <- function(
    chnlSettings,
    .data,
    indBatchList,
    pathProject) {
  if (!.profileEnabled()) {
    return(.profileOriginalGateInit(
      chnlSettings = chnlSettings,
      .data = .data,
      indBatchList = indBatchList,
      pathProject = pathProject
    ))
  }

  .profileInit(pathProject)
  .profileWithContext(
    .profileTime(
      .profileOriginalGateInit(
        chnlSettings = chnlSettings,
        .data = .data,
        indBatchList = indBatchList,
        pathProject = pathProject
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
    gateName = NULL,
    calcCytPos = TRUE,
    stage,
    pathProject) {
  if (!.profileEnabled()) {
    return(.profileOriginalGateCytPos(
      chnlSettings = chnlSettings,
      indBatchList = indBatchList,
      .data = .data,
      gateName = gateName,
      calcCytPos = calcCytPos,
      stage = stage,
      pathProject = pathProject
    ))
  }

  .profileWithContext(
    .profileTime(
      .profileOriginalGateCytPos(
        chnlSettings = chnlSettings,
        indBatchList = indBatchList,
        .data = .data,
        gateName = gateName,
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
    pathProject) {
  if (!.profileEnabled()) {
    return(.profileOriginalGateStats(
      .data = .data,
      gateTbl = gateTbl,
      calcCytPosGates = calcCytPosGates,
      chnlSettings = chnlSettings,
      indBatchList = indBatchList,
      pathProject = pathProject
    ))
  }

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
  out
}

.gateChnl <- function(
    .data,
    indBatchList,
    chnlSettings,
    gateTbl = NULL,
    tolGateSingle = NULL,
    calcCytPosGates,
    pathProject,
    stage) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGateChnl(
      .data = .data,
      indBatchList = indBatchList,
      chnlSettings = chnlSettings,
      gateTbl = gateTbl,
      tolGateSingle = tolGateSingle,
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
        gateTbl = gateTbl,
        tolGateSingle = tolGateSingle,
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
    pathProject) {
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
    pathProject) {
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
    gateStatsTbl,
    gateTblCtrl,
    chnlSettings,
    stage,
    pathProject,
    control = list(),
    filterOtherCytPos,
    calcCytPosGates,
    indBatchList) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGetCpCluster(
      .data = .data,
      gateTbl = gateTbl,
      gateStatsTbl = gateStatsTbl,
      gateTblCtrl = gateTblCtrl,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject,
      control = control,
      filterOtherCytPos = filterOtherCytPos,
      calcCytPosGates = calcCytPosGates,
      indBatchList = indBatchList
    ))
  }

  .profileTime(
    .profileOriginalGetCpCluster(
      .data = .data,
      gateTbl = gateTbl,
      gateStatsTbl = gateStatsTbl,
      gateTblCtrl = gateTblCtrl,
      chnlSettings = chnlSettings,
      stage = stage,
      pathProject = pathProject,
      control = control,
      filterOtherCytPos = filterOtherCytPos,
      calcCytPosGates = calcCytPosGates,
      indBatchList = indBatchList
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
    plot = TRUE,
    probMin = 0.1,
    bias,
    pathProject,
    stage) {
  if (!.profileEnabled() || !identical(stage, "init")) {
    return(.profileOriginalGetCpUnsLocCondition(
      exTblUnsBias = exTblUnsBias,
      exTblStimNoMin = exTblStimNoMin,
      chnlSettings = chnlSettings,
      exTblStimOrig = exTblStimOrig,
      exTblUnsOrig = exTblUnsOrig,
      plot = plot,
      probMin = probMin,
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
        plot = plot,
        probMin = probMin,
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
    chnlSettings) {
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
    chnlSettings) {
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

.getCpUnsLocAntimodeDensity <- function(
    expr,
    chnlSettings,
    originalBw = NULL,
    mtd = c("taut_string", "kde")) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocAntimodeDensity(
      expr = expr,
      chnlSettings = chnlSettings,
      originalBw = originalBw,
      mtd = mtd
    ))
  }

  .profileTime(
    .profileOriginalGetCpUnsLocAntimodeDensity(
      expr = expr,
      chnlSettings = chnlSettings,
      originalBw = originalBw,
      mtd = mtd
    ),
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
    chnlSettings = list()) {
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
    shapeLowerBoundX = NA_real_) {
  if (!.profileEnabled() || !.profileInitialSampleActive()) {
    return(.profileOriginalGetCpUnsLocFilterMarginal(
      dataMod = dataMod,
      chnlSettings = chnlSettings,
      probCol = probCol,
      antimodeX = antimodeX,
      threshold = threshold,
      dominance = dominance,
      globalLowerBoundX = globalLowerBoundX,
      shapeLowerBoundX = shapeLowerBoundX
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
      shapeLowerBoundX = shapeLowerBoundX
    ),
    level = "sample_detail",
    major = "initial_gating",
    minor = "local_fdr",
    operation = "marginal_filtering"
  )
}
