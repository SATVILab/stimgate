.simBandwidthBsFreq <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  probExact = FALSE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  samplePerturbationSd,
  conditionPerturbationSd,
  clusterPerturbationSd,
  backgroundRelativeToResponse,
  ncellUnsRelativeToStim,
  covEvMin = 1,
  covEvMax = 2,
  tolClust = 1e-7
) {
  purrr::map_df(seq_len(nIter), function(iterNum) {
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simMiscGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet") # Fixed case
    gs <- flowWorkspace::GatingSet(fs)
    labelsList <- outListExperiment[["labelsList"]] # Fixed case

    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-bw",
      paste0("iter-", iterNum, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    )
    on.exit({
      if (dir.exists(pathProject)) {
        unlink(pathProject, recursive = TRUE)
      }
    })
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE)

    batchList <- lapply(seq(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    # gate
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      bw = bw,
      tolClust = tolClust
    ))

    stopifnot(file.exists(file.path(pathProject, "gateStats.rds")))

    propBsTblTruth <- purrr::map_df(
      seq_len(nSample),
      function(sampleCurr) {
        labelsVec <- (sampleCurr - 1) * nCondition + seq(1, nCondition)
        purrr::map_df(labelsVec, function(i) {
          labelVec <- labelsList[[i]]
          stimCondition <- ifelse(i %% nCondition == 1, "unstim", "stim")
          f1p <- sum(grepl("^gp$", labelVec)) / length(labelVec)
          tibble::tibble(
            ind = sampleCurr |> as.character(),
            stim = stimCondition,
            chnl = "F1",
            F1 = f1p
          )
        })
      }
    ) |>
      tidyr::pivot_wider(
        names_from = stim,
        values_from = F1
      ) |>
      dplyr::mutate(
        propRespTruth = stim - unstim
      ) |>
      dplyr::rename(
        sample = ind
      ) |>
      dplyr::select(
        all_of(c("chnl", "sample", "propRespTruth"))
      )

    pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
    chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
    chnl <- chnlVec[[1]]
    ind <- 2

    propBsTblEstSmooth <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- seq(2, nSample * nCondition, by = 2) |>
        as.character()
      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", ind)
        probSmooth <- file.path(pathInd, "dataMod.rds") |>
          readRDS()
        if (!inherits(probSmooth, "data.frame")) {
          # it's zero if it doesn't inherit anything
          return(tibble::tibble(
            sample = as.character(as.numeric(ind) / nCondition),
            ind = ind,
            chnl = chnl,
            ncellCondition = nCellStim,
            ncellRespSmooth = 0,
            ncellRespPred = 0,
            propRespSmooth = 0,
            propRespPred = 0
          ))
        }
        ncellRespSmooth <- if (nrow(probSmooth) > 0L) {
          sum(probSmooth$probSmooth)
        } else {
          0
        }
        ncellRespPred <- if (nrow(probSmooth) > 0L) {
          sum(probSmooth$pred)
        } else {
          0
        }
        tibble::tibble(
          sample = as.character(as.numeric(ind) / nCondition),
          ind = ind,
          chnl = chnl,
          ncellCondition = nCellStim,
          ncellRespSmooth = ncellRespSmooth,
          ncellRespPred = ncellRespPred,
          propRespSmooth = ncellRespSmooth / ncellCondition,
          propRespPred = ncellRespPred / ncellCondition
        )
      })
    }) |>
      dplyr::select(
        sample,
        chnl,
        propRespSmooth,
        propRespPred
      ) |>
      tidyr::pivot_longer(
        cols = c(propRespSmooth, propRespPred),
        names_to = "method",
        values_to = "propRespEst"
      ) |>
      dplyr::mutate(
        threshold = NA_real_,
        thresholdOrigin = NA_character_,
        locGenerated = NA,
        locGeneratedDirect = NA,
        locSource = NA_character_,
        locReason = NA_character_,
        detailLevel = NA_character_
      )

    propBsTblDetailed <- try(
      getStimGatesDetailed(pathProject),
      silent = TRUE
    )
    propBsTblDetailed <- if (
      inherits(propBsTblDetailed, "try-error") ||
        !is.data.frame(propBsTblDetailed) ||
        nrow(propBsTblDetailed) == 0L
    ) {
      tibble::tibble()
    } else {
      propBsTblDetailed |>
        dplyr::filter(
          .data$detailLevel %in% c("condition", "sample", "cluster_final")
        ) |>
        dplyr::filter(is.finite(.data$propBs)) |>
        dplyr::mutate(
          sample = as.character(as.numeric(.data$ind) / nCondition),
          method = paste0("loc_", .data$detailLevel),
          propRespEst = .data$propBs
        ) |>
        dplyr::select(
          sample,
          chnl,
          method,
          propRespEst,
          threshold,
          thresholdOrigin,
          locGenerated,
          locGeneratedDirect,
          locSource,
          locReason,
          detailLevel
        )
    }

    propBsTblEst <- dplyr::bind_rows(
      propBsTblEstSmooth,
      propBsTblDetailed
    )

    comparisonTbl <- propBsTblTruth |>
      dplyr::left_join(
        propBsTblEst,
        by = c("sample", "chnl")
      )
    comparisonTbl |>
      dplyr::mutate(
        iter = iterNum
      ) |>
      dplyr::select(iter, chnl, dplyr::everything())
  })
}

.simBandwidthEstBw <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  probExact = FALSE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  samplePerturbationSd,
  conditionPerturbationSd,
  clusterPerturbationSd,
  backgroundRelativeToResponse,
  ncellUnsRelativeToStim,
  covEvMin = 1,
  covEvMax = 2,
  tolClust = 1e-7
) {
  purrr::map_df(seq_len(nIter), function(iterNum) {
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simMiscGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet") # Fixed case
    gs <- flowWorkspace::GatingSet(fs)
    labelsList <- outListExperiment[["labelsList"]] # Fixed case

    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-bw",
      paste0("iter-", iterNum, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    )
    on.exit({
      if (dir.exists(pathProject)) {
        unlink(pathProject, recursive = TRUE)
      }
    })
    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE)

    batchList <- lapply(seq(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    # gate
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      bwMtd = bwMtd,
      tolClust = tolClust
    ))

    stopifnot(file.exists(file.path(pathProject, "gateStats.rds")))

    pathDirIntInit <- file.path(pathProject, "intermediateData", "init")
    chnlVec <- list.dirs(pathDirIntInit, full.names = FALSE, recursive = FALSE)
    chnl <- chnlVec[[1]]
    ind <- 2

    bwTbl <- purrr::map_df(chnlVec, function(chnl) {
      indVecStim <- seq(2, nSample * nCondition, by = 2) |>
        as.character()
      purrr::map_df(indVecStim, function(ind) {
        pathInd <- file.path(pathDirIntInit, chnl, "ind", ind)
        bwEst <- file.path(pathInd, "bwCpUnsLoc.rds") |>
          readRDS()
        if (!is.numeric(bwEst)) {
          # it's zero if it doesn't inherit anything
          return(tibble::tibble(
            sample = as.character(as.numeric(ind) / nCondition),
            ind = ind,
            chnl = chnl,
            bw = bwEst
          ))
        }

        tibble::tibble(
          sample = as.character(as.numeric(ind) / nCondition),
          ind = ind,
          chnl = chnl,
          bw = bwEst
        )
      })
    }) |>
      dplyr::mutate(
        iter = iterNum
      ) |>
      dplyr::select(iter, chnl, dplyr::everything())
  })
}

.simPlotGate <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nIter,
  biasUns,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = NULL,
  bwMax = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL,
  probExact = FALSE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  samplePerturbationSd,
  conditionPerturbationSd,
  clusterPerturbationSd,
  backgroundRelativeToResponse,
  ncellUnsRelativeToStim,
  covEvMin = 1,
  covEvMax = 2,
  tolClust = 1e-7,
  markerToPlot = "MarkerF1" # Specify univariate marker here
) {
  allPlots <- lapply(seq_len(nIter), function(iterNum) {
    # 1. Setup Simulation Parameters
    nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
    nCellByCondition <- c(nCellUns, nCellStim)
    transformationFunc <- .simMiscGetTrans(transformation)
    meanExprMat <- matrix(
      c(0, meanPos),
      byrow = TRUE,
      ncol = 1
    )
    clusterLabelVec <- c("gn", "gp")
    probResponseUns <- probResponse * backgroundRelativeToResponse
    probVecUns <- c(1 - probResponseUns, probResponseUns)
    probResponseVecByStimCondition <- list(c(-probResponse, probResponse))

    # 2. Simulate Cytometry Experiment
    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      samplePerturbationSd = samplePerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    fs <- as(outListExperiment[["flowFrameList"]], "flowSet")
    gs <- flowWorkspace::GatingSet(fs)

    # 3. Setup Project Directory for gating
    pathProject <- file.path(
      tempdir(),
      "stimgate",
      "sim-plot",
      paste0("iter-", iterNum, "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    )

    on.exit(
      {
        if (dir.exists(pathProject)) {
          unlink(pathProject, recursive = TRUE)
        }
      },
      add = TRUE
    )

    if (dir.exists(pathProject)) {
      unlink(pathProject, recursive = TRUE)
    }
    dir.create(pathProject, recursive = TRUE)

    batchList <- lapply(seq(nSample), function(i) {
      seq((i - 1) * nCondition + 1, i * nCondition)
    })

    # 4. Execute Gating
    Sys.setenv("STIMGATE_INTERMEDIATE" = "TRUE")
    invisible(gateStim(
      .data = gs,
      pathProject = pathProject,
      popGate = "root",
      batchList = batchList,
      marker = paste0("MarkerF", seq_len(nMarker)),
      bw = bw,
      biasUns = biasUns,
      bwMtd = bwMtd,
      bwMin = bwMin,
      bwMax = bwMax,
      bwAdj = bwAdj,
      bwNcellMin = bwNcellMin,
      bwNcellMax = bwNcellMax,
      bwCluster = bwCluster,
      tolClust = tolClust
    ))

    # 5. Extract plots for each sample
    samplePlots <- lapply(seq_along(batchList), function(sampleIdx) {
      # batchList[[sampleIdx]] grabs both unstim and stim conditions
      indCurr <- batchList[[sampleIdx]]

      # Assuming condition 1 is unstim and condition 2 is stim based on setup
      conditionLabels <- c("Unstimulated", "Stimulated")
      names(conditionLabels) <- indCurr

      # Generate the overlaid plot
      pList <- plotStim(
        ind = indCurr,
        .data = gs,
        pathProject = pathProject,
        marker = markerToPlot,
        indLab = conditionLabels,
        grid = FALSE, # Return raw list of ggplots instead of rendering a cowplot grid
        showGate = TRUE, # Overlay gates
        excMin = TRUE
      )

      # plotStim returns a list (length 1 if univariate). Extract it.
      if (!is.null(pList) && length(pList) > 0) {
        p <- pList[[1]]
        # Update title to track iteration and sample
        p <- p +
          ggplot2::ggtitle(
            sprintf(
              "Iter: %d | Sample: %d | %s",
              iterNum,
              sampleIdx,
              markerToPlot
            )
          )
        return(p)
      } else {
        return(NULL)
      }
    })

    return(samplePlots)
  })

  # Flatten the nested lists so you get a single continuous list of ggplot objects
  do.call(c, allPlots)
}

#' Estimate bandwidths directly from simulated data, without running gateStim
#'
#' This mirrors the bandwidth-estimation part of cp_uns_loc:
#'   1. simulate unstim/stim data
#'   2. optionally exclude the minimum values, as gateStim does by default
#'   3. optionally cap values at the same max-density x used by cp_uns_loc
#'   4. estimate bw separately for stim and unstim
#'   5. return min(bw_stim, bw_uns)
#'
#' @keywords internal
.simBandwidthEstBwDirect <- function(
  nSample = 10L,
  nMarker = 1L,
  nCondition = 2L,
  nCluster = 2L,
  nIter = 10L,
  biasUns = 0.05,
  bw = NULL,
  bwMtd = "hpi1",
  bwMin = 1e-10,
  bwMax = 1e10,
  bwFallback = NULL,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  bwCluster = NULL, # retained only for signature compatibility
  tolClust = NULL, # retained only for signature compatibility
  probExact = TRUE,
  nCellStim,
  probResponse,
  meanPos,
  transformation,
  backgroundRelativeToResponse = 0.2,
  ncellUnsRelativeToStim = 1,
  covEvMin = 2,
  covEvMax = 2,
  excMin = TRUE,
  capStimRange = TRUE,
  summarise = TRUE
) {
  if (!identical(as.integer(nMarker), 1L)) {
    stop("This helper currently expects nMarker = 1.")
  }
  if (!identical(as.integer(nCondition), 2L)) {
    stop("This helper currently expects nCondition = 2.")
  }
  if (!identical(as.integer(nCluster), 2L)) {
    stop("This helper currently expects nCluster = 2.")
  }

  nCellUns <- round(nCellStim * ncellUnsRelativeToStim)
  nCellByCondition <- c(nCellUns, nCellStim)

  transformationFunc <- .simBandwidthGetTrans(transformation)

  meanExprMat <- matrix(
    c(0, meanPos),
    byrow = TRUE,
    ncol = 1
  )

  clusterLabelVec <- c("gn", "gp")

  probResponseUns <- probResponse * backgroundRelativeToResponse
  probVecUns <- c(1 - probResponseUns, probResponseUns)

  probResponseVecByStimCondition <- list(
    c(-probResponse, probResponse)
  )

  raw_tbl <- purrr::map_dfr(seq_len(nIter), function(iterNum) {
    outListExperiment <- simCytExperiment(
      nSample = nSample,
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = "gaussianOnly",
      meanExprMat = meanExprMat,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      samplePerturbationSd = 0,
      conditionPerturbationSd = 0,
      clusterPerturbationSd = 0,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    flowFrameList <- outListExperiment[["flowFrameList"]]

    purrr::map_dfr(seq_len(nSample), function(sampleCurr) {
      indUns <- (sampleCurr - 1L) * nCondition + 1L
      indStim <- indUns + 1L

      x_uns <- as.numeric(flowCore::exprs(flowFrameList[[indUns]])[, "F1"])
      x_stim <- as.numeric(flowCore::exprs(flowFrameList[[indStim]])[, "F1"])

      # cp_uns_loc applies bias to the unstim expression before density work.
      x_uns <- x_uns + (biasUns %||% 0)

      if (excMin) {
        x_uns <- .simBandwidthExcMin(x_uns)
        x_stim <- .simBandwidthExcMin(x_stim)
      }

      x_cap <- .simBandwidthCapForCpUnsLoc(
        x_stim = x_stim,
        x_uns = x_uns,
        capStimRange = capStimRange
      )

      bw_stim <- .simBandwidthBwOne(
        x = x_cap$x_stim,
        bwMtd = bwMtd,
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwFallback = bwFallback
      )

      bw_stim <- if (bw_stim == bwFallback) {
        NA_real_
      } else {
        bw_stim
      }

      bw_uns <- .simBandwidthBwOne(
        x = x_cap$x_uns,
        bwMtd = bwMtd,
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        bwFallback = bwFallback
      )

      bw_uns <- if (bw_uns == bwFallback) {
        NA_real_
      } else {
        bw_uns
      }

      bw_final <- if (!is.null(bw)) {
        bw
      } else {
        bw_vec <- c(bw_stim, bw_uns)
        bw_vec <- bw_vec[!is.na(bw_vec)]
        if (length(bw_vec) == 0) {
          NA_real_
        } else {
          min(bw_vec)
        }
      }

      tibble::tibble(
        transformation = transformation,
        prob_response = probResponse,
        n_cell = nCellStim,
        mean_pos = meanPos,
        bw_mtd = bwMtd,
        iter = iterNum,
        sample = as.character(sampleCurr),
        ind = as.character(indStim),
        chnl = "F1",
        n_cell_uns = nCellUns,
        n_cell_stim = nCellStim,
        n_uns_bw = length(x_cap$x_uns),
        n_stim_bw = length(x_cap$x_stim),
        max_dens_x = x_cap$max_dens_x,
        bw_uns = bw_uns,
        bw_stim = bw_stim,
        bw = bw_final,
        bw_source = dplyr::case_when(
          !is.null(bw) ~ "fixed",
          is.finite(bw_stim) & is.finite(bw_uns) & bw_stim <= bw_uns ~ "stim",
          is.finite(bw_stim) & is.finite(bw_uns) ~ "unstim",
          TRUE ~ NA_character_
        )
      )
    })
  })

  if (!summarise) {
    return(raw_tbl)
  }

  .simBandwidthSummariseBw(raw_tbl)
}


#' Run the direct bandwidth estimator over a Simulation-Bandwidth sim_grid
#'
#' @keywords internal
.simBandwidthEstBwDirectGrid <- function(
  sim_grid,
  nSample = 10L,
  nIter = 10L,
  biasUns = 0.05,
  bwMin = 1e-10,
  bwMax = 1e10,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  probExact = TRUE,
  backgroundRelativeToResponse = 0.2,
  ncellUnsRelativeToStim = 1,
  covEvMin = 2,
  covEvMax = 2,
  excMin = TRUE,
  capStimRange = TRUE,
  summarise = TRUE
) {
  raw_tbl <- purrr::map_dfr(seq_len(nrow(sim_grid)), function(i) {
    row <- sim_grid[i, , drop = FALSE]

    out <- tryCatch(
      .simBandwidthEstBwDirect(
        nSample = nSample,
        nIter = nIter,
        biasUns = biasUns,
        bwMtd = row$bw_mtd[[1]],
        bwMin = bwMin,
        bwMax = bwMax,
        bwAdj = bwAdj,
        bwNcellMin = bwNcellMin,
        bwNcellMax = bwNcellMax,
        probExact = probExact,
        nCellStim = row$n_cell[[1]],
        probResponse = row$prob_response[[1]],
        meanPos = row$mean_pos[[1]],
        transformation = row$transformation[[1]],
        backgroundRelativeToResponse = backgroundRelativeToResponse,
        ncellUnsRelativeToStim = ncellUnsRelativeToStim,
        covEvMin = covEvMin,
        covEvMax = covEvMax,
        excMin = excMin,
        capStimRange = capStimRange,
        summarise = FALSE
      ),
      error = function(e) {
        tibble::tibble(
          iter = NA_integer_,
          sample = NA_character_,
          ind = NA_character_,
          chnl = "F1",
          n_cell_uns = NA_real_,
          n_cell_stim = row$n_cell[[1]],
          n_uns_bw = NA_integer_,
          n_stim_bw = NA_integer_,
          max_dens_x = NA_real_,
          bw_uns = NA_real_,
          bw_stim = NA_real_,
          bw = NA_real_,
          bw_source = NA_character_,
          error = e$message
        )
      }
    )

    # Avoid duplicating scenario columns if the direct helper already added them.
    out <- out |>
      dplyr::select(
        -dplyr::any_of(c(
          "transformation",
          "prob_response",
          "n_cell",
          "mean_pos",
          "bw_mtd"
        ))
      )

    dplyr::bind_cols(
      row[rep(1L, nrow(out)), , drop = FALSE],
      out
    )
  })

  if (!summarise) {
    return(raw_tbl)
  }

  .simBandwidthSummariseBw(raw_tbl)
}


#' Summarise raw bandwidth estimates by simulation scenario
#'
#' @keywords internal
.simBandwidthSummariseBw <- function(.data) {
  group_vars <- intersect(
    c(
      "transformation",
      "prob_response",
      "n_cell",
      "mean_pos",
      "bw_mtd",
      "chnl"
    ),
    names(.data)
  )

  .data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_est = sum(is.finite(.data$bw)),
      n_error = sum(!is.na(.data$error %||% NA_character_)),
      bw_mean = mean(.data$bw, na.rm = TRUE),
      bw_median = stats::median(.data$bw, na.rm = TRUE),
      bw_q05 = stats::quantile(.data$bw, 0.05, na.rm = TRUE),
      bw_q25 = stats::quantile(.data$bw, 0.25, na.rm = TRUE),
      bw_q75 = stats::quantile(.data$bw, 0.75, na.rm = TRUE),
      bw_q95 = stats::quantile(.data$bw, 0.95, na.rm = TRUE),
      bw_min = min(.data$bw, na.rm = TRUE),
      bw_max = max(.data$bw, na.rm = TRUE),
      bw_uns_median = stats::median(.data$bw_uns, na.rm = TRUE),
      bw_stim_median = stats::median(.data$bw_stim, na.rm = TRUE),
      n_source_uns = sum(.data$bw_source == "unstim", na.rm = TRUE),
      n_source_stim = sum(.data$bw_source == "stim", na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$n_cell,
      dplyr::desc(.data$prob_response),
      .data$transformation,
      .data$mean_pos,
      .data$bw_mtd
    )
}


#' Bandwidth estimate for one vector
#'
#' This is the vector-only equivalent of
#' .getCpUnsLocGetDensRawDensitiesBwInit(), but with nrd0 and SJ returning
#' the numeric density bandwidth via $bw.
#'
#' @keywords internal
.simBandwidthBwOne <- function(
  x,
  bwMtd,
  bwMin,
  bwMax,
  bwAdj,
  bwNcellMin,
  bwNcellMax,
  bwFallback
) {
  x <- x[is.finite(x)]

  if (length(x) < 2L || length(unique(x)) < 2L) {
    return(bwMin %||% NA_real_)
  }

  if (!is.null(bwNcellMin) && is.finite(bwNcellMin) && length(x) < bwNcellMin) {
    iqrX <- diff(stats::quantile(x, c(0.75, 0.25), na.rm = TRUE))
    sdX <- abs(iqrX) / 1.5

    if (!is.finite(sdX) || sdX <= 0) {
      sdX <- stats::sd(x, na.rm = TRUE)
    }
    if (!is.finite(sdX) || sdX <= 0) {
      sdX <- .Machine$double.eps
    }

    x <- sample(x, replace = TRUE, size = bwNcellMin) +
      stats::rnorm(bwNcellMin, mean = 0, sd = sdX / 10)
  }

  if (!is.null(bwNcellMax) && is.finite(bwNcellMax) && length(x) > bwNcellMax) {
    x <- sample(x, size = bwNcellMax, replace = FALSE)
  }

  bw_calc <- tryCatch(
    switch(
      bwMtd,
      "nrd0" = stats::density(x, bw = "nrd0")$bw,
      "sj" = stats::density(x, bw = "SJ")$bw,
      "hpi0" = ks::hpi(x = x, deriv.order = 0),
      "hpi1" = ks::hpi(x = x, deriv.order = 1),
      "hpi2" = ks::hpi(x = x, deriv.order = 2),
      "hpi3" = ks::hpi(x = x, deriv.order = 3),
      stop("Unrecognised bwMtd: ", bwMtd)
    ),
    error = function(e) NA_real_
  )

  bw_calc <- as.numeric(bw_calc)[1]

  if (!is.finite(bw_calc)) {
    if (is.null(bwFallback)) {
      stop("Bandwidth calculation failed and no bwFallback provided.")
    }
    return(bwFallback %||% NA_real_)
  }

  bw_calc <- bw_calc * bwAdj

  if (!is.null(bwMin)) {
    bw_calc <- max(bwMin, bw_calc)
  }
  if (!is.null(bwMax)) {
    bw_calc <- min(bwMax, bw_calc)
  }

  bw_calc
}


#' @keywords internal
.simBandwidthExcMin <- function(x) {
  x <- x[is.finite(x)]
  x[x > min(x, na.rm = TRUE)]
}


#' @keywords internal
.simBandwidthCapForCpUnsLoc <- function(
  x_stim,
  x_uns,
  capStimRange
) {
  if (!capStimRange || length(x_stim) < 2L) {
    return(list(
      x_stim = x_stim,
      x_uns = x_uns,
      max_dens_x = NA_real_
    ))
  }

  range_stim <- range(x_stim, na.rm = TRUE)
  max_dens_x <- max(x_stim, na.rm = TRUE) - 0.05 * diff(range_stim)

  list(
    x_stim = pmin(x_stim, max_dens_x),
    x_uns = pmin(x_uns, max_dens_x),
    max_dens_x = max_dens_x
  )
}


#' @keywords internal
.simBandwidthGetTrans <- function(transformation) {
  if (exists(".simMiscGetTrans", mode = "function")) {
    return(.simMiscGetTrans(transformation))
  }

  switch(
    transformation,
    "gamma" = calc_gamma,
    "gaussian" = calc_gaussian,
    "skew" = calc_skew,
    stop("Transformation not recognised: ", transformation)
  )
}
