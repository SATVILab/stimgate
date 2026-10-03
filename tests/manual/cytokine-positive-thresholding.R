devtools::load_all()

baselineGates <- tibble::tribble(
  ~chnl, ~marker, ~gate,
  "IFNG", "IFN-gamma", 5,
  "IL2", "IL-2", 1
)

makeCytPosExample <- function(
  conditionalShape = c("bimodal", "unimodal"),
  nOtherCytPos = 120L
) {
  conditionalShape <- match.arg(conditionalShape)
  stopifnot(nOtherCytPos >= 2L)

  nBackground <- 180L
  nLeft <- floor(nOtherCytPos / 2)

  ifngConditional <- if (conditionalShape == "bimodal") {
    c(
      seq(1, 2, length.out = nLeft),
      seq(5.5, 6.5, length.out = nOtherCytPos - nLeft)
    )
  } else {
    seq(5.5, 6.5, length.out = nOtherCytPos)
  }

  tibble::tibble(
    IFNG = c(
      seq(-2, 0, length.out = nBackground),
      ifngConditional
    ),
    IL2 = c(
      seq(0, 0.8, length.out = nBackground),
      seq(1.2, 2, length.out = nOtherCytPos)
    )
  )
}

inspectCytPosCase <- function(
  ex,
  label,
  gateTbl = baselineGates,
  bwMin = 0.05,
  ind = 2L
) {
  stage <- "cytPos"
  pathProject <- file.path(
    tempdir(),
    paste0("stimgate-manual-cytpos-", label)
  )
  unlink(pathProject, recursive = TRUE, force = TRUE)
  dir.create(pathProject, recursive = TRUE, showWarnings = FALSE)

  oldIntermediate <- Sys.getenv(
    "STIMGATE_INTERMEDIATE",
    unset = NA_character_
  )
  on.exit(
    if (is.na(oldIntermediate)) {
      Sys.unsetenv("STIMGATE_INTERMEDIATE")
    } else {
      Sys.setenv(STIMGATE_INTERMEDIATE = oldIntermediate)
    },
    add = TRUE
  )
  Sys.setenv(STIMGATE_INTERMEDIATE = as.character(ind))

  basePos <- .getCytPosBasePos(
    ex = ex,
    gateTblInd = gateTbl
  )

  markerResults <- lapply(gateTbl$chnl, function(chnlCurr) {
    gateCyt <- .getCpPosGatesChnl(
      chnlCurr = chnlCurr,
      ex = ex,
      gateTblInd = gateTbl,
      basePos = basePos,
      bwMin = bwMin,
      ind = ind,
      stage = stage,
      pathProject = pathProject
    )

    saved <- function(suffix) {
      readRDS(
        .intSavePathSave(
          pathProject = pathProject,
          stage = stage,
          ind = ind,
          name = paste0(chnlCurr, "_", suffix)
        )
      )
    }

    incVec <- saved("incVec")
    cpOrig <- saved("cpOrig")
    shapeRef <- saved("shapeReference")
    cpTaut <- saved("cpTaut")
    tautInfo <- saved("tautStringInfo")
    cpFinal <- saved("cpCytPosFinal")

    conditionalValues <- ex[[chnlCurr]][incVec]
    tautDensity <- .getCpUnsLocAntimodeDensity(conditionalValues)

    thresholdTbl <- tibble::tibble(
      threshold = c(
        shapeRef$lowerX,
        cpOrig,
        cpTaut,
        cpFinal
      ),
      type = c(
        "marginal lower boundary",
        "original gate",
        "candidate aggressive gate",
        "final gateCyt"
      )
    ) |>
      dplyr::filter(is.finite(.data$threshold))

    distributionTbl <- dplyr::bind_rows(
      tibble::tibble(
        value = ex[[chnlCurr]],
        distribution = "full stimulated marginal"
      ),
      tibble::tibble(
        value = conditionalValues,
        distribution = "positive for another cytokine"
      )
    )

    distributionPlot <- ggplot2::ggplot(
      distributionTbl,
      ggplot2::aes(x = .data$value)
    ) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$distribution),
        ncol = 1,
        scales = "free_y"
      ) +
      ggplot2::geom_vline(
        data = thresholdTbl,
        ggplot2::aes(
          xintercept = .data$threshold,
          linetype = .data$type
        ),
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        title = paste(label, chnlCurr),
        x = "Expression",
        y = "Cell count",
        linetype = NULL
      )

    tautPlot <- if (is.null(tautDensity)) {
      NULL
    } else {
      ggplot2::ggplot(
        tibble::tibble(
          x = tautDensity$x,
          y = tautDensity$y
        ),
        ggplot2::aes(x = .data$x, y = .data$y)
      ) +
        ggplot2::geom_step() +
        ggplot2::geom_vline(
          data = thresholdTbl,
          ggplot2::aes(
            xintercept = .data$threshold,
            linetype = .data$type
          ),
          inherit.aes = FALSE
        ) +
        ggplot2::labs(
          title = paste(label, chnlCurr, "taut-string density"),
          x = "Expression",
          y = "Density",
          linetype = NULL
        )
    }

    marker <- gateTbl$marker[match(chnlCurr, gateTbl$chnl)]
    cells <- tibble::as_tibble(ex) |>
      dplyr::mutate(
        cell = dplyr::row_number(),
        currentMarker = marker,
        currentPositive = basePos$pos[[chnlCurr]],
        nOtherPositive = basePos$nPos -
          as.integer(basePos$pos[[chnlCurr]]),
        included = incVec
      ) |>
      dplyr::select(
        .data$cell,
        .data$currentMarker,
        .data$currentPositive,
        .data$nOtherPositive,
        .data$included,
        dplyr::everything()
      )

    list(
      summary = tibble::tibble(
        case = label,
        chnl = chnlCurr,
        marker = marker,
        nCells = nrow(ex),
        nOtherCytPos = tautInfo$nOtherCytPos,
        shapeReason = shapeRef$reason,
        peakX = shapeRef$peakX,
        windowWidth = shapeRef$windowWidth,
        lowerX = shapeRef$lowerX,
        gateOriginal = cpOrig,
        candidateAggressive = cpTaut,
        candidateReason = tautInfo$reason,
        gateCyt = gateCyt
      ),
      cells = cells,
      shapeReference = shapeRef,
      tautStringInfo = tautInfo,
      tautDensity = tautDensity,
      distributionPlot = distributionPlot,
      tautPlot = tautPlot
    )
  })

  names(markerResults) <- gateTbl$chnl

  list(
    pathProject = pathProject,
    baselineGates = gateTbl,
    summary = dplyr::bind_rows(
      lapply(markerResults, function(x) x$summary)
    ),
    markers = markerResults
  )
}

accepted <- inspectCytPosCase(
  makeCytPosExample("bimodal"),
  "accepted-bimodal"
)
fallback <- inspectCytPosCase(
  makeCytPosExample("unimodal"),
  "fallback-unimodal"
)

summaryTbl <- dplyr::bind_rows(
  accepted$summary,
  fallback$summary
)
print(summaryTbl)

# Cell identities used for the IFNG conditional distributions.
print(
  accepted$markers$IFNG$cells |>
    dplyr::filter(.data$included)
)
print(
  fallback$markers$IFNG$cells |>
    dplyr::filter(.data$included)
)

# Compare the full marginal and conditional distributions with the original,
# candidate and final thresholds.
print(accepted$markers$IFNG$distributionPlot)
print(accepted$markers$IFNG$tautPlot)
print(fallback$markers$IFNG$distributionPlot)
print(fallback$markers$IFNG$tautPlot)

# Useful variants:
# inspectCytPosCase(makeCytPosExample("bimodal", nOtherCytPos = 8L), "few-cells")
# Edit the component ranges in makeCytPosExample() to move the antimode closer
# to or farther below the original IFNG gate.
