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
    tautInfo <- saved("tautStringInfo")
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

    conditional <- ex[[chnlCurr]][incVec]

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
        candidateAggressive = tautInfo$threshold,
        candidateReason = tautInfo$reason,
        gateCyt = gateCyt
      ),
      cells = cells,
      marginal = ex[[chnlCurr]],
      conditional = conditional,
      shapeReference = shapeRef,
      tautStringInfo = tautInfo,
      tautDensity = .getCpUnsLocAntimodeDensity(conditional)
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

plotCytPosMarker <- function(result) {
  thresholds <- c(
    lower = result$shapeReference$lowerX,
    original = result$summary$gateOriginal,
    candidate = result$summary$candidateAggressive,
    final = result$summary$gateCyt
  )
  thresholds <- thresholds[is.finite(thresholds)]
  lineType <- seq_along(thresholds)

  addThresholds <- function() {
    graphics::abline(v = thresholds, lty = lineType)
    graphics::legend(
      "topright",
      legend = names(thresholds),
      lty = lineType,
      bty = "n"
    )
  }

  oldPar <- graphics::par(mfrow = c(3, 1))
  on.exit(graphics::par(oldPar), add = TRUE)

  graphics::hist(
    result$marginal,
    breaks = 30,
    main = "Full stimulated marginal",
    xlab = "Expression"
  )
  addThresholds()

  graphics::hist(
    result$conditional,
    breaks = 30,
    main = "Positive for at least one other cytokine",
    xlab = "Expression"
  )
  addThresholds()

  if (is.null(result$tautDensity)) {
    graphics::plot.new()
    graphics::title("Taut-string density unavailable")
  } else {
    graphics::plot(
      result$tautDensity$x,
      result$tautDensity$y,
      type = "s",
      main = "Conditional taut-string density",
      xlab = "Expression",
      ylab = "Density"
    )
    addThresholds()
  }
}

accepted <- inspectCytPosCase(
  makeCytPosExample("bimodal"),
  "accepted-bimodal"
)
fallback <- inspectCytPosCase(
  makeCytPosExample("unimodal"),
  "fallback-unimodal"
)

stopifnot(
  accepted$markers$IFNG$summary$gateCyt <
    accepted$markers$IFNG$summary$gateOriginal,
  fallback$markers$IFNG$summary$gateCyt ==
    fallback$markers$IFNG$summary$gateOriginal
)

summaryTbl <- dplyr::bind_rows(
  accepted$summary,
  fallback$summary
)
print(baselineGates)
print(summaryTbl)

# Exact cell identities entering the IFNG conditional distributions.
print(
  accepted$markers$IFNG$cells |>
    dplyr::filter(.data$included)
)
print(
  fallback$markers$IFNG$cells |>
    dplyr::filter(.data$included)
)

plotCytPosMarker(accepted$markers$IFNG)
plotCytPosMarker(fallback$markers$IFNG)

# Useful variants:
# inspectCytPosCase(makeCytPosExample("bimodal", nOtherCytPos = 8L), "few-cells")
# Edit the component ranges in makeCytPosExample() to move the antimode closer
# to or farther below the original IFNG gate.
