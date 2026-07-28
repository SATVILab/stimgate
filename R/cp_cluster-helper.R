# Defaults for joint-density clustering and direct-threshold transfer.
# The values may be overridden through the internal `control` list supplied to
# `.getCpCluster()`.
#' @keywords internal
.getCpClusterControlUpdate <- function(control) {
  defaults <- list(
    locClusterLeftFraction = 0.8,
    locClusterLeftQuantile = 0.1,
    locClusterDensityN = 128L,
    locClusterWinsorLower = 0.15,
    locClusterImputeQuantile = 0.60,
    locClusterWinsorUpper = 0.85,
    locClusterMinDirect = 5L,
    locClusterMaxClusters = 6L,
    locClusterNstart = 25L,
    locClusterGapB = 50L,
    locClusterSeed = 1L
  )

  for (name in names(defaults)) {
    if (!name %in% names(control) || is.null(control[[name]])) {
      control[[name]] <- defaults[[name]]
    }
  }

  control$locClusterDensityN <- max(
    16L,
    as.integer(control$locClusterDensityN)[1]
  )
  control$locClusterMinDirect <- max(
    1L,
    as.integer(control$locClusterMinDirect)[1]
  )
  control$locClusterMaxClusters <- max(
    1L,
    as.integer(control$locClusterMaxClusters)[1]
  )
  control$locClusterNstart <- max(
    1L,
    as.integer(control$locClusterNstart)[1]
  )
  control$locClusterGapB <- max(
    1L,
    as.integer(control$locClusterGapB)[1]
  )
  control$locClusterSeed <- as.integer(control$locClusterSeed)[1]

  probs <- c(
    control$locClusterWinsorLower,
    control$locClusterImputeQuantile,
    control$locClusterWinsorUpper
  )
  if (
    any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) ||
      !(
        probs[1] <= probs[2] &&
          probs[2] <= probs[3]
      )
  ) {
    stop(
      "Cluster quantiles must satisfy 0 <= lower <= impute <= upper <= 1."
    )
  }

  control
}
