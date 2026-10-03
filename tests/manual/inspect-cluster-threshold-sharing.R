# Manual inspection of paired-density clustering and threshold sharing.
#
# Run from the repository root:
#   source("tests/manual/inspect-cluster-threshold-sharing.R")
#
# The objects created below are intentionally left in the calling environment so
# that densities, cluster assignments and threshold changes can be inspected.

devtools::load_all()

control <- stimgate:::.getCpClusterControlUpdate(list(
  nGrid = 64L,
  gapBootstraps = 20L,
  kmeansNstart = 5L,
  seed = 1L
))

chnl_settings <- list(
  chnlCut = "expr",
  bwMtd = "nrd0",
  bwAdj = 1,
  bwMin = 0.01,
  bwMax = 1,
  bwFallback = 0.08
)

make_ex <- function(values) {
  ex <- data.frame(expr = as.numeric(values))
  attr(ex, "chnlCut") <- "expr"
  ex
}

base_expr <- c(
  seq(0.03, 0.20, length.out = 200),
  seq(0.45, 0.75, length.out = 200)
)

make_pair <- function(ind, stim_shift, uns_shift) {
  list(
    ind = ind,
    batch = paste0("batch_", ind),
    stim = make_ex(base_expr + stim_shift),
    uns = make_ex(base_expr + uns_shift)
  )
}

similar_ids <- sprintf("similar_%02d", seq_len(8L))
similar_pairs <- stats::setNames(
  lapply(
    similar_ids,
    make_pair,
    stim_shift = 0.12,
    uns_shift = 0.03
  ),
  similar_ids
)

ex_lookup <- c(
  similar_pairs,
  list(
    distinct_01 = make_pair(
      "distinct_01",
      stim_shift = 0.58,
      uns_shift = 0.04
    ),
    distinct_02 = make_pair(
      "distinct_02",
      stim_shift = 0.10,
      uns_shift = 0.62
    )
  )
)

# The first seven similar samples are direct local-FDR donors. The first and
# seventh are deliberately low/high so winsorising can be inspected. The eighth
# similar sample is non-direct and should borrow from its cluster. The two
# clearly different samples are also non-direct and should retain their original
# thresholds when their clusters contain no direct donor.
gate_tbl <- tibble::tibble(
  ind = c(similar_ids, "distinct_01", "distinct_02"),
  gate = c(1.00, 1.30, 1.35, 1.40, 1.45, 1.50, 3.00, 2.00, 2.20, 2.40),
  locGenerated = c(rep(TRUE, 7L), rep(FALSE, 3L)),
  locGeneratedDirect = c(rep(TRUE, 7L), rep(FALSE, 3L)),
  locSource = c(rep("local_fdr", 7L), rep("high_fallback", 3L)),
  locReason = c(
    rep("direct_local_fdr_threshold", 7L),
    rep("no_direct_local_fdr_threshold", 3L)
  )
) |>
  stimgate:::.getCpClusterLocGateTblPrepare()

direct <- gate_tbl$locGeneratedDirect %in% TRUE &
  is.finite(gate_tbl$gate)

common_bw <- stimgate:::.getCpClusterLocCommonBw(
  indDirect = gate_tbl$ind[direct],
  exLookup = ex_lookup,
  chnlSettings = chnl_settings
)

expr_range <- stimgate:::.getCpClusterLocExprRange(ex_lookup)
left_upper_x <- control$leftThresholdFrac * stats::quantile(
  gate_tbl$gate[direct],
  probs = control$leftThresholdQuantile,
  na.rm = TRUE
)[[1]]

density_grid <- stimgate:::.getCpClusterLocDensityGrid(
  exprMin = expr_range[["min"]],
  leftUpperX = left_upper_x,
  nGrid = control$nGrid
)

feature_tbl <- stimgate:::.getCpClusterLocJointFeatureTbl(
  exLookup = ex_lookup,
  densityGrid = density_grid,
  bw = common_bw
)

cluster_obj <- stimgate:::.getCpClusterLocClusters(
  featureTbl = feature_tbl,
  control = control
)
cluster_tbl <- cluster_obj$clusterTbl

loc_tbl <- gate_tbl |>
  dplyr::left_join(cluster_tbl, by = "ind")

cluster_threshold_tbl <- stimgate:::.getCpClusterLocApplyQuantiles(
  locTbl = loc_tbl,
  commonBw = common_bw,
  control = control,
  nInitialClusters = cluster_obj$nInitialClusters
) |>
  dplyr::arrange(.data$grp, .data$ind)

pair_tbl <- purrr::map_df(ex_lookup, function(ex_pair) {
  stim <- stimgate:::.getCut(ex_pair$stim)
  uns <- stimgate:::.getCut(ex_pair$uns)
  tibble::tibble(
    ind = ex_pair$ind,
    batch = ex_pair$batch,
    stim_min = min(stim),
    stim_max = max(stim),
    uns_min = min(uns),
    uns_max = max(uns)
  )
}) |>
  dplyr::left_join(cluster_tbl, by = "ind")

inspection_tbl <- cluster_threshold_tbl |>
  dplyr::transmute(
    grp = .data$grp,
    ind = .data$ind,
    original = .data$cpOrigQuantMin,
    direct_donor = .data$locGeneratedDirect,
    cluster_n_direct = .data$locClusterNDirect,
    q15 = .data$locClusterQ15,
    q60 = .data$locClusterQ60,
    q85 = .data$locClusterQ85,
    final = .data$cpJoinTgOrig,
    action = .data$locClusterAction,
    adjusted = .data$locClusterAdjusted,
    source = .data$locSource,
    reason = .data$locReason
  )

direct_thresholds <- loc_tbl |>
  dplyr::filter(
    .data$locGeneratedDirect %in% TRUE,
    is.finite(.data$gate)
  ) |>
  dplyr::group_by(.data$grp) |>
  dplyr::summarise(
    direct_thresholds = paste(sort(.data$gate), collapse = ", "),
    .groups = "drop"
  )

cluster_summary_tbl <- inspection_tbl |>
  dplyr::distinct(
    .data$grp,
    .data$cluster_n_direct,
    .data$q15,
    .data$q60,
    .data$q85
  ) |>
  dplyr::left_join(direct_thresholds, by = "grp") |>
  dplyr::arrange(.data$grp)

feature_cols <- stimgate:::.getCpClusterLocFeatureCols(feature_tbl)
feature_long <- feature_tbl |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(feature_cols),
    names_to = "feature",
    values_to = "density"
  ) |>
  dplyr::mutate(
    condition = dplyr::if_else(
      startsWith(.data$feature, "uns_"),
      "unstimulated",
      "stimulated"
    ),
    grid_index = as.integer(sub("^.*x", "", .data$feature)),
    expression = density_grid[.data$grid_index]
  ) |>
  dplyr::left_join(cluster_tbl, by = "ind")

feature_plot <- ggplot2::ggplot(
  feature_long,
  ggplot2::aes(
    x = .data$expression,
    y = .data$density,
    colour = .data$condition
  )
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::facet_wrap(
    ggplot2::vars(.data$grp, .data$ind),
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "Paired density features used for clustering",
    subtitle = sprintf(
      "Common bandwidth %.4f on a %d-point grid",
      common_bw,
      length(density_grid)
    ),
    x = "Expression on common grid",
    y = "Normalised density feature",
    colour = NULL
  ) +
  ggplot2::theme_minimal()

threshold_long <- inspection_tbl |>
  dplyr::select(
    .data$grp,
    .data$ind,
    .data$direct_donor,
    original = .data$original,
    final = .data$final
  ) |>
  tidyr::pivot_longer(
    cols = c("original", "final"),
    names_to = "threshold_stage",
    values_to = "threshold"
  )

threshold_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = inspection_tbl,
    ggplot2::aes(
      x = .data$ind,
      xend = .data$ind,
      y = .data$original,
      yend = .data$final
    ),
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    data = threshold_long,
    ggplot2::aes(
      x = .data$ind,
      y = .data$threshold,
      colour = .data$threshold_stage,
      shape = .data$direct_donor
    ),
    size = 2.5
  ) +
  ggplot2::facet_wrap(ggplot2::vars(.data$grp), scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Original and cluster-adjusted thresholds",
    x = NULL,
    y = "Threshold",
    colour = NULL,
    shape = "Direct donor"
  ) +
  ggplot2::theme_minimal()

cat("\nInput stim/unstim pairs and cluster assignments:\n")
print(pair_tbl, n = Inf)

cat("\nCommon bandwidth:\n")
print(common_bw)

cat("\nCommon density grid summary:\n")
print(tibble::tibble(
  n = length(density_grid),
  min = min(density_grid),
  max = max(density_grid),
  step = density_grid[[2]] - density_grid[[1]]
))

cat("\nWithin-cluster direct-threshold distribution and transfer quantiles:\n")
print(cluster_summary_tbl, n = Inf)

cat("\nOriginal/final thresholds with provenance:\n")
print(inspection_tbl, n = Inf)

if (interactive()) {
  print(feature_plot)
  print(threshold_plot)
}
