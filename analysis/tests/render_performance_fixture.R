# Bounded visual fixture, without gateStim() or simulation. Run from checkout root:
# Rscript --no-init-file analysis/tests/render_performance_fixture.R /tmp/stimgate-performance-figures
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else file.path(tempdir(), "stimgate-performance-figures")
devtools::load_all(quiet = TRUE)
for (file in c("analysis-runtime.R", "analysis-plot-style.R", "analysis-mcse.R",
  "sim-bandwidth-analysis-plot.R", "sim-compare-freq_bs.R", "sim-compare-performance-plot.R")) {
  source(file.path("scripts", "r", file))
}
raw <- tidyr::expand_grid(n_cell = c(1000L, 10000L), iter = 1:20,
  sample = as.character(1:20), method = c("stimgate", "fbeta")) |>
  dplyr::mutate(
    sim_seed = 12345L, transformation = "gaussian", prob_response = 0.01,
    propRespTruth = 0.01, error = NA_character_,
    rel_error = dplyr::case_when(
      iter <= 10L & sample == "1" ~ ifelse(method == "stimgate", 0.4, 0.8),
      iter <= 5L & sample == "2" ~ -0.2,
      TRUE ~ 0),
    propRespEst = propRespTruth * (1 + rel_error)
  )
keys <- c("n_cell", "method", "transformation", "prob_response")
maxima <- .simCompareDatasetMaxSummary(raw, keys, expected_datasets = 20L, mcse = TRUE)
main <- .simCompareSignedErrorSummary(raw, keys, mcse = TRUE)
plots <- list(
  pooled_signed_percentiles = .simComparePlotSignedError(main, mcse = TRUE),
  dataset_max_occurrence = .simCompareDatasetMaxPlot(maxima, "occurrence", mcse = TRUE),
  dataset_max_signed_severity = .simCompareDatasetMaxPlot(maxima, "severity", mcse = TRUE)
)
plots$pooled_ratio_percentiles <- .simBandwidthRatioPlot(plots$pooled_signed_percentiles)
plots$dataset_max_ratio_severity <- .simBandwidthRatioPlot(plots$dataset_max_signed_severity)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (name in names(plots)) {
  .analysis_save_fig(plots[[name]], file.path(out, paste0(name, ".png")),
    width = 23, height = 18, mcse_mode = "both")
}
utils::write.csv(.simCompareDatasetMaxCoverage(maxima),
  file.path(out, "dataset_max_coverage.csv"), row.names = FALSE)
message("Fixture figures: ", normalizePath(out))
