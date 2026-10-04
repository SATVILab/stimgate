# Changelog

## stimgate 0.99.15

### Performance

- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  computes combination gate statistics in a single pass over each tube,
  rather than once per marker combination, and streams expression by
  channel instead of loading a whole batch.

### Bug fixes

- In combination statistics, stimulated samples without any gates now
  report `NA` counts (with cell counts retained) rather than counts of
  zero.

## stimgate 0.99.14

### New features

- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  accepts sample names in `indBatchList`, resolved to indices as in
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md).

## stimgate 0.99.13

### Breaking changes

- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  validates `batchList`: each batch needs its unstimulated sample first
  and at least one stimulated sample; a shared unstimulated sample must
  be first in every batch, and a stimulated sample may belong to only
  one batch.
- [`getBatchList()`](https://satvilab.github.io/stimgate/reference/getBatchList.md)
  errors when a group has more than one unstimulated sample.

### Bug fixes

- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  matches unstimulated gates to batches by stimulated-sample membership,
  so indices of different digit widths (e.g. 9 and 10) and stimulated
  samples without gate rows no longer fail.
- [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md)
  reads saved local-FDR bandwidths from stimulated samples, rather than
  always falling back to `"nrd0"`.

## stimgate 0.99.12

### New features

- Cytometry entry points now accept FCS paths/directories, flowSets,
  cytosets, individual frames, lists of numeric matrices/data frames,
  and long data frames with a `sample` column. Non-GatingSet inputs
  expose only the root population.
- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  accepts sample names in `batchList` and saves resolved indices.

## stimgate 0.99.8

### New features

- `gateStim(bwScope = )` chooses which samples share the scalar
  local-FDR bandwidth: `"cytokine"` (new default; one trimmed-mean
  bandwidth per channel from about 100 tubes), `"cluster"` (one
  bandwidth per cluster of similarly shaped tubes) or `"sample"`
  (previous per-sample behaviour). The chosen values are reported and
  saved in the channel settings, so they can be inspected and fixed via
  `bw`.
- `bwCluster` is no longer estimated automatically. When `NULL`,
  cluster-based threshold sharing uses the shared local-FDR bandwidth; a
  supplied `bwCluster` now takes precedence rather than being a
  fallback.

## stimgate 0.99.5

### Breaking changes

- Removed the no-op
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  arguments `locLeftLowRel`, `locLeftLowAbs`, `locLeftCellFrac`,
  `locLeftLengthFrac` and `locTolRefPeak`.
- Dropped the `gtools` dependency.

### New features

- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  writes `manifest.csv` (one row per sample: `ind`, `batch`, `fileName`,
  `nCellPos`, `written`, `reason`) and returns it invisibly, so samples
  without an FCS file are explained.
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  attaches `attr(, "nCellPos")`, giving the positive-cell count for
  every requested sample, including those with none.
- `getStimGatesDetailed(pop = )` now filters by population.

### Bug fixes

- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  exports events from the chosen `pop`, passes inferred channels on, and
  no longer deletes existing output when gate preparation fails.
- Cluster bandwidths honour the `norm*` settings; known cytokine
  positivity is kept when another marker’s expression is NA.
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  with `ind = NULL` no longer reuses the first population’s samples for
  later populations.

## stimgate 0.99.4

### Breaking changes

- Removed the no-op
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  arguments `normExtraJitterFrac` and `normPeakFrac`.
- Removed the unused `gateUnsMethod` argument from
  [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  and
  [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md).
- `axisLimits()` is no longer exported.
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  returns zero rows, rather than one NA row, for samples with no
  positive cells.

### Bug fixes

- `plotStim(excMin = FALSE)` keeps univariate plots; bivariate gate
  lines are drawn on their own channel’s axis.
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  with `chnlGate`/`markerGate` no longer errors on marker-gated
  projects, and `bias = TRUE` reads the saved marker-keyed settings.
- Unstimulated cell counts are re-checked against `minCell` after
  cytokine-positive cells are removed.
- Per-channel settings are validated with exact names (partial matching
  wrongly rejected valid settings).

### Note

The snake_case API (`stimgate_gate()`, `stimgate_fcs_write()`, …) was
removed in 0.99.3; code written for it needs stimgate \< 0.99.3.

## stimgate 0.3.1.9013

### Major features

- Main gating functionality via `stimgate_gate()` to identify responding
  cells
- Statistics generation with `get_stats()`
- Visualization capabilities through `stimgate_plot()`
- FCS file writing with `stimgate_fcs_write()`
- Gate table extraction via `get_gate_tbl()`

### Minor improvements and bug fixes

- Improved documentation and examples
- Enhanced package structure for BioConductor submission
- Added proper dependency specifications
