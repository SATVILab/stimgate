# stimgate 0.99.28

## Breaking changes

- When `biasUns` is not supplied, `gateStim()` now sets it to `biasUnsFactor`
  times the bandwidth every sample of the marker uses: a fixed `bw`, or the
  shared bandwidth with `stimControl(bwScope = "cytokine")` (the default).
  Previously it used the fallback bandwidth, a separate estimate that could
  differ noticeably from the bandwidth actually used. Per-sample and
  per-cluster bandwidths still use the fallback.

## New features

- `stimControl(locThresholdMethod = "cap")` keeps the filtered-region gate
  unless the background-subtracted frequency above it is more than
  `locThresholdCap` (default 1.3) times the sum of fitted response
  probabilities. The gate then moves up to the lowest value where it is
  within that limit. This limits the large over-estimates the `"region"`
  method can give, while keeping its lower gates elsewhere.

# stimgate 0.99.27

## Bug fixes

- `gateStim()` now finds the unstimulated sample's main peak even when it lies
  below every stimulated cell, for example when the stimulated negatives are
  shifted upwards by more than `biasUns`. Previously a small bump among the
  responding cells could be taken as the peak, which excluded most genuine
  responders and set the gate far too high. Gates change only for samples
  whose unstimulated cells extend beyond the stimulated cells' range.

# stimgate 0.99.26

## Breaking changes

- `gateStim()` now sets each local-FDR gate at the lower edge of the region
  kept by its filtering steps (`stimControl(locThresholdMethod = "region")`,
  the new default). Previously the gate was raised further, until the
  background-subtracted frequency above it matched the sum of the fitted
  response probabilities. Because those probabilities are deliberately
  conservative, that extra step could exclude genuine responders, so gates
  are now generally lower and frequencies higher. Use
  `stimControl(locThresholdMethod = "match")` to reproduce the previous gates.
  Threshold sharing (`clusterGates`), combining gates within a batch and
  cytokine-positive refinement still apply afterwards, and
  `getStimStats()` counts cells above the gate actually applied.

## New features

- `getStimGatesDetailed()` reports the threshold method
  (`locThresholdMethod`) and the filtered-region boundary (`locRegionX`) for
  each condition-level gate, and `stimgateMetaReadSettingsChnls()` records
  the method used for each marker. The probability-sum estimate stays in
  `propBsEst` as a diagnostic; under `"region"` it is not the frequency
  above the gate.

# stimgate 0.99.25

## Bug fixes

- `getStimGates()` and `getStimStats()` name clustered gates `loc_minClust`,
  matching the `loc_min` form of unclustered gates (was `locminClust`).

# stimgate 0.99.24

## Breaking changes

- `stimControl()` now defaults to `bwMtd = "nrd0"` (was `"hpi1"`),
  `bwNcellMax = 10000` (was 100000) and `bwNcellMin = bwNcellMax` (was 100),
  so every tube is downsampled or upsampled to 10,000 cells before its
  bandwidth is chosen.
- Shared bandwidths (`bwScope = "cytokine"` or `"cluster"`) now prefer tubes
  with at least `bwNcellMax` cells. If there are too few, they add tubes
  chosen at random from the next sizes down (9,000-10,000 cells first, then
  8,000-9,000, and so on to 5,000 by default), and estimate every bandwidth on
  `bwNcellMax` cells, upsampling smaller tubes.
- With `biasUns = NULL`, `gateStim()` now sets `biasUns` to `biasUnsFactor`
  times the fallback bandwidth, rather than a quarter of that. With the
  default `biasUnsFactor = 1`, the automatic `biasUns` is therefore four times
  larger than before.

# stimgate 0.99.23

## Bug fixes

- `gateStim()` now places each local-FDR gate just below the cell value it
  selects, rather than on it. Gates count cells strictly above them, so the
  selected cell and any ties were previously left out in both tubes; with
  very rare responses this could reduce the estimate to zero.
  The gate now sits in the empty gap below that cell, at most twice the
  density bandwidth at that cell below it (including adaptive bandwidths),
  so it counts exactly the cells the threshold was selected for.

# stimgate 0.99.21

## Bug fixes

- `getStimGates()` preserves threshold-generation provenance for clustered
  gates, including cluster-derived thresholds and retained finite fallback gates.

# stimgate 0.99.20

## Bug fixes

- `gateStim()` now excludes cells exactly at the gate from sample-level
  frequencies, matching the strict gate applied to cells.

## Analysis

- Reject stale caches and unsafe partial HPC outputs. Preserve comparator
  errors, threshold provenance and finite-estimate coverage in summaries.
- Report method differences and uncertainty across paired simulated datasets,
  and retain negative background-subtracted estimates in signed-error plots.

# stimgate 0.99.19

## Breaking changes

- `stimControl()` defaults to `bwScope = "cytokine"` again, so `gateStim()`
  shares one local-FDR bandwidth per marker. Per-sample estimation remains
  available with `bwScope = "sample"`.

# stimgate 0.99.16

## Breaking changes

- `stimControl()` now defaults to `bwScope = "sample"`, so `gateStim()`
  estimates the local-FDR bandwidth separately for each sample again. Shared
  per-marker bandwidths remain available with `bwScope = "cytokine"`.

# stimgate 0.99.15

## Performance

- `gateStim()` computes combination gate statistics in a single pass over each
  tube, rather than once per marker combination, and streams expression by
  channel instead of loading a whole batch.

## Bug fixes

- In combination statistics, stimulated samples without any gates now report
  `NA` counts (with cell counts retained) rather than counts of zero.

# stimgate 0.99.14

## New features

- `writeStimFCS()` accepts sample names in `indBatchList`, resolved to indices
  as in `gateStim()`.

# stimgate 0.99.13

## Breaking changes

- `gateStim()` validates `batchList`: each batch needs its unstimulated sample
  first and at least one stimulated sample; a shared unstimulated sample must be
  first in every batch, and a stimulated sample may belong to only one batch.
- `getBatchList()` errors when a group has more than one unstimulated sample.

## Bug fixes

- `writeStimFCS()` matches unstimulated gates to batches by stimulated-sample
  membership, so indices of different digit widths (e.g. 9 and 10) and stimulated
  samples without gate rows no longer fail.
- `plotStim()` reads saved local-FDR bandwidths from stimulated samples, rather
  than always falling back to `"nrd0"`.

# stimgate 0.99.12

## New features

- Cytometry entry points now accept FCS paths/directories, flowSets, cytosets,
  individual frames, lists of numeric matrices/data frames, and long data frames
  with a `sample` column. Non-GatingSet inputs expose only the root population.
- `gateStim()` accepts sample names in `batchList` and saves resolved indices.

# stimgate 0.99.8

## New features

* `gateStim(bwScope = )` chooses which samples share the scalar local-FDR bandwidth: `"cytokine"` (new default; one trimmed-mean bandwidth per channel from about 100 tubes), `"cluster"` (one bandwidth per cluster of similarly shaped tubes) or `"sample"` (previous per-sample behaviour). The chosen values are reported and saved in the channel settings, so they can be inspected and fixed via `bw`.
* `bwCluster` is no longer estimated automatically. When `NULL`, cluster-based threshold sharing uses the shared local-FDR bandwidth; a supplied `bwCluster` now takes precedence rather than being a fallback.

# stimgate 0.99.5

## Breaking changes

* Removed the no-op `gateStim()` arguments `locLeftLowRel`, `locLeftLowAbs`, `locLeftCellFrac`, `locLeftLengthFrac` and `locTolRefPeak`.
* Dropped the `gtools` dependency.

## New features

* `writeStimFCS()` writes `manifest.csv` (one row per sample: `ind`, `batch`, `fileName`, `nCellPos`, `written`, `reason`) and returns it invisibly, so samples without an FCS file are explained.
* `getStimExpr()` attaches `attr(, "nCellPos")`, giving the positive-cell count for every requested sample, including those with none.
* `getStimGatesDetailed(pop = )` now filters by population.

## Bug fixes

* `writeStimFCS()` exports events from the chosen `pop`, passes inferred channels on, and no longer deletes existing output when gate preparation fails.
* Cluster bandwidths honour the `norm*` settings; known cytokine positivity is kept when another marker's expression is NA.
* `getStimExpr()` with `ind = NULL` no longer reuses the first population's samples for later populations.

# stimgate 0.99.4

## Breaking changes

* Removed the no-op `gateStim()` arguments `normExtraJitterFrac` and `normPeakFrac`.
* Removed the unused `gateUnsMethod` argument from `getStimExpr()` and `plotStim()`.
* `axisLimits()` is no longer exported.
* `getStimExpr()` returns zero rows, rather than one NA row, for samples with no positive cells.

## Bug fixes

* `plotStim(excMin = FALSE)` keeps univariate plots; bivariate gate lines are drawn on their own channel's axis.
* `getStimExpr()` with `chnlGate`/`markerGate` no longer errors on marker-gated projects, and `bias = TRUE` reads the saved marker-keyed settings.
* Unstimulated cell counts are re-checked against `minCell` after cytokine-positive cells are removed.
* Per-channel settings are validated with exact names (partial matching wrongly rejected valid settings).

## Note

The snake_case API (`stimgate_gate()`, `stimgate_fcs_write()`, ...) was removed in 0.99.3; code written for it needs stimgate < 0.99.3.

# stimgate 0.3.1.9013

## Major features

* Main gating functionality via `stimgate_gate()` to identify responding cells
* Statistics generation with `get_stats()` 
* Visualization capabilities through `stimgate_plot()`
* FCS file writing with `stimgate_fcs_write()`
* Gate table extraction via `get_gate_tbl()`

## Minor improvements and bug fixes

* Improved documentation and examples
* Enhanced package structure for BioConductor submission
* Added proper dependency specifications