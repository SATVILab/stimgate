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