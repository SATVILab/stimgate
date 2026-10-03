# Changelog

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
