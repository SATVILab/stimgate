stimgate
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/SATVILab/stimgate/graph/badge.svg)](https://app.codecov.io/gh/SATVILab/stimgate)
<!-- badges: end -->

`stimgate` is an R package that finds cells that may have responded to
immune stimulation in flow cytometry data. For each donor, it compares
the stimulated samples with an unstimulated sample and sets a gate for
each marker based on that donor's background, instead of using one fixed
cut-off for everyone.

You can give it a `flowWorkspace::GatingSet`, other
`flowCore`/`flowWorkspace` objects, FCS files or expression matrices.
`gateStim()` does the gating and saves the gates, expression values and
response statistics in a folder. The other functions read from that
folder, so you can plot results or export cells without gating again.

## Method overview

For each marker, StimGate compares how strongly cells express it in the
stimulated and unstimulated samples. It places the gate where cells
start to be more common in the stimulated sample than the unstimulated
one would explain. It then makes gates more consistent across similar
samples, and can adjust a marker's gate using cells that are already
positive for another cytokine. The results give, for each sample, the
percentage of cells positive for each marker and each combination of
markers, minus the unstimulated background.

The method is still being developed. See the function help pages for its
settings.

## Installation

The development version can be installed from GitHub. `stimgate` depends
on Bioconductor flow-cytometry packages, including `flowCore` and
`flowWorkspace`.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("flowCore", "flowWorkspace"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("SATVILab/stimgate")
```

The package requires R \>= 4.4.0.

## Minimal example

`getExampleData()` loads a small example dataset that comes with the
package, so you can try StimGate without your own FCS files.

``` r
library(stimgate)

example_data <- getExampleData()
gs <- flowWorkspace::load_gs(example_data$pathGs)

path_project <- file.path(tempdir(), "stimgate-example")

path_project <- gateStim(
  pathProject = path_project,
  .data = gs,
  popGate = "root",
  batchList = example_data$batchList,
  marker = example_data$marker
)

# Gates and response statistics saved by gateStim()
gates <- getStimGates(path_project)
gate_stats <- getStimStats(path_project)

# Plot the samples in the first matched batch
plotStim(
  ind = example_data$batchList[[1]],
  .data = gs,
  pathProject = path_project,
  marker = example_data$marker,
  grid = TRUE
)
```

For real data, `gateStim()` needs your cytometry data, the markers or
channels to gate, and a `batchList` that groups each donor's samples
with the unstimulated sample first. The defaults suit most data. To tune
the method, change `biasUns` or `bw`, pass other settings with `control
= stimControl(...)`, or set them for single markers with
`markerControl`.

## Main functions

- `gateStim()` finds the gates and saves the results.
- `getStimGates()` and `getStimStats()` read the saved gates and
  response statistics.
- `plotStim()` plots expression for one or two markers, with the gates.
- `getStimExpr()` reads the saved expression values, optionally only for
  positive cells.
- `getStimGatesDetailed()` reads extra detail on how each gate was
  chosen, if it was saved.
- `writeStimFCS()` saves the positive cells as FCS files.
- `getBatchList()` builds `batchList` from a table describing your
  samples.
- `getExampleData()` loads the packaged example dataset.

## Repository structure

The package code is in `R/` and its tests are in `tests/testthat/`. The
`analysis/` folder holds research analyses (simulations and method
comparisons); you don't need them to use the package.

Guidance for working on the code is in `AGENTS.md`. `BUILDLOG.md`
records past runs of the research analyses.

## Licence

`stimgate` is released under GPL (\>= 3).
