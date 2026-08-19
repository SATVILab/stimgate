# stimgate

`stimgate` is an R package for identifying cells that may have responded
to immune stimulation in flow cytometry data. It compares stimulated
samples with matched unstimulated controls and estimates marker-specific
expression gates relative to the unstimulated background, rather than
applying one fixed positivity threshold across samples.

The package works with
[`flowWorkspace::GatingSet`](https://rdrr.io/pkg/flowWorkspace/man/GatingSet-class.html)
objects. Its main entry point,
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
performs the gating workflow and writes the resulting gates, cached
expression data and response statistics to a project directory.
Plotting, expression extraction and FCS export functions can then use
those saved results without rerunning the full gating procedure.

## Method overview

For each selected marker and gated population, StimGate compares
stimulated and matched unstimulated expression distributions. The
current workflow uses density-based local response-probability/local-FDR
thresholding, followed by cross-sample refinement of the resulting
gates. It can also refine gates using the target-marker distribution
among cells positive for other cytokines. The final project output
includes marker-level gates and statistics for marker-positive and
marker-combination populations.

The statistical method is still under active development, so the README
should give only this high-level overview. Detailed algorithmic settings
belong in the function documentation and research analyses rather than
here.

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

The canonical packaged example data are provided by
[`getExampleData()`](https://satvilab.github.io/stimgate/reference/getExampleData.md),
so the basic workflow can be tried without external FCS files.

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

# Gating statistics written by gateStim()
gate_stats <- readRDS(file.path(path_project, "gateStats.rds"))

# Plot the samples in the first matched batch
plotStim(
  ind = example_data$batchList[[1]],
  .data = gs,
  pathProject = path_project,
  marker = example_data$marker,
  grid = TRUE
)
```

For real data,
[`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
requires a compatible flow-cytometry object, the population to gate, the
marker or channel names to analyse, and a `batchList` describing the
matched samples. See the function documentation for the full set of
method and bandwidth controls.

## Main functions

- [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md)
  runs the StimGate gating workflow and saves the project results.
- [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md)
  plots univariate or bivariate expression together with the fitted
  gates.
- [`getStimExpr()`](https://satvilab.github.io/stimgate/reference/getStimExpr.md)
  reads expression data saved by a StimGate project and can optionally
  apply the saved gates.
- [`getStimGatesDetailed()`](https://satvilab.github.io/stimgate/reference/getStimGatesDetailed.md)
  reads detailed saved threshold diagnostics when these are available.
- [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md)
  exports marker-positive cells to FCS files using fitted gates.
- [`getBatchList()`](https://satvilab.github.io/stimgate/reference/getBatchList.md)
  constructs matched sample batches from sample metadata.
- [`getExampleData()`](https://satvilab.github.io/stimgate/reference/getExampleData.md)
  loads the packaged canonical dataset for examples and package tests.

## Repository structure

The package implementation is in `R/`, with tests in `tests/testthat/`.
Method-development, simulation and comparison analyses are kept
separately under `analysis/`; these are research/development analyses
and are not required to use the package.

Developer and coding-agent setup instructions are in `AGENTS.md`. Build
history for the research project is recorded in `BUILDLOG.md`.

## Licence

`stimgate` is released under GPL (\>= 3).
