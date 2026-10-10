# Assemble response probabilities with peak-density negative widths

Widths and their per-tube diagnostics use the same density grids as the
main peaks; the probability comparison retains its original grid.

## Usage

``` r
.getCpUnsLocGetProbTbl(
  densTblRaw,
  stage,
  cpMin,
  exVecStimThreshold,
  exVecUnsThreshold,
  shiftedPeak = NULL
)
```
