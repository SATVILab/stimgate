# Use a high fitted-probability point when the derivative reference fails

The first point reaching a fixed fraction of the maximum fitted
probability supplies a conservative right-hand reference when no
derivative peak is identifiable. It is calculated from whichever
ordinary or shape-restricted probability fit the user selected.

## Usage

``` r
.getCpUnsLocHighProbabilityReference(
  dataMod,
  probCol,
  fraction = 0.85,
  derivativeInfo = NULL,
  shapeApplied = FALSE
)
```
