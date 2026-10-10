# Build candidate thresholds using a discounted response-probability sum

All threshold methods share the bounded leading-run discount; the
estimate sums candidate rows' weighted fitted probabilities over the
original stimulated cell count.

## Usage

``` r
.getCpUnsLocGetCpDataThreshold(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsOrig,
  pathProject,
  stage,
  densityBw = NULL
)
```

## Arguments

- densityBw:

  numeric, list or NULL Local-FDR density bandwidth; adaptive objects
  use the shared bandwidth at the lowest estimate value. Default: NULL.
