# Find the marginal reference threshold and scan bins to its left

Find the marginal reference threshold and scan bins to its left

## Usage

``` r
.getCpUnsLocFilterMarginal(
  dataMod,
  chnlSettings,
  probCol,
  antimodeX = NULL,
  threshold = NULL,
  dominance = NULL,
  globalLowerBoundX = NA_real_,
  shapeLowerBoundX = NA_real_,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL
)
```

## Arguments

- exTblStimOrig, exTblUnsOrig:

  data.frame or NULL Original sample expression, without unstimulated
  bias, for marginal span frequencies. Default: NULL (skip trimming when
  original expression is unavailable).
