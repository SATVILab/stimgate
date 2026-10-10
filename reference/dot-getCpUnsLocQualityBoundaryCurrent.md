# Obtain the quality-based lower boundary starting at x_clear

Obtain the quality-based lower boundary starting at x_clear

## Usage

``` r
.getCpUnsLocQualityBoundaryCurrent(
  dataMod,
  chnlSettings,
  probCol,
  xClear,
  lowerBoundX = NA_real_,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL
)
```

## Arguments

- exTblStimOrig, exTblUnsOrig:

  data.frame or NULL Original sample expression, without unstimulated
  bias, for marginal span frequencies. Default: NULL (skip trimming when
  original expression is unavailable).
