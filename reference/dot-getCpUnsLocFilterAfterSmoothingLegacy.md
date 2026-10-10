# Apply all filtering steps after smoothing

Apply all filtering steps after smoothing

## Usage

``` r
.getCpUnsLocFilterAfterSmoothingLegacy(
  dataMod,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  chnlSettings,
  exTblStimOrig = NULL,
  exTblUnsOrig = NULL
)
```

## Arguments

- exTblStimOrig, exTblUnsOrig:

  data.frame or NULL Original sample expression, without unstimulated
  bias, for marginal span frequencies. Default: NULL (skip trimming when
  original expression is unavailable).
