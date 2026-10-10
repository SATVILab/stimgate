# Apply current post-smoothing filtering for the ordinary local-FDR route

Apply current post-smoothing filtering for the ordinary local-FDR route

## Usage

``` r
.getCpUnsLocFilterAfterSmoothing(
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
