# Apply the marginal-bin scan and trim unsupported boundary extensions

Empty bins are pending: they neither move the cut nor interrupt
consecutive non-empty rejections. An accepted non-empty bin pulls in
intervening empty bins and up to two rejected non-empty bins. Each
acceptance records a span \[new cut, previous cut). From the last span
backwards, undo steps whose raw stimulated fraction minus raw
unstimulated fraction is not positive, stopping at the first positive
contribution. Denominators include all original cells. scanTbl records
pending, accepted and finally retained bins; trimTbl records each step's
counts, contribution and final retention, with nStepsTrimmed and
trimReason explaining the trim outcome. finalStartX is the final cut.

## Usage

``` r
.getCpUnsLocFilterMarginalBins(
  dataMod,
  chnlSettings,
  probCol,
  startX,
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
