# Return the full successful smoothing result

predVec and meanAbsError may be supplied when they were already
calculated during fit validation. This avoids repeating the full-data
prediction.

## Usage

``` r
.getCpUnsLocGetProbSmoothActualResponseSuccess(
  fit,
  dataMod,
  chnlSettings = list(),
  method = NA_character_,
  predVec = NULL,
  meanAbsError = NULL
)
```
