# Use the region boundary unless its frequency exceeds the capped estimate

Keeps the region boundary (`regionX`, as under "region") when the
background-subtracted frequency strictly above it is at most `cap` times
the probability-sum estimate. Otherwise the gate is placed below the
lowest candidate cell value at or above the boundary whose frequency
(cells at or above it) is within the cap, as matching places its gate
below the selected cell. When no candidate is within the cap, the
matching choice is used.

## Usage

``` r
.getCpUnsLocGetCpCap(
  dataThreshold,
  regionX,
  cap,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage,
  exTblStimOrig,
  exTblUnsOrig,
  densityBw = NULL
)
```

## Arguments

- cap:

  numeric Largest allowed ratio of frequency to estimate.

## Value

list from `.getCpUnsLocConditionOut()`. A selected cell is attribute
`cpSelected`; `locCapExceededAbove` records whether any higher
candidate's frequency exceeds the cap again.
