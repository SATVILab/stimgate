# Use the lower boundary of the filtered region as the gate

The gate is the filtering boundary itself (`xSum`); cells strictly above
it are positive. The probability-sum frequency estimate is not matched,
no empirical cell value is selected and the gate is not moved below a
cell. The no-response and non-finite cases fall back exactly as in
matching.

## Usage

``` r
.getCpUnsLocGetCpRegion(
  dataThreshold,
  regionX,
  exTblStimNoMin,
  exTblUnsBias,
  cpMin,
  stage
)
```

## Arguments

- regionX:

  numeric Lower boundary of the filtered region.

## Value

list from `.getCpUnsLocConditionOut()`.
