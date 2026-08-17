# Apply the upper-population ratio correction

The target ratio is measured on the pre-transformation scale as the
distance from the upper mean to the lower reference mean, divided by the
combined spread of the lower and upper populations. The achieved ratio
is measured after the transformation in the same way. Only the upper
population is adjusted: its mean distance is multiplied by c and its SD
is multiplied by 1 / c.

## Usage

``` r
.simCytRatioAdjustUpper(
  yUpper,
  xUpperRaw,
  xLowerRaw,
  yLower,
  lowerMeanRawReference,
  lowerMeanTransReference,
  eps = sqrt(.Machine$double.eps)
)
```
