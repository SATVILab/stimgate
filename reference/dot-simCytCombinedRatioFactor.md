# Solve the correction factor with combined lower and upper variance

We adjust only the upper population. If c is the mean-distance
multiplier, the upper SD multiplier is 1 / c. c is chosen so that

## Usage

``` r
.simCytCombinedRatioFactor(
  delta,
  sdLower,
  sdUpper,
  targetRatio,
  eps = sqrt(.Machine$double.eps)
)
```

## Details

c \* delta / sqrt(sd_lower^2 + (sd_upper / c)^2) == target_ratio.
