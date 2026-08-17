# Calculate bandwidth using ordinary or background-normalised methods

Calculate bandwidth using ordinary or background-normalised methods

## Usage

``` r
.bwCalcOne(
  x,
  bwMtd,
  bwAdj = 1,
  bwNcellMin = NULL,
  bwNcellMax = NULL,
  normPeakFrac = 0.1,
  normPeakMinRel = 0.75,
  normExtraFrac = 0.2,
  normExtraMax = Inf,
  normExtraJitterFrac = 0.25,
  normLambda = seq(-2, 2, length.out = 81),
  normDensityN = 512L,
  normExcessBwMtd = "hpi3",
  normExcessNcell = 10000L,
  normAdaptiveNcell = 2500L,
  bwAdaptiveCore = NULL,
  bwAdaptiveExtra = NULL,
  bwAdaptiveCrossover = NULL,
  bwAdaptiveTransitionWidth = 0,
  normMtd = "moments",
  adaptive = FALSE
)
```
