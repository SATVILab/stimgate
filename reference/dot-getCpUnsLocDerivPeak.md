# Select the left-most derivative peak meeting alpha

Flat-topped peaks are represented by the left-most point of the plateau.
With `minRiseProb`, a peak counts only if its own rise lifts the
probability to at least `minRiseProb`: the probability where the
derivative first falls back to `riseFrac` times the peak height (the top
of that rise).

## Usage

``` r
.getCpUnsLocDerivPeak(
  x,
  prob,
  deriv,
  alpha = 0.75,
  leftRiseFrac = 0.15,
  minRiseProb = NA_real_,
  riseFrac = 0.2
)
```
