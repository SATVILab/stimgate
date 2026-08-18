# Bound the marginal scan by the stimulated peak's descending shoulder

This does not itself define a filtering threshold. When no antimode was
identified, it prevents the marginal scan from considering bins below
the point where the negative stimulated-density derivative has flattened
to the requested fraction of its most negative value on the peak's right
shoulder.

## Usage

``` r
.getCpUnsLocMarginalDensityLowerBound(density, peakX, fraction = 1/200)
```
