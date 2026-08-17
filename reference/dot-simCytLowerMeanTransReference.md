# Lower reference mean after transformation

For gamma this is the standard gamma transform of the raw lower
reference. For raw skew, calc_skew() is stochastic and depends on the
full vector, so the transformed lower reference should be estimated from
the realised lower population for that condition and marker.

## Usage

``` r
.simCytLowerMeanTransReference(
  transformationFunc,
  lowerMeanRawReference,
  yLower = NULL
)
```
