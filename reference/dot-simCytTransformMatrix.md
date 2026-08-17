# Apply a transformation to a simulated cluster before ratio correction

Ratio correction is done at the condition level, because it needs the
lower component spread as well as the upper component spread.

## Usage

``` r
.simCytTransformMatrix(simData, transformationFunc)
```
