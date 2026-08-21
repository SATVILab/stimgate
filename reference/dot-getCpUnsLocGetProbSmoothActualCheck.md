# Check whether a fitted smoother is acceptable

Retained as a separate helper for existing internal callers and tests.
Unlike the previous implementation, this does not calculate a derivative
table just to decide whether the fit should be retained.

## Usage

``` r
.getCpUnsLocGetProbSmoothActualCheck(fit, dataMod)
```
