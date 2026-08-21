# Predict one fitted smoother over the full model data

This evaluates the expensive full-data prediction once and calculates
the fit diagnostic used to decide whether the smoother is acceptable.
The derivative is deliberately not calculated here because a rejected
fit does not need one.

## Usage

``` r
.getCpUnsLocGetProbSmoothFitEval(fit, dataMod)
```
