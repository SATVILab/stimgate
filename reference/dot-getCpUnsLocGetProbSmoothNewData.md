# Construct the minimal prediction data required by the smoother

The fitted SCAM has only the expression values, named `x`, as a
predictor (see
[`.fitScam()`](https://satvilab.github.io/stimgate/reference/dot-fitScam.md)),
so prediction does not require copying dataMod.

## Usage

``` r
.getCpUnsLocGetProbSmoothNewData(x)
```
