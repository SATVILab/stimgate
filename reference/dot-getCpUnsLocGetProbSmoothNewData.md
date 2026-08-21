# Construct the minimal prediction data required by the smoother

The fitted SCAM contains only the expression channel as a predictor, so
prediction does not require copying every column of dataMod.

## Usage

``` r
.getCpUnsLocGetProbSmoothNewData(chnl, x)
```
