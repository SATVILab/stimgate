# Fit one response-probability curve

The ordinary fit uses the existing preliminary probability filter. The
shape-restricted fit instead uses every density-grid point remaining
after the shape threshold, because that threshold has already defined
the admissible modelling region.

## Usage

``` r
.getCpUnsLocGetProbFit(
  exTblStimNoMin,
  exTblStimThreshold,
  exTblUnsThreshold,
  exTblUnsBias,
  bias,
  exTblUnsOrig,
  stage,
  pathProject,
  chnlSettings,
  applyPreliminaryFilter = TRUE,
  peakX = NULL,
  windowWidth = NULL
)
```
