# Fit one response-probability curve

The ordinary fit uses the existing preliminary probability filter. The
shape-restricted fit instead uses every density-grid point remaining
after the shape threshold, because that threshold has already defined
the admissible modelling region. It retains the ordinary
negative-population width and its diagnostics rather than measuring the
truncated distribution.

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
  windowWidth = NULL,
  shiftedPeakRef = NULL,
  windowWidthInfo = NULL
)
```

## Arguments

- windowWidthInfo:

  list or NULL Ordinary per-tube width diagnostics.
