# Filter and select the sample's local-FDR threshold

Marginal boundary extensions use exTblStimOrig and unbiased
exTblUnsOrig, with complete original-tube denominators, before threshold
selection.

## Usage

``` r
.getCpUnsLocGetCp(
  dataMod,
  exTblStimOrig,
  exTblStimNoMin,
  exTblUnsOrig,
  exTblUnsBias,
  bias,
  cpMin,
  stage,
  pathProject,
  chnlSettings = list()
)
```
