# Get detailed gate diagnostics

Read detailed local-FDR threshold diagnostics saved during gating,
including condition-level, sample-level and final cluster-level
thresholds with the corresponding background-subtracted frequencies.

## Usage

``` r
getStimGatesDetailed(
  pathProject,
  pop = NULL,
  marker = NULL,
  chnl = NULL,
  save = FALSE,
  pathSave = NULL
)
```

## Arguments

- pathProject:

  character. Path to the project directory.

- pop:

  character. Optional population name(s) to retain. Population is
  currently recorded as `NA` for intermediate diagnostics.

- marker:

  character. Optional marker name(s) to retain.

- chnl:

  character. Optional channel name(s) to retain.

- save:

  logical. If TRUE, save the detailed table as an RDS file.

- pathSave:

  character. Optional path for the saved RDS file. Defaults to
  `file.path(pathProject, "gatesDetailed.rds")`.

## Value

A tibble with one row per saved threshold diagnostic.
