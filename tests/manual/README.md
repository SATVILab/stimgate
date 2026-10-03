# Manual inspection workflows

This directory contains small, reproducible scripts for inspecting current StimGate behaviour interactively. They are exploratory aids, not automated regression tests, and deliberately sit outside `tests/testthat/`.

Run scripts from the repository root. For example:

```r
source("tests/manual/inspect-cluster-threshold-sharing.R")
```

or, for console/table output without interactive plots:

```bash
Rscript tests/manual/inspect-cluster-threshold-sharing.R
```

Each script calls `devtools::load_all()` so it uses the current checkout. Keep inputs small and deterministic, and call current package functions or helpers rather than copying package algorithms. If an exploratory check becomes a stable behavioural contract, promote that invariant to an automated test under `tests/testthat/`.
