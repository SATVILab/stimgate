# Manual inspection workflows

These scripts are for reproducible, exploratory inspection of the current StimGate implementation. They sit outside `tests/testthat/`, so `devtools::test()` does not run them.

Run them from the repository root. Each script loads the current checkout with `devtools::load_all()` and then uses package functions and saved intermediates rather than copying the algorithm into the script.

For example:

```r
source("tests/manual/cytokine-positive-thresholding.R")
```

The scripts may print small tables and plots, and may leave temporary intermediate files under the session `tempdir()` so they can be inspected after the run.
