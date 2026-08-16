# stimgate

This is an R package for flow cytometry analysis, intended for eventual submission to BioConductor.
The package identifies cells that have possibly responded to stimulation by comparing the unstimulated and stimulated tubes from the same sample using outlier-based gating methods.

## Package Overview

`stimgate` applies statistical methods to flow cytometry data from the `flowCore` and `flowWorkspace` ecosystem to identify stimulation-specific cell responses. The package is particularly useful for analyzing immune cell activation assays where stimulated and unstimulated conditions are compared.

## Code standards

### Required before each commit

- Run `devtools::document()` to update documentation
- Run `styler::style_pkg()` to ensure consistent code formatting
- Run `devtools::test()` to ensure all tests pass
- Run `lintr::lint_package()` to check for linting violations

### Development flow

- Build and reload: Use `devtools::load_all()` to reload the package in your R session
- Test: Run `devtools::test()` to execute all unit tests
- Documentation: Use `devtools::document()` to update package documentation
- Coverage: Run `covr::report()` to check test coverage
- Full check: Run `devtools::check()` to perform a comprehensive package check

### Environment setup

This package uses `renv` for dependency management:
- The `.Rprofile` file automatically activates `renv` when you start R in the project directory
- Always ensure your working directory is the project root (containing `DESCRIPTION`) when running R code
- Dependencies are locked in `renv.lock` and managed through `renv::restore()`
- If dependencies are missing, run `renv::restore()` to install them

### Test data files

- Test data files (`.rds` files in `tests/testthat/`) may need to be regenerated when major package dependencies are updated, especially when upgrading to new major versions of packages like `ggplot2`
- If tests fail with objects behaving unexpectedly after a dependency upgrade, check if the test data was created with an older version and needs regeneration
- Regenerate test data files using the current package versions to ensure compatibility

## Repository structure

- `R/`: Core R source code
  - `UtilsCytoRSV-chnl_lab.R`: Channel label utilities (get markers/channels from cytometry objects)
  - `UtilsCytoRSV-plot_cyto.R`: Cytometry plotting utilities
  - `UtilsGGSV-axis_limits.R`: ggplot2 axis limit helpers
  - `bw_norm_helpers.R`: Shared bandwidth helpers for standard and normalised bandwidth methods
  - `check.R`: Input validation helpers
  - `chnl_settings.R`: Complete channel parameter list with all required settings
  - `cp-sub.R`: Auxiliary functions for getting clusters
  - `cp_cluster-helper.R`: Helper functions for grouping thresholds from like distributions
  - `cp_cluster.R`: Gets thresholds by grouping thresholds from like distributions
  - `cp_uns_loc.R`: Gets thresholds by comparing stim and unstim distributions (main local-FDR entry point)
  - `cp_uns_loc_density.R`: Local-FDR density and raw probability estimation
  - `cp_uns_loc_derivative.R`: Appendix derivative thresholds for local-FDR filtering
  - `cp_uns_loc_filtering.R`: Local-FDR filtering after probability smoothing
  - `cp_uns_loc_output.R`: Local-FDR diagnostics, metadata, and output assembly
  - `cp_uns_loc_smoothing.R`: Local-FDR probability smoothing
  - `cp_uns_loc_threshold.R`: Local-FDR response estimate and final threshold
  - `cyt_pos_gates-helper.R`: Helper functions for cytokine-positive cell gates
  - `cyt_pos_gates.R`: Functions for more aggressive gates applied to cytokine-positive cells
  - `cytoUtils-cytokine_cutpoint.R`: Cytokine cutpoint utilities — part of the StimGate package algorithm (used by `cp-sub.R`); not comparison-only despite the filename
  - `debug.R`: Debugging utilities (`.debug()`) and global variable declarations
  - `ex.R`: Extract expression matrices from GatingSets
  - `example_data.R`: Creates example GatingSet for examples and testing
  - `fcs_write.R`: Write FCS files of cytokine-positive cells (`writeStimFCS`)
  - `functionsForBenchmarking-Cyt.R`: Benchmarking helpers for cytokine simulation (used internally by `example_data.R` to create the example GatingSet)
  - `gate.R`: Main entry point for gating (`gateStim`)
  - `gate_batch-helper.R`: Helper functions for gating batches of samples
  - `gate_batch.R`: Gate batches of samples
  - `gate_chnl-helper.R`: Helper functions for gating individual channels
  - `gate_chnl.R`: Gate individual channels
  - `gates.R`: Extract the identified gates/thresholds (`getStimGates`, `getStimGatesDetailed`)
  - `ind_batch.R`: Get the list of indices grouped by batch
  - `openCyto-find_peaks_and_valleys.R`: Peak and valley finding utilities (from openCyto) — part of the StimGate package algorithm (used by `.cytokineCutpoint` in `cytoUtils-cytokine_cutpoint.R`); not comparison-only despite the filename
  - `peaks_and_troughs.R`: Peak and trough detection helpers
  - `pipe.R`: Pipe operator and related utilities
  - `plot_gate.R`: Plot the identified gates (`plotStim`)
  - `pos_ind.R`: Identify the indices of the cytokine-positive cells
  - `stats-helper-overall.R`: Helper functions for overall statistics
  - `stats-helper.R`: Helper functions for statistics
  - `stats.R`: Get statistics for the identified gates
  - `verify.R`: Input verification helpers
- `scripts/`:
  - `python/`: Python helper scripts used by analysis (not part of the R package)
    - `fbeta.py`: Richards F-beta thresholding implementation (comparison method)
  - `r/`: R helper scripts used by analysis only — not loaded by `devtools::load_all()` and not part of the stimgate namespace
    - `functionsForBenchmarking-Pheno.R`: Benchmarking helpers for phenotype simulation
    - `sim-bandwidth.R`: Simulation bandwidth utilities (sourced explicitly by `analysis/2-*.qmd`, `3-*.qmd`, `4-*.qmd`, `5-*.qmd`, `6-*.qmd` and their split variants)
    - `sim-bw-adaptive.R`: Adaptive bandwidth simulation helpers
    - `sim-compare-freq_bs.R`: Bootstrap frequency comparison for simulation (comparison analysis only; sourced explicitly by `analysis/7-sim-compare-freq_bs.qmd`)
    - `sim-misc.R`: Miscellaneous simulation utilities (sourced explicitly by all simulation analysis files)
    - `sim-trans.R`: Simulation transformation utilities (sourced explicitly by `analysis/1-sim-trans.qmd`)
- `_reference/`:
  - Reference material for the package, which at this stage is just old scripts that I didn't want to completely delete at the time.
- `.github/`: GitHub associated files
- `data-raw/`: Raw data files used for testing and examples
- `man/`: Automatically generated documentation files
- `renv/`: R package environment management files
- `tests/`: Unit tests for the package
- `_dependencies.R`: Explicitly listed dependencies for `renv` to pick up via `implicit` dependencies
- `.gitignore`: Git ignore file to exclude unnecessary files from version control
- `.Rbuildignore`: R build ignore file to exclude unnecessary files from package builds
- `DESCRIPTION`: Package metadata file
- `LICENSE`: License file for the package
- `LICENSE.md`: License file in Markdown format
- `README.md`: Package overview and usage documentation
- `renv.lock`: Lock file for the `renv` package management system, capturing the state of the R package environment

## Key Guidelines

1. Begin each internal function with a `.`
2. Use the `.debug()` function for debug messages, which takes a message and an optional value. Debug output is written to a temp file and controlled by the `STIMGATE_DEBUG` environment variable (set to "true", "yes", "y", or "1" to enable). Never add a debug flag or parameter to function signatures.
3. Use the `stage` parameter to track algorithm stages ("init", "cytPos", or "single"). Pass the stage parameter through function calls to enable intermediate data saving via `.intSave()` or `.intSaveNm()` functions. Intermediate saving is controlled by the `STIMGATE_INTERMEDIATE` environment variable.
4. Add unit tests using `testthat` for all new functionality
5. Validate inputs and provide meaningful error messages
6. Explicitly refer to all packages used, rather than using `@import` or `@importFrom`, with the exception of `ggplot2` functions and the `flowCore::exprs` function
7. Use `@export` for functions that should be available to users
8. When running code from the project, you must always have as your working directory the root of the project, i.e. the directory containing the `DESCRIPTION` file. This is especially important when the project uses `renv`, as otherwise the `.Rprofile` will not be sourced and the package environment will not be set up correctly
9. Never update `.Rd` files manually; use `devtools::document()` to regenerate them
10. Documentation in `.R` files for parameters must always have the format `@param param_name <type_of_input> <info>`
  - For example, `@param pathProject character Path to project.`
  - Of course, the information can be longer, and where there are multiple options, it can be mentioned:
    - For example, `logical or character`, or `"always", "never" or "automatic"`
  - If there is a default specified, it should be stated at the end in the documentation of that parameter, e.g. `Default is "automatic".`
  - Where longer parameter descriptions are required, usually you should describe overall what the parameter is for, and then list the options with a short description of each
  - If there are really nitty-gritty details, then you can use `@details` to provide more information and refer to it in the parameter description
11. Always make sure there is no unnecessary trailing whitespace, including for blank lines
12. Always make sure that any edited files have a final newline at the end of the file
13. Tests written should not simply test whether it errors out correctly, but should also test that the output is as expected
14. Never use `return` for the last line of a function, but only when you want to return early from a function
15. For testing, never add source commands at the top to source files in the `R` folder, as these are automatically sourced (effectively) by the `testthat` package
16. The tests need to pass on Mac, Windows and Ubuntu, so be aware of that (e.g. file path availability and specification (forward or backward slashes, root directories, etc.), case sensitivity, etc.). You do not run the tests on all three platforms, but you should think about them passing on all three platforms
17. Tests should verify observable behavior and outputs, not internal implementation details. Avoid testing for the existence of internal functions or variables (those starting with `.`). Focus on testing that functions produce the correct output, handle errors appropriately, and have the expected side effects
18. Always update this copilot-instructions.md file when you identify new coding patterns, guidelines, or best practices during issue resolution. This ensures future agents benefit from your learnings and the instructions stay current with the evolving codebase
19. **pkgdown website maintenance**: Whenever functions are added, removed, or have their export status changed (via `@export`), update the `_pkgdown.yml` file to ensure the reference section accurately reflects all exported functions. The `_pkgdown.yml` file organizes functions into logical categories for the package website. After making changes, verify the pkgdown configuration is valid by running `pkgdown::check_pkgdown()` if pkgdown is available
20. **Taut-string density**: The piecewise-constant taut-string density used for antimode detection is provided by the internal helper `.tautStringPmden()` (in `cp_uns_loc_filtering.R`), which wraps `cytoUtils::tautstring()`. Do not re-introduce `ftnonpar` or vendor C++ source for this purpose. `cytoUtils` is listed under `Remotes: RGLab/cytoUtils` in DESCRIPTION. Note: the `cytoUtils` DESCRIPTION declares GPL-2, but some of its FAUST-derived C++ source files carry GPL-3-or-later notices — this is an upstream discrepancy; stimgate depends on the compiled package and does not vendor those source files.
21. **Comparison code vs package code**: `R/` contains only StimGate implementation code. Simulation helpers, benchmarking utilities, and other analysis-only code belong in `scripts/r/` and must be sourced explicitly by the relevant `analysis/*.qmd` or targets files. `analysis/7-sim-compare-freq_bs.qmd` is the sole consumer of `scripts/r/sim-compare-freq_bs.R`. The simulation bandwidth and misc helpers (`sim-bandwidth.R`, `sim-misc.R`, `sim-trans.R`, etc.) are sourced by the relevant analysis files via explicit `source()` calls after the `devtools::load_all()` block. Note that `R/functionsForBenchmarking-Cyt.R` is an exception: it stays in `R/` because `getExampleData()` (an exported function) calls `simCytExperiment()` internally.
22. **Follow-up task — remove cp-sub.R tailgate dependency**: `R/cp-sub.R` calls `.cytokineCutpoint()` (from `R/cytoUtils-cytokine_cutpoint.R`, which in turn uses `R/openCyto-find_peaks_and_valleys.R`). StimGate already has equivalent logic in `R/bw_norm_helpers.R` (`.bwNormFindRightFlat`) and `R/peaks_and_troughs.R`. Replacing the openCyto tailgate call in `cp-sub.R` with StimGate's own implementation is a separate follow-up task that requires careful testing before the openCyto-derived files can be removed from `R/`.
23. **Audit status for `.getCpTg()` migration (issue #157/#158)**: remaining call sites are catalogued by `.get_cp_tg_call_audit()` and summarised by `.get_cp_tg_migration_note_157()`. Current default behaviour still constructs `tgClust` control gates in `.gateBatchAll()`, but the current local-FDR cluster quantile implementation does not consume `gateTblCtrl`, so this branch is dead plumbing for current outputs. The other `.getCpTg()` callers are in legacy single-positive paths that are optional/backwards-compatible only.

## Testing Best Practices

When writing tests with `testthat`, follow these important practices:

1. **Avoid library() calls in test files**: Never use `library(testthat)` or similar calls at the top of test files. The testthat package is automatically loaded when tests are run.

2. **Keep tests independent**: Each test should be self-contained and not rely on global state created outside of test blocks. If you need setup data:
   - Create it within each test that needs it, OR
   - Use shared setup only when all tests in the file require it, AND
   - Ensure cleanup happens within the test or use `withr::defer()` for automatic cleanup

3. **Variable scope in tests**: When a test creates its own test data and variables (like `example_data`, `gs`, `pathProject_2`), always use those local variables throughout that test. Never mix local and global variables from different scopes.

4. **Avoid code outside test blocks**: Do not place cleanup code (like `unlink()`) or other operations outside of `test_that()` blocks. Such code:
   - Executes in an unpredictable order relative to tests
   - Can interfere with tests that run after it
   - Makes the test file harder to understand and maintain
   - May cause race conditions in parallel test execution

5. **Cleanup within tests**: Each test should clean up its own temporary files and directories:
   ```r
   test_that("example test", {
     temp_dir <- file.path(tempdir(), "test_output")
     # ... test code ...
     unlink(temp_dir, recursive = TRUE)
   })
   ```

6. **Shared test fixtures**: If multiple tests need the same setup data:
   - Either create the data in each test (preferred for independence)
   - OR create it once at the top of the file, but document that all tests depend on it
   - Never delete shared fixtures until the very end of the file (if at all)

7. **Test file structure**: A well-structured test file should have:
   - Optional shared setup code at the top (clearly documented)
   - All test cases in `test_that()` blocks
   - Each test responsible for its own cleanup
   - Optional final cleanup at the very end (only for shared resources)
