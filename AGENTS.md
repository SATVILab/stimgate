# AGENTS.md — Configuration for AI Coding Agents

This file is the **canonical source of truth** for all AI coding agents
(e.g. Google Jules, GitHub Copilot) working on the `stimgate`
repository.

> **Instruction update rule for AI agents:** Always update `AGENTS.md`
> when you identify new coding patterns, guidelines, or best practices
> during issue resolution. This ensures all agents benefit from
> learnings and instructions remain unified and current.

------------------------------------------------------------------------

## 1. Core Philosophy / Project Context

`stimgate` is an R package for flow cytometry analysis, intended for
eventual submission to [Bioconductor](https://bioconductor.org/). It
identifies cells that have possibly responded to immune stimulation by
applying outlier-based gating to flow cytometry data. The core idea is
to compare an *unstimulated* tube with a *stimulated* tube from the same
donor sample and flag cells whose marker expression in the stimulated
condition is unusually high relative to the unstimulated background.

**Key architectural patterns:**

- The main user entry point is
  [`gateStim()`](https://satvilab.github.io/stimgate/reference/gateStim.md),
  which writes intermediate results to a `pathProject` directory on disk
  and returns that path.
- Downstream helpers
  ([`getStimGates()`](https://satvilab.github.io/stimgate/reference/getStimGates.md),
  [`getStimGatesDetailed()`](https://satvilab.github.io/stimgate/reference/getStimGatesDetailed.md),
  [`plotStim()`](https://satvilab.github.io/stimgate/reference/plotStim.md),
  [`writeStimFCS()`](https://satvilab.github.io/stimgate/reference/writeStimFCS.md))
  all accept `pathProject` and read from that directory.
- Internal functions are prefixed with `.` and are not exported.
- The package integrates tightly with the `flowCore` / `flowWorkspace`
  Bioconductor ecosystem; input data are `GatingSet` objects.
- `renv` is used for reproducible dependency management with two
  profiles: `bioc_container` (Bioconductor Docker) and
  `non_bioc_container` (standard R).

------------------------------------------------------------------------

## 2. Tech Stack & Tooling

| Layer | Tool / Package |
|----|----|
| Language | R (≥ 4.4.0) |
| Package framework | `devtools`, `roxygen2`, `testthat` (≥ 3.0.0) |
| Documentation | `roxygen2`, `pkgdown` |
| Flow cytometry I/O | `flowCore`, `flowWorkspace` |
| Data manipulation | `dplyr`, `purrr`, `tidyr`, `tibble`, `stringr`, `rlang` |
| Plotting | `ggplot2`, `cowplot` |
| Statistical modelling | `scam`, `mgcv` |
| Clustering | `cluster`, `gtools` |
| Dependency management | `renv` |
| CI | GitHub Actions (R-CMD-check, pkgdown, Codecov) |

**Do NOT:**

- Use `@import` or `@importFrom` directives in roxygen comments;
  explicitly qualify all package calls with `pkg::fun()`. The only
  exceptions are `ggplot2` functions and
  [`flowCore::exprs`](https://rdrr.io/pkg/Biobase/man/exprs.html), which
  may be called without a namespace qualifier and do not require
  `@importFrom` tags.
- Modify `.Rd` files manually; regenerate them with
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
- Use [`return()`](https://rdrr.io/r/base/function.html) as the last
  line of a function; use it only for early returns.
- Add [`library()`](https://rdrr.io/r/base/library.html) calls inside
  test files.

------------------------------------------------------------------------

## 3. Environment Setup & Dependency Management

Run package commands from the **repository root** (the directory
containing `DESCRIPTION`).

### Local development / ordinary agent environments

Outside CI and the GitHub Copilot cloud agent, `.Rprofile` selects the
appropriate `renv` profile and activates it automatically when R starts
in the project directory. If the project library needs to be restored,
run:

``` r

# 1. Install renv (if not already installed)
install.packages("renv")

# 2. Restore the active renv profile
renv::restore()

# 3. Load the package in your R session
devtools::load_all()
```

For a clean environment on a new machine you may also need to install
Bioconductor dependencies explicitly:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("flowCore", "flowWorkspace"))
```

### GitHub Copilot cloud agent

The Copilot cloud agent is set up by
`.github/workflows/copilot-setup-steps.yml`. That workflow installs R,
package and system dependencies, and the development tools needed for
package work. The repository `.Rprofile` deliberately treats the Copilot
agent as CI and does **not** activate `renv` there.

Do not run `renv::restore()` or require an activated `renv` project in
the Copilot cloud agent. Use the pre-installed development library
directly:

``` r

devtools::load_all()
devtools::test()
```

If these packages are unexpectedly unavailable in Copilot, treat that as
an environment-setup failure rather than switching to `renv`.

------------------------------------------------------------------------

## 4. Build, Test & Quality Instructions

Run the following commands from the **repository root** in an R session:

``` r

# Load and reload package
devtools::load_all()

# Run all unit tests
devtools::test()

# Regenerate documentation from roxygen comments
devtools::document()

# Style package code formatting
styler::style_pkg()

# Check for linting violations
lintr::lint_package()

# Full R CMD check (mimics CRAN / Bioconductor checks)
devtools::check()

# Check test coverage
covr::report()
```

Equivalent shell commands (via `Rscript`):

``` bash
Rscript -e "devtools::test()"
Rscript -e "devtools::check()"
Rscript -e "devtools::document()"
Rscript -e "styler::style_pkg()"
Rscript -e "lintr::lint_package()"
```

### Checklist before each commit / opening a PR

1.  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
2.  `styler::style_pkg()`
3.  `lintr::lint_package()`
4.  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)

### Website Maintenance (`pkgdown`)

Whenever functions are added, removed, or have their export status
changed (via `@export`), update `_pkgdown.yml` to ensure the reference
section accurately reflects all exported functions. Verify the site
configuration with:

``` r

pkgdown::check_pkgdown()
```

------------------------------------------------------------------------

## 5. Repository Structure

- `R/`: Core R source code for the package.
  - `UtilsCytoRSV-chnl_lab.R`: Channel label utilities (get
    markers/channels from cytometry objects).
  - `UtilsCytoRSV-plot_cyto.R`: Cytometry plotting utilities.
  - `UtilsGGSV-axisLimits.R`: `ggplot2` axis limit helpers.
  - `bw_norm_helpers.R`: Shared bandwidth helpers for standard and
    normalised bandwidth methods.
  - `check.R`: Input validation helpers.
  - `chnl_settings.R`: Complete channel parameter list with all required
    settings.
  - `cp-sub.R`: Auxiliary functions for getting clusters.
  - `cp_cluster-helper.R`: Helper functions for grouping thresholds from
    like distributions.
  - `cp_cluster.R`: Gets thresholds by grouping thresholds from like
    distributions.
  - `cp_uns_loc.R`: Gets thresholds by comparing stim and unstim
    distributions (main local-FDR entry point).
  - `cp_uns_loc_density.R`: Local-FDR density and raw probability
    estimation.
  - `cp_uns_loc_derivative.R`: Appendix derivative thresholds for
    local-FDR filtering.
  - `cp_uns_loc_filtering.R`: Local-FDR filtering after probability
    smoothing.
  - `cp_uns_loc_output.R`: Local-FDR diagnostics, metadata, and output
    assembly.
  - `cp_uns_loc_smoothing.R`: Local-FDR probability smoothing.
  - `cp_uns_loc_threshold.R`: Local-FDR response estimate and final
    threshold.
  - `cpp11.R`: Automatically generated C++ wrapper functions via
    `cpp11`.
  - `cyt_pos_gates-helper.R`: Helper functions for cytokine-positive
    cell gates.
  - `cyt_pos_gates.R`: Functions for more aggressive gates applied to
    cytokine-positive cells.
  - `debug.R`: Debugging utilities
    ([`.debug()`](https://satvilab.github.io/stimgate/reference/dot-debug.md))
    and global variable declarations.
  - `ex.R`: Extract expression matrices from `GatingSet` objects.
  - `example_data.R`: Creates example `GatingSet` for examples and
    testing
    ([`getExampleData()`](https://satvilab.github.io/stimgate/reference/getExampleData.md)).
  - `fcs_write.R`: Write FCS files of cytokine-positive cells
    (`writeStimFCS`).
  - `functionsForBenchmarking-Cyt.R`: Benchmarking helpers for cytokine
    simulation (used internally by `example_data.R`).
  - `gate.R`: Main entry point for gating (`gateStim`).
  - `gate_batch-helper.R`: Helper functions for gating batches of
    samples.
  - `gate_batch.R`: Gate batches of samples.
  - `gate_chnl-helper.R`: Helper functions for gating individual
    channels.
  - `gate_chnl.R`: Gate individual channels.
  - `gates.R`: Extract the identified gates/thresholds (`getStimGates`,
    `getStimGatesDetailed`).
  - `getCpTg_audit.R`: Audit helpers for `.getCpTg()` migration
    tracking.
  - `ind_batch.R`: Get the list of indices grouped by batch.
  - `peaks_and_troughs.R`: Peak and trough detection helpers.
  - `pipe.R`: Pipe operator and related utilities.
  - `plot_gate.R`: Plot the identified gates (`plotStim`).
  - `pos_ind.R`: Identify the indices of the cytokine-positive cells.
  - `stats-helper-overall.R`: Helper functions for overall statistics.
  - `stats-helper.R`: Helper functions for statistics.
  - `stats.R`: Get statistics for the identified gates.
  - `verify.R`: Input verification helpers.
- `scripts/`:
  - Shell scripts (`dev.sh`, `install.sh`, `patch.sh`, `minor.sh`,
    `major.sh`, `dev-*.sh`) for workflow, benchmarking, and version
    bumping.
  - `python/`: Python helper scripts used by analysis (not part of the R
    package).
    - `fbeta.py`: Richards F-beta thresholding implementation
      (comparison method).
  - `r/`: R helper scripts used by analysis only — not loaded by
    [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
    and not part of the `stimgate` namespace.
    - `functionsForBenchmarking-Pheno.R`: Benchmarking helpers for
      phenotype simulation.
    - `sim-bandwidth.R`: Simulation bandwidth utilities.
    - `sim-bw-adaptive.R`: Adaptive bandwidth simulation helpers.
    - `sim-compare-freq_bs.R`: Bootstrap frequency comparison for
      simulation.
    - `sim-misc.R`: Miscellaneous simulation utilities.
    - `sim-trans.R`: Simulation transformation utilities.
- `src/`: C++ source code compiled into the package via `cpp11`
  (`cpPmden.cpp`, `stimgate_cppmden.cpp`, `tautstring.cpp`, etc.).
- `analysis/`: Quarto (`.qmd`) documents for research, simulation, and
  benchmarking analysis.
- `vignettes/`: Package vignettes (`stimgate.Rmd`).
- `inst/`: Installed package material (e.g. `COPYRIGHTS`).
- `.github/`: GitHub CI workflows and Copilot setup.
- `data-raw/`: Raw data files used for testing and examples.
- `man/`: Automatically generated documentation files.
- `renv/`: R package environment management files (`renv.lock`,
  `.Rprofile`).
- `tests/`: Unit tests for the package (`tests/testthat/`).
- `_dependencies.R`: Explicitly listed dependencies for `renv`.
- `DESCRIPTION`: Package metadata file.

------------------------------------------------------------------------

## 6. Coding Style & Conventions

### Naming conventions

- **Exported functions**: `camelCase`, no leading dot (e.g. `gateStim`,
  `plotStim`).
- **Internal functions**: `camelCase` with a leading dot
  (e.g. `.getThreshold`). Begin each internal function with a `.`.

### Debugging & Intermediate Data Saving

- Use `.debug(msg, val)` inside internal functions for debug messages.
  Debug output is written to a temp file and controlled by the
  `STIMGATE_DEBUG` environment variable (set to `"true"`, `"yes"`,
  `"y"`, or `"1"` to enable). Do **NOT** add a debug flag or parameter
  to function signatures.
- Use the `stage` parameter to track algorithm stages (`"init"`,
  `"cytPos"`, or `"single"`). Pass `stage` through function calls to
  enable intermediate data saving via `.intSave()` or `.intSaveNm()`
  functions. Intermediate saving is controlled by the
  `STIMGATE_INTERMEDIATE` environment variable.

### Function Signatures & Returns

- Validate inputs and provide meaningful error messages.
- Never use [`return()`](https://rdrr.io/r/base/function.html) as the
  last line of a function; use it only for early returns.

### Documentation (roxygen2)

- Every exported function must have `@param`, `@return`, and `@export`
  tags.
- Parameter docs must follow the exact format:
  `@param param_name <type> <description>. Default: <value>.`
  - Example: `@param pathProject character Path to project.`
  - Multiple types/options example: `logical or character`, or
    `"always", "never" or "automatic"`.
  - Stating default: `Default: "automatic".` or
    `Default is "automatic".`
- Use `@details` for complex explanations rather than overloading
  `@param`.
- Whitespace: Ensure no trailing whitespace on lines, and ensure files
  end with a final newline.

### Package Namespace

- Reference all external functions explicitly as `pkg::fun()`.
- Exceptions: `ggplot2` functions and
  [`flowCore::exprs`](https://rdrr.io/pkg/Biobase/man/exprs.html) may be
  called without a namespace qualifier and do not require `@importFrom`
  tags.

------------------------------------------------------------------------

## 7. Specific Package Policies & Design Notes

1.  **Taut-string density**: The piecewise-constant taut-string density
    used for antimode detection is provided by the internal helper
    [`.tautStringPmden()`](https://satvilab.github.io/stimgate/reference/dot-tautStringPmden.md)
    (in `cp_uns_loc_filtering.R`), which wraps the native FAUST-derived
    C++ implementation `stimgate_cpPmden()` compiled via `cpp11`
    (`src/stimgate_cppmden.cpp` and `src/cpPmden.cpp`).
2.  **Comparison code vs. package code**: `R/` contains only StimGate
    implementation code. Benchmark comparisons against the legacy
    `openCyto` tailgate call `openCyto:::.cytokineCutpoint()` from the
    `openCyto` package directly in `scripts/r/sim-compare-freq_bs.R` and
    `analysis/7-sim-compare-freq_bs.qmd`. Note that
    `R/functionsForBenchmarking-Cyt.R` is an exception: it stays in `R/`
    because
    [`getExampleData()`](https://satvilab.github.io/stimgate/reference/getExampleData.md)
    (an exported function) calls
    [`simCytExperiment()`](https://satvilab.github.io/stimgate/reference/simCytExperiment.md)
    internally.
3.  **Legacy comparator policy**: Tailgate comparator functions are
    invoked directly from the `openCyto` package via
    `openCyto:::.cytokineCutpoint()`. Do not reintroduce vendored legacy
    tailgate helpers under `scripts/r/` or `R/`.
4.  **Audit status for `.getCpTg()` migration (issue \#157/#158)**:
    Remaining call sites are catalogued by `.get_cp_tg_call_audit()` and
    summarised by `.get_cp_tg_migration_note_157()`. Current default
    behaviour still constructs `tgClust` control gates in
    `.gateBatchAll()`, but the current local-FDR cluster quantile
    implementation does not consume `gateTblCtrl`, so this branch is
    dead plumbing for current outputs. Single-positive gating branches
    have been removed per issue \#196.

------------------------------------------------------------------------

## 8. Testing Best Practices & Guidelines

Place all unit tests in `tests/testthat/` as `test-<topic>.R` files.

1.  **Avoid [`library()`](https://rdrr.io/r/base/library.html) calls in
    test files**: Never use
    [`library(testthat)`](https://testthat.r-lib.org) or similar calls
    at the top of test files. The `testthat` package is automatically
    loaded when tests are run.
2.  **Keep tests independent**: Each test should be self-contained and
    not rely on global state created outside of test blocks.
3.  **Variable scope in tests**: When a test creates its own test data
    and variables (e.g. `example_data`, `gs`), always use those local
    variables throughout that test. Never mix local and global variables
    from different scopes.
4.  **Avoid code outside test blocks**: Do not place code outside
    `test_that()` blocks (except shared setup data required by every
    test in the file). Operations like
    [`unlink()`](https://rdrr.io/r/base/unlink.html) outside blocks
    execute in unpredictable order, cause race conditions, and interfere
    with parallel execution.
5.  **Cleanup within tests**: Each test must clean up its own temporary
    files/directories created during execution (e.g.,
    `unlink(tmp_dir, recursive = TRUE)` or
    [`withr::defer()`](https://withr.r-lib.org/reference/defer.html)).
6.  **Shared test fixtures**: If multiple tests need the same setup
    data, create it within each test or create it once at the top with
    clear documentation. Never delete shared fixtures mid-file.
7.  **Test data files compatibility**: Test data files (`.rds` in
    `tests/testthat/`) may need regeneration when major dependencies
    (e.g., `ggplot2`) upgrade. Regenerate test data files using current
    package versions if objects behave unexpectedly after dependency
    updates.
8.  **Test observable behaviour**: Tests must verify observable outputs
    and behavior, not internal implementation details or the existence
    of internal (`.`-prefixed) functions.
9.  **Cross-platform compatibility**: Tests must pass on macOS, Windows,
    and Ubuntu. Use
    [`file.path()`](https://rdrr.io/r/base/file.path.html) (never
    hard-coded `/` or `\\` separators) and avoid platform-specific
    paths.
