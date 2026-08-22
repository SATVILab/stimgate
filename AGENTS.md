# AGENTS.md — Configuration for AI Coding Agents

This file is the **canonical source of truth** for all AI coding agents
(e.g. Google Jules, GitHub Copilot) working on the `stimgate` repository.

> **Instruction update rule for AI agents:**
> Always update `AGENTS.md` when you identify new coding patterns, guidelines,
> or best practices during issue resolution. This ensures all agents benefit from
> learnings and instructions remain unified and current.

---

## 1. Core Philosophy / Project Context

`stimgate` is an R package for flow cytometry analysis, intended for eventual
submission to [Bioconductor](https://bioconductor.org/). It identifies cells that
have possibly responded to immune stimulation by applying outlier-based gating
to flow cytometry data. The core idea is to compare an *unstimulated* tube with
a *stimulated* tube from the same donor sample and flag cells whose marker
expression in the stimulated condition is unusually high relative to the
unstimulated background.

**Key architectural patterns:**

- The main user entry point is `gateStim()`, which writes intermediate
  results to a `pathProject` directory on disk and returns that path.
- Downstream helpers (`getStimGates()`, `getStimGatesDetailed()`,
  `plotStim()`, `writeStimFCS()`) all accept `pathProject` and read from
  that directory.
- Internal functions are prefixed with `.` and are not exported.
- The package integrates tightly with the `flowCore` / `flowWorkspace`
  Bioconductor ecosystem; input data are `GatingSet` objects.
- `renv` is used for reproducible dependency management with two profiles:
  `bioc_container` (Bioconductor Docker) and `non_bioc_container` (standard R).

---

## 2. Tech Stack & Tooling

| Layer | Tool / Package |
|---|---|
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

- Use `@import` or `@importFrom` directives in roxygen comments; explicitly
  qualify all package calls with `pkg::fun()`. The only exceptions are
  `ggplot2` (imported wholesale via `#' @import ggplot2` in `R/misc.R`) and
  `flowCore::exprs`, which may be called without a namespace qualifier and do
  not require `@importFrom` tags.
- Modify `.Rd` files manually; regenerate them with `devtools::document()`.
- Use `return()` as the last line of a function; use it only for early returns.
- Add `library()` calls inside test files.

---

## 3. Environment Setup & Dependency Management

Run package commands from the **repository root** (the directory containing
`DESCRIPTION`).

### Local development / ordinary agent environments

Outside CI and the GitHub Copilot cloud agent, `.Rprofile` selects the
appropriate `renv` profile and activates it automatically when R starts in the
project directory. If the project library needs to be restored, run:

```r
# 1. Install renv (if not already installed)
install.packages("renv")

# 2. Restore the active renv profile
renv::restore()

# 3. Load the package in your R session
devtools::load_all()
```

For a clean environment on a new machine you may also need to install
Bioconductor dependencies explicitly:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("flowCore", "flowWorkspace"))
```

### GitHub Copilot cloud agent

The Copilot cloud agent is set up by
`.github/workflows/copilot-setup-steps.yml`. That workflow installs R, package
and system dependencies, and the development tools needed for package work.
The repository `.Rprofile` deliberately treats the Copilot agent as CI and
does **not** activate `renv` there.

Do not run `renv::restore()` or require an activated `renv` project in the
Copilot cloud agent. Use the pre-installed development library directly:

```r
devtools::load_all()
devtools::test()
```

If these packages are unexpectedly unavailable in Copilot, treat that as an
environment-setup failure rather than switching to `renv`.

---

## 4. Build, Test & Quality Instructions

Run the following commands from the **repository root** in an R session:

```r
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

```bash
Rscript -e "devtools::test()"
Rscript analysis/tests/run_analysis_tests.R
Rscript -e "devtools::check()"
Rscript -e "devtools::document()"
Rscript -e "styler::style_pkg()"
Rscript -e "lintr::lint_package()"
```

### Checklist before each commit / opening a PR

1. `devtools::document()`
2. `styler::style_pkg()`
3. `lintr::lint_package()`
4. `devtools::test()`
5. If `analysis/` or `scripts/r/` changed, `Rscript analysis/tests/run_analysis_tests.R`

### Analysis / Repository Integration Tests

In addition to the package test suite (`tests/testthat/`), there is a
separate analysis integration test suite in `analysis/tests/testthat/`.

**When to use which suite:**

| Test type | Location | Run with |
|---|---|---|
| Package unit/integration tests | `tests/testthat/` | `devtools::test()` |
| Analysis helper / scripts/r / QMD drift checks | `analysis/tests/testthat/` | see below |

**What belongs in the analysis test suite (`analysis/tests/testthat/`):**

- Checks that `scripts/r/` helpers source cleanly in dependency order.
- Checks that QMD files do not call `scripts/r` helpers through `stimgate:::`.
- Checks that analysis wrapper parameters forwarded to `gateStim()` still exist
  in the current package API.
- Checks that removed arguments (e.g. `calcSinglePosGates`) are not reintroduced.
- Smoke calls for representative `.simBandwidth*()` / comparison-wrapper functions.
- Numerical agreement between `.simBandwidthBwOne()` and `stimgate:::.bwCalcOne()`.

**What belongs in the package test suite (`tests/testthat/`):**

- Tests of exported package functions.
- Tests of internal package functions where drift would break the package.
- Do **not** place tests whose subject is solely `scripts/r/` or `analysis/` code here.

**Running the analysis test suite:**

The analysis suite loads the package from the current checkout via
`devtools::load_all()` before running tests, so it always tests the current
source rather than any previously installed version. The GitHub Actions job in
`.github/workflows/analysis-integration.yaml` runs this suite independently of
`R CMD check`/`devtools::test()`. Run from the repository root:

```r
# Shell
# Rscript analysis/tests/run_analysis_tests.R

# R session (after devtools::load_all())
devtools::load_all()
testthat::test_dir("analysis/tests/testthat")
```

Or via the runner script:

```r
source("analysis/tests/run_analysis_tests.R")
```

Equivalent shell command:

```bash
Rscript analysis/tests/run_analysis_tests.R
```

### Website Maintenance (`pkgdown`)

Whenever functions are added, removed, or have their export status changed (via
`@export`), update `_pkgdown.yml` to ensure the reference section accurately
reflects all exported functions. Verify the site configuration with:

```r
pkgdown::check_pkgdown()
```

---

## 5. Repository Structure

- `R/`: Core R source code for the installed package.
  - `UtilsCytoRSV-chnl_lab.R`: Channel label utilities (get markers/channels from cytometry objects).
  - `UtilsCytoRSV-plot_cyto.R`: Cytometry plotting utilities.
  - `UtilsGGSV-axisLimits.R`: `ggplot2` axis limit helpers.
  - `bw_norm_helpers.R`: Shared bandwidth helpers for standard and normalised bandwidth methods.
  - `check.R`: Input validation helpers.
  - `chnl_settings.R`: Complete channel parameter list with all required settings.
  - `cp-sub.R`: Auxiliary functions for getting clusters.
  - `cp_cluster-helper.R`: Helper functions for grouping thresholds from like distributions.
  - `cp_cluster.R`: Gets thresholds by grouping thresholds from like distributions.
  - `cp_uns_loc.R`: Gets thresholds by comparing stim and unstim distributions (main local-FDR entry point).
  - `cp_uns_loc_density.R`: Local-FDR density and raw probability estimation.
  - `cp_uns_loc_derivative.R`: Appendix derivative thresholds for local-FDR filtering.
  - `cp_uns_loc_filtering.R`: Local-FDR filtering after probability smoothing.
  - `cp_uns_loc_output.R`: Local-FDR diagnostics, metadata, and output assembly.
  - `cp_uns_loc_smoothing.R`: Local-FDR probability smoothing.
  - `cp_uns_loc_threshold.R`: Local-FDR response estimate and final threshold.
  - `cpp11.R`: Automatically generated C++ wrapper functions via `cpp11`.
  - `cyt_pos_gates-helper.R`: Helper functions for cytokine-positive cell gates.
  - `cyt_pos_gates.R`: Functions for more aggressive gates applied to cytokine-positive cells.
  - `debug.R`: Debugging utilities (`.debug()`) and global variable declarations.
  - `ex.R`: Extract expression matrices from `GatingSet` objects.
  - `example_data.R`: Loads the canonical shipped example dataset via `getExampleData()` and does not simulate new data.
  - `fcs_write.R`: Write FCS files of cytokine-positive cells (`writeStimFCS`).
  - `gate.R`: Main entry point for gating (`gateStim`).
  - `gate_batch-helper.R`: Helper functions for gating batches of samples.
  - `gate_batch.R`: Gate batches of samples.
  - `gate_chnl-helper.R`: Helper functions for gating individual channels.
  - `gate_chnl.R`: Gate individual channels.
  - `gates.R`: Extract the identified gates/thresholds (`getStimGates`, `getStimGatesDetailed`).
  - `getCpTg_audit.R`: Audit helpers for `.getCpTg()` migration tracking.
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
  - Shell scripts (`dev.sh`, `install.sh`, `patch.sh`, `minor.sh`, `major.sh`, `dev-*.sh`) for workflow, benchmarking, and version bumping.
  - `python/`: Python helper scripts used by analysis (not part of the R package).
    - `fbeta.py`: Richards F-beta thresholding implementation (comparison method).
  - `r/`: Developer-side R analysis/simulation helpers used for research, benchmarking, and fixture regeneration. These are not loaded by `devtools::load_all()` and are not part of the installed package.
  - `analysis-runtime.R`: Shared QMD execution/runtime plumbing for parameter lookup, env overrides, chunk validation and atomic RDS output.
  - `functionsForBenchmarking-Cyt.R`: Cytokine simulation utilities.
  - `functionsForBenchmarking-Pheno.R`: Benchmarking helpers for phenotype simulation.
  - `sim-bandwidth.R`: Simulation bandwidth utilities.
  - `sim-compare-freq_bs.R`: Bootstrap frequency comparison for simulation.
  - `sim-misc.R`: Miscellaneous simulation utilities.
  - `sim-trans.R`: Simulation transformation utilities.
- `src/`: C++ source code compiled into the package via `cpp11` (`cpPmden.cpp`, `stimgate_cppmden.cpp`, `tautstring.cpp`, etc.).
- `analysis/`: Quarto (`.qmd`) documents for research, simulation, and benchmarking analysis.
- `vignettes/`: Package vignettes (`stimgate.Rmd`).
- `inst/extdata/`: Canonical saved example datasets consumed by `getExampleData()` and by package examples/tests.
- `inst/`: Installed package material (e.g. `COPYRIGHTS`).
- `.github/`: GitHub CI workflows and Copilot setup.
- `data-raw/`: Raw data files used for testing and examples.
- `man/`: Automatically generated documentation files.
- `renv/`: R package environment management files (`renv.lock`, `.Rprofile`).
- `tests/`: Unit tests for the package (`tests/testthat/`).
- `_dependencies.R`: Explicitly listed dependencies for `renv`.
- `DESCRIPTION`: Package metadata file.

### Analysis code layering

For new or moved analysis code, use this layering:

1. `R/`: installed StimGate package implementation only.
2. `scripts/r/sim-*.R`: developer-side simulation/domain helpers and wrappers around the current package API.
3. Generic analysis runtime helpers under `scripts/r/`: reusable QMD execution plumbing such as parameter/environment handling, chunk validation and atomic output writing.
4. Analysis-specific helpers under `scripts/r/`: substantial orchestration, restart/collation, IO and plotting helpers that should not live inline in QMDs.
5. `analysis/*.qmd`: scientific settings, analysis calls, result-specific transformations and presentation.

Source analysis helper files explicitly in dependency order. Do not move analysis-only
helpers into `R/` unless they have genuinely become part of the installed package
implementation or API. Keep large domain helper files such as `sim-bandwidth.R`
focused on their domain rather than using them as catch-all locations for generic
QMD runtime or unrelated plotting/orchestration code.

---

## 6. Coding Style & Conventions

### Naming conventions

- **Exported functions**: `camelCase`, no leading dot (e.g. `gateStim`, `plotStim`).
- **Internal functions**: `camelCase` with a leading dot (e.g. `.getThreshold`). Begin each internal function with a `.`.

### Debugging & Intermediate Data Saving

- Use `.debug(msg, val)` inside internal functions for debug messages. Debug output
  is written to a temp file and controlled by the `STIMGATE_DEBUG` environment
  variable (set to `"true"`, `"yes"`, `"y"`, or `"1"` to enable).
  Do **NOT** add a debug flag or parameter to function signatures.
- Use the `stage` parameter to track algorithm stages (`"init"`, `"cytPos"`, or
  `"single"`). Pass `stage` through function calls to enable intermediate data
  saving via `.intSave()` or `.intSaveNm()` functions. Intermediate saving is
  controlled by the `STIMGATE_INTERMEDIATE` environment variable.

### Function Signatures & Returns

- Validate inputs and provide meaningful error messages.
- Never use `return()` as the last line of a function; use it only for early returns.

### Documentation (roxygen2)

- Every exported function must have `@param`, `@return`, and `@export` tags.
- Parameter docs must follow the exact format:
  `@param param_name <type> <description>. Default: <value>.`
  - Example: `@param pathProject character Path to project.`
  - Multiple types/options example: `logical or character`, or `"always", "never" or "automatic"`.
  - Stating default: `Default: "automatic".` or `Default is "automatic".`
- Use `@details` for complex explanations rather than overloading `@param`.
- Whitespace: Ensure no trailing whitespace on lines, and ensure files end with a final newline.

### Package Namespace

- Reference all external functions explicitly as `pkg::fun()`.
- Exceptions: `ggplot2` is imported wholesale via `#' @import ggplot2` in
  `R/misc.R`, so `ggplot2` functions and `flowCore::exprs` may be called without
  a namespace qualifier and do not require `@importFrom` tags.

---

## 7. Specific Package Policies & Design Notes

1. **Taut-string density**:
   The piecewise-constant taut-string density used for antimode detection is
   provided by the internal helper `.tautStringPmden()` (in
   `cp_uns_loc_filtering.R`), which wraps the native FAUST-derived C++ implementation
   `stimgate_cpPmden()` compiled via `cpp11` (`src/stimgate_cppmden.cpp` and `src/cpPmden.cpp`).
2. **Comparison code vs. package code**:
   `R/` contains only StimGate implementation code. Benchmark comparisons against
   the tailgate method call `cytoUtils:::.cytokine_cutpoint()` from the
   `cytoUtils` package are implemented directly in `scripts/r/sim-compare-freq_bs.R`
   and `analysis/7-sim-compare-freq_bs.qmd`. Cytokine simulation logic remains in
   `scripts/r/` and is not installed with the package.
3. **Legacy comparator policy**:
   Tailgate comparator functions are invoked directly from the `cytoUtils` package via
   `cytoUtils:::.cytokine_cutpoint()`. Do not reintroduce vendored legacy tailgate helpers
   under `scripts/r/` or `R/`.
4. **Temporary migration status for `.getCpTg()` (issues #157/#158)**:
   This is a current-state note rather than a permanent design rule. Verify it against
   the current implementation and relevant issues before relying on it in later work.
   At the time of this update, remaining call sites are catalogued by
   `.get_cp_tg_call_audit()` and summarised by `.get_cp_tg_migration_note_157()`.
   Current default behaviour still constructs `tgClust` control gates in
   `.gateBatchAll()`, but the current local-FDR cluster quantile implementation does
   not consume `gateTblCtrl`, so this branch is dead plumbing for current outputs.
   Single-positive gating branches have been removed per issue #196.
5. **Simulation engine migration to `simcyto` (issue #288 / umbrella #271)**:
   Generic cytometry simulations and post-simulation transformations are progressively
   migrating to the exported `simcyto` package API (e.g. `simcyto::simCytExperiment()`,
   `simcyto::simCytTransform*()`). `analysis/3-sim-bw-est-base.qmd` uses `simcyto` and does
   not source `functionsForBenchmarking-Cyt.R`. StimGate scientific scenario calculations
   and bandwidth/gating orchestration remain StimGate-side under `scripts/r/`.

---

## 8. Testing Best Practices & Guidelines

Package unit/integration tests belong in `tests/testthat/`. Tests whose subject is
analysis code, `scripts/r/` helpers or QMD/package-API drift belong in
`analysis/tests/testthat/`. Use `test-<topic>.R` filenames in both suites.

1. **Avoid `library()` calls in test files**:
   Never use `library(testthat)` or similar calls at the top of test files.
   The `testthat` package is automatically loaded when tests are run.
2. **Keep tests independent**:
   Each test should be self-contained and not rely on global state created outside
   of test blocks.
3. **Variable scope in tests**:
   When a test creates its own test data and variables (e.g. `example_data`, `gs`),
   always use those local variables throughout that test. Never mix local and global
   variables from different scopes.
4. **Avoid code outside test blocks**:
   Do not place code outside `test_that()` blocks (except shared setup data required
   by every test in the file). Operations like `unlink()` outside blocks execute in
   unpredictable order, cause race conditions, and interfere with parallel execution.
5. **Cleanup within tests**:
   Each test must clean up its own temporary files/directories created during execution
   (e.g., `unlink(tmp_dir, recursive = TRUE)` or `withr::defer()`).
6. **Shared test fixtures**:
   If multiple tests need the same setup data, create it within each test or create it
   once at the top with clear documentation. Never delete shared fixtures mid-file.
7. **Test data files compatibility**:
   Test data files (`.rds` in `tests/testthat/`) may need regeneration when major
   dependencies (e.g., `ggplot2`) upgrade. Regenerate test data files using current
   package versions if objects behave unexpectedly after dependency updates.
8. **Test observable behaviour and explicit integration contracts**:
   Package tests should verify observable outputs and behaviour rather than merely
   asserting implementation details or the existence of internal (`.`-prefixed)
   functions. Analysis integration tests may directly check helper/API contracts when
   the purpose is to catch drift between `scripts/r/`, QMDs and the installed package.
9. **Cross-platform compatibility**:
   Tests must pass on macOS, Windows, and Ubuntu. Use `file.path()` (never hard-coded
   `/` or `\\` separators) and avoid platform-specific paths.
10. **Use the package-shipped example data for routine tests and examples**:
    The package ships one canonical deterministic cytometry example dataset in
    `inst/extdata/stimgate_example_data/` (2 samples × 2 conditions × 2 markers ×
    ~10,000 cells per condition, seed 42). Load it with:
    ```r
    exampleData <- getExampleData()
    ```
    Package tests and examples should use `getExampleData()` rather than sourcing
    developer-side simulation utilities. Keep simulation code in `scripts/r/` for
    deliberate analysis/fixture-generation work only. To regenerate the dataset
    after intentional changes to its structure, run
    `source("data-raw/create_test_fixture.R")` from the repository root in a
    clean R session (no `devtools::load_all()` required).
