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
| Clustering | `cluster` |
| Dependency management | `renv` |
| CI | GitHub Actions (R-CMD-check, pkgdown, Codecov) |

**Do NOT:**

- Use `@import` or `@importFrom` directives in roxygen comments; explicitly
  qualify all package calls with `pkg::fun()`. The only exceptions are
  `ggplot2` (imported wholesale via `#' @import ggplot2` in `R/stimgate-package.R`) and
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

### GitHub Project administration

For GitHub issue or Project-administration work, use the shared
`github-projects` skill from
`MiguelRodo/github-projects-skill/skills/github-projects/` and read
`.projects/project.md`. Keep its lightweight Project environment separate from
the R package-development setup above.

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

### Checklist before opening a PR

1. `devtools::document()`
2. `styler::style_pkg()`
3. `lintr::lint_package()`
4. `devtools::test()`
5. If `analysis/` or `scripts/r/` changed, `Rscript analysis/tests/run_analysis_tests.R`

### Test runs while iterating

The full suites are slow. While iterating, run only the test files that cover
the code you changed, e.g. `devtools::test(filter = "cp_uns_loc|pos_ind")` or
`testthat::test_file()` for analysis tests. Run the full suite once, on the
finished change, before opening the PR; CI runs it again.

When several agents work in parallel (subagents, separate worktrees), each
agent runs targeted tests only and the coordinating agent runs the full suite
once on the combined result. Worktrees share one `git stash`, so parallel
agents must not use it; use a patch file or a temporary commit instead.

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
- Checks that QMD setup code does not overwrite sourced helper-function names
  with computed values; use a distinct variable name for helper return values.
- Checks that chunked simulation QMDs use per-scenario deterministic seeds and
  validate complete cross-chunk collation before promoting canonical results.
- Persist each per-scenario output atomically before writing its completed/error
  marker, and pass required run/chunk paths explicitly to progress helpers. This
  keeps restart markers consistent with durable output files.
- Estimator-comparison simulations should use the same simulated dataset for
  rows that differ only by estimator or estimator-tuning settings. Derive the
  data-generation seed from the biological scenario, not from the estimator,
  cap, chunk index or worker scheduling.
- Keep requested estimator settings distinct from the estimator path actually
  used after fallbacks. Preserve and summarise fallback provenance rather than
  labelling fallback rows as though they used the requested estimator.
- Active simulation chunks must collate only their own chunk outputs; canonical
  cross-chunk reads happen after promotion. A render with simulations disabled
  must use the read-only current-results context and must not create staging state.
- `run_plots = FALSE` must stop before optional plot/report chunks; multi-chunk
  simulation renders should not write shared plot files concurrently.
- Comparison analyses must fail before simulation when a required competitor
  dependency is unavailable. Do not let a missing package/script be converted
  into an algorithmic fallback and then score that fallback as a real method result.
- When plotting a summary over a simulation grid, every varying scenario
  dimension must be filtered, faceted or included in the plot grouping. Do not
  connect or aggregate distinct scenario settings into one line implicitly.
- Controlled mismatch/degradation simulations should use common random numbers
  within each baseline biological scenario when the mismatch itself is
  deterministic, so curve differences are not driven by different simulated draws.
- Compatibility wrappers for optional upstream features must apply an effect
  exactly once. If upstream supports the feature, pass it through without also
  applying a local fallback; otherwise neutralise the upstream global effect
  before applying the local selective fallback.
- Comparator exceptions in benchmarking analyses must remain explicit runtime
  errors. A numerical fallback may be retained for diagnostics, but the
  exception must not be silently promoted or scored as a valid prediction.
- Transactional simulation/collation chunks must not use Quarto
  `error: true`; validation and promotion errors must fail the render/job.
- For adaptive normalised bandwidth estimation, `normAdaptiveNcell` controls
  the fixed-size synthetic core/extra samples. Do not vary `bwNcellMax` as if
  it controlled that adaptive branch unless the implementation changes.
- When an estimator can legitimately fail to return a finite scientific
  estimate, retain that failure as analysis data (for example with
  `n_*_finite` / `prop_*_finite`) rather than hiding it behind a magic
  numeric fallback or averaging only successful estimates without reporting
  coverage. Distinguish estimator failure from infrastructure/runtime errors.
- Simulation wrappers that claim to mirror a current package calculation must
  use the same preprocessing as the package implementation. If a wrapper keeps
  a legacy preprocessing option for other analyses, set the current behaviour
  explicitly in the QMD rather than relying on the wrapper default.
- For end-to-end background-subtracted-frequency performance, score the final
  sample-level `loc_sample` `propRespEst` against `propRespTruth`.
  `propBsEst` is an internal local-FDR diagnostic used during threshold
  selection and must not silently replace the final frequency estimand.
- Preserve StimGate threshold provenance in method-comparison outputs. A finite
  high-value fallback is still a fallback: use `locGenerated`,
  `locGeneratedDirect`, `locSource` and `locReason` from the final gate
  table rather than inferring success solely from `is.finite(threshold)`.
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

### Continuous integration (GitHub Actions)

Pull-request CI is deliberately minimal and Windows-based, because Windows
installs CRAN and Bioconductor binaries while Ubuntu compiles the
`flowWorkspace` stack from source.

| Workflow | Pull requests | Otherwise |
|---|---|---|
| `R-CMD-check.yaml` | `windows-latest` (release) only | Same on pushes to `master`; full OS/R matrix on published releases, manual runs (`full` input) and PRs labelled `full-check` |
| `analysis-integration.yaml` | `windows-latest`, path-filtered | Same on pushes to `master` |
| `pkgdown.yaml` | Not run | Builds and deploys on `master`, releases and manual runs (Ubuntu) |
| `test-coverage.yaml` | Not run | `master` and manual runs (Windows) |
| `document.yaml` | Pushes touching `R/` (Windows) | Also runnable manually |

- Add the `full-check` label to a PR, or run R-CMD-check manually, when a
  change needs Linux, macOS or older-R coverage, and before releases.
- `analysis-integration.yaml` installs Python, `numpy` and `reticulate` and
  sets `RETICULATE_PYTHON`, because the F-beta comparator tests call
  `scripts/python/fbeta.py`.
- In CI, `.Rprofile` must keep preferring the `RSPM` repository URL exported by
  `r-lib/actions/setup-r`. pak resolves packages in a subprocess that skips the
  site profile but sources `.Rprofile`; falling back to the source-only
  `https://packagemanager.posit.co/cran/latest` there makes every package
  build from source.
- Never put comments inside a `setup-r-dependencies` `extra-packages: |`
  block: YAML keeps them as text and pak treats them as package names.
- `setup-r-dependencies` uses `cache: always` where a failing job should still
  save its package library for later runs.
- Ubuntu jobs depend on that cache: an uncached run compiles the Bioconductor
  `flowWorkspace` stack (RProtoBufLib, Rhdf5lib, cytolib, ...) from source in
  about 10 minutes, while a cached run installs dependencies in 1-4 minutes.
  Pull-request branches restore `master`'s caches.

---

## 5. Repository Structure

- `R/`: Core R source code for the installed package.
  - `UtilsCytoRSV-chnl_lab.R`: Channel label utilities (get markers/channels from cytometry objects).
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
  - `profile.R`: Structured debug timing, incremental persistence and final collation helpers.
  - `zz_profile_instrumentation.R`: Debug-only profiling wrappers around selected workflow boundaries.
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

Plot-construction helpers under `scripts/r/` should return plot objects without
creating directories or writing files. Keep filesystem side effects in the
corresponding save/orchestration helper or QMD.

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
  is written to `pathProject/debug/debug.txt` and controlled by the `STIMGATE_DEBUG`
  environment variable (set to `"true"`, `"yes"`, `"y"`, or `"1"` to enable).
  Do **NOT** add a debug flag or parameter to function signatures.
- Debug mode also enables structured timing profiles under `pathProject/profile/`.
  The whole-run timer and directory lifecycle are managed from `gateStim()` via
  an `on.exit()` handler. Each completed timing unit is saved immediately as a
  small RDS record under `profile/raw/`; successful and failed runs both collate
  all readable records into `profile/profile.rds` and `profile/profile.csv`
  (recording `"completed"` or `"failed"` status). A new debug `gateStim()` run
  removes previous `debug/` and `profile/` directories first. Profiling and debug
  logging failures must never fail the gating run.
- Keep profiling selective. Record broad workflow stages generally, but detailed
  per-sample timings are currently reserved for initial local-FDR gating. The
  selected sample-level detail is the probability model, density/bandwidth work,
  antimode work, filtering/threshold work and marginal filtering. Do not add a
  timer for every debug statement or small helper without evidence that it is
  useful for diagnosing run time.
- `zz_profile_instrumentation.R` deliberately loads after the implementation files
  and wraps selected internal functions without changing their arguments. Preserve
  wrapped function signatures when profiling boundaries change.
- Use the `stage` parameter to track algorithm stages (`"init"`, `"cytPos"`, or
  `"single"`). Pass `stage` through function calls to enable intermediate data
  saving via `.intSave()` or `.intSaveNm()` functions. Intermediate saving is
  controlled by the `STIMGATE_INTERMEDIATE` environment variable.

### Saved expression and stimulation gates

Saved expression includes unstimulated samples, while final stimulation gate
tables omit them. Positivity helpers must treat channels with no gate for the
current sample as all-FALSE, preserving one logical value per cell.
Completed `chnlSettings.rds` settings are keyed by marker labels, although saved
expression columns use channel names. Resolve that mapping before applying the
saved `biasUns`; channels without a saved bias use zero.

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
  `R/stimgate-package.R`, so `ggplot2` functions and `flowCore::exprs` may be called without
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
4. **F-beta comparator provenance**:
   `scripts/python/fbeta.py` is adapted from the Richards et al. (2014) positivity
   threshold implementation. Preserve its F-beta scoring, standard parameters,
   automatic bin count and moving-average smoothing. The deliberate StimGate-side
   adaptation is that common histogram edges span both the stimulated and
   unstimulated distributions; the published code derives them from the negative
   distribution alone. Document any further deviation explicitly in the relevant
   analyses and comparison issue. Reticulate-backed Python environments must be
   created inside the R process that uses them and kept local to that run; never
   store them in a global R cache that can be serialised to multisession workers.
   For the ACS analysis, requesting a Tailgate/F-beta run removes both prior
   comparator `result.rds` files and recomputes them. Existing results are read only
   when comparator execution is disabled.
5. **Temporary migration status for `.getCpTg()` (issues #157/#158)**:
   This is a current-state note rather than a permanent design rule. Verify it against
   the current implementation and relevant issues before relying on it in later work.
   At the time of this update, remaining call sites are catalogued by
   `.get_cp_tg_call_audit()` and summarised by `.get_cp_tg_migration_note_157()`.
   Current default behaviour still constructs `tgClust` control gates in
   `.gateBatchAll()`, but the current local-FDR cluster quantile implementation does
   not consume `gateTblCtrl`, so this branch is dead plumbing for current outputs.
   Single-positive gating branches have been removed per issue #196.
6. **Simulation engine migration to `simcyto` (issues #288/#289/#291/#295 / umbrella #271)**:
   Generic cytometry simulations, post-simulation transformations, and condition-mismatch
   controls are progressively migrating to the exported `simcyto` package API (e.g.
   `simcyto::simCytExperiment()`, `simcyto::simCytTransform*()`). `analysis/2-sim-bw-freq_bs-global.qmd`,
   `analysis/3-sim-bw-est-base.qmd`, `analysis/7-sim-compare-freq_bs.qmd`, and
   `analysis/8-sim-compare-freq_bs-batch.qmd` use `simcyto` and do not source
   `functionsForBenchmarking-Cyt.R`. StimGate scientific scenario calculations, downstream
   comparison orchestration, and method evaluations remain StimGate-side under `scripts/r/`.
7. **Standardised simulation and plotting controls across analysis QMDs (issue #299)**:
   All analysis QMDs follow a unified execution control pattern sourced from `scripts/r/analysis-runtime.R`:
   - YAML headers declare `params: run_simulations: true, run_plots: false` (along with any chunking parameters).
   - Setup chunks initialise `run_simulations` (defaulting to `FALSE` in interactive execution) and `run_plots` (defaulting to `TRUE` in interactive execution) via `.as_flag(.get_qmd_param_env(...))`.
   - Environment variables `RUN_SIMULATIONS` and `RUN_PLOTS` override YAML parameters and interactive default values.
   - Expensive simulation chunks are guarded with `if (isTRUE(run_simulations))`, and plotting chunks are guarded with `if (isTRUE(run_plots))`.
   - Collation chunks read cached output RDS files unconditionally so downstream summaries and diagnostics work whether simulations just ran or were loaded from cache.
8. **Run-scoped staging, progress and promotion for expensive analysis simulations (issue #304)**:
   Expensive simulation analyses that support resumable per-scenario/per-chunk outputs must use shared run-management helpers from `scripts/r/analysis-runtime.R`:
   - Treat each logical run as a unique run ID (`analysis_run_id` QMD param or `ANALYSIS_RUN_ID` env var; auto-generated when absent).
   - Write run outputs to `cache/sim/<analysis-key>/staging/<YYYY-MM-DD>/<run-id>/...` and keep canonical outputs in `cache/sim/<analysis-key>/current/`.
   - Write run progress/state to `cache/log/analysis/<analysis-key>/<YYYY-MM-DD>/<run-id>/` (`progress.txt`, `manifest.rds`, `status.rds`, chunk/job subdirs).
   - For external chunking, all chunks of one logical run must use the same run ID and write under the same staged run directory, separated by chunk labels.
   - Slurm chunk launchers render the current top-level QMD and pass chunk controls through environment variables. Do not create or launch physical split-QMD copies; all chunks of one submission must receive the same `ANALYSIS_RUN_ID`.
   - Never promote on partial/incomplete runs. Promote only after required chunks are complete and collated outputs validate.
   - Promotion updates `current/` only after a complete staged run is available; failed/interrupted staged runs remain inspectable and resumable.
   - Read canonical outputs through `.analysis_current_file()`, which requires a `COMPLETE` marker, a readable manifest for the requested analysis key, and any analysis-specific semantic version required by the caller.
   - To read canonical results without running the simulation chunk (so no `run_ctx` exists), collation chunks fall back to `.analysis_results_context()`, a read-only stand-in whose staging paths point at `current/` and which creates no run state. Guard all writes, chunk marking and promotion with `if (!isTRUE(run_ctx$read_only))`.
   - Record scientific and semantic settings in the run manifest. Reusing an explicit run ID must match those settings; only operational controls such as plotting, simulation execution and the current chunk index may differ across invocations.

9. **Versioning before the first Bioconductor release**:
   Keep `Version` in `DESCRIPTION` at `0.99.z` (three components, no `-n`
   suffix) until stimgate's first Bioconductor release, bumping `z` for each
   change worth marking. Do not move to `0.100.0` or higher; Bioconductor sets
   the release version itself.

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
   `/` or `\\` separators) and avoid platform-specific paths. Pull-request CI runs on
   Windows, so in particular:
   - Pass only a file-name prefix to `tempfile()`; a full path as the pattern is
     prepended with `tempdir()` again, which is invalid on Windows.
   - Compare paths after `normalizePath(path, winslash = "/", mustWork = FALSE)`,
     since equivalent paths may differ in separator style.
   - Embed only forward-slash paths in R code run through `Rscript -e`;
     Windows backslashes are escape sequences there.
   - Use `skip_on_os("windows")`, with a comment giving the reason, for checks of
     Unix-only process or signal behaviour.
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

<!-- github-projects:start -->
## 9. Ponytail

For coding, refactoring, bug-fixing, review and implementation design, read
`.agents/skills/ponytail/SKILL.md` and apply Ponytail in **full** mode by
default. Do not use **ultra** unless the operator explicitly requests it.
Provenance is in `.agents/skills/ponytail/README.md`.

Ponytail is subordinate to settled behaviour. Precedence is:

1. the current explicit operator instruction or issue acceptance criteria;
2. the conventions and policies in this file;
3. the exported API and supported analysis/QMD contracts;
4. the smallest implementation.

Do not simplify away input validation, `.debug()`/profiling/intermediate-save
plumbing, run-scoped staging and promotion for analysis simulations, comparator
provenance, or roxygen documentation of exported functions. The generic
Ponytail "one runnable check" suggestion is not a test cap here: add focused
`testthat` coverage in the appropriate suite (Section 8) and run the
pre-commit checklist (Section 4).

## GitHub issues and Projects

For GitHub issue or Project administration, use
`.agents/skills/github-projects/SKILL.md` and read
`.projects/project.md` before acting.
<!-- github-projects:end -->
