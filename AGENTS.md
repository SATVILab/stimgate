# AGENTS.md — Configuration for AI Coding Agents

This file provides guidance for autonomous coding agents (e.g. Google
Jules, GitHub Copilot) working on the `stimgate` repository.

------------------------------------------------------------------------

## 1. Core Philosophy / Project Context

`stimgate` is an R package intended for submission to
[Bioconductor](https://bioconductor.org/). It identifies cells that have
possibly responded to immune stimulation by applying outlier-based
gating to flow cytometry data. The core idea is to compare an
*unstimulated* tube with a *stimulated* tube from the same donor sample
and flag cells whose marker expression in the stimulated condition is
unusually high relative to the unstimulated background.

**Key architectural patterns:**

- The main user entry point is
  [`stimgate_gate()`](https://satvilab.github.io/stimgate/reference/stimgate_gate.md),
  which writes intermediate results to a `path_project` directory on
  disk and returns that path.
- Downstream helpers
  ([`get_stats()`](https://satvilab.github.io/stimgate/reference/get_stats.md),
  `get_gate_tbl()`,
  [`stimgate_plot()`](https://satvilab.github.io/stimgate/reference/stimgate_plot.md),
  [`stimgate_fcs_write()`](https://satvilab.github.io/stimgate/reference/stimgate_fcs_write.md))
  all accept `path_project` and read from that directory.
- Internal functions are prefixed with `.` and are not exported.
- The package integrates tightly with the `flowCore` / `flowWorkspace`
  Bioconductor ecosystem; input data are `GatingSet` objects.
- `renv` is used for reproducible dependency management with two
  profiles: `bioc_container` (Bioconductor Docker) and
  `non_bioc_container` (standard R).

------------------------------------------------------------------------

## 2. Tech Stack & Tooling

| Layer                 | Tool / Package                                          |
|-----------------------|---------------------------------------------------------|
| Language              | R (≥ 4.4.0)                                             |
| Package framework     | `devtools`, `roxygen2`, `testthat` (≥ 3.0.0)            |
| Documentation         | `roxygen2`, `pkgdown`                                   |
| Flow cytometry I/O    | `flowCore`, `flowWorkspace`                             |
| Data manipulation     | `dplyr`, `purrr`, `tidyr`, `tibble`, `stringr`, `rlang` |
| Plotting              | `ggplot2`, `cowplot`                                    |
| Statistical modelling | `scam`, `mgcv`                                          |
| Clustering            | `cluster`, `gtools`                                     |
| Dependency management | `renv`                                                  |
| CI                    | GitHub Actions (R-CMD-check, pkgdown, Codecov)          |

**Do NOT:**

- Use `@import` or `@importFrom` directives in roxygen comments;
  explicitly qualify all package calls with `pkg::fun()`. The only
  exceptions are `ggplot2` functions and
  [`flowCore::exprs`](https://rdrr.io/pkg/Biobase/man/exprs.html), which
  may be called without a namespace qualifier and do not require
  `@importFrom` tags.
- Modify `.Rd` files manually; regenerate them with
  `devtools::document()`.
- Use [`return()`](https://rdrr.io/r/base/function.html) as the last
  line of a function; use it only for early returns.
- Add [`library()`](https://rdrr.io/r/base/library.html) calls inside
  test files.

------------------------------------------------------------------------

## 3. Setup Commands

All commands below must be run from the **repository root** (the
directory containing `DESCRIPTION`) so that `.Rprofile` is sourced and
`renv` activates correctly.

``` r
# 1. Install renv (if not already installed)
install.packages("renv")

# 2. Restore the package library from renv.lock
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
BiocManager::install(c("flowCore", "flowWorkspace", "HDCytoData"))
```

------------------------------------------------------------------------

## 4. Build & Test Instructions

Run the following from the **repository root** in an R session:

``` r
# Run all unit tests
devtools::test()

# Full R CMD check (mimics CRAN / Bioconductor checks)
devtools::check()

# Regenerate documentation from roxygen comments
devtools::document()

# Check test coverage
covr::report()
```

Equivalent shell commands (via `Rscript`):

``` bash
Rscript -e "devtools::test()"
Rscript -e "devtools::check()"
Rscript -e "devtools::document()"
```

------------------------------------------------------------------------

## 5. Coding Style & Conventions

### Formatting

- Run `styler::style_pkg()` before every commit to apply consistent
  formatting.
- Run `lintr::lint_package()` to catch style violations before opening a
  PR.

``` r
styler::style_pkg()
lintr::lint_package()
```

### Naming conventions

- **Exported functions**: `snake_case`, no leading dot
  (e.g. `stimgate_gate`).
- **Internal functions**: `snake_case` with a leading dot
  (e.g. `.get_threshold`).
- **Debug parameter**: always named `.debug` in internal functions.
- **Debug messages**: use `.debug_msg(flag, message, optional_value)`.

### Documentation (roxygen2)

- Every exported function must have a `@param`, `@return`, and `@export`
  tag.
- Parameter docs follow the format:
  `@param param_name <type> <description>. Default: <value>.`
- Use `@details` for complex explanations rather than overloading
  `@param`.

### Package namespace

- Reference all external functions as `pkg::fun()`.
- Exceptions: `ggplot2` functions and
  [`flowCore::exprs`](https://rdrr.io/pkg/Biobase/man/exprs.html) may be
  called without the namespace qualifier and do not require
  `@importFrom` tags.

### Tests

- Place tests in `tests/testthat/` as `test-<topic>.R` files.
- Tests use `testthat` 3rd-edition conventions (`expect_*` helpers).
- Do **not** place any code outside `test_that()` blocks except shared
  setup data that every test in the file requires.
- Each test cleans up its own temporary files
  (e.g. `unlink(tmp, recursive = TRUE)`).
- Tests must pass on macOS, Windows, and Ubuntu; use
  [`file.path()`](https://rdrr.io/r/base/file.path.html) (not hard-coded
  `/` or `\` separators) and avoid platform-specific paths.
- Test observable behaviour and outputs, not internal implementation
  details or the existence of internal (`.`-prefixed) functions.

### Commits & PRs

- Keep commits focused and atomic.
- Commit message format: `<verb> <what>` in the imperative mood
  (e.g. `Add threshold helper for cluster method`).
- Before opening a PR:
  1.  `styler::style_pkg()`
  2.  `lintr::lint_package()`
  3.  `devtools::document()`
  4.  `devtools::test()`
- Update `_pkgdown.yml` whenever exported functions are added, removed,
  or renamed, and verify the site config with
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html).
