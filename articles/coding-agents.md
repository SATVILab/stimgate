# Setting up coding agents for StimGate

## Overview

StimGate can be developed with remote coding agents as long as the agent
environment contains the same core R, Bioconductor, system and analysis
dependencies used by the repository.

This guide currently covers:

- Google Jules
- OpenAI Codex Cloud
- GitHub Copilot coding agent

The repository-level instructions for coding agents live in `AGENTS.md`.
Agents should read that file before changing code. The environment
scripts described below are kept under `scripts/agents/` where
appropriate so that the setup is version-controlled alongside the
package rather than copied independently between agent products.

The recommended cloud-agent environment is deliberately different from
ordinary local development. In local development, StimGate may activate
an `renv` profile through `.Rprofile`. Cloud agents instead use a
directly installed development library and should run with `CI=true`,
which tells the repository not to activate `renv`.

The current setup targets:

- R 4.6
- Bioconductor 3.23
- the current `SATVILab/simcyto` default branch
- Quarto
- the StimGate `DESCRIPTION` dependencies and Suggests
- development tools such as `devtools`, `rcmdcheck`, `roxygen2`,
  `styler`, `lintr`, `covr` and `pkgdown`
- analysis dependencies used by the simulation and comparison workflows,
  including `future`, `furrr`, `progressr`, `reticulate` and `cytoUtils`

`simcyto` is intentionally not pinned to a commit. StimGate development
normally targets the latest version in `SATVILab/simcyto`. Codex
maintenance explicitly refreshes it, and the Copilot and analysis GitHub
Actions workflows install the unpinned repository package.

## Shared principles

### Prefer repository scripts or workflows over duplicated UI scripts

The recommended approach is to make the agent platform use configuration
stored in the repository whenever possible. This keeps the authoritative
setup in Git, makes changes reviewable, and avoids having different
Jules, Codex and documentation copies drift apart.

For Codex and Jules, repository setup scripts live under
`scripts/agents/`. GitHub Copilot uses the repository workflow
`.github/workflows/copilot-setup-steps.yml` directly.

If an agent product only accepts pasted setup text, copy the current
contents of the corresponding repository file rather than maintaining a
separate unmanaged version.

### Do not restore `renv` in the cloud-agent environment

The cloud environment installs the packages it needs directly and then
runs commands such as:

``` r

#| eval: false
devtools::load_all()
devtools::test()
```

Do not run `renv::restore()` merely because the project contains `renv`
files. For these cloud environments, a missing package should normally
be treated as an environment-setup problem.

### Keep `simcyto` current

Use:

``` text
SATVILab/simcyto
```

rather than a commit-pinned package specification unless a particular
historical analysis explicitly requires an old API.

This rule also applies to CI and coding-agent setup. A stale `simcyto`
pin can make the agent test StimGate against an API that is no longer
the one used by current development.

### Validate the environment, not just the installer exit code

The setup environments finish by checking that the important packages
are actually available and that:

``` r

#| eval: false
devtools::load_all()
```

succeeds. Intermediate package-install failures may be recoverable,
especially for Bioconductor dependencies, so the Codex and Jules scripts
use a pak-first, Bioconductor-repair, pak-again strategy before applying
the final hard checks.

## OpenAI Codex Cloud

### Recommended Codex settings

Create a Codex Cloud environment for the StimGate repository and use the
universal container image.

The preinstalled language versions in the universal image can normally
be left at their defaults. Python is the only preinstalled language that
is directly relevant to current StimGate analysis code; the comparison
workflow calls the historical F-beta Python implementation through
`reticulate`. The setup script separately installs NumPy and the R-side
Python integration requirements.

Add the following persistent environment variable:

``` text
CI=true
```

This setting is important. The setup script itself runs in a separate
shell from later agent commands, so exporting `CI=true` only inside the
setup script is not sufficient. StimGate’s `.Rprofile` uses `CI=true` to
avoid activating `renv` in the cloud-agent session.

Enable container caching. The initial R/Bioconductor installation is
substantial, and caching avoids rebuilding the environment for every
task.

### Setup script

Choose **Manual** setup and use:

``` bash
bash scripts/agents/setup-codex.sh
```

The full script is stored at:

``` text
scripts/agents/setup-codex.sh
```

It installs R and Bioconductor, the system libraries required by the
package, Quarto, StimGate’s development dependencies, the latest
`SATVILab/simcyto`, and the additional analysis tooling used by the
repository.

### Maintenance script

Use:

``` bash
bash scripts/agents/maintenance-codex.sh
```

The maintenance script is intentionally lighter than the initial setup.
It:

1.  installs any newly declared StimGate development dependencies;
2.  explicitly refreshes `SATVILab/simcyto` from its current default
    branch;
3.  ensures the repository-wide tools that are not declared in
    `DESCRIPTION` remain available; and
4.  checks that StimGate still loads.

The full script is stored at:

``` text
scripts/agents/maintenance-codex.sh
```

### Internet access

Agent internet access is useful for StimGate because an agent may need
to inspect current package documentation or fetch a GitHub dependency
after the cached environment was created.

A practical configuration is:

- Agent internet access: **On**
- Domain allowlist: **Common dependencies**, plus the domains needed for
  GitHub, R, Posit Package Manager and Bioconductor
- HTTP methods: **All methods** if that is the only option that permits
  Git-over-HTTPS operations

Useful additional domains include:

``` text
github.com
api.github.com
codeload.github.com
raw.githubusercontent.com
objects.githubusercontent.com
r-lib.github.io
cloud.r-project.org
cran.r-project.org
r-project.org
packagemanager.posit.co
p3m.dev
bioconductor.org
```

Keep the allowlist narrower than unrestricted internet access. The setup
script itself should contain the normal software-installation endpoints,
while task-time network access should be limited to domains the agent is
likely to need.

### After setup

A normal coding task can use the current checkout directly. The standard
package test sequence is documented in `AGENTS.md`, but common commands
are:

``` bash
Rscript -e "devtools::load_all(); devtools::test()"
Rscript analysis/tests/run_analysis_tests.R
```

Run the analysis integration suite when changes touch `analysis/` or
`scripts/r/`, rather than for every unrelated package change.

## Google Jules

### Setup

For Jules, use the repository setup script:

``` bash
bash scripts/agents/setup-jules.sh
```

The full script is stored at:

``` text
scripts/agents/setup-jules.sh
```

The Jules version follows the same dependency strategy as the Codex
setup, but retains the Jules checkout convention used by the working
StimGate environment.

It sets `CI=true` during setup so that R processes launched by the
installer do not activate the project `renv` environment. If Jules
provides a persistent environment-variable setting for a project, also
set `CI=true` there so later agent R sessions use the same direct
development library.

### Network and repository access

Jules needs access to the StimGate repository and to the package sources
used during environment setup. If network controls are available, allow
the same R, Bioconductor, Posit and GitHub domains listed for Codex
rather than enabling unrelated destinations.

### Keeping the environment current

When the Jules environment is rebuilt, `setup-jules.sh` installs the
latest `SATVILab/simcyto` rather than a pinned revision.

If a long-lived Jules environment becomes stale, rerun the setup or
explicitly refresh the GitHub package from R:

``` r

#| eval: false
pak::pkg_install("SATVILab/simcyto", upgrade = TRUE, ask = FALSE)
```

## GitHub Copilot coding agent

### How StimGate configures Copilot

The GitHub Copilot coding agent does not use the Codex or Jules setup
scripts. Its environment is configured by the repository workflow:

``` text
.github/workflows/copilot-setup-steps.yml
```

That file is the source of truth for Copilot’s cloud environment. Keep
Copilot environment changes there so that every Copilot task gets the
same setup and the configuration remains reviewable in Git.

The workflow currently:

- installs the Linux system libraries needed by StimGate;
- checks out the repository;
- installs the current R release with Posit Package Manager enabled;
- installs StimGate package and development dependencies;
- installs `devtools`, `rcmdcheck` and `decor`;
- installs the latest `SATVILab/simcyto` from its default branch;
- verifies the important R packages and
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html);
  and
- installs Quarto.

The repository `.Rprofile` treats the Copilot agent as CI, so Copilot
should use the directly installed development library rather than
running `renv::restore()`.

### `simcyto` policy for Copilot

Do not pin `simcyto` to a commit in
`.github/workflows/copilot-setup-steps.yml`.

The expected package specification is:

``` text
SATVILab/simcyto
```

This ensures newly started Copilot tasks are built against the latest
`simcyto` development API.

The analysis integration workflow follows the same rule.
`.github/workflows/analysis-integration.yaml` also installs unpinned
`SATVILab/simcyto`, so Copilot, analysis CI and current StimGate
development test the same moving dependency rather than different
historical revisions.

Pinning a SHA is only appropriate for a deliberately historical
reproduction or version-regression task. If that is required, make the
pin task-specific rather than changing the default Copilot setup.

### Using Copilot on StimGate

Before changing code, the agent should follow `AGENTS.md`, which defines
the repository architecture, coding style, testing rules and which test
suite to run for each kind of change.

Typical package work can validate with:

``` bash
Rscript -e "devtools::load_all(); devtools::test()"
```

Changes involving `analysis/` or `scripts/r/` should also run:

``` bash
Rscript analysis/tests/run_analysis_tests.R
```

The Copilot setup workflow is intentionally closer to the GitHub Actions
environment than to the heavier Codex/Jules scripts. That makes it fast
enough for repeated cloud-agent tasks while still covering the normal
package-development path.

### Updating the Copilot environment

When StimGate gains a new dependency:

1.  add ordinary package dependencies to `DESCRIPTION` where they
    belong;
2.  add agent-only or repository-wide development tools to
    `copilot-setup-steps.yml` only when they are genuinely needed for
    Copilot tasks;
3.  keep `SATVILab/simcyto` unpinned; and
4.  extend the verification step if the dependency is important enough
    that a missing installation should make the environment fail
    immediately.

Avoid duplicating a large independent shell installer for Copilot unless
the GitHub Actions setup becomes unable to express the required
environment cleanly.

## Alternative environment strategies

### Use the GitHub Actions environment as the minimal baseline

The repository’s GitHub Actions workflows are useful references for the
minimum environment needed for package and analysis tests. They are
deliberately more compact than the Codex and Jules cloud-agent scripts.
This is a good option when an agent only needs package tests and not the
full range of research analyses.

### Use `renv` for local or long-lived interactive development

The direct cloud-agent installation is not intended to replace the
project’s `renv` profiles for ordinary local development. `renv` remains
useful when an exact reproducible local library is more important than
fast cloud-agent provisioning.

### Pin dependencies only for reproducibility exercises

The agent environment normally follows the current `simcyto` repository.
Pinning a GitHub SHA is appropriate when reproducing a historical result
or debugging a version-specific regression, but it should be an explicit
task-level decision rather than the default development environment.

### Reduce the environment for narrow tasks

An agent working only on documentation or a small package function may
not need Quarto, the comparison stack, or all analysis dependencies. A
smaller environment can start faster. The full Codex and Jules scripts
are intentionally broad because those agents are expected to handle
package code, tests and current analyses without frequent environment
rebuilds.

## Troubleshooting

### R starts by activating `renv`

Check that `CI=true` is present in the agent’s task-time environment,
not only exported during setup. For Copilot, also check that the
repository `.Rprofile` still recognises the Copilot/CI environment.

### `simcyto` is installed but appears stale

In Codex, rerun the maintenance script:

``` bash
bash scripts/agents/maintenance-codex.sh
```

For another environment, run:

``` r

#| eval: false
pak::pkg_install("SATVILab/simcyto", upgrade = TRUE, ask = FALSE)
```

For Copilot, inspect `.github/workflows/copilot-setup-steps.yml` and
confirm that `SATVILab/simcyto` is unpinned. A newly provisioned Copilot
environment should then install the current default branch.

### A Bioconductor package fails during installation

The Codex and Jules setup scripts deliberately continue through a repair
cycle rather than treating the first pak failure as final. If the final
validation still fails, inspect the missing package and its Linux system
requirements before changing the package code. For Copilot, diagnose the
`setup-r-dependencies` job and system-library step in the setup
workflow.

### StimGate installs but analysis tests fail

Remember that the analysis code has dependencies beyond the installed
package `DESCRIPTION`, including `simcyto` and the comparison tooling.
Use the full Codex/Jules setup when the task involves `analysis/` or
`scripts/r/`. For Copilot, add a missing analysis dependency to its
workflow only when Copilot is expected to execute that analysis path.

## Adding another coding agent

When adding another remote coding agent, prefer this pattern:

1.  add a small agent-specific script under `scripts/agents/`, or use
    the platform’s repository-native workflow if it has one;
2.  reuse the same R/Bioconductor and package strategy unless the
    platform supplies a better native equivalent;
3.  keep `AGENTS.md` as the coding-agent source of truth for repository
    conventions;
4.  add the new platform to this guide; and
5.  avoid copying version-pinned dependency lists into several unrelated
    configuration files unless reproducibility requires it.
