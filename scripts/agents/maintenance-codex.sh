set -euxo pipefail

cd "$(git rev-parse --show-toplevel)"

export DEBIAN_FRONTEND=noninteractive
export PKG_SYSREQS=true
export PKG_SYSREQS_VERBOSE=true

if [ "${CI:-}" != "true" ]; then
    echo "ERROR: Configure CI=true as a persistent Codex environment variable."
    exit 1
fi

R_MINOR="4.6"
BIOC_VERSION="3.23"
POSIT_CRAN="https://packagemanager.posit.co/cran/latest/bin/linux/noble-x86_64/${R_MINOR}"

cat > /tmp/stimgate-agent-repos.R <<EOF
options(
    repos = {
        repos <- BiocManager::repositories(version = "$BIOC_VERSION")
        repos <- repos[!names(repos) %in% c("CRAN", "RSPM")]
        c(repos, CRAN = "$POSIT_CRAN")
    }
)
EOF

if ! Rscript --vanilla -e 'quit(status = if (requireNamespace("pak", quietly = TRUE)) 0 else 1)'; then
    sudo env CI=true Rscript --vanilla -e '
    install.packages(
        "pak",
        repos = sprintf(
            "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
            .Platform$pkgType,
            R.Version()$os,
            R.Version()$arch
        )
    )
    '
fi

# Pick up any newly declared StimGate dependencies without broadly upgrading
# the whole cached environment.
sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
    Rscript --vanilla -e '
source("/tmp/stimgate-agent-repos.R")
pak::local_install_dev_deps(".", upgrade = FALSE, ask = FALSE)
'

# simcyto deliberately tracks the current SATVILab default branch. It is not
# pinned because StimGate development normally targets the latest simcyto API.
sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
    Rscript --vanilla -e '
source("/tmp/stimgate-agent-repos.R")
pak::pkg_install("SATVILab/simcyto", upgrade = TRUE, ask = FALSE)
'

# Ensure repository-wide tools that are intentionally outside DESCRIPTION are
# still available. Do not force upgrades for these on every maintenance run.
sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
    Rscript --vanilla -e '
source("/tmp/stimgate-agent-repos.R")
pak::pkg_install(
    c(
        "devtools", "rcmdcheck", "decor", "roxygen2", "styler", "lintr", "covr", "pkgdown",
        "SATVILab/projr",
        "future", "furrr", "progressr",
        "reticulate",
        "RGLab/cytoUtils"
    ),
    upgrade = FALSE,
    ask = FALSE
)
'

Rscript --vanilla -e '
required <- c(
    "flowCore", "flowWorkspace", "Biobase", "testthat", "devtools",
    "projr", "simcyto", "future", "furrr", "progressr", "reticulate", "cytoUtils"
)
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Codex maintenance left packages missing: ", paste(required[!ok], collapse = ", "))
stopifnot(identical(Sys.getenv("CI"), "true"))
devtools::load_all()
sim_desc <- utils::packageDescription("simcyto")
cat("StimGate Codex maintenance complete\n")
cat("simcyto: ", as.character(sim_desc[["Version"]]), "\n", sep = "")
if (!is.null(sim_desc[["RemoteSha"]])) cat("simcyto GitHub SHA: ", sim_desc[["RemoteSha"]], "\n", sep = "")
cat("stimgate loaded successfully\n")
'
