set -euxo pipefail

# Jules currently checks the repository out under /app.
cd /app

export CI=true

export DEBIAN_FRONTEND=noninteractive
export PKG_SYSREQS=true
export PKG_SYSREQS_VERBOSE=true

R_MINOR="4.6"
BIOC_VERSION="3.23"

# Use Posit Package Manager's pak-compatible Linux binary URL. Do not replace
# this with the newer /bin/linux/noble-x86_64/<R> form: pak currently
# classifies that repository as source and rebuilds ordinary CRAN packages.
POSIT_CRAN="https://packagemanager.posit.co/cran/__linux__/noble/latest"

sudo apt-get update -qq
sudo apt-get install --no-install-recommends -y \
    ca-certificates curl dirmngr gnupg lsb-release software-properties-common wget git \
    libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev libgit2-dev \
    libjpeg-dev libpng-dev libx11-dev libssl-dev libxml2-dev pandoc \
    python3 python3-numpy

wget -qO /tmp/cran_ubuntu_key.asc \
    https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc
sudo install -m 0644 /tmp/cran_ubuntu_key.asc /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository -y \
    "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update -qq
sudo apt-get install --no-install-recommends -y r-base r-base-dev

Rscript --vanilla -e "
expected <- '${R_MINOR}'
actual <- paste(R.version\$major, sub('\\\\..*$', '', R.version\$minor), sep = '.')
if (!identical(actual, expected)) stop('Expected R ', expected, '.x but found ', R.version.string)
cat('Using ', R.version.string, '\n', sep = '')
"

sudo env CI=true Rscript --vanilla -e "
options(repos = c(CRAN = '${POSIT_CRAN}'))
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)
cat('Bioconductor version: ', as.character(BiocManager::version()), '\n', sep = '')
"

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

cat > /tmp/stimgate-agent-repos.R <<EOF
options(
    repos = {
        repos <- BiocManager::repositories(version = "$BIOC_VERSION")
        repos <- repos[!names(repos) %in% c("CRAN", "RSPM")]
        c(repos, CRAN = "$POSIT_CRAN")
    }
)
EOF

sudo env CI=true Rscript --vanilla -e '
source("/tmp/stimgate-agent-repos.R")
cat("pak repository status:\n")
print(pak::repo_status())
'

# Install the latest stable Quarto release.
QUARTO_VERSION="$(curl -fsSL https://api.github.com/repos/quarto-dev/quarto-cli/releases/latest \
    | python3 -c 'import json,sys; print(json.load(sys.stdin)["tag_name"].lstrip("v"))')"
wget -qO /tmp/quarto.deb \
    "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb"
sudo dpkg -i /tmp/quarto.deb
rm -f /tmp/quarto.deb
quarto --version

run_pak_bioc_core() {
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::pkg_install(
        c("flowCore", "flowWorkspace", "Biobase"),
        upgrade = FALSE,
        ask = FALSE
    )
    '
}

run_bioc_repair() {
    sudo env CI=true Rscript --vanilla -e "
    source('/tmp/stimgate-agent-repos.R')
    BiocManager::install(
        c('flowCore', 'flowWorkspace', 'Biobase'),
        version = '${BIOC_VERSION}',
        ask = FALSE,
        update = FALSE
    )
    "
}

run_pak_local_deps() {
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::local_install_dev_deps(".", upgrade = FALSE, ask = FALSE)
    '
}

run_pak_dev_tools() {
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::pkg_install(
        c("devtools", "rcmdcheck", "decor", "roxygen2", "styler", "lintr", "covr", "pkgdown"),
        upgrade = FALSE,
        ask = FALSE
    )
    '
}

run_pak_analysis_deps() {
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::pkg_install(
        c("future", "furrr", "progressr", "reticulate"),
        upgrade = FALSE,
        ask = FALSE
    )
    '
}

run_pak_github_deps() {
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::pkg_install(
        c("SATVILab/projr", "SATVILab/simcyto", "RGLab/cytoUtils"),
        upgrade = FALSE,
        ask = FALSE
    )
    '
}

run_stage() {
    local label="$1"
    local fn="$2"

    echo
    echo "################################################################"
    echo "# ${label}"
    echo "################################################################"

    if "$fn"; then
        return 0
    fi

    echo "${label} failed."
    return 1
}

# Keep the initial setup in bounded transactions. In particular, install the
# Bioconductor runtime stack before asking pak to resolve all local dev deps.
# This avoids one large mixed CRAN/Bioconductor transaction and makes a stalled
# download or source build much easier to diagnose.
BIOC_RC=0
LOCAL_RC=0
TOOLS_RC=0
ANALYSIS_RC=0
GITHUB_RC=0

run_stage "STAGE 1: core Bioconductor packages" run_pak_bioc_core || BIOC_RC=$?

if [ "$BIOC_RC" -ne 0 ]; then
    echo "Trying BiocManager repair before retrying the core Bioconductor stage."
    run_bioc_repair || true
    BIOC_RC=0
    run_stage "STAGE 1 RETRY: core Bioconductor packages" run_pak_bioc_core || BIOC_RC=$?
fi

run_stage "STAGE 2: StimGate DESCRIPTION development dependencies" run_pak_local_deps || LOCAL_RC=$?
run_stage "STAGE 3: package-development tools" run_pak_dev_tools || TOOLS_RC=$?
run_stage "STAGE 4: analysis CRAN dependencies" run_pak_analysis_deps || ANALYSIS_RC=$?
run_stage "STAGE 5: GitHub dependencies" run_pak_github_deps || GITHUB_RC=$?

if [ "$LOCAL_RC" -ne 0 ] || [ "$TOOLS_RC" -ne 0 ] || [ "$ANALYSIS_RC" -ne 0 ] || [ "$GITHUB_RC" -ne 0 ]; then
    echo
    echo "One or more pak stages failed. Repairing installed system requirements and retrying only failed stages."

    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    try(pak::sysreqs_fix_installed(), silent = FALSE)
    ' || true

    if [ "$LOCAL_RC" -ne 0 ]; then
        LOCAL_RC=0
        run_stage "STAGE 2 RETRY: StimGate DESCRIPTION development dependencies" run_pak_local_deps || LOCAL_RC=$?
    fi
    if [ "$TOOLS_RC" -ne 0 ]; then
        TOOLS_RC=0
        run_stage "STAGE 3 RETRY: package-development tools" run_pak_dev_tools || TOOLS_RC=$?
    fi
    if [ "$ANALYSIS_RC" -ne 0 ]; then
        ANALYSIS_RC=0
        run_stage "STAGE 4 RETRY: analysis CRAN dependencies" run_pak_analysis_deps || ANALYSIS_RC=$?
    fi
    if [ "$GITHUB_RC" -ne 0 ]; then
        GITHUB_RC=0
        run_stage "STAGE 5 RETRY: GitHub dependencies" run_pak_github_deps || GITHUB_RC=$?
    fi
fi

Rscript --vanilla -e '
required <- c(
    "flowCore", "flowWorkspace", "Biobase",
    "testthat", "devtools", "rcmdcheck", "decor", "roxygen2", "styler", "lintr", "covr", "pkgdown",
    "projr", "simcyto", "future", "furrr", "progressr", "reticulate", "cytoUtils"
)
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Jules R environment is incomplete. Missing: ", paste(required[!ok], collapse = ", "))
devtools::load_all()
cat("StimGate Jules environment ready\n")
cat("R: ", R.version.string, "\n", sep = "")
cat("Bioconductor: ", as.character(BiocManager::version()), "\n", sep = "")
sim_desc <- utils::packageDescription("simcyto")
cat("simcyto: ", as.character(sim_desc[["Version"]]), "\n", sep = "")
if (!is.null(sim_desc[["RemoteSha"]])) cat("simcyto GitHub SHA: ", sim_desc[["RemoteSha"]], "\n", sep = "")
'

quarto check
