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

# Install the latest stable Quarto release. Avoid pinning the cloud environment
# to an obsolete Quarto release while still using a deterministic release asset.
QUARTO_VERSION="$(curl -fsSL https://api.github.com/repos/quarto-dev/quarto-cli/releases/latest \
    | python3 -c 'import json,sys; print(json.load(sys.stdin)["tag_name"].lstrip("v"))')"
wget -qO /tmp/quarto.deb \
    "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb"
sudo dpkg -i /tmp/quarto.deb
rm -f /tmp/quarto.deb
quarto --version

run_pak_pass() {
    local rc=0

    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::local_install_dev_deps(".", upgrade = FALSE, ask = FALSE)
    ' || rc=1

    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    pak::pkg_install(
        c(
            "devtools", "rcmdcheck", "decor", "roxygen2", "styler", "lintr", "covr", "pkgdown",
            "SATVILab/projr",
            "SATVILab/simcyto",
            "future", "furrr", "progressr",
            "reticulate",
            "RGLab/cytoUtils"
        ),
        upgrade = FALSE,
        ask = FALSE
    )
    ' || rc=1

    return "$rc"
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

if run_pak_pass; then
    PAK_PASS_1=0
else
    PAK_PASS_1=$?
fi

if [ "$PAK_PASS_1" -ne 0 ]; then
    run_bioc_repair || true
fi

if run_pak_pass; then
    PAK_PASS_2=0
else
    PAK_PASS_2=$?
fi

if [ "$PAK_PASS_2" -ne 0 ]; then
    sudo env CI=true PKG_SYSREQS=true PKG_SYSREQS_VERBOSE=true \
        Rscript --vanilla -e '
    source("/tmp/stimgate-agent-repos.R")
    try(pak::sysreqs_fix_installed(), silent = FALSE)
    ' || true
    run_bioc_repair || true
    run_pak_pass || true
fi

Rscript --vanilla -e '
required <- c(
    "flowCore", "flowWorkspace", "Biobase",
    "testthat", "devtools", "rcmdcheck", "decor", "roxygen2", "styler", "lintr", "covr", "pkgdown",
    "projr", "simcyto", "future", "furrr", "progressr", "reticulate", "cytoUtils"
)
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Codex R environment is incomplete. Missing: ", paste(required[!ok], collapse = ", "))
stopifnot(identical(Sys.getenv("CI"), "true"))
devtools::load_all()
cat("StimGate Codex environment ready\n")
cat("R: ", R.version.string, "\n", sep = "")
cat("Bioconductor: ", as.character(BiocManager::version()), "\n", sep = "")
sim_desc <- utils::packageDescription("simcyto")
cat("simcyto: ", as.character(sim_desc[["Version"]]), "\n", sep = "")
if (!is.null(sim_desc[["RemoteSha"]])) cat("simcyto GitHub SHA: ", sim_desc[["RemoteSha"]], "\n", sep = "")
'

quarto check
