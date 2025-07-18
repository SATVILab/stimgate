.set_renv_profile <- function() {
  # set the renv profile based on the environment
  # if using the Bioconductor Docker image,
  # then use the bioc_container profile.
  # Otherwise use the non_bioc_container profile.
  # To use the non_bioc_container profile,
  # we have to override the system-wide
  # .Renviron settings used in the devcontainer feature config-r
  # to enable global cache accessible after image build.
  if (.is_bioc_container()) {
    Sys.setenv("RENV_PROFILE" = "bioc_container")
  } else {
    Sys.setenv("BIOCONDUCTOR_USE_CONTAINER_REPOSITORY" = "false")
    Sys.setenv("RENV_PROFILE" = "non_bioc_container")
    if (.is_devcontainer_config_r_feature()) {
      Sys.unsetenv("RENV_PATHS_CACHE")
      Sys.unsetenv("RENV_PATHS_LIBRARY_ROOT")
      Sys.unsetenv("RENV_PATHS_ROOT")
    }
  }
}

.is_bioc_container <- function() {
  Sys.getenv("BIOCONDUCTOR_NAME") == "bioconductor_docker"
}

.is_devcontainer_config_r_feature <- function() {
  # undo system-wide .Renviron settings
  # used in the devcontainer feature config-r
  # to enable global cache accessible after
  # image build.
  # just use renv defaults now.
  renv_paths_cache_match <- Sys.getenv("RENV_PATHS_CACHE") == "/renv/cache"
  renv_paths_library_root_match <- Sys.getenv("RENV_PATHS_LIBRARY_ROOT") ==
    "/workspaces/.local/lib/R/library"
  renv_paths_root_match <- Sys.getenv("RENV_PATHS_ROOT") == "/renv/local"
  renv_paths_cache_match &&
    renv_paths_library_root_match &&
    renv_paths_root_match
}

.is_ci <- function() {
  is_gha <- Sys.getenv("GITHUB_ACTIONS") == "true"
  is_ci <- Sys.getenv("CI") == "true"
  is_gha || is_ci
}


if (!.is_ci()) {
  .set_renv_profile()
  source("renv/activate.R")
} else {
  try({
    bioc_vec <- if (requireNamespace("BiocManager", quietly = TRUE)) {
      suppressWarnings(BiocManager::repositories())
    } else NULL
    repos_vec <- c(bioc_vec, getOption("repos"))
    if (isFALSE("CRAN" %in% names(repos_vec))) {
      repos_vec <- c(repos_vec, CRAN = "https://cloud.r-project.org/")
    }
    options(repos = repos_vec)
  }, silent = TRUE)
}

if (length(getOption("repos")) == 0L){
  options(repos = c("CRAN" = "https://cloud.r-project.org/"))
}

if (is.null(names(getOption("repos")))) {
  options(
    repos = c(getOption("repos"), "CRAN" = "https://cloud.r-project.org/")
  )
}

try(suppressWarnings(rm(
  .set_renv_profile, .is_ci,
  .is_bioc_container, .is_devcontainer_config_r_feature
)), silent = TRUE)
