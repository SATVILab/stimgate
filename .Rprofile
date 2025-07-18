.set_renv_profile <- function() {
  if (Sys.getenv("BIOCONDUCTOR_NAME") == "bioconductor_docker") {
    Sys.setenv("RENV_PROFILE" = "bioc_container")
  } else {
    Sys.setenv("BIOCONDUCTOR_USE_CONTAINER_REPOSITORY" = "false")
    Sys.setenv("RENV_PROFILE" = "non_bioc_container")
    # undo system-wide .Renviron settings
    # used in the devcontainer feature config-r
    # to enable global cache accessible after
    # image build.
    # just use renv defaults now.
    renv_paths_cache_match <- Sys.getenv("RENV_PATHS_CACHE") == "/renv/cache"
    renv_paths_library_root_match <- Sys.getenv("RENV_PATHS_LIBRARY_ROOT") ==
      "/workspaces/.local/lib/R/library"
    renv_paths_root_match <- Sys.getenv("RENV_PATHS_ROOT") == "/renv/local"
    renv_paths_match <- renv_paths_cache_match &&
      renv_paths_library_root_match && renv_paths_root_match
    if (renv_paths_match) {
      Sys.unsetenv("RENV_PATHS_CACHE")
      Sys.unsetenv("RENV_PATHS_LIBRARY_ROOT")
      Sys.unsetenv("RENV_PATHS_ROOT")
    }
    # if not using container cache, then likely fine
    # to use pak.
    # may be an issue in RAM-limited environments,
    # such as GitHub Actions, but we'll see
    # (building from source may be an issue
    # whether or not we use pak)
    # options("renv.config.pak.enabled" = TRUE)
  }
}

if (Sys.getenv("GITHUB_ACTIONS") != "true") {
  # only use renv when not in GitHub Actions
  .set_renv_profile()
  source("renv/activate.R")
}

