if (Sys.getenv("BIOCONDUCTOR_NAME") == "bioconductor_docker") {
  Sys.setenv("RENV_PROFILE" = "bioc_container")
} else {
  Sys.setenv("BIOCONDUCTOR_USE_CONTAINER_REPOSITORY" = "false")
  Sys.setenv("RENV_PROFILE" = "non_bioc_container")
}
source("renv/activate.R")
