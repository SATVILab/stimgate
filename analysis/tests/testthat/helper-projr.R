# Temporary projr project for tests that write analysis outputs. projr resolves
# paths from the project root and needs a VERSION file for "output", which it
# then places in its cache at _tmp/projr/v<version>/output.
.local_projr_root <- function(env = parent.frame()) {
  root <- normalizePath(withr::local_tempdir(.local_envir = env), winslash = "/")
  writeLines("0.0.1", file.path(root, "VERSION"))
  root
}

.projr_output_path <- function(root, ...) {
  file.path(root, "_tmp", "projr", "v0.0.1", "output", ...)
}
