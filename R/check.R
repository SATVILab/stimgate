.is_string <- function(x) {
  is.character(x) && length(x) == 1 && nzchar(x)
}

.assert_string <- function(x) {
  if (!.is_string(x)) {
    stop(paste0("Expected a non-empty string, got: ", class(x)))
  }
}

.is_string_vector <- function(x) {
  is.character(x) && length(x) >= 1 && all(nzchar(x))
}

.assert_string_vector <- function(x) {
  if (!.is_string_vector(x)) {
    stop(paste0("Expected a non-empty character vector, got class: ", class(x)))
  }
}