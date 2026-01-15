#' Check if object is a non-empty string
#'
#' @param x object Object to check
#' @return logical TRUE if x is a character vector of length 1 with non-zero characters
#' @keywords internal
.is_string <- function(x) {
  is.character(x) && length(x) == 1 && nzchar(x)
}

#' Assert that object is a non-empty string
#'
#' @param x object Object to check
#' @return invisible NULL if check passes, otherwise throws an error
#' @keywords internal
.assert_string <- function(x) {
  if (!.is_string(x)) {
    stop(paste0("Expected a non-empty string, got: ", class(x)))
  }
}

#' Check if object is a non-empty character vector
#'
#' @param x object Object to check
#' @return logical TRUE if x is a character vector of length >= 1 with all non-zero character elements
#' @keywords internal
.is_string_vector <- function(x) {
  is.character(x) && length(x) >= 1 && all(nzchar(x))
}

#' Assert that object is a non-empty character vector
#'
#' @param x object Object to check
#' @return invisible NULL if check passes, otherwise throws an error
#' @keywords internal
.assert_string_vector <- function(x) {
  if (!.is_string_vector(x)) {
    stop(paste0("Expected a non-empty character vector, got class: ", class(x)))
  }
}
