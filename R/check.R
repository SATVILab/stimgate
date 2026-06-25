#' Check if object is a non-empty string
#'
#' @param x object Object to check
#' @return logical TRUE if x is a character vector of length 1 with non-zero characters
#' @keywords internal
.isString <- function(x) {
  is.character(x) && length(x) == 1 && nzchar(x)
}

#' Assert that object is a non-empty string
#'
#' @param x object Object to check
#' @return invisible NULL if check passes, otherwise throws an error
#' @keywords internal
.assertString <- function(x) {
  if (!.isString(x)) {
    stop(paste0("Expected a non-empty string, got: ", class(x)))
  }
}

#' Check if object is a non-empty character vector
#'
#' @param x object Object to check
#' @return logical TRUE if x is a character vector of length >= 1 with all non-zero character elements
#' @keywords internal
.isStringVector <- function(x) {
  is.character(x) && length(x) >= 1 && all(nzchar(x))
}

#' Assert that object is a non-empty character vector
#'
#' @param x object Object to check
#' @return invisible NULL if check passes, otherwise throws an error
#' @keywords internal
.assertStringVector <- function(x) {
  if (!.isStringVector(x)) {
    stop(paste0("Expected a non-empty character vector, got class: ", class(x)))
  }
}
