.debug <- function(debug, msg, val = NULL) {
  if (!is.null(val)) {
    msg <- paste0(msg, ": ", val)
  }
  message(msg)
}
