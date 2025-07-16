.debug_msg <- function(.debug, msg, val = NULL) {
  if (!debug) {
    return(invisible(FALSE))
  }
  if (!is.null(val)) {
    msg <- paste0(msg, ": ", val)
  }
  message(msg)
  invisible(TRUE)
}
