test_that("analysis figures have no title or subtitle plumbing", {
  root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
  files <- c(
    list.files(file.path(root, "scripts", "r"), pattern = "\\.R$", full.names = TRUE),
    list.files(file.path(root, "analysis"), pattern = "\\.qmd$", full.names = TRUE)
  )
  files <- files[basename(files) != "sim-debug-loc.R"]
  check_labels <- function(expr, file) {
    if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return(invisible(NULL))
    if (is.call(expr)) {
      fn <- paste(deparse(expr[[1]]), collapse = "")
      expect_false(grepl("(^|::)ggtitle$", fn), info = file)
      if (grepl("(^|::)labs$", fn)) {
        expect_false(any(names(as.list(expr)[-1]) %in% c("title", "subtitle")), info = file)
      }
    }
    parts <- as.list(expr)
    for (i in seq_along(parts)) {
      # Required formals and omitted call arguments use the missing symbol.
      if (!identical(parts[[i]], quote(expr = ))) check_labels(parts[[i]], file)
    }
    invisible(NULL)
  }
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    if (endsWith(file, ".qmd")) {
      in_r <- FALSE
      code <- character()
      for (line in lines) {
        if (grepl("^```\\{r([ ,}]|$)", line)) {
          in_r <- TRUE
        } else if (in_r && grepl("^```\\s*$", line)) {
          in_r <- FALSE
          code <- c(code, "")
        } else if (in_r) code <- c(code, line)
      }
    } else code <- lines
    expect_false(any(grepl("\\bsubtitle\\s*=", code)), info = file)
    check_labels(parse(text = code), file)
  }
})
