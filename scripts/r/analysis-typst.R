# Typst source for tables destined for the thesis (Manuscript22RodoThesis,
# stanford-thesis Typst template). Analyses show the source in a `typst` code
# block (copy button) and save it in the projr `output-thesis` folder only for
# tables the thesis uses. Build the cells from the same formatted values as the
# HTML table so both show the same numbers.

# Escape text for Typst markup (not for cells that already hold Typst markup,
# such as `$beta = 0.8$`).
.analysis_typst_escape <- function(x) {
  gsub("([\\\\#*_$@<>\\[\\]`~])", "\\\\\\1", x, perl = TRUE)
}

# A three-rule `#figure(table(...))`. `cells` is a character matrix or data
# frame of Typst markup, `header` its column headings (`\\` breaks a line) and
# `bold` marks rows set in bold. Text is 9 pt and unjustified, so narrow
# columns do not hyphenate at the thesis's 12 pt body size.
.analysis_typst_table <- function(cells, header, caption, label,
                                  align = rep("left", length(header)),
                                  bold = rep(FALSE, nrow(cells))) {
  cells <- as.matrix(cells)
  if (ncol(cells) != length(header) || length(align) != length(header) ||
      length(bold) != nrow(cells)) {
    stop("Typst table cells, header, align and bold do not match.")
  }
  rows <- vapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    if (isTRUE(bold[[i]])) cell <- paste0("*", cell, "*")
    paste0("    ", paste0("[", cell, "]", collapse = ", "), ",")
  }, character(1))
  paste(c(
    "#figure(",
    "  {",
    "  set text(size: 9pt)",
    "  set par(justify: false)",
    "  table(",
    paste0("    columns: ", length(header), ","),
    paste0("    align: (", paste(align, collapse = ", "), "),"),
    "    stroke: none,",
    "    inset: (x: 4pt, y: 3pt),",
    "    table.hline(),",
    paste0("    table.header(", paste0("[", header, "]", collapse = ", "), "),"),
    "    table.hline(stroke: 0.5pt),",
    rows,
    "    table.hline(),",
    "  )",
    "  },",
    paste0("  caption: [", caption, "],"),
    paste0(") <", label, ">")
  ), collapse = "\n")
}

# `1e-05` style numbers as Typst math, e.g. `$1 times 10^(-5)$`.
.analysis_typst_sci <- function(x, digits = 2) {
  exponent <- floor(log10(abs(x)))
  mantissa <- vapply(signif(x / 10^exponent, digits), format, character(1))
  paste0("$", mantissa, " times 10^(", exponent, ")$")
}

# Show Typst source in a copyable code block (`results: asis` chunk; the QMD
# sets `code-copy: true`) and save it to `path` when given.
.analysis_typst_show <- function(source, path = NULL) {
  if (!is.null(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(source, path)
  }
  cat("\n```typst\n", source, "\n```\n\n", sep = "")
  invisible(source)
}
