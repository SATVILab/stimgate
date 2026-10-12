# Tables destined for the thesis (Manuscript22RodoThesis), written the way the
# thesis writes its own: a Quarto pipe table followed by a caption line
# `: **Title**. Caption. {#tbl-label}`, so Quarto numbers it and `@tbl-label`
# cross-references it. Maths uses LaTeX-style `$...$`, which Quarto converts.
# Analyses render the table in the HTML, show its source in a copyable block
# and save it (as .md) in the projr `output-thesis` folder only for tables the
# thesis uses. Build the cells from formatted values, so the HTML and the
# saved source show the same numbers.

# Escape characters that would break a pipe-table cell.
.analysis_md_escape <- function(x) {
  gsub("|", "\\|", x, fixed = TRUE)
}

# `cells`: character matrix or data frame of Markdown; `header`: column
# headings; `align`: "left" or "right" per column; `bold`: rows set in bold.
.analysis_md_table <- function(cells, header, caption, label,
                               align = rep("left", length(header)),
                               bold = rep(FALSE, nrow(cells))) {
  cells <- as.matrix(cells)
  if (ncol(cells) != length(header) || length(align) != length(header) ||
      length(bold) != nrow(cells)) {
    stop("Table cells, header, align and bold do not match.")
  }
  row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
  rule <- ifelse(align == "right", "---:", ":---")
  body <- vapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    if (isTRUE(bold[[i]])) cell <- ifelse(nzchar(cell), paste0("**", cell, "**"), cell)
    row(cell)
  }, character(1))
  paste(c(row(header), row(rule), body, "", paste0(": ", caption, " {#", label, "}")),
    collapse = "\n")
}

# Render the table in the HTML (`results: asis` chunk), then show its source in
# a copyable code block (the QMD sets `code-copy: true`); save it to `path`
# when given.
.analysis_md_show <- function(source, path = NULL) {
  if (!is.null(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(source, path)
  }
  cat("\n\n", source, "\n\n", sep = "")
  cat("Source for the thesis (copy with the button at the top right):\n\n")
  cat("````markdown\n", source, "\n````\n\n", sep = "")
  invisible(source)
}
