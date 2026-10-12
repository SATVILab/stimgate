# Tables destined for the thesis (Manuscript22RodoThesis, Quarto
# stanford-thesis-typst format). See AGENTS.md, "Thesis tables". Each table is a
# Quarto pipe table with a caption and a stable `{#tbl-...}` identifier, so
# Quarto numbers it, cross-references it and lists it, wrapped in a locally
# scoped Typst block that draws minimal academic rules: one below the
# headings, one at the bottom and thin pale-grey lines between row groups. The
# Typst block is dropped from HTML renders, so the HTML shows the plain table.
# Build cells from formatted values only (round for display, never the
# results). Analyses render the table, show its source in a copyable block and
# save it (as .md) in the projr `output-thesis` folder for tables the thesis uses.

# Escape characters that would break a pipe-table cell.
.analysis_md_escape <- function(x) {
  gsub("|", "\\|", x, fixed = TRUE)
}

# `cells`: character matrix of Markdown; `header`: short column headings;
# `align`: "left" or "right" per column; `bold`: logical matrix of cells set
# in bold (e.g. the selected setting's parameter cell); `caption`: sentences,
# one per source line; `groups`: row-group sizes (e.g. rows per method), for
# the separators. Returns the complete source: Typst opener, table, closer.
.analysis_md_table <- function(cells, header, caption, label,
                               align = rep("left", length(header)),
                               bold = NULL, groups = nrow(cells)) {
  cells <- as.matrix(cells)
  if (is.null(bold)) bold <- matrix(FALSE, nrow(cells), ncol(cells))
  if (ncol(cells) != length(header) || length(align) != length(header) ||
      !identical(dim(bold), dim(cells)) || sum(groups) != nrow(cells)) {
    stop("Table cells, header, align, bold and groups do not match.")
  }
  if (any(grepl("\n", cells, fixed = TRUE))) stop("Table cells must be single lines.")
  cells[bold & nzchar(cells)] <- paste0("**", cells[bold & nzchar(cells)], "**")
  mdRow <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
  rule <- ifelse(align == "right", "---:", ":---")
  body <- vapply(seq_len(nrow(cells)), function(i) mdRow(cells[i, ]), character(1))
  captionLines <- c(paste0(": ", caption[[1]]), caption[-1])
  last <- length(captionLines)
  captionLines[[last]] <- paste0(captionLines[[last]], " {#", label, "}")
  paste(c(
    .analysis_typst_rules_open(nrow(cells), groups), "",
    mdRow(header), mdRow(rule), body, "", captionLines, "",
    .analysis_typst_rules_close()
  ), collapse = "\n")
}

# Locally scoped Typst rules for a table with `n_rows` body rows. Typst
# numbers the header row 0, so the bottom rule sits under row `n_rows` and a
# separator under the last row of each group but the last.
.analysis_typst_rules_open <- function(n_rows, groups = n_rows) {
  separators <- utils::head(cumsum(groups), -1L)
  separator <- if (length(separators)) {
    paste0(
      "    } else if ", paste0("y == ", separators, collapse = " or "), " {\n",
      "      0.3pt + rgb(\"#D6D6D6\")\n"
    )
  } else {
    ""
  }
  paste0(
    "```{=typst}\n",
    "#[\n",
    "#set table(\n",
    "  stroke: (x, y) => (\n",
    "    bottom: if y == 0 or y == ", n_rows, " {\n",
    "      0.65pt + black\n",
    separator,
    "    } else {\n",
    "      none\n",
    "    },\n",
    "    rest: none,\n",
    "  ),\n",
    ")\n",
    "```"
  )
}

.analysis_typst_rules_close <- function() {
  "```{=typst}\n]\n```"
}

# Tolerances spanning orders of magnitude as powers of ten, e.g. 10^{-2} and
# 10^{-2.5} (a grid of half steps on the log10 scale).
.analysis_md_pow10 <- function(x, digits = 1) {
  exponent <- round(log10(x), digits)
  paste0("10^{", format(exponent, trim = TRUE, drop0trailing = TRUE), "}")
}

# Fixed decimals with leading zeros, the same number in a column; NA stays
# visible as "NA" rather than becoming zero.
.analysis_md_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

# Render the table in the HTML (`results: asis` chunk; Typst blocks are not
# shown there), then show its source in a copyable code block (the QMD sets
# `code-copy: true`); save it to `path` when given.
.analysis_md_show <- function(source, path = NULL) {
  if (!is.null(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(source, path)
  }
  cat("\n\n", source, "\n\n", sep = "")
  cat("Source for the thesis, with its Typst rule styling (copy with the button at the top right):\n\n")
  cat("`````markdown\n", source, "\n`````\n\n", sep = "")
  invisible(source)
}
