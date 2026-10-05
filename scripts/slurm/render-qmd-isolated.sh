#!/usr/bin/env bash
# Render a QMD from a temporary, job-specific copy in the same folder.
#
# Quarto names its working files (e.g. <name>.knit.md, <name>_files/) after
# the QMD and keeps them next to it. Chunk jobs of one analysis render the
# same QMD at the same time, so without this they collide over those files and
# all but the last job to finish fail at Quarto's final step. The copy is
# identical to the QMD, sits in the same folder (so the QMD still finds the
# checkout root) and is removed with its outputs when the job exits.
#
# Usage (from the project root): render_qmd_isolated <qmd_file> <tag>

render_qmd_isolated() {
  local qmd_file="$1"
  local tag="$2"
  local qmd_dir qmd_stem render_stem
  qmd_dir=$(dirname -- "$qmd_file")
  qmd_stem=$(basename -- "$qmd_file" .qmd)
  render_stem="${qmd_dir}/${qmd_stem}--${tag}"
  RENDER_QMD_ISOLATED_STEM="$render_stem"
  cp -- "$qmd_file" "${render_stem}.qmd"
  trap 'rm -rf -- "${RENDER_QMD_ISOLATED_STEM}.qmd" "${RENDER_QMD_ISOLATED_STEM}.html" "${RENDER_QMD_ISOLATED_STEM}.knit.md" "${RENDER_QMD_ISOLATED_STEM}_files"' EXIT

  local r_expr="qmd_file <- '${render_stem}.qmd'; if (requireNamespace('quarto', quietly = TRUE)) { quarto::quarto_render(input = qmd_file) } else { status <- system2('quarto', c('render', qmd_file)); if (!identical(status, 0L)) quit(status = status) }"
  apptainer-rscript -f stimgate -- "$r_expr"
}
