## 🧬 `stimgate`: Identify Responding Cells via Outlier Gating


<!-- badges: start -->
[![R-CMD-check](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/SATVILab/stimgate/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SATVILab/stimgate)
[![pkgdown](https://github.com/SATVILab/stimgate/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/SATVILab/stimgate/actions/workflows/pkgdown.yaml)
<!-- badges: end -->


`stimgate` is an R package for identifying cell populations that have responded to stimulation. It applies outlier gating functions to flow cytometry data, comparing unstimulated and stimulated tubes from the same sample to isolate stimulation-specific responses.

---

## 📦 Installation

```r
# Install BiocManager and devtools if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install all dependencies using Bioconductor and CRAN repos
devtools::install_github(
  "SATVILab/stimgate",
  repos = BiocManager::repositories(), dependencies = TRUE
)
```

---

## 🚀 Quick Start

```r
library(stimgate)

# Get example dataset
example_data <- get_example_data()
gs <- flowWorkspace::load_gs(example_data$path_gs)

# Run the stimgate pipeline
path_project <- stimgate_gate(
  path_project = file.path(tempdir(), "stimgate_example"),
  .data = gs,
  batch_list = example_data$batch_list,
  marker = example_data$marker,
  pop_gate = "root"
)

# Get statistics for the identified gates
stats <- get_stats(path_project)

# Extract gate information
gates <- get_gate_tbl(path_project)

# Visualize results
stimgate_plot(
  ind = seq_len(2),
  .data = gs,
  path_project = path_project,
  marker = example_data$marker,
)
```

---

## 🔑 Features

* **Outlier-based gating**: Identifies responding cells by comparing stimulated vs unstimulated samples
* **Batch processing**: Handles multiple donors or experimental batches simultaneously  
* **Background subtraction**: Reduces false positives from baseline cytokine activity
* **Flexible gating**: Supports various gating parameters and combinations
* **Comprehensive output**: Generates statistics, gate tables, and visualization plots
* **FCS export**: Write cytokine-positive cells to new FCS files for downstream analysis
* **Bioconductor compatibility**: Designed for integration with `flowCore`, `flowWorkspace`, and `CytoML`

---

## 📖 Documentation

* 📘 **Vignette**: `vignette("stimgate")` - Getting started guide
* 🌐 **Website**: [satvilab.github.io/stimgate](https://satvilab.github.io/stimgate/)
* 🔗 **GitHub**: [github.com/SATVILab/stimgate](https://github.com/SATVILab/stimgate)
* 📦 **Bioconductor**: *URL to be added upon release*
* 🔧 **Function reference**: `help(package = "stimgate")`

### Key Functions

- `stimgate_gate()`: Main gating function to identify cytokine-positive cells
- `get_stats()`: Generate comprehensive statistics from gating results  
- `stimgate_plot()`: Create bivariate hex and univariate density plots with gate overlays
- `get_gate_tbl()`: Extract gate thresholds and parameters
- `stimgate_fcs_write()`: Export cytokine-positive cells as FCS files

---

## 📌 Citation

```r
citation("stimgate")
```

> Rodo M. (2024). *stimgate: Identify responding cells as outliers*. R package version 0.3.1. https://github.com/SATVILab/stimgate

---

## 🤝 Contributing

We welcome bug reports, feature requests, and pull requests via [GitHub Issues](https://github.com/SATVILab/stimgate/issues).

### Development Guidelines

Before submitting changes, please ensure:
- Run `devtools::test()` to verify all tests pass
- Run `devtools::document()` to update documentation  
- Run `styler::style_pkg()` for consistent code formatting
- Run `lintr::lint_package()` to check for style violations

---

## 📄 License

This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.
