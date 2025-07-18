## 🧬 `stimgate`: Identify Responding Cells via Outlier Gating


<!-- badges: start -->
[![R-CMD-check](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SATVILab/stimgate/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/SATVILab/stimgate/graph/badge.svg)](https://app.codecov.io/gh/SATVILab/stimgate)
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
devtools::install_deps(repos = BiocManager::repositories())

# Install stimgate from GitHub
devtools::install_github("SATVILab/stimgate")
```

---

## 🚀 Quick Start

```r
library(stimgate)

# Load your GatingSet (flow cytometry data)
gs <- flowWorkspace::load_gs("path/to/gatingset")

# Define batch structure and markers to gate
batch_list <- list(
  batch1 = seq(1, 10),   # Sample indices for batch 1
  batch2 = seq(11, 20)   # Sample indices for batch 2
)

# Channel names to gate on
marker_channels <- c("IL2", "TNFa", "IFNg")

# Run the stimgate pipeline
result <- stimgate_gate(
  path_project = "/path/to/project",
  .data = gs,
  batch_list = batch_list,
  marker = marker_channels,
  pop_gate = "root"
)

# Get statistics for the identified gates
stats <- get_stats("/path/to/project")

# Extract gate information
gates <- get_gate_tbl("/path/to/project")

# Visualize results
plots <- plot_gate(
  ind = seq_len(3),
  .data = gs,
  path_project = "/path/to/project",
  marker = c("IL2", "TNFa")
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
* 🔗 **GitHub**: [github.com/SATVILab/stimgate](https://github.com/SATVILab/stimgate)
* 📦 **Bioconductor**: *URL to be added upon release*
* 🔧 **Function reference**: `help(package = "stimgate")`

### Key Functions

- `stimgate_gate()`: Main gating function to identify cytokine-positive cells
- `get_stats()`: Generate comprehensive statistics from gating results  
- `plot_gate()`: Create bivariate hex and univariate density plots with gate overlays
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
