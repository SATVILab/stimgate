% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate-gate_cp_uns_loc.R
\name{.get_cp_uns_loc_sample}
\alias{.get_cp_uns_loc_sample}
\title{Get cutpoint for a range of samples given the q-value and fdr}
\usage{
.get_cp_uns_loc_sample(
  cut_stim,
  cut_uns,
  cp_min,
  bw_min,
  ind_uns,
  ind_gate,
  min_cell,
  params,
  plot = TRUE,
  bias,
  path_project,
  debug = FALSE
)
}
\arguments{
\item{cut_stim}{list. List where each element are the marker expression readings
of the marker to be cut on for the cells in a sample. Note that the i-th element
in \code{cut_stim} must correspond to the i-th element in \code{q_list}, i.e.
must be related to the same marker in same cell population from the same blood sample and stimulation.}

\item{fdr}{numeric. A value between 0 and 1 specifying the false discovery rate
the sample should be cut at.}
}
\value{
Numeric vector. A cutpoint for each sample.
}
\description{
Calculate the cutpoint for each sample in a batch at a given FDR.
}
