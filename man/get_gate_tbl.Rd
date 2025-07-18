% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gates.R
\name{get_gate_tbl}
\alias{get_gate_tbl}
\title{Get gates}
\usage{
get_gate_tbl(path_project)
}
\arguments{
\item{path_project}{character. Path to the project directory.}
}
\value{
Gate table with gates for each sample for each marker.
}
\description{
Get all the gates for each of the markers gated.
}
\examples{
\donttest{
# Get gate table from project directory
gate_tbl <- get_gate_tbl("/path/to/project")
}
}
