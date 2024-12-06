% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate-plot_gate--sub.R
\name{.plot_logit}
\alias{.plot_logit}
\title{Plot proportion positive with fitted log reg line}
\usage{
.plot_logit(high_ind_tbl, cp, plot_prop)
}
\arguments{
\item{plot_prop}{ggplot2 object. Plot returned by \code{.plot_prop}.}
}
\value{
A ggplot2 object.
}
\description{
This is the same as \code{plot_prop}, except that instead of a
smoothed line of proportion positive, the fitted values
from the logistic regression used to obtain the mid-probability cut
is displayed instead.
}
