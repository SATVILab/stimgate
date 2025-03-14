% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate-get_pos_ind.R
\name{.get_pos_ind_mult}
\alias{.get_pos_ind_mult}
\title{Identify cells that are positive for at least two cytokines.}
\usage{
.get_pos_ind_mult(
  ex,
  gate_tbl,
  chnl = NULL,
  chnl_alt = NULL,
  gate_type_cyt_pos
)
}
\arguments{
\item{ex}{dataframe. Expression data.}

\item{gate_tbl}{dataframe.
Contains gates for specific
sample for specific gate method only (but for
all markers of interest).}

\item{chnl}{character vector. Cytokines that a cell must be positive for
in any combination (exclusively, i.e. a cell need not be positive for both
simultaneously if both are placed here). If \code{NULL},
then set to \code{unique(gate_tbl$chnl)}. Default is \code{NULL}.}

\item{chnl_alt}{character vector. Cytokine(s) that a cell must be positive
for in conjunction with at least one cytokine in \code{chnl}. If \code{NULL},
then set to whatever \code{chnl} is, which is either specified
by the parameter or set to  \code{setdiff(unique(gate_tbl$chnl), <chnl_crr>)} if \code{chnl} parameter is NULL,
where \code{<chnl_curr>} is one of the channels in \code{chnl} parameter. Default is \code{NULL}.}

\item{gate_type_cyt_pos}{"base" or "cyt". If \code{"base"}, then the initial thresholds are
used even for cells that are cytokine positive for at least one other cytokine.
If \code{"cyt"}, then the cytokine+ thresholds are used for cells that are cytokine
positive for at least one other cytokine.}
}
\value{
A logical vector with value \code{TRUE} if a cell is positive for
at least one cytokine in \code{chnl} and any cytokine in \code{chnl_alt}, and \code{FALSE}
otherwise.
}
\description{
Identify cells that express multiple cytokines. May
force each cytokine combination to contain at least one of a set of
cytokines (using \code{chnl} parameter) and specify which cytokines a cell
must express at least one of (using \code{chnl_alt} parameter). The
default is to consider a cell multi-functional if it expresses any cytokine
in conjunction with any other cytokine.
}
