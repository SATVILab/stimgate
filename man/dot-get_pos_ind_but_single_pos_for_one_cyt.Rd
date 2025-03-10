% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate-get_pos_ind.R
\name{.get_pos_ind_but_single_pos_for_one_cyt}
\alias{.get_pos_ind_but_single_pos_for_one_cyt}
\title{Identify cells that are positive for every cytokine except one}
\usage{
.get_pos_ind_but_single_pos_for_one_cyt(
  ex,
  gate_tbl,
  chnl_single_exc,
  chnl = NULL,
  gate_type_cyt_pos,
  gate_type_single_pos
)
}
\arguments{
\item{ex}{dataframe. Expression data.}

\item{gate_tbl}{dataframe.
Contains gates for specific
sample for specific gate method only (but for
all markers of interest).}

\item{chnl_single_exc}{character.
Channel that is the only one that a cell
may be left positive for and single-positive (i.e.
a cell cannot be positive for any other channel
on its own or this channel with any other channel).}

\item{chnl}{character.
Channels for which the cells
may be positive for. Default is \code{NULL},
in which case every chnl in unique(gate_tbl$chnl)
other than \code{chnl_single_exc} will be considered a channel that
a cell can be positive for.}
}
\value{
A logical vector, with TRUE for every cell that
is negative for every other cytokine and for the cytokine
in question with
any other cytokine. FALSE otherwise.
}
\description{
Identify cells that are positive for every
cytokine except one,
for which they may be single-positive.
}
