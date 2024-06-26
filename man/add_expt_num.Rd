% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.R
\name{add_expt_num}
\alias{add_expt_num}
\title{Set a distinct number for each experiment}
\usage{
add_expt_num(tb)
}
\arguments{
\item{tb}{A tibble of dose effect data for many experiments with a cell_line column. In this column, each experiment has a cell line name in the first row followed by blank rows.}
}
\value{
\item{A tibble}{It has a distinct number for each expt.}
}
\description{
Each experiment has a cell line names followed by blank rows. Use this information to set a distinct number for each experiment.
}
\examples{
library(decurve)
tb <- add_expt_num(tb)
}
