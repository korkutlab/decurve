% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.R
\name{effect2GR}
\alias{effect2GR}
\title{Calculate growth rates for a single experiment}
\usage{
effect2GR(tb, assay_duration, doubling_time)
}
\arguments{
\item{tb}{A tibble of dose effect data.

Columns: cell_line (string), drug1 (string), drug2 (string), dose1 (double), dose2 (double), effect (double)}

\item{assay_duration}{Assay duration.}

\item{doubling_time}{Doubling time of untreated cells.}
}
\value{
\item{A tibble}{Columns: cell_line, drug1, drug2, dose1, dose2, GR}
}
\description{
Given assay duration and doubling time, convert effect to GR for each dose.
}
\details{
GR(c) = 2^( (log2( x(c)/x(0) ) + T/Td) / (T/Td) ) - 1
= 2^( (log2( effect(c) ) + T/Td) / (T/Td) ) - 1

x's are the number of cells at the end of the assay
in an untreated (or vehicle-treated) control well (x(0)) and
in a drug-treated well (x(c)).

T is assay duration.

Td is doubling time of untreated cells.
}
\examples{
library(decurve)
tb <- effect2gr(tb, doubling_time = 44, assay_duration = 120)
}
