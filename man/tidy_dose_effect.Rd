% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.R
\name{tidy_dose_effect}
\alias{tidy_dose_effect}
\title{Tidy dose effect data of a single experiment}
\usage{
tidy_dose_effect(tb)
}
\arguments{
\item{tb}{A tibble of dose effect data. Each element is of string type.

Columns: cell_line, drug1, drug2, effect1_..., effect12_..., effect2_....

Rows:

first row: cell line name, drug1 name, drug2 name, backgrounds.

other rows: doses and effects.}
}
\value{
\item{A tidied tibble}{Columns: cell_line, drug1, drug2, dose1, dose2, effect}
}
\description{
tidy data to have attributes: cell_line, drug1, drug2, dose1, dose2, effect
}
\examples{
library(decurve)
tb <- tidy_dose_effect(tb)
}
