% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{pnorm_outliers}
\alias{pnorm_outliers}
\title{Detect outliers in a named vector using normalized z-score.}
\usage{
pnorm_outliers(x)
}
\arguments{
\item{x}{A named vector of values.}
}
\value{
A table of all the results including intermediate tables.
}
\description{
Function to detect outliers in a named vector of values.
It uses CDF of pnorm to get p-values
}
\examples{
pnorm_outliers(c(2, 4, 7, 8,9))
}
