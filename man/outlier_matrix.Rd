% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{outlier_matrix}
\alias{outlier_matrix}
\title{Detect outliers in a matrix.}
\usage{
outlier_matrix(m, type = "pnorm", name = "eRNA")
}
\arguments{
\item{m}{matrix of expression values gene/eRNA in column and sample in rows.
It expects the row.names to set to gene IDs.}

\item{type}{String indicating the p-value estimation method.
"pnorm" or "permutation"}

\item{name}{The name to use for the new gene column.
This option is to differentiate eRNA/mRNA runs.}
}
\value{
Table of results based on the output of the outlier functions.
}
\description{
Function to apply the two outlier function to a matrix.
This is specific to eRNAkit usage
}
\examples{
outlier_matrix(matrix(1:9, nrow = 3, ncol = 3), "pnorm")
}
