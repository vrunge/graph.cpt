% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ts_to_dists}
\alias{ts_to_dists}
\title{ts_to_dists}
\usage{
ts_to_dists(data, states, window_length)
}
\arguments{
\item{data}{matrix of data, size p x n for p dimensions and n data points}

\item{states}{list of reference matrices of size d}

\item{window_length}{the number of data points to take for computing each covariance matrix}
}
\value{
A matrix of distances between SPD matrices of size d x n
}
\description{
function combining the functions ts_to_SPDts and SPDts_to_dists
}
