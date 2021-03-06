% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mvrnormArma}
\alias{mvrnormArma}
\title{Sample from multivariate normal distribution with C++}
\usage{
mvrnormArma(n, mu, sigma, chol)
}
\arguments{
\item{n}{number of simulations.}

\item{mu}{mean vector.}

\item{sigma}{covariance matrix or Cholesky decomposition of the matrix (see chol).}

\item{chol}{integer, if 0 sigma is a covariance matrix, otherwise it is the Cholesky decomposition of the matrix.}
}
\value{
A matrix of size \eqn{d x n} containing the samples.
}
\description{
Simulates realizations from a multivariate normal with mean mu and covariance matrix sigma.
}

