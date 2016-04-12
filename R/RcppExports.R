# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Sample from multivariate normal distribution with C++
#' @description Simulates realizations from a multivariate normal with mean mu and covariance matrix sigma.
#' @param n number of simulations.
#' @param mu mean vector.
#' @param sigma covariance matrix or Cholesky decomposition of the matrix (see chol).
#' @param chol integer, if 0 sigma is a covariance matrix, otherwise it is the Cholesky decomposition of the matrix.
#' @return A matrix containing of size dxn containing the samples.
#' @export
mvrnormArma <- function(n, mu, sigma, chol) {
    .Call('ConservativeEstimates_mvrnormArma', PACKAGE = 'ConservativeEstimates', n, mu, sigma, chol)
}

#' @title Sample from truncated multivariate normal distribution with C++
#' @description Simulates realizations from a truncated multivariate normal with mean mu, covariance matrix sigma in the bounds lower upper.
#' @param n number of simulations.
#' @param mu mean vector.
#' @param sigma covariance matrix.
#' @param lower vector of lower bounds.
#' @param upper vector of upper bounds.
#' @param verb level of verbosity: if lower than 3 nothing, 3 minimal, 4 extended.
#' @return A matrix containing of size dxn containing the samples.
#' @export
trmvrnorm_rej_cpp <- function(n, mu, sigma, lower, upper, verb) {
    .Call('ConservativeEstimates_trmvrnorm_rej_cpp', PACKAGE = 'ConservativeEstimates', n, mu, sigma, lower, upper, verb)
}
