#' @title ConservativeEstimates package
#' @description Conservative estimation of excursion sets under Gaussian models.
#' @details Package: ConservativeEstimates \cr
#' Type: Package \cr
#' Version: 0.1.1 \cr
#' Date: 2016-05-03
#'
#' @author Dario Azzimonti (dario.azzimonti@@gmail.com) . Thanks to David Ginsbourger for the fruitful discussions and his help in testing and improving the package.
#' @docType package
#' @name ConservativeEstimates
#' @import microbenchmark mvtnorm
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dist lm pnorm quantile var
#' @useDynLib ConservativeEstimates
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Bolin, D. and Lindgren, F. (2015). Excursion and contour uncertainty regions for latent Gaussian models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1):85--106.
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Dickmann, F. and Schweizer, N. (2014). Faster comparison of stopping times by nested conditional Monte Carlo. arXiv preprint arXiv:1402.0243.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#'
#' Genz, A. and Bretz, F. (2009). Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics 195. Springer-Verlag.
#'
#' Horrace, W. C. (2005). Some results on the multivariate truncated normal distribution. Journal of Multivariate Analysis, 94(1):209--221.
#'
#' Robert, C. P. (1995). Simulation of truncated normal variables. Statistics and Computing, 5(2):121--125.
NULL
