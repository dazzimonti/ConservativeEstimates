#include <RcppArmadillo.h>
#include <algorithm>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//' @title Sample from multivariate normal distribution with C++
//' @description Simulates realizations from a multivariate normal with mean mu and covariance matrix sigma.
//' @param n number of simulations.
//' @param mu mean vector.
//' @param sigma covariance matrix or Cholesky decomposition of the matrix (see chol).
//' @param chol integer, if 0 sigma is a covariance matrix, otherwise it is the Cholesky decomposition of the matrix.
//' @return A matrix containing of size dxn containing the samples.
//' @export
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu,arma::mat sigma,int chol) {
  int ncols = mu.n_elem;
  arma::mat Y = arma::randn(ncols, n);
  if(chol){
    return arma::repmat(mu, 1, n) + sigma*Y;
  }else{
    return arma::repmat(mu, 1, n) + arma::chol(sigma)*Y;
  }
  return arma::repmat(mu, 1, n) + arma::chol(sigma)*Y;
}

//' @title Sample from truncated multivariate normal distribution with C++
//' @description Simulates realizations from a truncated multivariate normal with mean mu, covariance matrix sigma in the bounds lower upper.
//' @param n number of simulations.
//' @param mu mean vector.
//' @param sigma covariance matrix.
//' @param lower vector of lower bounds.
//' @param upper vector of upper bounds.
//' @param verb level of verbosity: if lower than 3 nothing, 3 minimal, 4 extended.
//' @return A matrix containing of size dxn containing the samples.
//' @export
// [[Rcpp::export]]
arma::mat trmvrnorm_rej_cpp(int n, arma::vec mu,arma::mat sigma, arma::vec lower, arma::vec upper,int verb) {
  int ncols = sigma.n_cols;
  arma::mat Y(ncols, n);
  arma::mat cholSigma(ncols,ncols);
  int samplesRemaining = n;
  int totalAccepted = 0;
  Environment mvtnorm("package:mvtnorm");
  Function mypmvn =  mvtnorm["pmvnorm"];
  double alpha = as<double>(mypmvn(NumericVector(lower.begin(),lower.end()),NumericVector(upper.begin(),upper.end()),NumericVector(mu.begin(),mu.end()),R_NilValue,sigma));
  cholSigma = arma::chol(sigma);
  if(verb>=3){
    Rcout << "Acceptance rate: " << alpha << "\n";
  }
//    if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
//    estimatedAlpha <- TRUE
int totalSamplesRun=0;
while(samplesRemaining>0)
    {
    int nPossibleSamples = ((samplesRemaining/alpha > 1000000) ? (samplesRemaining) : (std::max( (int) (ceil(samplesRemaining/alpha)),10)));
    arma::mat X = mvrnormArma(nPossibleSamples,mu,cholSigma,1);
    totalSamplesRun = totalSamplesRun+nPossibleSamples;
    arma::uvec goodIndices(nPossibleSamples);
    goodIndices.fill(-1);
    int j=0;
    for(int i=0; i<nPossibleSamples; i++){
      if(all(X.col(i) >= lower ) && all(X.col(i) <= upper )){
        goodIndices[j]=i;
        j++;
      }
    }
    goodIndices.resize(j);

    if(j==0){
      continue;
    }
    alpha = ((double) j) / nPossibleSamples;
    if(verb>=4){
      Rcout << "Current acceptance rate: " << alpha << "\n";
    }
//      if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
    int addedNow = std::min(j, samplesRemaining);
    Y.submat(0,totalAccepted,(ncols-1),(totalAccepted+addedNow-1)) = X.cols(goodIndices.head(addedNow));
    totalAccepted += addedNow;
    samplesRemaining -= addedNow;
    }
  if(verb>=3){
    Rcout << "\nTotal samples run " << totalSamplesRun << "\nTotal samples accepted " << totalAccepted << "\nRatio: " << (((double) totalAccepted)/totalSamplesRun) << "\nLast alpha: " << alpha;
  }
  return Y;//Rcpp::List::create(Rcpp::Named("mat")= X, Rcpp::Named("ind")=goodIndices,Rcpp::Named("matY")= Y);
}


