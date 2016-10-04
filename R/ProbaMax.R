##########
#' @title Probability of exceedance of maximum of Gaussian vector
#'
#'
#' @description Computes \eqn{P(max X > Thresh)}
#' with choice of algorithm between ANMC_Gauss and MC_Gauss.
#' The two most expensive parts are computed with the RCpp functions.
# [Version 1]
# INPUT
#' @param cBdg computational budget.
#' @param Thresh threshold.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param E discretization design for the field. If \code{NULL}, a simplex-lattice design {n,n} is used, with n=length(mu). In this case the choice of method=4,5 are not advised.
#' @param q number of active dimensions. Can be passed either as an integer or as numeric vector of length 2. The vector is the range where to search for the best number of active dimensions. If \code{NULL} q is selected as the best number of active dimensions in the feasible range.
#' @param pn coverage function vector.
#' @param lightReturn boolean, if true light return.
#' @param method method chosen to select the active dimensions.
#' @param verb level of verbosity (0-5), selects verbosity also for ANMC_Gauss (verb-1) and MC_Gauss (verb-1).
#' @param Algo choice of algorithm to compute the remainder Rq ("ANMC" or "MC").
#' @return A list containing
#' \itemize{
#'    \item{\code{probability}: }{The probability estimate}
#'    \item{\code{variance}: }{the variance of the probability estimate}
#'    \item{\code{q}:}{the number of selected active dimensions}
#' }
#' If \code{lightReturn=F} then the list also contains:
#' \itemize{
#'    \item{\code{aux_probabilities}: }{ a list with the probability estimates: \code{probability} the actual probability, \code{pq} the biased estimator \eqn{p_q}, \code{Rq} the conditional probability \eqn{R_q}}
#'    \item{\code{Eq}: }{the points of the design \eqn{E} selected for \eqn{p_q}}
#'    \item{\code{indQ}: }{the indices of the active dimensions chosen for \eqn{p_q}}
#'    \item{\code{resRq}: }{The list returned by the MC method used for \eqn{R_q}}
#' }
# the probability estimate ($probabilities) and its variance ($variance), if lightReturn=FALSE it also includes the active dimensions ($indQ), the list returned by the MC estimator for Rq ($resRq).
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Dickmann, F. and Schweizer, N. (2014). Faster comparison of stopping times by nested conditional Monte Carlo. arXiv preprint arXiv:1402.0243.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
ProbaMax = function(cBdg,Thresh,mu,Sigma,E=NULL,q=NULL,pn=NULL,lightReturn=T,method=4,verb=0,Algo="ANMC"){

  # initialize parameters
  n<-length(mu)

  if(is.null(pn))
    pn<-pnorm((Thresh-mu)/sqrt(diag(Sigma)))

  if(is.null(E)){
    E<-seq(0,1,length.out = length(mu))
  }

  if(!is.matrix(E))
    E<-as.matrix(E)

  # Select q and the active dimensions (if q is a number we select q active dims)
  # [pn = P(x in Gamma), interesting points here have high prob of not being in Gamma]
  indQ<-selectEq(q=q,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=(1-pn),method=method)

  # if q was given as a range here we reinitialize it as the number of active dims
  q<-length(indQ)

  if(verb>=1){
    cat("Computed indQ points, q = ",q,"\n")
  }

  Eq<-E[indQ,]

  # compute muEq
  muEq<-mu[indQ]

  # compute k(Eq,Eq)
  KEq<-Sigma[indQ,indQ]

  while(attr(chol(KEq,pivot=T),"rank")!=q){
    # ReSelect the appropriate q points
    indQ<-selectEq(q=q,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=(1-pn),method=method)
    #  Eq<-E[indQ]
    Eq<-E[indQ,]
    # compute muEq
    muEq<-mu[indQ]
    # compute k(Eq,Eq)
    KEq<-Sigma[indQ,indQ]
    if(verb>=2){
      cat("Covariance matrix non p.d., re-Computed indQ points\n")
    }
  }

  cholKeq<-chol(KEq)

  # Compute p'
  pPrime<- 1 - pmvnorm(upper = Thresh,mean = muEq,sigma = KEq)
  #  sSize<-N
  if(verb>=1){
    cat("Computed pPrime = ",pPrime,"\n")
  }
  if((1-pPrime)<attr(pPrime,"error")){
    if(verb>=2){
      cat("pPrime close to 1: pPrime=",pPrime,", error=",attr(pPrime,"error"),"\n")
    }
    if(lightReturn){
      res<-list(probabilities=list(probability=as.vector(pPrime),pPrime=pPrime,conditional=0),variance=(2/7*attr(pPrime,"error"))^2)
    }else{
      res<-list(probabilities=list(probability=as.vector(pPrime),pPrime=pPrime,conditional=0),variance=(2/7*attr(pPrime,"error"))^2,Uncond=NULL,Eq=Eq,indQ=indQ)
    }
    if(verb>=2){
      cat("Early return. \n")
    }
    return(res)
  }

  # problem = list defining the problem: with mandatory fields
  #         - muEq = mean vector of X^{q}
  #         - sigmaEq = covariance matrix of X^q
  #         - Thresh = threshold
  #         - muEmq = mean vector of X^{-q}
  #         - wwCondQ = ``weights'' for X^{-q} | X^q [Sigma^{−q,q}(Sigma^q)^{−1}]
  #         - sigmaCondQChol = Cholesky factorization of the conditional covariance matrix X^{-q} | X^q

  problem_Gauss<-list(muEq=muEq,sigmaEq=KEq,Thresh=Thresh,muEmq=mu[-indQ])
  invXX<-chol2inv(cholKeq)
  sigmaXY<-Sigma[indQ,-indQ]
  problem_Gauss$wwCondQ<-crossprod(sigmaXY,invXX)
  sigmaYcondX<-Sigma[-indQ,-indQ]-problem_Gauss$wwCondQ%*%sigmaXY
  problem_Gauss$sigmaCondQChol<-chol(sigmaYcondX)


  if(Algo=="ANMC"){
    if(verb>=1){
      cat("Starting ANMC. \n")
    }
    resMCQMC<-ANMC_Gauss(compBdg = cBdg,problem_Gauss,delta=0.45,typeReturn=2,verb=max(0,verb-1))
  }else{
    if(verb>=1){
      cat("Starting MC. \n")
    }
    resMCQMC<-MC_Gauss(compBdg = cBdg,problem_Gauss,delta=0.2,typeReturn=2,verb=max(0,verb-1))
  }


  proba<- pPrime + resMCQMC$estim*(1-pPrime)

  varpPrime<-(attr(pPrime,"error")/3)^2
  vars<- (1-resMCQMC$estim)^2*varpPrime+resMCQMC$varEst*(1-pPrime)^2 +resMCQMC$varEst*varpPrime

  if(lightReturn){
    res<-list(probability=as.vector(proba), variance=vars,q=q)
  }else{
    res<-list(probability=as.vector(proba), variance=vars, aux_probabilities=list(probability=as.vector(proba),pq=pPrime,Rq=resMCQMC$estim),Eq=Eq,indQ=indQ,resRq=resMCQMC)
  }

  return(res)
}

