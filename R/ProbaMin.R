##########
#' @title Probability of exceedance of minimum of Gaussian vector
#'
#'
#' @description Computes \eqn{P(min X \ge Thresh)}
#' with choice of algorithm between ANMC_Gauss and MC_Gauss.
#' The two most expensive parts are computed with the RCpp functions.
# [Version 1]
# INPUT
#' @param cBdg computational budget.
#' @param q number of active dimensions.
#' @param E discretization design for the field (somehow optional?).
#' @param Thresh threshold.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param pn coverage function vector.
#' @param lightReturn boolean, if true light return.
#' @param method method chosen to select the active dimensions.
#' @param verb level of verbosity (0-5), selects verbosity also for ANMC_Gauss (verb-1) and MC_Gauss (verb-1).
#' @param Algo choice of algorithm to compute the remainder Rq ("ANMC" or "MC").
#' @return A list containing
#' \itemize{
#'    \item{\code{probability}: }{The probability estimate}
#'    \item{\code{variance}: }{the variance of the probability estimate}
#' }
#' If \code{lightReturn=F} then the list also contains:
#' \itemize{
#'    \item{\code{aux_probabilities}: }{ a list with the probability estimates: \code{probability} the actual probability, \code{pq} the biased estimator \eqn{p_q}, \code{Rq} the conditional probability \eqn{R_q}}
#'    \item{\code{Eq}: }{the points of the design \eqn{E} selected for \eqn{p_q}}
#'    \item{\code{indQ}: }{the indices of the active dimensions chosen for \eqn{p_q}}
#'    \item{\code{resRq}: }{The list returned by the MC method used for \eqn{R_q}}
#'  }
# A list containing the probability estimate ($probabilities) and its variance ($variance), if lightReturn=FALSE it also includes the active dimensions ($indQ), the list returned by the MC estimator for Rq ($resRq).
#' @export
ProbaMin = function(cBdg,q,E,Thresh,mu,Sigma,pn=NULL,lightReturn=T,method=2,verb=0,Algo="ANMC"){

  # initialize parameters
  n<-length(mu)

  if(is.null(pn))
    pn<-pnorm((mu-Thresh)/sqrt(diag(Sigma)))

  if(!is.matrix(E))
    E<-as.matrix(E)

  # Select the appropriate q points
  indQ<-selectEq(q=q,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method)

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
    indQ<-selectEq(q=q,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method)
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
  pPrime<- 1 - pmvnorm(lower = Thresh,mean = muEq,sigma = KEq)
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
  #         - wwCondQ = ``weights'' for X^{-q} | X^q
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
    resMCQMC<-ANMC_Gauss(compBdg = cBdg,problem_Gauss,delta=0.45,type="m",typeReturn=2,verb=max(0,verb-1))
  }else{
    if(verb>=1){
      cat("Starting MC. \n")
    }
    resMCQMC<-MC_Gauss(compBdg = cBdg,problem = problem_Gauss,delta=0.2,type="m",typeReturn=2,verb=max(0,verb-1))
  }


  proba<- pPrime + resMCQMC$estim*(1-pPrime)

  varpPrime<-(attr(pPrime,"error")/3)^2
  vars<- (1-resMCQMC$estim)^2*varpPrime+resMCQMC$varEst*(1-pPrime)^2 +resMCQMC$varEst*varpPrime

  if(lightReturn){
    res<-list(probability=as.vector(proba), variance=vars)
  }else{
    res<- list(probability=as.vector(proba), variance=vars, aux_probabilities=list(probability=as.vector(proba),pq=pPrime,Rq=resMCQMC$estim),Eq=Eq,indQ=indQ,resRq=resMCQMC)
  }

  return(res)
}

