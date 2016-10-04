## Functions to select the points Eq
# Methods
# 0 => selects by taking equally spaced indexes in mu
# 1 => samples from pn
# 2 => samples from pn*(1-pn)
# 3 => samples from pn adjusting for the distance (tries to explore all modes)
# 4 => samples from pn*(1-pn) adjusting for the distance (tries to explore all modes)
# 5 => samples with equal probabilities
#' @title Select active dimensions for small dimensional estimate
#'
#' @description The function \code{selectEq} selects the active dimensions for the computation of \eqn{p_q} with an heuristic method.
#'
#' @param q either the fixed number of active dimensions or the range where the number of active dimensions is chosen with \code{selectQdims}. If \code{NULL} the function \code{selectQdims} is called.
#' @param E discretization design for the field.
#' @param Thresh threshold.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param pn coverage function (if NULL it is computed).
#' @param method integer chosen between \itemize{
#' \item 0  selects by taking equally spaced indexes in mu;
#' \item 1  samples from pn;
#' \item 2  samples from pn*(1-pn);
#' \item 3  samples from pn adjusting for the distance (tries to explore all modes);
#' \item 4  samples from pn*(1-pn) adjusting for the distance (tries to explore all modes);
#' \item 5  samples with uniform probabilities.
#' }
#' @param verb level of verbosity: 0 returns nothing, 1 returns minimal info
#' @return A vector of integers denoting the chosen active dimensions of the vector mu.
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
selectEq = function(q=NULL,E,Thresh,mu,Sigma,pn=NULL,method=1,verb=0){
  n<-length(mu)

  # If q is NULL we don't know how many q to select thus we use the sequential procedure
  if(is.null(q)){
    return(selectQdims(E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method,reducedReturn=T,verb=verb))
  }else if(length(q)==2){ # here we pass the vector of limits instead of a fixed q
    return(selectQdims(E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method,reducedReturn=T,verb=verb,limits = q))
  }

  if(method==0){
    indQ<-as.integer(seq(from = 1,to = n,length.out = q))
  }else if(method==1){
    if(is.null(pn)){
      pn<-pnorm((mu-Thresh)/sqrt(diag(Sigma)))
    }
    indQ<-sample.int(n,q,prob = pn)
  }else if(method==2){
    if(is.null(pn)){
      pn<-pnorm((mu-Thresh)/sqrt(diag(Sigma)))
    }
    indQ<-sample.int(n,q,prob = pn*(1-pn))
  }else if(method==3){
    distances<-as.matrix(dist(E))
    if(is.null(pn)){
      pn<-pnorm((mu-Thresh)/sqrt(diag(Sigma)))
    }

    indQ<-rep(-1,q)
    indQ[1]<-sample.int(n,1,prob = pn)
    dd<-1
    for(i in (2:q)){
    #  plot(pn,type='l')
    #  plot(dd,type='l')
      dd<-dd^0.8*distances[indQ[i-1],]
      dd<-dd/diff(range(dd))
    #  plot(dd,type='l',col=2)
  #      image(matrix(dd*pn,nrow=30),col=grey.colors(20))
  #      contour(matrix(dd*pn,nrow=30),add=T,nlevels=12)
   #     points(E[indQ[1:(i-1)],],pch=16)

    #  plot(dd*pn,type='l')  #image(as.image(dd*pn,q = nd),col=col)
      indQ[i]<-sample.int(n,1,prob = dd*pn)

#      points(t(E[indQ[i],]),pch=16,col=2)
    #  points(indQ[1:i],rep(0,i),col=2,pch=16)
 #     i=i+1
    }

  }else if(method==4){
    distances<-as.matrix(dist(E))
    if(is.null(pn)){
      pn<-pnorm((mu-Thresh)/sqrt(diag(Sigma)))
    }

    indQ<-rep(-1,q)
    indQ[1]<-sample.int(n,1,prob = pn*(1-pn))
    dd<-1
    for(i in (2:q)){
      #  plot(pn*(1-pn),type='l')
      #  plot(dd,type='l')
      dd<-dd^0.8*distances[indQ[i-1],]
      dd<-dd/diff(range(dd))
      #  plot(dd,type='l',col=2)

      #  plot(dd*pn*(1-pn),type='l')  #image(as.image(dd*pn,q = nd),col=col)
      indQ[i]<-sample.int(n,1,prob = dd*pn*(1-pn))
      #  points(indQ[1:i],rep(0,i),col=2,pch=16)
    }
  }else if(method==5){
    indQ<-sample.int(n,q)
  }else {
    indQ<-as.integer(seq(from = 1,to = n,length.out = q))
  }
  indQ<-sort(indQ)
  return(indQ)
}



#' @title Iteratively select active dimensions
#'
#' @description The function \code{selectQdims} iteratively selects the number of active dimensions and the dimensions themselves for the computation of \eqn{p_q}.
#' We keep increasing the number of dimensions until \eqn{p_{q}-p_{q-1}} is smaller than the error of the procedure.
#'
#' @param E discretization design for the field.
#' @param Thresh threshold.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param pn coverage function (if NULL it is computed).
#' @param method integer chosen between \itemize{
#' \item 0  selects by taking equally spaced indexes in mu;
#' \item 1  samples from pn;
#' \item 2  samples from pn*(1-pn);
#' \item 3  samples from pn adjusting for the distance (tries to explore all modes);
#' \item 4  samples from pn*(1-pn) adjusting for the distance (tries to explore all modes);
#' \item 5  samples with uniform probabilities.
#' }
#' @param reducedReturn boolean to select the type of return. See Value for further details.
#' @param verb level of verbosity: 0 returns nothing, 1 returns minimal info.
#' @param limits numeric vector of length 2 with q_min and q_max. If \code{NULL} initialized at c(10,300)
#' @return If \code{reducedReturn=F} returns a list containing
#' \itemize{
#'    \item{\code{indQ}: }{the indices of the active dimensions chosen for \eqn{p_q};}
#'    \item{\code{pq}: }{the biased estimator \eqn{p_q} with attribute \code{error}, the estimated absolute error;}
#'    \item{\code{Eq}: }{the points of the design \eqn{E} selected for \eqn{p_q};}
#'    \item{\code{muEq}: }{the subvector of \code{mu} selected for \eqn{p_q};}
#'    \item{\code{KEq}: }{the submatrix of \code{Sigma} composed by the indexes selected for \eqn{p_q}.}
#' }
#' Otherwise it returns only \code{indQ}.
#'
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
selectQdims = function(E,Thresh,mu,Sigma,pn=NULL,method=1,reducedReturn=T,verb=0,limits=NULL){
  if(verb>0)
    cat("\n selectQdims: find good q and select active dimensions \n")
  n<-length(mu)
  if(!is.null(limits)){
    q0=max(2,min(limits[1],n))
    qMax=min(limits[2],n,300)
  }else{
    q0=min(10,n)
    qMax=min(n,300)
  }
  if(verb>0)
    cat("\n Initial q:",q0)
  # compute the first pPrime
  indQ<-selectEq(q=q0,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method)
  Eq<-E[indQ,]
  # compute muEq
  muEq<-mu[indQ]
  # compute k(Eq,Eq)
  KEq<-Sigma[indQ,indQ]
  # Compute p_q
  pPrime<- 1 - pmvnorm(upper = Thresh,mean = muEq,sigma = KEq)
  err<-attr(pPrime,"error")
  deltaP<-1
  # q increment
  qIncrement<-min(10,ceiling(0.01*n))
  if(verb>0)
    cat(", Increments q:",qIncrement,"\n")
  flag=0
  q=q0
  while(flag<2){
    q<-min(q+qIncrement,qMax)
    indQ<-selectEq(q=q,E=E,Thresh=Thresh,mu=mu,Sigma=Sigma,pn=pn,method=method)
    Eq<-E[indQ,]
    # compute muEq
    muEq<-mu[indQ]
    # compute k(Eq,Eq)
    KEq<-Sigma[indQ,indQ]
    # Compute p_q
    temp<-1 - pmvnorm(upper = Thresh,mean = muEq,sigma = KEq,algorithm = GenzBretz(abseps = 0.01))
#    pPrime<-c(pPrime,temp)
    err<-attr(temp,"error")
    deltaP<-abs(temp-pPrime)/(temp+1)
    if(deltaP<=err)
      flag<-flag+1
    if(q==qMax)
      flag<-flag+2
    pPrime<-temp
  }
  if(verb>0)
    cat("Final q:",q,", deltaP: ",deltaP,", err:",err,"\n")

  res<-indQ
  if(!reducedReturn){
    pq<-1 - pmvnorm(upper = Thresh,mean = muEq,sigma = KEq)
    attr(pq,"error")<-attr(temp,"error")
    res<-list(indQ=indQ,pq=pq,Eq=Eq,muEq=muEq,KEq=KEq)
  }

  return(res)
}
