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
#' @param q number of active dimensions.
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
#' \item 5  samples with equal probabilities.
#' }
#' @return A vector of integers denoting the chosen active dimensions of the vector mu.
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
selectEq = function(q,E,Thresh,mu,Sigma,pn=NULL,method=1){
  n<-length(mu)
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
