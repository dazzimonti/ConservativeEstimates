# Version 0.3: optimized version for the Gaussian case
# changed minimum m0 and computation of nStar (now cBdg - time actually passed)
#' @title ANMC estimate for the remainder
#' @description Asymmetric nested Monte Carlo estimation of \eqn{P(max X^{-q} > Thresh | max X^{q} \le Thresh)} where X is a normal vector.
# Input:
#' @param compBdg  total computational budget in seconds.
#' @param problem list defining the problem with mandatory fields \itemize{
#'         \item muEq = mean vector of \eqn{X^{q}};
#'         \item sigmaEq = covariance matrix of \eqn{X^q};
#'         \item Thresh = threshold;
#'         \item muEmq = mean vector of \eqn{X^{-q}};
#'         \item wwCondQ = ``weights'' for \eqn{X^{-q} | X^q} [the vector \eqn{\Sigma^{-q,q}(\Sigma^q)^{-1}}];
#'         \item sigmaCondQChol = Cholesky factorization of the conditional covariance matrix \eqn{\Sigma^{-q | q}};
#'         }
#' @param delta total proportion of budget assigned to initial estimate (default 0.4), the actual proportion used might be smaller.
#' @param type type of excursion: "m", for minimum below threshold or "M", for maximum above threshold.
#' @param typeReturn integer chosen between \itemize{
#'          \item 0 a number with only the probability estimation;
#'          \item 1 light return: a list with the probability estimator, the variance of the estimator, the vectors of conditional quantities used to obtain m^* and the system dependent parameters;
#'          \item 2 heavy return: the same list as light return with also the computational times and additional intermediate parameters.
#' }
#' @param verb level of verbosity (0,1 for this function), also sets the verbosity of trmvrnorm_rej_cpp (to verb-1).
#'
#' @return A list containing the estimated probability of excursion, see typeReturn for details.
#' @references Azzimonti, D. and Ginsbourger, D. (2016). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Dickmann, F. and Schweizer, N. (2014). Faster comparison of stopping times by nested conditional Monte Carlo. arXiv preprint arXiv:1402.0243.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
ANMC_Gauss<-function(compBdg,problem,delta=0.4,type="M",typeReturn=0,verb=0){
  sizeX<-length(problem$muEq)
  sizeY<-length(problem$muEmq)

  if(type=="M"){
    upperTmvn = rep(problem$Thresh,sizeX)
    lowerTmvn = rep(-Inf,sizeX)
    gg = function(x){return(max(x)>problem$Thresh)}
  }else{
    upperTmvn = rep(Inf,sizeX)
    lowerTmvn = rep(problem$Thresh,sizeX)
    gg = function(x){return(min(x)<problem$Thresh)}
  }

  # Step 0: compute the cross-covariance and the necessary Cholesky decompositions
  # cholKeq<-chol(problem$sigmaEq)
  if(verb>0)
    cat("Starting ANMC... \n")

  # Step 1: estimate Cx0, beta0 to get
  # (possibly) reasonable n0,m0

  # estimate Cx
  timeInPart1<-get_nanotime()
  simsX<-trmvrnorm_rej_cpp(n = 1,mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1))
  time1SimX<-(get_nanotime()-timeInPart1)*1e-9

  nTests<-floor(compBdg*delta*0.4/time1SimX)

  ttX<-rep(0,nTests)
  ii<-floor(seq(from=1,to=120,by=120/nTests))
  for(i in seq(nTests)){
    timeIn<-get_nanotime()
    temp<-trmvrnorm_rej_cpp(n = ii[i],mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1))
    ttX[i]<-(get_nanotime()-timeIn)*1e-9*1.03
    simsX<-cbind(simsX,temp)
  }
  Cx0<-unname(lm(ttX~ii)$coefficients[1])
  Cx<-unname(lm(ttX~ii)$coefficients[2])


  # estimate alpha and beta with means
  timeIn<-get_nanotime()
  muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsX[,1]-problem$muEq)
  alpha<-(get_nanotime()-timeIn)*1e-9
  simsYcX<- mvrnormArma(n=1,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
  time1SimYcX<-(get_nanotime()-timeIn)*1e-9

  nBetas<-ceiling(compBdg*delta*0.1/(time1SimYcX-alpha))

  timeIn<-get_nanotime()
  simsYcX<- mvrnormArma(n=nBetas,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
  beta0<-(get_nanotime()-timeIn)*1e-9/nBetas

  # estimate time to evaluate function
  timeG<-microbenchmark(gg(1))
  tEvalG<-quantile(timeG$time,probs = 0.99,names = F)*1e-9

  C_adj<-compBdg*delta -time1SimX - sum(ttX) -alpha-beta0*nBetas -time1SimYcX - sum(timeG$time)*1e-9

  if(verb>1){
    cat("Time passed:",compBdg*delta-C_adj,"\n",
        "tEvalG:",tEvalG,"\n",
        "beta0: ",beta0,"\n",
        "Cx0: ",Cx0,"\n",
        "alpha: ",alpha,"\n")
  }


  n0<-ncol(simsX)

  m0<-max(30,floor(((C_adj-Cx0)/n0-Cx-alpha)/(beta0+tEvalG)))

  if(verb>1){
    cat("m0:",m0,"\n",
        "n0:",n0,"\n")
  }


  # generate simsX
  # simsX<-problem$simulatorX(n0)
  #cat("n0: ",n0,", m0: ",m0,"simsX: ",str(simsX),"\n")
  # create conditional covariance
  #  timeIn<-proc.time()
  simsYcondX<-matrix(0,nrow = sizeY,ncol = n0*m0)
  if(typeReturn>=2){ # not really needed, only for debug
    hatHatG0<-0
  }
  varCondX0<-rep(0,n0)
  expYcondX0<-rep(0,n0)

  gEval0<-rep(0,n0*m0)

  for(j in seq(n0)){

    muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsX[,j]-problem$muEq)
    simsYcondX[,((1:m0)+m0*(j-1))]<-mvrnormArma(n=m0,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)

    for(i in seq(m0)){
      gEval0[(i+m0*(j-1))]<-gg(simsYcondX[,(i+m0*(j-1))])
      #      expYcondX0[j]<-expYcondX0[j]+gEval[(i+m0*(j-1))]
    }
    expYcondX0[j]<-mean(gEval0[(1+m0*(j-1)):(m0*j)])
    #    for(i in seq(m0)){
    #  varCondX0[j]<-varCondX0[j] + (gEval[(i+m0*(j-1))]-expYcondX0[j])^2
    #}
    varCondX0[j]<- sum((gEval0[(1+m0*(j-1)):(m0*j)]-expYcondX0[j])^2)/(m0-1)
  }


  if(typeReturn>=2){
    hatHatG0<-mean(expYcondX0) # again not really needed
  }

  # Step 2: compute mStar
  ratioAAmB<-mean(apply(X = matrix(1:m0),MARGIN = 1,FUN = function(nn){var(gEval0[seq(nn,by=m0,length.out = n0)])*var(expYcondX0)/(cov(gEval0[seq(nn,by=m0,length.out = n0)],expYcondX0))^2-1}),na.rm = TRUE)
  if(!is.finite(ratioAAmB)){
    ratioAAmB<-0.25 # upper bound
  }
  mStar<-sqrt(((Cx0+Cx+3*alpha)*ratioAAmB)/(beta0+tEvalG))
  if(!is.finite(mStar)){
    if(verb>5){
      cat("mStar:",mStar,", replaced it with",m0,"\n")
    }
    mStar<-m0
  }
  if(verb>0){
    cat("Obtained mStar:",mStar,", ")
  }
  timePart1<-(get_nanotime()-timeInPart1)*1e-9

  # round it to the nearest integer
  epsMstar<-mStar-floor(mStar)
  if(epsMstar<0.5*(2*mStar+1 - sqrt(4*mStar^2+1))){
    mStar<-max(floor(mStar),2,na.rm = T)
  }else{
    mStar<-max(ceiling(mStar),2,na.rm = T)
  }

  # derive nStar
  nStar<-round((compBdg-timePart1*0.9-Cx0)/(Cx+alpha+(beta0)*mStar))

  if(verb>1){
    cat("Cx0: ",Cx0,", beta0: ",beta0, ", mean(varYcX0): ", mean(varCondX0), "var(expYx0): ",var(expYcondX0),", mStar: ",mStar,"nStar: ",nStar, "Time Part 1: ",timePart1,"(compBdg assigned: ",compBdg*delta,")\n")
  }
  nStar<-max(nStar,n0)
  # Step 3
  # we already have n0 sims for X and m0 sims of Y|X for each x sim
  # so we only need nStar-n0 simulations for X and mStar simulations
  # of Y conditional on the nStar-n0 xs, and mStar-m0 for the n0 Xs
  simsYcondXfull<-matrix(0,nrow = sizeY,ncol = nStar*mStar)

  # generate the missing nStar -n0 simulations of X
  # re-estimate also Cx for debug reasons
  #  timeIn<-proc.time()
  if(nStar>n0){
    simsXfull<-cbind(simsX,trmvrnorm_rej_cpp(n = (nStar-n0),mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1)))
  }else{
    simsXfull<-simsX
  }
  #  Cx<-max(unname((proc.time()-timeIn)[1]),0.001)/(nStar-n0)

  # generate all sims Y|X preserving the ones we already have.
  # compute estim and related quantities

  estim<-0
  expYcondXfull<-rep(0,nStar)
  if(typeReturn>=1){
    varCondXfull<-rep(0,nStar)
  }
  indM<-min(m0,mStar)

  ########
  # you need to save gEval in a vector and fit the ones you have already done in it!
  ########
  gEval<-rep(0,nStar*mStar)
  for(j in seq(nStar)){
    #    muYcondX<-problem$muY+wwYcondX%*%(simsXfull[,j]-problem$muX)

    if(j<=n0){  # we have already simulated X
      simsYcondXfull[,(1:indM+mStar*(j-1))]<-simsYcondX[,(1:indM+m0*(j-1))] # all i,j s.t. j<=n0 and i<=m0
      gEval[(1:indM+mStar*(j-1))]<-gEval0[(1:indM+m0*(j-1))]


      # here indM=mStar>m0 we need to simulate Y|X
      if(indM>m0){
        #        timeIn<-proc.time()
        muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsXfull[,j]-problem$muEq)
        simsYcondXfull[,((m0+1):mStar+mStar*(j-1))]<-mvrnormArma(n=(mStar-m0),mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
        for(i in seq(mStar)){
          # now we have the simulations, we can compute all estimates
          gEval[i+mStar*(j-1)]<-gg(simsYcondXfull[,(i+mStar*(j-1))])
        }
      }
    }else{ # here we haven't simulated X before so we need all Y|X
      #      timeIn<-proc.time()
      muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsXfull[,j]-problem$muEq)
      simsYcondXfull[,(1:mStar+mStar*(j-1))]<-mvrnormArma(n=mStar,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
      for(i in seq(mStar)){
        # now we have the simulations, we can compute all estimates
        gEval[i+mStar*(j-1)]<-gg(simsYcondXfull[,(i+mStar*(j-1))])
      }
    }


    #    expYcondXfull[j]<-expYcondXfull[j]/mStar
    expYcondXfull[j]<-sum(gEval[(1+mStar*(j-1)):(mStar*j)])/mStar
    if(typeReturn>=1){
      #      for(i in seq(mStar)){
      #        varCondXfull[j]<-varCondXfull[j] + (gEval[i+mStar*(j-1)]-expYcondXfull[j])^2
      #      }
      if(mStar>1){
        varCondXfull[j]<- sum((gEval[(1+mStar*(j-1)):(mStar*j)]-expYcondXfull[j])^2)/(mStar-1)
      }else{
        varCondXfull[j]<-0
      }
    }
  }
  #  betaFull<-max(mean(unname(betaFull)),0.001)
  estim<-mean(expYcondXfull)
  timeTot<-(get_nanotime()-timeInPart1)*1e-9

  if(verb>=1){
    cat("Total time: ",timeTot,"(compBdg: ",compBdg,")\n")
  }

  if(typeReturn==0){
    return(estim)
  }else{
    varG<-var(expYcondXfull)+mean(varCondXfull)
    varHatHatG<-varG/nStar-(mStar-1)/(nStar*mStar)*mean(varCondXfull)
    params<-list(n=nStar,m=mStar,Cx=Cx,Cx0=Cx0,alpha=alpha,beta=beta0,evalG=tEvalG)
    results<-list(estim=estim,varEst=varHatHatG,expYcondX=expYcondXfull,varYcondX=varCondXfull,params=params)
    if(typeReturn==1){
      return(results)
    }else{
      results$est1=mean(gEval)
      varG0<-var(expYcondX0)+mean(varCondX0)
      varEst0<-varG0/n0-(m0-1)/(n0*m0)*mean(varCondX0)
      results$params0<-list(n0=n0,m0=m0)
      results$values0<-list(estim0=hatHatG0,varEst0=varEst0,expYcondX0=expYcondX0,varYcondX0=varCondX0,ratio0=ratioAAmB)
      results$times<-list(part1=timePart1,total=timeTot)
      return(results)
    }
  }
  return("fail")
}
