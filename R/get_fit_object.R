#' Get fitted object from MCMC results
#' 
#' @param APS summed additive sender random effects
#' @param BPS summed additive sender random effects
#' @param UVPS summed multiplicative random effects
#' @param YPS summed Y posterior predictive values
#' @param BETA Matrix of draws for regression coefficient estimates
#' @param VC Matrix of draws for variance estimates
#' @param GOF Matrix of draws for goodness of fit calculations
#' @param Xlist List based version of design array
#' @param actorByYr List of actors by time point
#' @param startVals startVals for future model run
#' @param symmetric logical indicating whether model is symmetric
#' @param tryErrorChecks list with counts of MCMC errors
#' @param AIC Akaike Information Criterion (optional)
#' @param BIC Bayesian Information Criterion (optional)
#' @param model.name Name of the model (optional)
#' @return Fitted AME object
#' @author Peter Hoff, Shahryar Minhas
#' @export get_fit_object
get_fit_object <- function(
    APS, BPS, UVPS, YPS, 
    BETA, VC, GOF,
    Xlist, actorByYr, startVals,
    symmetric, tryErrorChecks,
    AIC=NA, BIC=NA, model.name=NULL
){
  
  # some labels and dims
  actors <- names(APS)
  pdLabs <- dimnames(YPS)[[3]]
  R <- ncol(startVals$U)
  N <- dim(YPS)[3]
  
  # posterior means 
  APM<-APS/nrow(VC) ; BPM<-BPS/nrow(VC)  
  UVPM<-UVPS/nrow(VC) ; YPM<-YPS/nrow(VC)
  EZ<-get_EZ_cpp(Xlist, 
                 apply(BETA,2,mean), outer(APM,BPM,"+"), 
                 UVPM, diag(nrow(UVPM)))
  
  # adding names
  rownames(UVPM)<-colnames(UVPM)<-actors
  dimnames(EZ)<-list(actors,actors,pdLabs)
  rownames(BETA)<-NULL
  
  # asymmetric uv
  if(!symmetric){
    UDV<-svd(UVPM)
    U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
    rownames(U)<-rownames(V)<-actors
  }
  
  # symmetric ul
  if(symmetric){
    ULUPM<-UVPM 
    eULU<-eigen(ULUPM) 
    eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
    U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
    L<-eULU$val[eR]   
    rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-actors
    for(t in 1:N){ 
      EZ[,,t]<-.5*(EZ[,,t]+t(EZ[,,t]))
      YPM[,,t]<-.5*(YPM[,,t]+t(YPM[,,t]))
    }
  }
  
  # reformat EZ and YPM as list objects
  EZ <- array_to_list(EZ, actorByYr, pdLabs)
  YPM <- array_to_list(YPM, actorByYr, pdLabs)
  
  # create fitted object
  if(symmetric){
    fit <- list(
      BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
      YPM=YPM,GOF=GOF, startVals=startVals, tryErrorChecks=tryErrorChecks,
      AIC=AIC, BIC=BIC, model.name=model.name)
  }
  if(!symmetric){
    fit <- list(
      BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
      YPM=YPM,GOF=GOF, startVals=startVals, tryErrorChecks=tryErrorChecks,
      AIC=AIC, BIC=BIC, model.name=model.name)
  }
  class(fit)<-"ame"
  return(fit)
}