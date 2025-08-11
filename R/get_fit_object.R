#' Get fitted object from MCMC results
#' 
#' @param APS summed additive sender random effects (or matrix for dynamic)
#' @param BPS summed additive receiver random effects (or matrix for dynamic)
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
#' @param U Latent sender positions (optional, for dynamic UV)
#' @param V Latent receiver positions (optional, for dynamic UV)
#' @param dynamic_uv logical indicating whether UV effects are dynamic
#' @param dynamic_ab logical indicating whether additive effects are dynamic
#' @return Fitted AME object
#' @author Peter Hoff, Shahryar Minhas
#' @export get_fit_object
get_fit_object <- function(
    APS, BPS, UVPS, YPS, 
    BETA, VC, GOF,
    Xlist, actorByYr, startVals,
    symmetric, tryErrorChecks,
    AIC=NA, BIC=NA, model.name=NULL,
    U=NULL, V=NULL, dynamic_uv=FALSE, dynamic_ab=FALSE
){
  
  # some labels and dims
  if(is.matrix(APS)) {
    actors <- rownames(APS)
  } else {
    actors <- names(APS)
  }
  pdLabs <- dimnames(YPS)[[3]]
  if(dynamic_uv && !is.null(U)) {
    if(length(dim(U)) == 3) {
      R <- dim(U)[2]
    } else {
      R <- ncol(U)
    }
  } else {
    R <- ncol(startVals$U)
  }
  N <- dim(YPS)[3]
  
  # posterior means 
  if(dynamic_ab && is.matrix(APS)) {
    # For dynamic effects, APS and BPS are already averaged matrices
    APM <- APS
    BPM <- BPS
    # Extract time-averaged values for EZ calculation
    APM_avg <- rowMeans(APM)
    BPM_avg <- rowMeans(BPM)
  } else {
    APM<-APS/nrow(VC) 
    BPM<-BPS/nrow(VC)
    APM_avg <- APM
    BPM_avg <- BPM
  }
  YPM<-YPS/nrow(VC)
  # Set dimension names for YPM
  dimnames(YPM) <- list(actors, actors, pdLabs)
  
  # Handle UVPM based on dynamic or static
  if(dynamic_uv && length(dim(UVPS)) == 3) {
    UVPM <- UVPS  # Already averaged, keep as 3D
    # For EZ calculation, use time-averaged UV
    UV_avg <- apply(UVPS, c(1,2), mean)
    EZ<-get_EZ_cpp(Xlist, 
                   apply(BETA,2,mean), outer(APM_avg,BPM_avg,"+"), 
                   UV_avg, diag(nrow(UV_avg)))
  } else {
    UVPM<-UVPS
    if(length(dim(UVPM)) == 2) {
      EZ<-get_EZ_cpp(Xlist, 
                     apply(BETA,2,mean), outer(APM_avg,BPM_avg,"+"), 
                     UVPM, diag(nrow(UVPM)))
    }
  }
  
  # Set dimension names for EZ
  dimnames(EZ) <- list(actors, actors, pdLabs)
  
  # adding names
  if(dynamic_uv && length(dim(UVPM)) == 3) {
    dimnames(UVPM) <- list(actors, actors, pdLabs)
  } else if(length(dim(UVPM)) == 2) {
    rownames(UVPM)<-colnames(UVPM)<-actors
  }
  rownames(BETA)<-NULL
  
  # asymmetric uv
  if(!symmetric){
    if(dynamic_uv && !is.null(U) && !is.null(V)) {
      # Use provided dynamic U and V
      if(length(dim(U)) == 3) {
        dimnames(U) <- list(actors, NULL, pdLabs)
        dimnames(V) <- list(actors, NULL, pdLabs)
      } else {
        rownames(U)<-rownames(V)<-actors
      }
    } else {
      # Standard SVD for static case
      UDV<-svd(UVPM)
      U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      rownames(U)<-rownames(V)<-actors
    }
  }
  
  # symmetric ul
  if(symmetric){
    if(dynamic_uv && !is.null(U)) {
      # Use provided dynamic U
      if(length(dim(U)) == 3) {
        dimnames(U) <- list(actors, NULL, pdLabs)
        ULUPM <- UVPM  # Keep 3D structure
        L <- NULL  # L not meaningful for dynamic case
      } else {
        rownames(U) <- actors
        ULUPM <- UVPM
        L <- NULL
      }
    } else {
      ULUPM<-UVPM 
      eULU<-eigen(ULUPM) 
      eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
      U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
      L<-eULU$val[eR]   
      rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-actors
    }
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
    # Add dynamic fields if applicable
    if(dynamic_ab && is.matrix(APM)) {
      fit$a_dynamic <- APM
      fit$APM <- APM_avg  # Overwrite with time-averaged for compatibility
    }
  }
  if(!symmetric){
    fit <- list(
      BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
      YPM=YPM,GOF=GOF, startVals=startVals, tryErrorChecks=tryErrorChecks,
      AIC=AIC, BIC=BIC, model.name=model.name)
    # Add dynamic fields if applicable
    if(dynamic_ab && is.matrix(APM)) {
      fit$a_dynamic <- APM
      fit$b_dynamic <- BPM
      fit$APM <- APM_avg  # Overwrite with time-averaged for compatibility
      fit$BPM <- BPM_avg
    }
  }
  class(fit)<-"ame"
  return(fit)
}