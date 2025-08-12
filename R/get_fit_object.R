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
#' @param start_vals start_vals for future model run
#' @param symmetric logical indicating whether model is symmetric
#' @param tryErrorChecks list with counts of MCMC errors
#' @param model.name Name of the model (optional)
#' @param U Latent sender positions (optional, for dynamic UV)
#' @param V Latent receiver positions (optional, for dynamic UV)
#' @param dynamic_uv logical indicating whether UV effects are dynamic
#' @param dynamic_ab logical indicating whether additive effects are dynamic
#' @param bip logical indicating whether the network is bipartite
#' @param rho_ab temporal correlation parameter for additive effects (optional)
#' @param rho_uv temporal correlation parameter for multiplicative effects (optional)
#' @param family character string specifying the model family (e.g., "binary", "normal", "poisson")
#' @param odmax vector of maximum ranks for ordinal or fixed rank nomination families
#' @param nA number of actors in first mode (for bipartite networks)
#' @param nB number of actors in second mode (for bipartite networks)
#' @param n_time number of time periods (for longitudinal models)
#' @return Fitted AME object
#' @author Shahryar Minhas
#' @export get_fit_object
get_fit_object <- function(
    APS, BPS, UVPS, YPS, 
    BETA, VC, GOF,
    Xlist, actorByYr, start_vals,
    symmetric, tryErrorChecks,
    model.name=NULL,
    U=NULL, V=NULL, dynamic_uv=FALSE, dynamic_ab=FALSE,
    bip=FALSE,
    rho_ab=NULL, rho_uv=NULL,
    family=NULL, odmax=NULL, nA=NULL, nB=NULL, n_time=NULL
){
  
  # some labels and dims
  if(is.matrix(APS)) {
    actors <- rownames(APS)
  } else {
    actors <- names(APS)
  }
  # If no actor names, create default ones
  if(is.null(actors)) {
    if(is.matrix(APS)) {
      actors <- paste0("Actor", 1:nrow(APS))
    } else {
      actors <- paste0("Actor", 1:length(APS))
    }
  }
  pdLabs <- dimnames(YPS)[[3]]
  if(dynamic_uv && !is.null(U)) {
    if(length(dim(U)) == 3) {
      R <- dim(U)[2]
    } else {
      R <- ncol(U)
    }
  } else {
    R <- ncol(start_vals$U)
  }
  N <- dim(YPS)[3]
  
  # posterior means 
  if(dynamic_ab && is.matrix(APS)) {
    # For dynamic effects, APS and BPS are already averaged matrices
    APM <- APS
    BPM <- BPS
    rownames(APM) <- rownames(BPM) <- actors
    # Extract time-averaged values for EZ calculation
    APM_avg <- rowMeans(APM)
    BPM_avg <- rowMeans(BPM)
    names(APM_avg) <- names(BPM_avg) <- actors
  } else {
    APM<-APS/nrow(VC) 
    BPM<-BPS/nrow(VC)
    names(APM) <- names(BPM) <- actors
    APM_avg <- APM
    BPM_avg <- BPM
  }
  YPM<-YPS/nrow(VC)
  # Set dimension names for YPM
  if(!bip) {
    dimnames(YPM) <- list(actors, actors, pdLabs)
  } else {
    # For bipartite, YPM is rectangular
    # Need to handle row and column actors separately
    # For now, skip dimnames for bipartite
  }
  
  # Handle UVPM based on dynamic or static
  if(!bip) {
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
  } else {
    # For bipartite, skip EZ computation for now
    UVPM <- UVPS
    EZ <- NULL
  }
  
  # adding names
  if(!bip) {
    if(dynamic_uv && length(dim(UVPM)) == 3) {
      dimnames(UVPM) <- list(actors, actors, pdLabs)
    } else if(length(dim(UVPM)) == 2) {
      rownames(UVPM)<-colnames(UVPM)<-actors
    }
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
      if(!bip) {
        UDV<-svd(UVPM)
        U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
        V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
        rownames(U)<-rownames(V)<-actors
      }
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
  if(!bip) {
    EZ <- array_to_list(EZ, actorByYr, pdLabs)
    YPM <- array_to_list(YPM, actorByYr, pdLabs)
  } else {
    # For bipartite, skip list conversion for now
    # EZ is already NULL, YPM needs special handling
  }
  
  # Reformat GOF from 3D array to named list if needed
  if(is.array(GOF) && length(dim(GOF)) == 3) {
    # GOF is [5 stats, N time points, nscan/odens+1 iterations]
    gof_names <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    GOF_list <- list()
    for(i in 1:5) {
      GOF_list[[gof_names[i]]] <- GOF[i,,]
    }
    GOF <- GOF_list
  }
  
  # create fitted object
  if(symmetric){
    fit <- list(
      BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
      YPM=YPM,GOF=GOF, start_vals=start_vals, tryErrorChecks=tryErrorChecks,
      model.name=model.name, family=family, symmetric=symmetric, odmax=odmax,
      mode=if(bip) "bipartite" else "unipartite")
    # Add dynamic fields if applicable
    if(dynamic_ab && is.matrix(APM)) {
      fit$a_dynamic <- APM
      fit$APM <- APM_avg  # Overwrite with time-averaged for compatibility
    }
    # Add rho values
    if(!is.null(rho_ab)) fit$rho_ab <- rho_ab
    if(!is.null(rho_uv)) fit$rho_uv <- rho_uv
  }
  if(!symmetric){
    fit <- list(
      BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
      YPM=YPM,GOF=GOF, start_vals=start_vals, tryErrorChecks=tryErrorChecks,
      model.name=model.name, family=family, symmetric=symmetric, odmax=odmax,
      mode=if(bip) "bipartite" else "unipartite")
    # Add dynamic fields if applicable
    if(dynamic_ab && is.matrix(APM)) {
      fit$a_dynamic <- APM
      fit$b_dynamic <- BPM
      fit$APM <- APM_avg  # Overwrite with time-averaged for compatibility
      fit$BPM <- BPM_avg
    }
    # Add rho values
    if(!is.null(rho_ab)) fit$rho_ab <- rho_ab
    if(!is.null(rho_uv)) fit$rho_uv <- rho_uv
  }
  # Add bipartite dimensions if applicable
  if(bip) {
    fit$nA <- nA
    fit$nB <- nB
  }
  # Add longitudinal info if applicable
  if(!is.null(n_time)) fit$n_time <- n_time
  # Add dynamic flags
  fit$dynamic_uv <- dynamic_uv
  fit$dynamic_ab <- dynamic_ab
  # Add Xlist if present
  if(!is.null(Xlist)) fit$Xlist <- Xlist
  class(fit)<-"ame"
  return(fit)
}