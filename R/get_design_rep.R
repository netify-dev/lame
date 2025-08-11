#' Create design array for replicate data
#' 
#' @usage get_design_rep(Y, Xdyad, Xrow, Xcol, actorSet, intercept, n, N, pr, pc, pd)
#' @param Y dependent variable in array format
#' @param Xdyad dyadic covariates in array format
#' @param Xrow sender covariates in array format
#' @param Xcol receiver covariates in array format
#' @param actorSet vector of actors
#' @param intercept logical indicating whether to include intercept
#' @param n number of actors
#' @param N number of replicates
#' @param pr number of receiver covariates
#' @param pc number of sender covariates
#' @param pd number of dyadic covariates
#' @return returns list of design array values necessary for ame_repL
#' @author Shahryar Minhas
#' 
#' @export get_design_rep

get_design_rep <- function(Y=Y, Xdyad, Xrow, Xcol, actorSet, intercept, n, N, pr, pc, pd){
  
  # construct design array
  X<-array(dim=c(n,n,pr+pc+pd+intercept,N)) 
  for (t in 1:N ){ 
    Xt<-design_array_listwisedel(Xrow[,,t],Xcol[,,t],Xdyad[,,,t],intercept,n) 
    
    # re-add intercept if it was removed
    if(dim(Xt)[3]<dim(X)[3] ){
      tmp<-array( dim=dim(Xt)+c(0,0,1),
                  dimnames=list(NULL,NULL,c('intercept',dimnames(Xt)[[3]])) ) 
      tmp[,,1]<-1 ; tmp[,,-1]<-Xt   
      Xt<-tmp
    }
    
    X[,,,t]<-Xt
  } 
  if( (pr+pc+pd+intercept)>0 ){
    dimnames(X) = list(actorSet,actorSet,dimnames(Xt)[[3]],dimnames(Y)[[3]])
  }
  dimnames(X)[[4]]<-dimnames(Y)[[3]]
  
  # add NAs from design array to Y
  if( (pr + pc + pd) > 0){
    for(t in 1:N){
      for(p in 1:dim(X)[3]){
        # find NAs
        naMat <- X[,,p,t] ; naMat[!is.na(naMat)] <- 1 ; naMat[is.nan(naMat)] <- NA
        # add NAs to Y
        Y[,,t] <- Y[,,t] * naMat
      }
    } ; rm(naMat)
  }
  
  # turn NAs in design array to zero
  X[is.na(X)]<-0
  
  # useful design array transformations to do up front
  ## these are used in calculation of rbeta
  Xlist <- lapply(1:N, function(t){ array(X[,,,t], dim=dim(X)[1:3], dimnames=dimnames(X)[1:3]) })
  XrLong <- apply(X, c(1,3,4), sum)                 # row sum
  XcLong <- apply(X, c(2,3,4), sum)                 # col sum
  mXLong <- apply(X, c(3,4), c)                     # design matrix
  mXtLong <- apply(aperm(X, c(2,1,3,4)), c(3,4), c) # dyad-transposed design matrix
  # regression sums of squares
  if( (pr+pc+pd+intercept)>0 ){
    xxLong <- array(apply(mXLong, 3, function(x){ t(x)%*%x }), dim=c(dim(X)[3],dim(X)[3],N))
    xxTLong <- array(unlist(lapply(1:N, function(t){ t(mXLong[,,t]) %*% mXtLong[,,t] })), dim=dim(xxLong))
  } else { xxLong<-NULL ; xxTLong<-NULL  }
  
  return(
    list(
      Y=Y, X=X, Xlist=Xlist, XrLong=XrLong, XcLong=XcLong,
      mXLong=mXLong, mXtLong=mXtLong, xxLong=xxLong, xxTLong=xxTLong
    )
  )
}