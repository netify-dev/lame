#' Get fitted object from MCMC results
#' 
#' @param startVals List object that is null or contains 
#' starting values
#' @param Y dependent variable in array format
#' @param family character vector (e.g. 'bin', 'nrm') specifying
#' family type
#' @param xP number of exogenous covariates
#' @param rvar logical indicating whether to include sender random
#' effects
#' @param cvar logical indicating whether to include receiver random
#' effects
#' @param R Number of dimensions for multiplicative effects
#' @param odmax vector of maximum ranks for cbin/frn families (optional)
#' @return List of starting values for MCMC
#' @author Peter Hoff, Shahryar Minhas
#' @export get_start_vals

get_start_vals <- function(startVals, Y, family, xP, rvar, cvar, R, odmax = NULL){
  
  # dims
  N <- dim(Y)[3]
  
  if(is.null(startVals)){
    # starting Z values
    Z<-array(dim=dim(Y))
    for (t in 1:N)
    {
      if(family=="normal"){Z[,,t]<-Y[,,t] }
      if(family=="tobit"){Z[,,t]<-Y[,,t] ; Z[,,t][Y[,,t]==0]<-min(Y[,,t][Y[,,t]>0],na.rm=TRUE)/2 }
      if(family=="ordinal"){Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),ncol(Y[,,t]))} 
      if(family=="rrl")
      {  
        Z[,,t]<-matrix(t(apply(Y[,,t],1,zscores)),nrow(Y[,,t]),ncol(Y[,,t])) 
      }  
      if(family=="binary" ){
        Z[,,t]<-matrix(zscores(Y[,,t]),nrow(Y[,,t]),nrow(Y[,,t]))
        # zyMax <- max(Z[,,t][Y[,,t]==0],na.rm=TRUE)
        zyMax <- if (sum(Y[,,t] == 0, na.rm = TRUE) > 0) {
          max(Z[,,t][Y[,,t] == 0], na.rm = TRUE)
        } else {
          min(Z[,,t][Y[,,t] == 1], na.rm = TRUE) - 1e-6
        }
        # zyMin <- min(Z[,,t][Y[,,t]==1],na.rm=TRUE)  
        zyMin <- if (sum(Y[,,t] == 1, na.rm = TRUE) > 0) {
          min(Z[,,t][Y[,,t] == 1], na.rm = TRUE)
        } else {
          max(Z[,,t][Y[,,t] == 0], na.rm = TRUE) + 1e-6
        }
        z01<-.5*(zyMax+zyMin ) 
        Z[,,t]<-Z[,,t] - z01
      } 
      if(family=="poisson"){
        Z[,,t]<-log(Y[,,t]+1) 
        diag(Z[,,t])<-0
      }
      
      if(is.element(family,c("cbin","frn")))
      {
        if(is.null(odmax)) {
          stop("odmax must be provided for family=cbin/frn")
        }
        Z[,,t]<-Y[,,t]
        for(i in 1:nrow(Y[,,t]))
        {
          yi<-Y[i,,t]
          zi<-zscores(yi)
          rnkd<-which( !is.na(yi) & yi>0 ) 
          if(length(rnkd)>0 && min(zi[rnkd])<0)
          { 
            zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
          }
          
          if(length(rnkd)<odmax[i]) 
          {
            urnkd<-which( !is.na(yi) & yi==0 ) 
            if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
          }
          
          Z[i,,t]<-zi
        } 
      }
    }
    
    # starting values for missing entries  
    ZA<-Z
    for (t in 1:N)
    { 
      mu<-mean(Z[,,t],na.rm=TRUE)
      a<-rowMeans(Z[,,t],na.rm=TRUE) ; b<-colMeans(Z[,,t],na.rm=TRUE)
      # a[is.na(a)] <- mean(a, na.rm=TRUE) ; b[is.na(b)] <- mean(b, na.rm=TRUE)
      a[is.na(a)] <- 0 ; b[is.na(b)] <- 0
      ZA[,,t]<-mu + outer(a,b,"+")
    }
    Z[is.na(Z)]<-ZA[is.na(Z)] 
    
    # other starting values
    beta<-rep(0,xP) 
    s2<-1 
    rho<-0
    Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
    U<-V<-matrix(0, nrow(Y[,,1]), R) 
  } # close of startVals condition    
  
  
  # unpack startVals list if applicable
  if(!is.null(startVals)){
    Z<-startVals$Z ; beta<-startVals$beta ; a<-startVals$a ; b<-startVals$b
    U<-startVals$U ; V<-startVals$V ; rho<-startVals$rho ; s2<-startVals$s2
    Sab<-startVals$Sab
  }
  
  return(
    list(
      Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab
    )
  )
  
}