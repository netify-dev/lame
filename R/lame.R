#' AME model fitting routine for replicated relational data
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to replicated relational data of
#' various types. 
#' 
#' This command provides posterior inference for parameters in AME models of
#' independent replicated relational data, assuming one of six possible data
#' types/models:
#' 
#' "nrm": A normal AME model.
#' 
#' "bin": A binary probit AME model.
#' 
#' "ord": An ordinal probit AME model. An intercept is not identifiable in this
#' model.
#' 
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' @usage lame(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl")
#' , cvar = TRUE, dcor = !symmetric, nvar=TRUE,  R = 0, family ="nrm",
#' intercept=!is.element(family,c("rrl","ord")),
#' symmetric=FALSE,
#' odmax=NULL, seed = 1,
#' nscan = 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE)
#' @param Y a T length list of n x n relational matrices, where T 
#' corresponds to the number of replicates (over time, for example). 
#' See family below for various data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)? 
#' @param dcor logical: fit a dyadic correlation (asymmetric case)?
#' @param nvar logical: fit nodal random effects (symmetric case)? 
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param family character: one of "nrm","bin","ord","cbin","frn","rrl" - see
#' the details below
#' @param intercept logical: fit model with an intercept?
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' families
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param startVals List from previous model run containing parameter
#' @param periodicSave logical: indicating whether to periodically save
#' MCMC results
#' @param saveInterval quantile interval indicating when to save during post
#' burn-in phase.
#' @param outFile character vector indicating name and path in which
#' file should be stored if periodicSave is selected. For example,
#' on an Apple OS outFile="~/Desktop/ameFit.rda".
#' starting values for new MCMC
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} 
#' \item{V}{posterior mean of multiplicative column effects v (asymmetric case)}
#' \item{UVPM}{posterior mean of UV}
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' \item{startVals}{Final parameter values from MCMC, can be used as the input
#' for a future model run.}
#' @author Peter Hoff, Yanjun He, Shahryar Minhas
#' @examples
#' 
#' data(YX_bin_list) 
#' fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="bin")
#' # you should run the Markov chain much longer than this
#' 
#' @export lame
lame <- function(
    Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL, 
    rvar = !(family=="rrl") , cvar = TRUE, dcor = !symmetric, 
    nvar = TRUE, 
    R = 0,
    family = "nrm",
    intercept = !is.element(family,c("rrl","ord")),
    symmetric = FALSE,
    odmax = NULL,
    seed = 1, nscan = 10000, burn = 500, odens = 25,
    plot = FALSE, print = FALSE, gof = TRUE, 
    startVals = NULL, periodicSave=FALSE, outFile=NULL,
    saveInterval=0.25
)
{
  #
  if( nscan %% odens !=0  ){ stop('"odens" must be a multiple of "nscan"')}
  
  # set random seed 
  set.seed(seed)
  
  # get actor info
  actorByYr <- lapply(Y, rownames)
  actorSet <- sort(unique(unlist( actorByYr ))) ; n <- length(actorSet)
  
  # reset odmax param
  odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), n )
  
  # calc savePoints
  savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
  savePoints <- savePoints[round(quantile(1:length(savePoints), probs=seq(saveInterval,1,saveInterval)))]
  
  # check formatting of input objects
  checkFormat(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol)
  
  # set diag to NA 
  N<-length(Y) ; pdLabs <- names(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})
  
  # convert into large array format
  arrayObj<-listToArray(actorSet, Y, Xdyad, Xrow, Xcol)
  Y<-arrayObj$Y ; Xdyad<-arrayObj$Xdyad ; Xrow<-arrayObj$Xrow
  Xcol<-arrayObj$Xcol ; rm(arrayObj)
  
  # force binary if binary family specified 
  if(is.element(family,c("bin","cbin"))) { Y<-1*(Y>0) } 
  
  # observed and max outdegrees 
  if(is.element(family,c("cbin","frn","rrl")) ){
    odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) } 
  }
  
  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }
  
  # construct design matrix    
  pr<-length(Xrow[,,1])/n
  pc<-length(Xcol[,,1])/n
  pd<-length(Xdyad[,,,1])/n^2
  designObj <- getDesignRep(
    Y=Y,Xdyad=Xdyad,Xrow=Xrow,Xcol=Xcol,actorSet=actorSet,
    intercept=intercept,n=n,N=N,pr=pr,pc=pc,pd=pd)
  Y<-designObj$Y ; X<-designObj$X ; Xlist<-designObj$Xlist
  XrLong<-designObj$XrLong ; XcLong<-designObj$XcLong
  mXLong<-designObj$mXLong ; mXtLong<-designObj$mXtLong
  xxLong<-designObj$xxLong ; xxTLong<-designObj$xxTLong ; rm(designObj)
  
  # design matrix warning for rrl
  if( family=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0)
      & !any( apply(X,c(3),function(x){var(c(x))})==0) )
  {
    cat("WARNING: row effects are not estimable using this procedure ","\n")
  }
  
  # design matrix warning for rrl and ord
  if( is.element(family,c("ord","rrl")) & 
      any( apply(X,c(3),function(x){var(c(x))})==0 ) )
  {
    cat("WARNING: an intercept is not estimable using this procedure ","\n")
  }
  
  # construct matrix of ranked nominations for frn, rrl   
  if(is.element(family,c("frn","rrl")))
  {
    ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
    YL<-list()
    for (t in 1:N) 
    {
      YL.t<-NULL
      warn<-FALSE
      for(i in 1:nrow(Y[,,1]))
      {
        yi<-Y[i,,t] ; rnkd<-which( !is.na(yi)&yi>0 )
        if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
        yi[rnkd]<-rank(yi[rnkd],ties.method="random")
        Y[i,,t]<-yi
        YL.t<-rbind(YL.t, match(1:ymx,yi))
      }
      YL[[t]]<-YL.t
      if(warn){cat("WARNING: Random reordering used to break ties in ranks\n")}
    }
  }
  
  # Get starting values for MCMC
  startValsObj <- getStartVals(startVals,Y,family,xP=dim(X)[3],rvar,cvar,R)
  Z<-startValsObj$Z ; beta<-startValsObj$beta ; a<-startValsObj$a
  b<-startValsObj$b ; U<-startValsObj$U ; V<-startValsObj$V
  rho<-startValsObj$rho ; s2<-startValsObj$s2 ; Sab<-startValsObj$Sab
  rm(list=c('startValsObj','startVals'))
  
  # helpful mcmc params
  symLoopIDs <- lapply(1:(nscan + burn), function(x){ rep(sample(1:nrow(U)),4) })  
  asymLoopIDs <- lapply(1:(nscan + burn), function(x){ sample(1:R) })  
  tryErrorChecks<-list(s2=0,betaAB=0,rho=0,UV=0)  
  iter <- 1
  
  # output items
  BETA <- matrix(nrow = nscan/odens, ncol = dim(X)[3] - pr*symmetric)
  VC<-matrix(nrow=nscan/odens,ncol=5-3*symmetric)
  UVPS <- U %*% t(V) * 0
  APS<-BPS<-rep(0,nrow(Y[,,1]))
  YPS<-array(0,dim=dim(Y),dimnames=dimnames(Y)) 
  # GOF <- matrix(NA, nrow=(nscan/odens)+1, ncol=5,
  #               dimnames=list(c('obs',1:(nscan/odens)),c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep", "trans.dep")))
  # GOF[1,] <- rowMeans(apply(Y,3,gofstats))
  # 5 corresponds to number of stats we want
  GOF <- array(NA, dim=c(5, N, (nscan/odens)+1))
  
  # add label to first dim
  dimnames(GOF)[[1]] <- c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep","trans.dep")
  
  # add row/col names
  if(!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
  
  # add time dimnames
  dimnames(GOF)[[3]] <- c('obs', 1:(nscan/odens))
  
  # fill in first entry with values from OBSERVED 
  GOF[,,1] <- apply(Y, 3, gofstats)
  
  names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
  
  # names of parameters, asymmetric case  
  if(!symmetric)
  { 
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
    colnames(BETA) <- dimnames(X)[[3]] 
  }   
  
  # names of parameters, symmetric case
  if(symmetric)
  { 
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
    bni<-bnames[1*intercept] 
    bnn<-gsub("row",bnames[rb],replacement="node") 
    bnd<-bnames[-c(1*intercept,rb,cb)]
    colnames(BETA)<-c(bni,bnn,bnd) 
  }    
  
  # MCMC
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE)) 
  
  if(burn!=0){
    pbBurn <- txtProgressBar(min=1,max=burn,style=3)
    cat('\nStarting burn-in period...\n')
  }
  if(!print){pbMain <- txtProgressBar(min=burn+1,max=nscan+burn,style=3)}
  for (s in 1:(nscan + burn)) 
  { 
    
    # update Z
    E.nrm<-array(dim=dim(Z))
    EZ <- get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
    for(t in 1:N ){
      if(family=="nrm")
      { 
        Z[,,t]<-rZ_nrm_fc(Z[,,t],EZ[,,t],rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ[,,t]
      }
      if(family=="bin"){ Z[,,t]<-rZ_bin_fc(Z[,,t],EZ[,,t],rho,Y[,,t]) }
      if(family=="ord"){ Z[,,t]<-rZ_ord_fc(Z[,,t],EZ[,,t],rho,Y[,,t]) }
      if(family=="cbin"){Z[,,t]<-rZ_cbin_fc(Z[,,t],EZ[,,t],rho,Y[,,t],odmax,odobs)}
      if(family=="frn")
      { 
        Z[,,t]<-rZ_frn_fc(Z[,,t],EZ[,,t],rho,Y[,,t],YL[[t]],odmax,odobs)
      }
      if(family=="rrl"){ Z[,,t]<-rZ_rrl_fc(Z[,,t],EZ[,,t],rho,Y[,,t],YL[[t]]) } 
    }
    
    # update s2
    if (family=="nrm"){
      s2New<-try(
        rs2_rep_fc_cpp(E.nrm,solve(matrix(c(1,rho,rho,1),2,2))), 
        silent=TRUE)
      if(class(s2New)!='try-error'){ s2 <- s2New } else { tryErrorChecks$s2<-tryErrorChecks$s2+1 }
    }
    
    # update beta, a b
    if( (pr+pc+pd+intercept)>0 ){
      iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
      tmp<-eigen(Sabs) ; k<-sum(zapsmall(tmp$val)>0 )
      G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
      betaABCalc <- try(
        rbeta_ab_rep_fc_cpp(
          ZT=sweep(Z,c(1,2),U%*%t(V)), Xr=XrLong, Xc=XcLong, mX=mXLong, mXt=mXtLong,
          XX=xxLong, XXt=xxTLong, iSe2=iSe2, Sabs=Sabs, k=k, G=G ),
        silent = TRUE)
    } else {
      betaABCalc <- try(
        rbeta_ab_rep_fc(sweep(Z,c(1,2),U%*%t(V)), Sab, rho, X, s2),
        silent = TRUE)
    }
    if(class(betaABCalc)!='try-error'){
      beta <- c(betaABCalc$beta)
      a <- c(betaABCalc$a) * rvar
      b <- c(betaABCalc$b) * cvar
      if(symmetric){ a<-b<-(a+b)/2 }
    } else { tryErrorChecks$betaAB<-tryErrorChecks$betaAB+1  }
    
    # update Sab - full SRM
    if(rvar & cvar & !symmetric)
    {
      Sab<-solve(rwish(solve(diag(2)+crossprod(cbind(a,b))),3+nrow(Z[,,1])))
    }
    
    # update Sab - rvar only
    if (rvar & !cvar & !symmetric)
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(a^2))/2)
    }
    
    # update Sab - cvar only 
    if (!rvar & cvar & !symmetric) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y[,,t]))/2, (1 + sum(b^2))/2)
    }
    
    # update Sab - symmetric case
    if(symmetric & nvar)
    {
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]   
    }
    
    # update rho
    if(dcor)
    {
      E.T <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
      rhoNew<-try( rrho_mh_rep_cpp(E.T, rho,s2), silent=TRUE )
      if(class(rhoNew)!='try-error'){ rho<-rhoNew } else { tryErrorChecks$rho<-tryErrorChecks$rho+1 }
    }
    
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }
    
    # update U,V
    if (R > 0)
    {
      E <- Z-get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U*0, V*0 )
      shrink<- (s>.5*burn)
      
      if(symmetric)
      { 
        EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
        UV<-try(
          rUV_sym_fc_cpp(EA, U, V, 
                         s2/dim(E)[3], shrink, symLoopIDs[[s]]-1), silent=TRUE )
        if(class(UV)=='try-error'){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
      }
      if(!symmetric){
        UV <- try(
          rUV_rep_fc_cpp(E, U, V, rho, s2,
                         mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)),
                         maxmargin=1e-6, shrink, asymLoopIDs[[s]]-1 ), silent = TRUE )
        if(class(UV)=='try-error'){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
      }      
      
      U<-UV$U ; V<-UV$V
    }
    
    # burn-in countdown
    if(burn!=0){setTxtProgressBar(pbBurn,s)}
    
    # store parameter values and monitor the MC
    if(s==burn+1&!print&burn!=0){cat('\nBurn-in period complete...');close(pbBurn)}
    if(s%%odens==0 & s>burn) 
    { 
      
      # store BETA and VC - symmetric case 
      if(symmetric){
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA[iter,]<-sbeta
        VC[iter,]<-c(Sab[1,1],s2)
      }
      
      # store BETA and VC - asymmetric case 
      if(!symmetric){
        BETA[iter,]<-beta
        VC[iter,]<- c(Sab[upper.tri(Sab, diag = T)], rho,s2)
      }
      
      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 
      
      # simulate from posterior predictive 
      EZ <- get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V );dimnames(EZ) <- dimnames(Y)
      Ys <- EZ*0
      for (t in 1:N)
      {
        if(symmetric){ EZ[,,t]<-(EZ[,,t]+t(EZ[,,t]))/2 }
        
        if(family=="bin"){ Ys[,,t]<-simY_bin(EZ[,,t],rho) }
        if(family=="cbin"){ Ys[,,t]<-1*(simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t])>0)}
        if(family=="frn"){ Ys[,,t]<-simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t]) }
        if(family=="rrl"){ Ys[,,t]<-simY_rrl(EZ[,,t],rho,odobs,YO=Y[,,t] ) }
        if(family=="nrm"){ Ys[,,t]<-simY_nrm(EZ[,,t],rho,s2) }
        if(family=="ord"){ Ys[,,t]<-simY_ord(EZ[,,t],rho,Y[,,t]) }
        
        if(symmetric)
        {  
          Yst<-Ys[,,t] ; Yst[lower.tri(Yst)]<-0 ; Ys[,,t]<-Yst+t(Yst)
        }
      } 
      
      # update posterior sum
      YPS<-YPS+Ys
      
      # save posterior predictive GOF stats
      #if(gof){Ys[is.na(Y)]<-NA ;GOF[(iter)+1,]<-rowMeans(apply(Ys,3,gofstats))}
      if(gof){
        
        Ys[is.na(Y)] <- NA
        GOF[,,(iter)+1] <- apply(Ys,3,gofstats)
      }
      
      # print MC progress 
      if(print)
      {
        cat('\n',s,
            round(apply(BETA[1:iter,,drop=FALSE],2,mean),2),":",
            round(apply(VC[1:iter,,drop=FALSE],2,mean),2),"\n")
        if (have_coda & nrow(VC[1:iter,,drop=FALSE]) > 3 & length(beta)>0) 
        {
          cat(round(coda::effectiveSize(BETA[1:iter,,drop=FALSE])), "\n")
        }
      }
      
      # periodic save
      if(periodicSave & s %in% savePoints & !is.null(outFile)){
        # save startVals for future model runs
        startVals <- list(Z=Z,beta=beta,a=a,b=b,U=U,V=V,rho=rho,s2=s2,Sab=Sab)
        fit <- getFitObject( APS=APS, BPS=BPS, UVPS=UVPS, YPS=YPS, 
                             BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr,
                             startVals=startVals, symmetric=symmetric, tryErrorChecks=tryErrorChecks)
        save(fit, file=outFile) ; rm(list=c('fit','startVals'))
      }
      
      # plot MC results
      if(plot & s==(burn+nscan))
      {
        # plot VC
        paramPlot(VC)
        
        # plot BETA
        if(length(beta)>0) 
        {
          betaIndices<-split(1:ncol(BETA), ceiling(seq_along(1:ncol(BETA))/5))
          for(bIndex in betaIndices){
            paramPlot( BETA[,bIndex,drop=FALSE] ) }
        }
        
        # plot GOF
        if(gof)
        {
          suppressMessages( print( gofPlot(GOF, symmetric=symmetric) ) )
        }
      } # plot code if applicable
      iter<-iter+1
    } # post burn-in
    if(!print){setTxtProgressBar(pbMain,s)}
  } # end MCMC  
  if(!print){close(pbMain)}
  
  # save startVals for future model runs
  startVals <- list( Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
  
  # output
  fit <- getFitObject( APS=APS, BPS=BPS, UVPS=UVPS, YPS=YPS, 
                       BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr, 
                       startVals=startVals, symmetric=symmetric, tryErrorChecks=tryErrorChecks)
  class(fit) <- "lame" # set class
  return(fit) # output object to workspace
  
}