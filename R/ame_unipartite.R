#' AME model fitting for unipartite (square) networks
#' 
#' Internal function that handles AME model fitting for unipartite networks.
#' This function contains the core logic extracted from the main ame() function
#' specifically for square adjacency matrices.
#' 
#' @inheritParams ame
#' @return An ame object with posterior samples and estimates
#' @keywords internal
#' @noRd
ame_unipartite <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
                          rvar = !(family=="rrl"), cvar = TRUE, dcor = !symmetric, 
                          nvar=TRUE, R = 0, family="normal",
                          intercept=!is.element(family,c("rrl","ordinal")), 
                          symmetric=FALSE,
                          odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
                          prior=list(), g=NA,
                          seed = 6886, nscan = 10000, burn = 500, odens = 25,
                          plot=TRUE, print = TRUE, gof=TRUE, 
                          start_vals=NULL, periodic_save=FALSE, out_file=NULL,
                          save_interval=0.25, model.name=NULL) {
  
  # set random seed
  set.seed(seed)
  
  # set diag to NA
  diag(Y) <- NA
  
  # force binary if binary family specified
  if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
  
  # handle tobit family
  if(family=="tobit") { Y[Y<0]<-0 } 
  
  # observed and max outdegrees 
  if(is.element(family,c("cbin","frn","rrl"))) {
    odobs<-apply(Y>0,1,sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y)) }
  }
  
  # Network dimensions
  n <- nrow(Y)
  
  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }
  
  # Set g-prior default
  if(is.na(g)) {
    if(family=="normal") { 
      g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)
    } else if(family=="tobit") {
      g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)*4
    } else {
      g <- sum(!is.na(Y))
    }
  }
  
  # process prior list and set defaults
  if(!is.list(prior)) { prior <- list() }
  if(is.null(prior$Sab0)) { prior$Sab0 <- diag(2) }
  if(is.null(prior$eta0)) { prior$eta0 <- round(4 + 3*n/100) }
  if(is.null(prior$etaab)) { prior$etaab <- round(4 + 3*n/100) }
  
  # construct design matrix
  pr<-length(Xrow)/n
  pc<-length(Xcol)/n
  pd<-length(Xdyad)/n^2
  
  X <- design_array(Xrow, Xcol, Xdyad, intercept, n)
  
  # design matrix warning for rrl 
  if( family=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0) 
      & !any( apply(X,3,function(x){var(c(x))})==0) ) {
    cli::cli_warn("Row effects are not estimable using this procedure")
  } 
  
  # design matrix warning for rrl and ord
  if( is.element(family,c("ord","rrl")) & 
      any( apply(X,3,function(x){var(c(x))})==0 ) ) {
    cli::cli_warn("An intercept is not estimable using this procedure")
  } 
  
  # construct matrix of ranked nominations for frn, rrl 
  if(is.element(family,c("frn","rrl"))) {
    ymx<-max(apply(1*(Y>0),1,sum,na.rm=TRUE))
    YL<-NULL
    warn<-FALSE
    for(i in 1:nrow(Y)) {
      yi<-Y[i,] ; rnkd<-which( !is.na(yi)&yi>0 )
      if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
      yi[rnkd]<-rank(yi[rnkd],ties.method="random")
      Y[i,]<-yi
      YL<-rbind(YL, match(1:ymx,yi))
    }
    if(warn){cli::cli_warn("Random reordering used to break ties in ranks")}
  }
  
  # starting Z values
  if(family=="normal") { Z<-Y }
  if(family=="tobit") { 
    Z<-Y 
    if(sum(Y>0, na.rm=TRUE) > 0) {
      Z[Y==0]<-min(Y[Y>0],na.rm=TRUE)/2 
    }
  }
  if(family=="ordinal") { Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) } 
  if(family=="rrl") { Z<-matrix(t(apply(Y,1,zscores)),nrow(Y),ncol(Y)) }  
  if(family=="binary") { 
    Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) 
    if(sum(Y==0, na.rm=TRUE) > 0 && sum(Y==1, na.rm=TRUE) > 0) {
      z01<- .5* ( max(Z[Y==0],na.rm=TRUE) + min(Z[Y==1],na.rm=TRUE) ) 
      Z<-Z - z01
    }
  } 
  if(family=="poisson") { 
    Z<-log(Y+1) 
    diag(Z)<-0
  }
  
  if(is.element(family,c("cbin","frn"))) {
    Z<-Y
    for(i in 1:nrow(Y)) {
      yi<-Y[i,]
      zi<-zscores(yi)
      rnkd<-which( !is.na(yi) & yi>0 ) 
      if(length(rnkd)>0 && min(zi[rnkd])<0) { 
        zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
      }
      
      if(length(rnkd)<odmax[i]) {
        urnkd<-which( !is.na(yi) & yi==0 ) 
        if(length(urnkd) > 0 && max(zi[urnkd])>0) { 
          zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 
        }
      }
      
      Z[i,]<-zi
    } 
  }
  
  # starting values for missing entries 
  mu<-mean(Z,na.rm=TRUE) 
  a<-rowMeans(Z,na.rm=TRUE) ; b<-colMeans(Z,na.rm=TRUE)  
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  ZA<-mu + outer(a,b,"+") 
  Z[is.na(Z)]<-ZA[is.na(Z)] 
  
  # other starting values (use start_vals if provided)
  if(!is.null(start_vals)) {
    beta <- start_vals$beta
    s2 <- start_vals$s2
    a <- start_vals$a
    b <- start_vals$b
    if(!is.null(start_vals$rho)) rho <- start_vals$rho else rho <- 0
    if(!is.null(start_vals$Sab)) Sab <- start_vals$Sab
    if(!is.null(start_vals$U)) U <- start_vals$U
    if(!is.null(start_vals$V)) V <- start_vals$V
    if(!is.null(start_vals$Z)) Z <- start_vals$Z
  } else {
    beta<-rep(0,dim(X)[3]) 
    s2<-1
  }
  
  # Initialize remaining values if not provided by start_vals
  if(is.null(start_vals)) {
    rho<-0
    Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
    U<-V<-matrix(0, nrow(Y), R)
  }
  
  # output items - pre-allocate for efficiency
  n_samples <- floor((nscan + burn * (burn > 0)) / odens) - floor(burn / odens)
  BETA <- matrix(NA_real_, nrow = n_samples, ncol = dim(X)[3] - pr*symmetric)
  VC <- matrix(NA_real_, nrow = n_samples, ncol = 5 - 3*symmetric)
  sample_idx <- 0  # Track current sample index
  UVPS <- U %*% t(V) * 0
  APS<-BPS<- rep(0,nrow(Y))
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
  
  if(gof) {
    GOF <- matrix(NA_real_, nrow = n_samples + 1, ncol = 5)
    GOF[1,] <- gof_stats(Y)
    gof_idx <- 1
    rownames(GOF) <- c("obs", rep("", n_samples))
  } else {
    GOF <- matrix(gof_stats(Y), 1, 5)
    rownames(GOF) <- "obs"
  }  
  colnames(GOF) <- c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep","trans.dep")
  
  names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y)
  
  # names of parameters, asymmetric case 
  if(!symmetric) {
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
    colnames(BETA) <- dimnames(X)[[3]] 
  }
  
  # names of parameters, symmetric case
  if(symmetric) {
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
    bni<-bnames[1*intercept] 
    bnn<-gsub("row",bnames[rb],replacement="node")       
    bnd<-bnames[-c(1*intercept,rb,cb)]
    colnames(BETA)<-c(bni,bnn,bnd) 
  }    
  
  # marginal means and regression sums of squares
  Xr<-apply(X,c(1,3),sum)            # row sum
  Xc<-apply(X,c(2,3),sum)            # col sum
  mX<- apply(X,3,c)                  # design matrix
  mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
  XX<-t(mX)%*%mX                     # regression sums of squares
  XXt<-t(mX)%*%mXt                   # crossproduct sums of squares
  
  # Pre-compute X with attributes ONCE (not in loop)
  X_precomp <- X
  attributes(X_precomp) <- c(attributes(X), list(Xr=Xr, Xc=Xc, mX=mX, 
                                                 mXt=mXt, XX=XX, XXt=XXt))
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  # Set up save points for periodic saving
  savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
  savePoints <- savePoints[round(quantile(1:length(savePoints), 
                                         probs=seq(save_interval,1,save_interval)))]
  
  # Initialize running statistics for efficient diagnostics
  beta_running_mean <- rep(0, dim(X)[3] - pr*symmetric)
  vc_running_mean <- rep(0, 5 - 3*symmetric)
  running_sample_count <- 0
  
  # Pre-select family-specific functions to avoid string comparisons in loop
  update_Z_fn <- switch(family,
    normal = rZ_nrm_fc,
    tobit = rZ_tob_fc,
    binary = rZ_bin_fc,
    ordinal = rZ_ord_fc,
    cbin = function(Z, EZ, rho, Y) rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs),
    frn = function(Z, EZ, rho, Y) rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs),
    rrl = function(Z, EZ, rho, Y) rZ_rrl_fc(Z, EZ, rho, Y, YL),
    poisson = rZ_pois_fc,
    stop(paste("Unknown family:", family))
  )
  
  simulate_Y_fn <- switch(family,
    normal = function(EZ, rho, s2) simY_nrm(EZ, rho, s2),
    tobit = function(EZ, rho, s2) simY_tob(EZ, rho, s2),
    binary = function(EZ, rho, s2) simY_bin(EZ, rho),
    ordinal = function(EZ, rho, s2) simY_ord(EZ, rho, Y),
    cbin = function(EZ, rho, s2) 1*(simY_frn(EZ, rho, odmax, YO=Y)>0),
    frn = function(EZ, rho, s2) simY_frn(EZ, rho, odmax, YO=Y),
    rrl = function(EZ, rho, s2) simY_rrl(EZ, rho, odobs, YO=Y),
    poisson = function(EZ, rho, s2) simY_pois(EZ),
    stop(paste("Unknown family:", family))
  )
  
  # Family-specific flags
  needs_s2 <- family %in% c("normal", "tobit", "poisson")
  
  # Initialize cached components for EZ computation
  Xbeta_cache <- if(dim(X)[3] > 0) Xbeta(X, beta) else matrix(0, n, n)
  ab_outer_cache <- outer(a, b, "+")
  UV_cache <- U %*% t(V)
  EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
  
  for (s in 1:(nscan + burn)) {
    
    # update Z (use cached EZ)
    EZ <- EZ_cache
    
    # Update Z using pre-selected function
    rho_eff <- rho
    Z <- if(family %in% c("normal", "tobit", "poisson")) {
      update_Z_fn(Z, EZ, rho_eff, s2, Y)
    } else {
      update_Z_fn(Z, EZ, rho_eff, Y)
    }
    
    # update s2
    if (needs_s2) {
      s2<-rs2_fc(Z,rho_eff,offset=EZ)  
    }
    
    # update beta, a b with g-prior
    tmp <- rbeta_ab_fc(Z, Sab, rho, X_precomp, s2, offset=UV_cache, g=g)
    beta_new <- tmp$beta
    a_new <- tmp$a * rvar
    b_new <- tmp$b * cvar
    
    if(symmetric){ a_new<-b_new<-(a_new+b_new)/2 }
    
    # Update cache if parameters changed
    if(dim(X)[3] > 0 && !isTRUE(all.equal(beta, beta_new))) {
      Xbeta_cache <- Xbeta(X, beta_new)
    }
    if(!isTRUE(all.equal(a, a_new)) || !isTRUE(all.equal(b, b_new))) {
      ab_outer_cache <- outer(a_new, b_new, "+")
    }
    
    beta <- beta_new
    a <- a_new
    b <- b_new
    
    # Update EZ cache with new components
    EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
    
    # update Sab using unified function
    if(is.element(family,c("normal","tobit","ordinal"))) { 
      Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                     rvar=rvar, cvar=cvar, symmetric=symmetric)
    }
    
    # special updates for discrete families
    if(family=="binary") {
      if(rvar & cvar & !symmetric) {
        tmp<-raSab_bin_fc(Z,Y,a,b,Sab,Sab0=prior$Sab0/prior$etaab,eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    if(family=="cbin") {
      if(rvar & cvar & !symmetric) {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,Sab0=prior$Sab0/prior$etaab,
                          eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    if(family=="frn") { 
      if(rvar & cvar & !symmetric) {
        tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs,Sab0=prior$Sab0/prior$etaab,
                         eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    # update rho
    if(dcor) {
      rho<-rrho_mh(Z, rho, s2, offset=EZ_cache)
    }
    
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }
    
    # update U,V
    if(R > 0) {
      shrink<- (s>.5*burn)
      offset <- Xbeta_cache + ab_outer_cache
      
      if(symmetric){ 
        E<-Z-offset ; E<-.5*(E+t(E))
        UV<-rUV_sym_fc(E, U, V, s2,shrink) 
      }
      if(!symmetric){
        # Compute Suv for multiplicative effects
        if(shrink) {
          Suv <- rSuv_fc(U, V, Suv0=diag(2*R)/(n+R+2), kappa0=n+R+2)
        } else {
          Suv <- diag(n, nrow=2*R)
        }
        UV<-rUV_fc(Z, U, V, Suv, rho, s2, offset) 
      }
      
      U<-UV$U ; V<-UV$V
      # Update UV cache
      UV_cache <- U %*% t(V)
      EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
    }    
    
    # burn-in countdown
    if(print & s%%odens==0 & s<=burn){
      message(paste0(round(100*s/burn,1), "% burn-in complete"))
    }
    
    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) {  
      sample_idx <- sample_idx + 1
      
      # save BETA and VC - symmetric case 
      if(symmetric) {
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2 
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA[sample_idx,] <- sbeta
        
        VC[sample_idx,] <- c(Sab[1,1],s2)
        
        # Update running means
        beta_running_mean <- beta_running_mean + (sbeta - beta_running_mean) / sample_idx
        vc_running_mean <- vc_running_mean + (c(Sab[1,1],s2) - vc_running_mean) / sample_idx
      }
      
      # save BETA and VC - asymmetric case 
      if(!symmetric) {
        BETA[sample_idx,] <- beta
        VC[sample_idx,] <- c(Sab[upper.tri(Sab, diag = T)], rho, s2)
        
        # Update running means
        beta_running_mean <- beta_running_mean + (beta - beta_running_mean) / sample_idx
        vc_running_mean <- vc_running_mean + (c(Sab[upper.tri(Sab, diag = T)], rho, s2) - vc_running_mean) / sample_idx
      }
      
      running_sample_count <- sample_idx
      
      # update posterior sums of random effects
      UVPS <- UVPS + UV_cache  # Use cached value
      APS <- APS + a
      BPS <- BPS + b 
      
      # simulate from posterior predictive (use cached EZ)
      EZ <- EZ_cache
      if(symmetric){ EZ<-(EZ+t(EZ))/2 }
      
      # Simulate Y for posterior predictive (use pre-selected function)
      rho_eff <- rho
      Ys <- simulate_Y_fn(EZ, rho_eff, s2)
      
      if(symmetric){ Ys[lower.tri(Ys)]<-0 ; Ys<-Ys+t(Ys)  }
      
      # update posterior sum
      YPS<-YPS+Ys
      
      # save posterior predictive GOF stats
      if(gof){ 
        Ys[is.na(Y)]<-NA
        gof_idx <- gof_idx + 1
        GOF[gof_idx,] <- gof_stats(Ys)
      }
      
      # print MC progress (use running means for efficiency)
      if (print) {
        beta_means <- round(beta_running_mean, 2)
        vc_means <- round(vc_running_mean, 2)
        cli::cli_text("Iteration {.val {s}}: beta = [{.field {paste(beta_means, collapse=', ')}}], VC = [{.field {paste(vc_means, collapse=', ')}}]")
        
        # Only compute expensive diagnostics occasionally
        if (have_coda & running_sample_count > 100 & running_sample_count %% 100 == 0 & length(beta) > 0) {
          # Compute ESS only every 100 samples after first 100
          eff_sizes <- round(coda::effectiveSize(BETA[1:sample_idx,, drop=FALSE]))
          cli::cli_text("  Effective sizes: [{.emph {paste(eff_sizes, collapse=', ')}}]")
        }
      }
      
    }
    
    # periodic save
    if(periodic_save & s %in% savePoints & !is.null(out_file)){
      # save start_vals for future model runs
      start_vals <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
      
      fit <- list(APM=APS/(sum(s>burn & s%%odens==0)), 
                 BPM=BPS/(sum(s>burn & s%%odens==0)),
                 BETA=BETA, VC=VC, GOF=GOF, 
                 start_vals=start_vals, symmetric=symmetric,
                 model.name=model.name)
      save(fit, file=out_file)
      rm(list=c('fit','start_vals'))
    }
    
  } # end MCMC   
  
  # output 
  
  # Trim matrices to actual sample size
  if(exists("sample_idx") && sample_idx > 0) {
    BETA <- BETA[1:sample_idx, , drop=FALSE]
    VC <- VC[1:sample_idx, , drop=FALSE]
    if(gof && exists("gof_idx")) {
      GOF <- GOF[1:gof_idx, , drop=FALSE]
    }
  }
  
  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)
  # Add names
  if(is.null(names(APM)) && !is.null(rownames(Y))) {
    names(APM) <- rownames(Y)
  }
  if(is.null(names(BPM)) && !is.null(rownames(Y))) {
    names(BPM) <- rownames(Y)
  }
  # If still no names, create default ones
  if(is.null(names(APM))) {
    names(APM) <- paste0("Actor", 1:length(APM))
  }
  if(is.null(names(BPM))) {
    names(BPM) <- paste0("Actor", 1:length(BPM))
  }
  UVPM<-UVPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  
  EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(APM,BPM,"+")+UVPM
  
  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  
  # Transform EZ to count scale for Poisson family
  if(family == "poisson") {
    EZ <- exp(EZ)
    EZ[EZ > 1e6] <- 1e6
  }
  
  dimnames(YPM)<-dimnames(EZ)<-dimnames(Y)
  rownames(BETA)<-NULL
  
  # asymmetric output 
  if(!symmetric) {
    if(R > 0) {
      UDV<-svd(UVPM)
      U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      rownames(U)<-rownames(V)<-rownames(Y)
    } else {
      # No multiplicative effects
      U <- NULL
      V <- NULL
    }
    # save start_vals for future model runs
    start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
    fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
                YPM=YPM,GOF=GOF,start_vals=start_vals_final,model.name=model.name,
                mode="unipartite",family=family,symmetric=symmetric,odmax=odmax)
  }
  
  # symmetric output
  if(symmetric) {
    ULUPM<-UVPM
    if(R > 0) {
      eULU<-eigen(ULUPM)
      eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
      U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
      L<-eULU$val[eR]
      rownames(U)<-rownames(ULUPM)
    } else {
      U <- NULL
      L <- NULL
    }
    rownames(ULUPM)<-colnames(ULUPM)<-rownames(Y)
    EZ<-.5*(EZ+t(EZ)) ; YPM<-.5*(YPM+t(YPM)) 
    # save start_vals for future model runs
    start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=U, rho=rho, s2=s2, Sab=Sab)
    fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
              YPM=YPM,GOF=GOF,start_vals=start_vals_final,model.name=model.name,
              mode="unipartite",family=family,symmetric=symmetric,odmax=odmax)
  } 
  
  # Ensure scalars are properly typed
  if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
  if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
  
  class(fit) <- "ame"
  fit
}