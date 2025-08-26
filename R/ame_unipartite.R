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
ame_unipartite <- function(
  Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
  rvar = !(family=="rrl"), cvar = TRUE, dcor = !symmetric, 
  nvar=TRUE, R = 0, family="normal",
  intercept=!is.element(family,c("rrl","ordinal")), 
  symmetric=FALSE,
  odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
  prior=list(), g=NA,
  seed = 6886, nscan = 10000, burn = 500, odens = 25,
  print = TRUE, gof=TRUE, custom_gof=NULL,
  start_vals=NULL, periodic_save=FALSE, out_file=NULL,
  save_interval=0.25, model.name=NULL,
  posterior_opts = NULL, use_sparse_matrices = FALSE
  ){
  
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
  
  # Initialize posterior sample storage if requested
  U_samples <- NULL ; V_samples <- NULL
  a_samples <- NULL ; b_samples <- NULL
  
  if(!is.null(posterior_opts)) {
    if(posterior_opts$save_UV && R > 0) {
      # Calculate thinned sample size
      n_UV_samples <- ceiling(n_samples / posterior_opts$thin_UV)
      U_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
      if(!symmetric) {
        V_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
      }
    }
    if(posterior_opts$save_ab && (rvar || cvar)) {
      # Calculate thinned sample size
      n_ab_samples <- ceiling(n_samples / posterior_opts$thin_ab)
      if(rvar) a_samples <- matrix(NA_real_, n, n_ab_samples)
      if(cvar && !symmetric) b_samples <- matrix(NA_real_, n, n_ab_samples)
    }
  }
  
  # Track sample indices for thinned storage
  UV_sample_idx <- 0
  ab_sample_idx <- 0
  
  UVPS <- U %*% t(V) * 0
  APS<-BPS<- rep(0,nrow(Y))
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
  
  if(gof) {
    # Determine number of GOF statistics
    n_base_stats <- 5  # For unipartite: sd.rowmean, sd.colmean, dyad.dep, cycle.dep, trans.dep
    
    # Process custom GOF functions
    custom_gof_names <- NULL
    n_custom_stats <- 0
    if(!is.null(custom_gof)) {
      if(is.function(custom_gof)) {
        # Single function - get example output to determine size
        tryCatch({
          test_output <- custom_gof(Y)
          n_custom_stats <- length(test_output)
          custom_gof_names <- names(test_output)
          if(is.null(custom_gof_names)) {
            custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
          }
        }, error = function(e) {
          warning("Custom GOF function failed on initial test. It will be skipped.")
          custom_gof <- NULL
        })
      } else if(is.list(custom_gof)) {
        # List of functions
        n_custom_stats <- length(custom_gof)
        custom_gof_names <- names(custom_gof)
        if(is.null(custom_gof_names)) {
          custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
        }
      }
    }
    
    # Initialize GOF matrix
    n_total_stats <- n_base_stats + n_custom_stats
    GOF <- matrix(NA_real_, nrow = n_samples + 1, ncol = n_total_stats)
    
    # Compute initial GOF for observed data
    base_gof <- gof_stats(Y)
    GOF[1, 1:n_base_stats] <- base_gof
    
    # Compute initial custom GOF if provided
    if(!is.null(custom_gof)) {
      if(is.function(custom_gof)) {
        tryCatch({
          custom_vals <- custom_gof(Y)
          GOF[1, (n_base_stats + 1):(n_base_stats + length(custom_vals))] <- custom_vals
        }, error = function(e) {})
      } else if(is.list(custom_gof)) {
        for(i in seq_along(custom_gof)) {
          if(is.function(custom_gof[[i]])) {
            tryCatch({
              GOF[1, n_base_stats + i] <- custom_gof[[i]](Y)
            }, error = function(e) {})
          }
        }
      }
    }
    
    gof_idx <- 1
    rownames(GOF) <- c("obs", rep("", n_samples))
    
    # Set column names
    base_names <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    all_names <- c(base_names, custom_gof_names)
    colnames(GOF) <- all_names
  } else {
    GOF <- matrix(gof_stats(Y), 1, 5)
    rownames(GOF) <- "obs"
    colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    custom_gof <- NULL  # Don't compute custom GOF if gof=FALSE
  }
  
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
  attributes(X_precomp) <- c(
    attributes(X), list(Xr=Xr, Xc=Xc, mX=mX, mXt=mXt, XX=XX, XXt=XXt))
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  # Set up save points for periodic saving
  savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
  savePoints <- savePoints[
    round(quantile(1:length(savePoints), 
    probs=seq(save_interval,1,save_interval))) ]
  
  # Initialize running statistics for efficient diagnostics
  beta_running_mean <- rep(0, dim(X)[3] - pr*symmetric)
  vc_running_mean <- rep(0, 5 - 3*symmetric)
  running_sample_count <- 0
  
  # Initialize posterior sample storage if requested
  UV_sample_idx <- 0
  ab_sample_idx <- 0
  if(!is.null(posterior_opts)) {
    if(posterior_opts$save_UV && R > 0) {
      n_UV_samples <- floor((nscan) / (odens * posterior_opts$thin_UV))
      if(n_UV_samples > 0) {
        U_samples <- array(NA, dim = c(n, R, n_UV_samples))
        V_samples <- array(NA, dim = c(n, R, n_UV_samples))
      }
    }
    if(posterior_opts$save_ab && (rvar || cvar)) {
      n_ab_samples <- floor((nscan) / (odens * posterior_opts$thin_ab))
      if(n_ab_samples > 0) {
        a_samples <- matrix(NA, n, n_ab_samples)
        b_samples <- matrix(NA, n, n_ab_samples)
      }
    }
  }
  
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
  
  # Initialize progress bars if printing
  if(print) {
    cli::cli_h2("Running MCMC for unipartite network")
    cli::cli_text("Dimensions: {.val {n}} x {.val {n}} nodes")
    cli::cli_text("R = {.val {R}}, symmetric = {.val {symmetric}}")
    cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")
    
    # Estimate time impact of GOF
    if(gof) {
      n_gof_stats <- 5
      if(!is.null(custom_gof)) {
        if(is.function(custom_gof)) {
          tryCatch({
            n_gof_stats <- n_gof_stats + length(custom_gof(Y))
          }, error = function(e) {})
        } else if(is.list(custom_gof)) {
          n_gof_stats <- n_gof_stats + length(custom_gof)
        }
      }
    }
    
    # Start timing
    start_time <- Sys.time()
    
    if(burn > 0) {
      cli::cli_progress_bar("Burn-in", total = burn, clear = TRUE)
    }
  }
  
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
    
    # Update progress bars
    if(print) {
      if(s <= burn) {
        cli::cli_progress_update(id = NULL)
        if(s == burn) {
          cli::cli_progress_done()
          cli::cli_progress_bar("Sampling", total = nscan, clear = TRUE)
        }
      } else if(s > burn) {
        cli::cli_progress_update(id = NULL)
      }
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
      
      # Save posterior samples if requested
      if(!is.null(posterior_opts)) {
        # Save U/V samples (thinned)
        if(posterior_opts$save_UV && R > 0) {
          if(sample_idx %% posterior_opts$thin_UV == 0) {
            UV_sample_idx <- UV_sample_idx + 1
            if(UV_sample_idx <= dim(U_samples)[3]) {
              U_samples[,,UV_sample_idx] <- U
              if(!symmetric && !is.null(V_samples)) {
                V_samples[,,UV_sample_idx] <- V
              }
            }
          }
        }
        
        # Save a/b samples (thinned)
        if(posterior_opts$save_ab) {
          if(sample_idx %% posterior_opts$thin_ab == 0) {
            ab_sample_idx <- ab_sample_idx + 1
            if(rvar && !is.null(a_samples) && ab_sample_idx <= ncol(a_samples)) {
              a_samples[,ab_sample_idx] <- a
            }
            if(cvar && !symmetric && !is.null(b_samples) && ab_sample_idx <= ncol(b_samples)) {
              b_samples[,ab_sample_idx] <- b
            }
          }
        }
      }
      
      running_sample_count <- sample_idx
      
      # update posterior sums of random effects
      UVPS <- UVPS + UV_cache  # Use cached value
      APS <- APS + a
      BPS <- BPS + b 
      
      # Save posterior samples if requested
      if(!is.null(posterior_opts)) {
        if(posterior_opts$save_UV && R > 0 && sample_idx %% posterior_opts$thin_UV == 0) {
          UV_sample_idx <- UV_sample_idx + 1
          if(UV_sample_idx <= dim(U_samples)[3]) {
            U_samples[,,UV_sample_idx] <- U
            V_samples[,,UV_sample_idx] <- V
          }
        }
        if(posterior_opts$save_ab && (rvar || cvar) && sample_idx %% posterior_opts$thin_ab == 0) {
          ab_sample_idx <- ab_sample_idx + 1
          if(ab_sample_idx <= ncol(a_samples)) {
            a_samples[,ab_sample_idx] <- a
            b_samples[,ab_sample_idx] <- b
          }
        }
      } 
      
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
        # Compute base GOF statistics
        GOF[gof_idx, 1:5] <- gof_stats(Ys)
        
        # Compute custom GOF statistics if provided
        if(!is.null(custom_gof)) {
          if(is.function(custom_gof)) {
            tryCatch({
              custom_vals <- custom_gof(Ys)
              GOF[gof_idx, (6):(5 + length(custom_vals))] <- custom_vals
            }, error = function(e) {})
          } else if(is.list(custom_gof)) {
            for(i in seq_along(custom_gof)) {
              if(is.function(custom_gof[[i]])) {
                tryCatch({
                  GOF[gof_idx, 5 + i] <- custom_gof[[i]](Ys)
                }, error = function(e) {})
              }
            }
          }
        }
      }
      
      # Update progress occasionally (no verbose output)
      
    }
    
    # periodic save
    if(periodic_save & s %in% savePoints & !is.null(out_file)){
      # save start_vals for future model runs
      start_vals <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
      
      fit <- list(APM=APS/(sum(s>burn & s%%odens==0)), 
                 BPM=BPS/(sum(s>burn & s%%odens==0)),
                 BETA=BETA, VC=VC, GOF=GOF, X=X, Y=Y,
                 start_vals=start_vals, symmetric=symmetric,
                 model.name=model.name)
      save(fit, file=out_file)
      rm(list=c('fit','start_vals'))
    }
    
  } # end MCMC   
  
  # Close progress bar if printing
  if(print) {
    cli::cli_progress_done()
    
    # Report timing
    if(exists("start_time")) {
      end_time <- Sys.time()
      elapsed <- difftime(end_time, start_time, units = "secs")
      cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
    }
  }
  
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
  UVPM<-UVPS/nrow(VC)  # Needed temporarily for computing U,V via SVD
  YPM<-YPS/nrow(VC) 
  
  # note to self, EZ can be reconstructed as: 
  # Xbeta(X,colMeans(BETA)) + outer(APM,BPM,"+") + U%*%t(V)
  # leave otu to save memory
  
  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  
  dimnames(YPM)<-dimnames(Y)
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
    fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,L=NULL,
                YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
                mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R)
    
    # Add posterior samples if saved
    if(!is.null(posterior_opts)) {
      if(posterior_opts$save_UV && !is.null(U_samples) && UV_sample_idx > 0) {
        actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
        if(actual_UV_samples > 0) {
          fit$U_samples <- U_samples[,,1:actual_UV_samples, drop=FALSE]
          if(!is.null(V_samples)) {
            fit$V_samples <- V_samples[,,1:actual_UV_samples, drop=FALSE]
          }
        }
      }
      if(posterior_opts$save_ab && ab_sample_idx > 0) {
        if(rvar && !is.null(a_samples)) {
          actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
          if(actual_ab_samples > 0) {
            fit$a_samples <- a_samples[,1:actual_ab_samples, drop=FALSE]
          }
        }
        if(cvar && !is.null(b_samples)) {
          actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
          if(actual_ab_samples > 0) {
            fit$b_samples <- b_samples[,1:actual_ab_samples, drop=FALSE]
          }
        }
      }
    }
  }
  
  # symmetric output
  if(symmetric) {
    ULUPM<-UVPM
    if(R > 0) {
      eULU<-eigen(ULUPM)
      eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
      U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
      L<-diag(eULU$val[eR], nrow=R)  # Create diagonal matrix from eigenvalues
      rownames(U)<-rownames(ULUPM)
    } else {
      U <- NULL
      L <- NULL
    }
    YPM<-.5*(YPM+t(YPM))  # Symmetrize YPM
    # save start_vals for future model runs
    start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=U, rho=rho, s2=s2, Sab=Sab)
    # note, ULUPM not stored: can be reconstructed as U %*% L %*% t(U)
    # BPM is NULL for symmetric (same as APM)
    fit<-list(BETA=BETA,VC=VC,APM=APM,BPM=NULL,U=U,V=NULL,L=L,
              YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
              mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R)  # For symmetric, V = U
    
    # Add posterior samples if saved
    if(!is.null(posterior_opts)) {
      if(posterior_opts$save_UV && !is.null(U_samples) && UV_sample_idx > 0) {
        actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
        if(actual_UV_samples > 0) {
          fit$U_samples <- U_samples[,,1:actual_UV_samples, drop=FALSE]
          if(!is.null(V_samples)) {
            fit$V_samples <- V_samples[,,1:actual_UV_samples, drop=FALSE]
          }
        }
      }
      if(posterior_opts$save_ab && ab_sample_idx > 0) {
        if(rvar && !is.null(a_samples)) {
          actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
          if(actual_ab_samples > 0) {
            fit$a_samples <- a_samples[,1:actual_ab_samples, drop=FALSE]
          }
        }
        if(cvar && !is.null(b_samples)) {
          actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
          if(actual_ab_samples > 0) {
            fit$b_samples <- b_samples[,1:actual_ab_samples, drop=FALSE]
          }
        }
      }
    }
  } 
  
  # Ensure scalars are properly typed
  if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
  if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
  
  class(fit) <- "ame"
  
  # Apply sparse matrix optimization if requested
  if(use_sparse_matrices) {
    fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
  }
  
  fit
}