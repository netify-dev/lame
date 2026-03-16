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
	verbose = TRUE, gof=TRUE, custom_gof=NULL,
	start_vals=NULL, periodic_save=FALSE, out_file=NULL,
	save_interval=0.25, model.name=NULL,
	posterior_opts = NULL, use_sparse_matrices = FALSE
	){
	
	#### parameter setup ####
	set.seed(seed)
	diag(Y) <- NA
	if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
	if(family=="tobit") { Y[Y<0]<-0 }

	if(is.element(family,c("cbin","frn","rrl"))) {
		odobs<-apply(Y>0,1,sum,na.rm=TRUE) 
		if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y)) }
	}
	
	n <- nrow(Y)
	if(is.null(rownames(Y))) {
		rownames(Y) <- paste0("node", 1:n)
	}
	if(is.null(colnames(Y))) {
		colnames(Y) <- paste0("node", 1:n)
	}
	
	if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }
	
	if(is.na(g)) {
		if(family=="normal") { 
			g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)
		} else if(family=="tobit") {
			g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)*4
		} else {
			g <- sum(!is.na(Y))
		}
	}
	
	if(!is.list(prior)) { prior <- list() }
	vscale <- 1
	vdfmlt <- 1
	if(!(family %in% c("binary", "ordinal", "cbin", "frn", "rrl"))) {
		if(is.null(prior$Sab0)) {
			prior$Sab0 <- diag(2)
		}
		if(is.null(prior$Suv0) && R > 0) {
			prior$Suv0 <- diag(2*R)
		}
		if(is.null(prior$kappa0) && R > 0) {
			prior$kappa0 <- 2*R + 2
		}
		if(is.null(prior$eta0)) {
			prior$eta0 <- round(4 + 3*n/100)
		}
	}
	
	#### design matrix ####
	pr<-length(Xrow)/n
	pc<-length(Xcol)/n
	pd<-length(Xdyad)/n^2

	X <- design_array(Xrow, Xcol, Xdyad, intercept, n)
	
	Xnames <- dimnames(X)[[3]]
	if(is.null(Xnames) || length(Xnames) == 0) {
		if(dim(X)[3] > 0) {
			Xnames <- paste0("X", 1:dim(X)[3])
		} else {
			Xnames <- character(0)
		}
	}
	
	if( family=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0) 
			& !any( apply(X,3,function(x){var(c(x))})==0) ) {
		cli::cli_warn("Row effects are not estimable using this procedure")
	} 
	
	if( is.element(family,c("ord","rrl")) & 
			any( apply(X,3,function(x){var(c(x))})==0 ) ) {
		cli::cli_warn("An intercept is not estimable using this procedure")
	} 
	
	#### Z initialization ####
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
	
	#### MCMC initialization ####
	mu<-mean(Z,na.rm=TRUE)
	a<-if(rvar) rowMeans(Z,na.rm=TRUE) else rep(0, nrow(Z))
	b<-if(cvar) colMeans(Z,na.rm=TRUE) else rep(0, ncol(Z))
	a[is.na(a)]<-0 ; b[is.na(b)]<-0 
	ZA<-mu + outer(a,b,"+") 
	Z[is.na(Z)]<-ZA[is.na(Z)] 
	
	if(!is.null(start_vals)) {
		if(!is.null(start_vals$beta)) beta <- start_vals$beta else beta <- rep(0,dim(X)[3])
		if(!is.null(start_vals$s2)) s2 <- start_vals$s2 else s2 <- 1
		if(!is.null(start_vals$a)) a <- start_vals$a 
		if(!is.null(start_vals$b)) b <- start_vals$b
		if(!is.null(start_vals$rho)) rho <- start_vals$rho else rho <- 0
		if(!is.null(start_vals$Sab)) Sab <- start_vals$Sab
		if(!is.null(start_vals$U)) U <- start_vals$U else U <- matrix(0, nrow(Y), R)
		if(!is.null(start_vals$V)) V <- start_vals$V else V <- matrix(0, nrow(Y), R)
		if(!is.null(start_vals$Z)) Z <- start_vals$Z
	} else {
		beta<-rep(0,dim(X)[3])
		if(dim(X)[3]>0 && !all(is.na(Z))) {
			if(is.null(attributes(X)$XX)) {
				X <- precomputeX(X)
			}
			tryCatch({
				beta <- solve(attributes(X)$XX + diag(dim(X)[3])) %*% 
								crossprod(attributes(X)$mX, c(Z))
				beta <- as.vector(beta)
			}, error = function(e) {
				beta <- rep(0, dim(X)[3])
			})
		}
		
		s2<-1
		
		if(dim(X)[3] > 0) {
			E_tmp <- Z - Xbeta(X, beta)
			a <- rowMeans(E_tmp, na.rm=TRUE) * rvar
			b <- colMeans(E_tmp, na.rm=TRUE) * cvar
			a[is.na(a)] <- 0
			b[is.na(b)] <- 0
			
			if(family %in% c("binary", "ordinal", "cbin", "frn", "rrl")) {
				if(family == "binary") {
					ydist <- table(Y)
					ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
					YB <- 1 * (Y != ymode)
					ybar <- mean(YB, na.rm=TRUE)
					mu_bin <- qnorm(ybar)
					E <- (YB - ybar) / dnorm(qnorm(ybar))
					diag(E) <- 0
					a_tmp <- rowMeans(E, na.rm=TRUE)
					b_tmp <- colMeans(E, na.rm=TRUE)
					a_tmp[is.na(a_tmp)] <- 0
					b_tmp[is.na(b_tmp)] <- 0
					vscale <- mean(diag(cov(cbind(a_tmp, b_tmp))))
					PHAT <- pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
					vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
				} else {
					vscale <- mean(diag(cov(cbind(a, b))))
					PHAT <- pnorm(mu + outer(a, b, "+"))
					vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
				}
				
				if(is.null(prior$Sab0)) {
					prior$Sab0 <- diag(2) * vscale
				}
				if(is.null(prior$Suv0) && R > 0) {
					prior$Suv0 <- diag(2*R) * vscale
				}
				if(is.null(prior$eta0)) {
					prior$eta0 <- round(4 * vdfmlt)
				}
				if(is.null(prior$kappa0) && R > 0) {
					prior$kappa0 <- round((2*R + 2) * vdfmlt)
				}
			}
		}
		
		E_init <- Z - Xbeta(X, beta) - mu - outer(a, b, "+")
		if(dcor && !all(is.na(E_init))) {
			rho <- cor(c(E_init[upper.tri(E_init)]), c(t(E_init)[upper.tri(E_init)]), 
								 use="complete.obs")
			if(is.na(rho)) rho <- 0
		} else {
			rho <- 0
		}
		U<-V<-matrix(0, nrow(Y), R)
		if(R > 0) {
			E_init <- Z - mu - outer(a, b, "+")
			if(!all(is.na(E_init)) && any(E_init != 0, na.rm=TRUE)) {
				E_init[is.na(E_init)] <- 0
				tryCatch({
					svd_E <- svd(E_init)
					if(length(svd_E$d) >= R) {
						scaling_factor <- 0.5
						U <- svd_E$u[, 1:R, drop=FALSE] %*% diag(sqrt(svd_E$d[1:R] * scaling_factor), nrow=R)
						V <- svd_E$v[, 1:R, drop=FALSE] %*% diag(sqrt(svd_E$d[1:R] * scaling_factor), nrow=R)
					}
				}, error = function(e) {
					U <- matrix(0, nrow(Y), R)
					V <- matrix(0, nrow(Y), R)
				})
			}
		}
	}
	
	if(!exists("Sab") || is.null(Sab)) {
		if(rvar && cvar) {
			Sab<-tryCatch({
				S <- cov(cbind(a,b)) * tcrossprod(c(rvar, cvar))
				if(any(is.na(S)) || det(S) <= 0) {
					diag(c(max(0.1, var(a, na.rm=TRUE)), max(0.1, var(b, na.rm=TRUE))))
				} else {
					S
				}
			}, error = function(e) {
				diag(c(0.1, 0.1))
			})
		} else if(rvar && !cvar) {
			Sab<-matrix(c(max(0.1, var(a)), 0, 0, 0.1), 2, 2)
		} else if(!rvar && cvar) {
			Sab<-matrix(c(0.1, 0, 0, max(0.1, var(b))), 2, 2)
		} else {
			Sab<-diag(c(0.1, 0.1))
		}
	}
	
	#### storage allocation ####
	n_samples <- floor((nscan + burn * (burn > 0)) / odens) - floor(burn / odens)
	BETA <- matrix(NA_real_, nrow = n_samples, ncol = dim(X)[3] - pr*symmetric)
	VC <- matrix(NA_real_, nrow = n_samples, ncol = 5 - 3*symmetric)
	sample_idx <- 0
	U_samples <- NULL ; V_samples <- NULL
	a_samples <- NULL ; b_samples <- NULL
	
	if(!is.null(posterior_opts)) {
		if(posterior_opts$save_UV && R > 0) {
			n_UV_samples <- ceiling(n_samples / posterior_opts$thin_UV)
			U_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
			if(!symmetric) {
				V_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
			}
		}
		if(posterior_opts$save_ab && (rvar || cvar)) {
			n_ab_samples <- ceiling(n_samples / posterior_opts$thin_ab)
			if(rvar) a_samples <- matrix(NA_real_, n, n_ab_samples)
			if(cvar && !symmetric) b_samples <- matrix(NA_real_, n, n_ab_samples)
		}
	}
	
	UV_sample_idx <- 0
	ab_sample_idx <- 0
	
	UVPS <- U %*% t(V) * 0
	APS<-BPS<- rep(0,nrow(Y))
	YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
	
	if(gof) {
		n_base_stats <- 5
		custom_gof_names <- NULL
		n_custom_stats <- 0
		if(!is.null(custom_gof)) {
			if(is.function(custom_gof)) {
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
				n_custom_stats <- length(custom_gof)
				custom_gof_names <- names(custom_gof)
				if(is.null(custom_gof_names)) {
					custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
				}
			}
		}

		n_total_stats <- n_base_stats + n_custom_stats
		GOF <- matrix(NA_real_, nrow = n_samples + 1, ncol = n_total_stats)
		
		base_gof <- gof_stats(Y, mode = "unipartite")
		GOF[1, 1:n_base_stats] <- base_gof
		
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
		base_names <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		all_names <- c(base_names, custom_gof_names)
		colnames(GOF) <- all_names
	} else {
		GOF <- matrix(gof_stats(Y, mode = "unipartite"), 1, 5)
		rownames(GOF) <- "obs"
		colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		custom_gof <- NULL
	}
	
	#### pre-MCMC setup ####
	names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y)

	if(!symmetric) {
		colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
		colnames(BETA) <- dimnames(X)[[3]] 
	}
	
	if(symmetric) {
		colnames(VC) <- c("va", "ve")  
		rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
		bnames<-dimnames(X)[[3]]
		bni<-bnames[1*intercept] 
		bnn<-gsub("row",bnames[rb],replacement="node")       
		bnd<-bnames[-c(1*intercept,rb,cb)]
		colnames(BETA)<-c(bni,bnn,bnd) 
	}    
	
	Xr<-apply(X,c(1,3),sum)
	Xc<-apply(X,c(2,3),sum)
	mX<- apply(X,3,c)
	mXt<-apply(aperm(X,c(2,1,3)),3,c)
	XX<-t(mX)%*%mX
	XXt<-t(mX)%*%mXt

	X_precomp <- X
	attributes(X_precomp) <- c(
		attributes(X), list(Xr=Xr, Xc=Xc, mX=mX, mXt=mXt, XX=XX, XXt=XXt))
	
	have_coda<-suppressWarnings(
		try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
	
	savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
	savePoints <- savePoints[
		round(quantile(1:length(savePoints), 
		probs=seq(save_interval,1,save_interval))) ]
	
	beta_running_mean <- rep(0, dim(X)[3] - pr*symmetric)
	vc_running_mean <- rep(0, 5 - 3*symmetric)
	running_sample_count <- 0
	
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
	
	needs_s2 <- family %in% c("normal", "tobit", "poisson")
	
	Xbeta_cache <- if(dim(X)[3] > 0) Xbeta(X, beta) else matrix(0, n, n)
	ab_outer_cache <- outer(a, b, "+")
	UV_cache <- U %*% t(V)
	EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
	
	if(verbose) {
		cli::cli_h2("Running MCMC for unipartite network")
		cli::cli_text("Dimensions: {.val {n}} x {.val {n}} nodes")
		cli::cli_text("R = {.val {R}}, symmetric = {.val {symmetric}}")
		cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")
		
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
		
		start_time <- Sys.time()
		
		if(burn > 0) {
			cli::cli_progress_bar("Burn-in", total = burn, clear = TRUE)
		}
	}
	
	#### MCMC loop ####
	for (s in 1:(nscan + burn)) {
		
		EZ <- EZ_cache
		rho_eff <- rho
		Z <- if(family %in% c("normal", "tobit", "poisson")) {
			update_Z_fn(Z, EZ, rho_eff, s2, Y)
		} else {
			update_Z_fn(Z, EZ, rho_eff, Y)
		}
		
		if (needs_s2) {
			s2<-rs2_fc(Z,rho_eff,offset=EZ)  
		}
		
		if(dcor) {
			rho<-rrho_mh(Z, rho, s2, offset=EZ)
		}
		
		if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }
		
		tmp <- rbeta_ab_fc(Z, Sab, rho, X_precomp, s2, offset=UV_cache, g=g)
		beta_new <- tmp$beta
		a_new <- tmp$a * rvar
		b_new <- tmp$b * cvar
		
		if(symmetric){ a_new<-b_new<-(a_new+b_new)/2 }
		
		if(dim(X)[3] > 0 && !isTRUE(all.equal(beta, beta_new))) {
			Xbeta_cache <- Xbeta(X, beta_new)
		}
		if(!isTRUE(all.equal(a, a_new)) || !isTRUE(all.equal(b, b_new))) {
			ab_outer_cache <- outer(a_new, b_new, "+")
		}
		
		beta <- beta_new
		a <- a_new
		b <- b_new
		
		EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
		
		if(R > 0) {
			shrink<- (s>.5*burn)
			offset <- Xbeta_cache + ab_outer_cache
			
			if(symmetric){ 
				E<-Z-offset ; E<-.5*(E+t(E))
				UV<-rUV_sym_fc(E, U, V, s2,shrink) 
			}
			if(!symmetric){
				E <- Z - offset
				
				Suv0_use <- if(!is.null(prior$Suv0)) prior$Suv0 else diag(2*R)
				kappa0_use <- if(!is.null(prior$kappa0)) prior$kappa0 else (2*R + 2)
				Suv <- rSuv_fc(U, V, Suv0=Suv0_use/kappa0_use, kappa0=kappa0_use)
				
				UV<-rUV_fc(E, U, V, Suv, rho, s2, offset=0) 
			}
			
			U<-UV$U ; V<-UV$V
			UV_cache <- U %*% t(V)
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
		}
		
		if(is.element(family,c("normal","tobit","ordinal"))) { 
			Sab <- rSab_fc(a, b, Sab0=prior$Sab0, eta0=prior$eta0, 
										 rvar=rvar, cvar=cvar, symmetric=symmetric)
		}
		
		if(family=="binary") {
			if(rvar & cvar & !symmetric) {
				Sab <- rSab_fc(a, b, Sab0=prior$Sab0, eta0=prior$eta0)
				tmp<-raSab_bin_fc(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
				Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(a, b, Sab0=prior$Sab0, eta0=prior$eta0, 
											 rvar=rvar, cvar=cvar, symmetric=symmetric)
			}
		}
		
		if(family=="cbin") {
			if(rvar & cvar & !symmetric) {
				tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,Sab0=prior$Sab0,
													eta0=prior$eta0)
				Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(a, b, Sab0=prior$Sab0, eta0=prior$eta0, 
											 rvar=rvar, cvar=cvar, symmetric=symmetric)
			}
		}
		
		if(family=="frn") { 
			if(rvar & cvar & !symmetric) {
				tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs,Sab0=prior$Sab0,
												 eta0=prior$eta0)
				Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(a, b, Sab0=prior$Sab0, eta0=prior$eta0, 
											 rvar=rvar, cvar=cvar, symmetric=symmetric)
			}
		}
		
		
		if(verbose) {
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
		
		if(s%%odens==0 & s>burn) {  
			sample_idx <- sample_idx + 1
			
			if(symmetric) {
				br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2 
				sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
				BETA[sample_idx,] <- sbeta
				
				VC[sample_idx,] <- c(Sab[1,1],s2)
				
				beta_running_mean <- beta_running_mean + (sbeta - beta_running_mean) / sample_idx
				vc_running_mean <- vc_running_mean + (c(Sab[1,1],s2) - vc_running_mean) / sample_idx
			}

			if(!symmetric) {
				BETA[sample_idx,] <- beta
				VC[sample_idx,] <- c(Sab[upper.tri(Sab, diag = T)], rho, s2)
				
				beta_running_mean <- beta_running_mean + (beta - beta_running_mean) / sample_idx
				vc_running_mean <- vc_running_mean + (c(Sab[upper.tri(Sab, diag = T)], rho, s2) - vc_running_mean) / sample_idx
			}
			
			if(!is.null(posterior_opts)) {
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
			
			UVPS <- UVPS + UV_cache
			APS <- APS + a
			BPS <- BPS + b 
			
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
			
			EZ <- EZ_cache
			if(symmetric){ EZ<-(EZ+t(EZ))/2 }
			
			rho_eff <- rho
			Ys <- simulate_Y_fn(EZ, rho_eff, s2)
			
			if(symmetric){ Ys[lower.tri(Ys)]<-0 ; Ys<-Ys+t(Ys)  }
			
			YPS<-YPS+Ys
			
			if(gof){
				Ys[is.na(Y)]<-NA
				gof_idx <- gof_idx + 1
				GOF[gof_idx, 1:5] <- gof_stats(Ys, mode = "unipartite")
				
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
			
		}
		
		if(periodic_save & s %in% savePoints & !is.null(out_file)){
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
	
	#### output assembly ####
	if(verbose) {
		cli::cli_progress_done()

		if(exists("start_time")) {
			end_time <- Sys.time()
			elapsed <- difftime(end_time, start_time, units = "secs")
			cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
		}
	}
	
	# trim matrices to actual sample size
	if(exists("sample_idx") && sample_idx > 0) {
		BETA <- BETA[1:sample_idx, , drop=FALSE]
		VC <- VC[1:sample_idx, , drop=FALSE]
		if(gof && exists("gof_idx")) {
			GOF <- GOF[1:gof_idx, , drop=FALSE]
		}
	}
	
	APM<-APS/nrow(VC)
	BPM<-BPS/nrow(VC)
	if(!is.null(rownames(Y))) {
		names(APM) <- rownames(Y)
		names(BPM) <- rownames(Y)
	} else {
		names(APM) <- paste0("Actor", 1:length(APM))
		names(BPM) <- paste0("Actor", 1:length(BPM))
	}
	UVPM<-UVPS/nrow(VC)
	YPM<-YPS/nrow(VC)
	names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
	
	dimnames(YPM)<-dimnames(Y)
	rownames(BETA)<-NULL
	
	if(!symmetric) {
		if(R > 0) {
			UDV<-svd(UVPM)
			U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
			V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
			rownames(U)<-rownames(V)<-rownames(Y)
		} else {
			U <- NULL
			V <- NULL
		}
		start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
		fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,L=NULL,
								YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
								mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R,
								n=n,p=dim(X)[3],rvar=rvar,cvar=cvar,dcor=dcor,
								X_names=colnames(BETA),actor_names=rownames(Y))
		
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

	if(symmetric) {
		ULUPM<-UVPM
		if(R > 0) {
			eULU<-eigen(ULUPM)
			eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
			U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
			L<-diag(eULU$val[eR], nrow=R)
			rownames(U)<-rownames(ULUPM)
		} else {
			U <- NULL
			L <- NULL
		}
		YPM<-.5*(YPM+t(YPM))
		start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=U, rho=rho, s2=s2, Sab=Sab)
		fit<-list(BETA=BETA,VC=VC,APM=APM,BPM=NULL,U=U,V=NULL,L=L,
							YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
							mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R,
							call=match.call(),n=n,p=dim(X)[3],rvar=rvar,cvar=cvar,dcor=dcor,
							X_names=colnames(BETA),actor_names=rownames(Y))

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

	if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
	if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
	
	class(fit) <- "ame"
	
	if(use_sparse_matrices) {
		fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
	}
	
	fit
}