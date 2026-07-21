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
	rvar = TRUE, cvar = TRUE, dcor = !symmetric,
	nvar=TRUE, R = 0, family="normal",
	intercept=!(family == "ordinal"),
	symmetric=FALSE,
	odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
	prior=list(), g=NA,
	seed = 6886, nscan = 10000, burn = 500, odens = 25,
	verbose = TRUE, gof=TRUE, custom_gof=NULL,
	start_vals=NULL, periodic_save=FALSE, out_file=NULL,
	save_interval=0.25, model.name=NULL,
	posterior_opts = NULL, use_sparse_matrices = FALSE,
	ordinal_cutpoints = c("data_induced", "explicit")
	){
	
	#### parameter setup ####
	# seed locally: restore the global rng stream on exit so a downstream
	# random draw is not silently perturbed by having fit a model
	if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		.old_seed <- get(".Random.seed", envir = globalenv())
		on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
		        add = TRUE)
	}
	set.seed(seed)
	diag(Y) <- NA
	if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }

	if(is.element(family,c("cbin","frn"))) {
		odobs<-apply(Y>0,1,sum,na.rm=TRUE)
		# reject na / <=0 odmax with observed nominations
		if (!is.null(odmax)) {
			if (anyNA(odmax)) {
				cli::cli_abort(c(
					"{.arg odmax} contains {.val NA}.",
					"i" = "Supply a positive integer (scalar or per-row vector), or pass {.code odmax = NULL} to infer from {.arg Y}."))
			}
			if (any(odmax <= 0L) && any(Y > 0, na.rm = TRUE)) {
				cli::cli_abort(c(
					"{.arg odmax} <= 0 but {.arg Y} contains nominations.",
					"i" = "{.code odmax = 0} prohibits any nominations; inconsistent with the observed data."))
			}
		}
		# default odmax to the max observed outdegree when the user passed
		# NULL, matching lame() and the bipartite path (and the behaviour
		# the odmax = NA error message above advertises)
		if (is.null(odmax)) {
			odmax <- max(odobs, na.rm = TRUE)
			cli::cli_inform(c(
				"i" = "{.arg odmax} not supplied for family {.val {family}}; using the maximum observed outdegree ({odmax})."))
		}
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
	
	# `g` is the scalar zellner g-prior scale. a non-scalar g doesn't match
	# what rbeta_ab_fc() expects (xx/g with a length-p vector miscomputes
	# the inverse-prior matrix); reject it explicitly rather than silently
	# producing the wrong posterior.
	if (length(g) > 1L) {
		cli::cli_abort(c(
			"{.arg g} must be a single scalar (length {length(g)} supplied).",
			"i" = "Per-coefficient g-priors are not currently supported in the unipartite path."))
	}
	if (length(g) == 1L && is.na(g)) {
		if(family=="normal") {
			g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)
		} else {
			g <- sum(!is.na(Y))
		}
	}
	
	if(!is.list(prior)) {
		if (!is.null(prior)) {
			cli::cli_warn(c(
				"{.arg prior} must be a named list; got {.cls {class(prior)[1]}}.",
				"i" = "Ignoring it and fitting with the default priors."))
		}
		prior <- list()
	}
	vscale <- 1
	vdfmlt <- 1
	if(!(family %in% c("binary", "ordinal", "cbin", "frn"))) {
		# scale the additive/multiplicative-effects priors to the data,
		# matching the adaptive scaling the binary path already applies:
		# an absolute diag(2) prior biases va/vb toward 1 at small n
		# (~+50% at n=30 when the true variance is 0.3) and breaks scale
		# equivariance of the intercept/additive-effect split
		if(is.null(prior$Sab0) || (is.null(prior$Suv0) && R > 0)) {
			Ymom <- if (family == "poisson") log1p(pmax(Y, 0)) else Y
			Ec <- Ymom - mean(Ymom, na.rm = TRUE)
			a_tmp <- rowMeans(Ec, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
			b_tmp <- colMeans(Ec, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
			vscale_ab <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
			if (!is.finite(vscale_ab) || vscale_ab <= 0) {
				vscale_ab <- max(stats::var(c(Ymom), na.rm = TRUE), 1e-12)
			}
		}
		if(is.null(prior$Sab0)) {
			prior$Sab0 <- diag(2) * vscale_ab
		}
		if(is.null(prior$Suv0) && R > 0) {
			prior$Suv0 <- diag(2*R) * vscale_ab
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

	# propagate covariate missingness into y: a dyad with a missing covariate
	# is treated as an unobserved tie (data-augmented), instead of silently
	# imputing the covariate as 0 and biasing its coefficient. matches lame().
	if ((pr + pc + pd) > 0) {
		n_na_before <- sum(is.na(Y))
		Y <- .ame_propagate_cov_na(Y, Xrow, Xcol, Xdyad)
		n_cov_na <- sum(is.na(Y)) - n_na_before
		if (n_cov_na > 0 && getOption("lame.warn.na", TRUE)) {
			cli::cli_inform(c(
				"i" = "{n_cov_na} {.arg Y} cell{?s} had a missing covariate and {?is/are} treated as unobserved (data-augmented).",
				"i" = "The covariate coefficients are estimated from the complete dyads only."))
		}
	}

	X <- design_array(Xrow, Xcol, Xdyad, intercept, n, warn = FALSE)
	
	Xnames <- dimnames(X)[[3]]
	if(is.null(Xnames) || length(Xnames) == 0) {
		if(dim(X)[3] > 0) {
			Xnames <- paste0("X", 1:dim(X)[3])
		} else {
			Xnames <- character(0)
		}
	}
	
	if( family=="ordinal" &&
			isTRUE(any( apply(X,3,function(x){var(c(x))})==0 )) ) {
		cli::cli_warn("An intercept is not estimable using this procedure")
	}

	#### z initialization ####
	if(family=="frn") {
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
	if(family=="ordinal") { Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) }
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
	
	#### mcmc initialization ####
	mu<-mean(Z,na.rm=TRUE)
	a<-if(rvar) rowMeans(Z,na.rm=TRUE) else rep(0, nrow(Z))
	b<-if(cvar) colMeans(Z,na.rm=TRUE) else rep(0, ncol(Z))
	a[is.na(a)]<-0 ; b[is.na(b)]<-0 
	ZA<-mu + outer(a,b,"+") 
	Z[is.na(Z)]<-ZA[is.na(Z)] 
	
	# accumulator for in-loop sampler failures (parity with lame()). exposed on
	# the fit object as fit$tryerrorchecks so the user can audit silent fallbacks.
	tryErrorChecks <- list(beta = 0L, Sab = 0L, s2 = 0L, rho = 0L, UV = 0L)
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
			
			if(family %in% c("binary", "ordinal", "cbin", "frn")) {
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
			# an exactly symmetric Y gives rho = 1, which makes the dyadic
			# covariance singular and crashes the samplers; keep the start
			# strictly inside (-1, 1)
			rho <- min(max(rho, -0.9), 0.9)
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
	# normalize a partially-specified posterior_opts so bare
	# posterior_opts$save_UV / $save_ab reads below cannot evaluate to NULL
	# (which crashes `NULL && ...` in the allocation and storage guards).
	if(!is.null(posterior_opts)) {
		posterior_opts$save_UV <- isTRUE(posterior_opts$save_UV)
		posterior_opts$save_ab <- isTRUE(posterior_opts$save_ab)
		if(is.null(posterior_opts$thin_UV)) posterior_opts$thin_UV <- 1
		if(is.null(posterior_opts$thin_ab)) posterior_opts$thin_ab <- 1
	}
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
			# symmetric fits store V = U L per draw as well, so the latent
			# similarity U L U' is reconstructable draw-by-draw
			V_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
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
		# tracker for silent custom_gof failures: count + first message,
		# surfaced as one warning at the end of mcmc (not per iteration).
		.cg_tracker <- .new_custom_gof_err_tracker()
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
					cli::cli_warn(c(
						"!" = "{.arg custom_gof} failed on the initial test call; the function will be skipped for this fit.",
						"i" = "Error: {conditionMessage(e)}",
						"i" = "Check that {.arg custom_gof} accepts a single matrix {.arg Y} and returns a numeric scalar / vector."))
					custom_gof <<- NULL
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
				}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
			} else if(is.list(custom_gof)) {
				for(i in seq_along(custom_gof)) {
					if(is.function(custom_gof[[i]])) {
						tryCatch({
							GOF[1, n_base_stats + i] <- custom_gof[[i]](Y)
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
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
	
	#### pre-mcmc setup ####
	names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y)

	if(!symmetric) {
		colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
		colnames(BETA) <- dimnames(X)[[3]] 
	}
	
	if(symmetric) {
		colnames(VC) <- c("va", "ve")
		rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
		bnames<-dimnames(X)[[3]]
		# guard against `-c(0, integer(0))` short-circuit: when intercept = false
		# and pr = 0, the negative index becomes `-0` and silently empties the
		# vector. build the exclude index explicitly and skip the [-excl] step
		# when there is nothing to drop.
		excl_idx <- c(if (intercept) 1L else integer(0), rb, cb)
		bni <- if (intercept) bnames[1L] else character(0)
		bnn <- if (length(rb) > 0L) gsub("row", bnames[rb], replacement = "node") else character(0)
		bnd <- if (length(excl_idx) > 0L) bnames[-excl_idx] else bnames
		colnames(BETA) <- c(bni, bnn, bnd)
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
	# sample storage is allocated once above (the U_samples/V_samples/
	# a_samples/b_samples block that guards on symmetric and per-effect
	# NULLs, with UV_sample_idx/ab_sample_idx zeroed there).

	# explicit-cutpoint state init (ordinal_cutpoints = "explicit"). when
	# active we override the ordinal z sampler to the explicit-alpha path
	# (rz_ord_(sym_)explicit_fc) and run a cowles (1996) mh on alpha after
	# each z sweep.
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	use_explicit_cutpoints <- (family == "ordinal" &&
	                            identical(ordinal_cutpoints, "explicit"))
	if (use_explicit_cutpoints) {
		ord_lvls <- sort(unique(c(Y[is.finite(Y)])))
		K_ord <- length(ord_lvls)
		if (K_ord < 2L) {
			cli::cli_abort(c(
				"{.arg ordinal_cutpoints} = {.val explicit} requires at least 2 distinct observed categories.",
				"i" = "Found {K_ord}; check {.arg Y} or use a different family."))
		}
		Y_int <- matrix(match(Y, ord_lvls), nrow(Y), ncol(Y),
		                dimnames = dimnames(Y))
		alpha <- .init_alpha_from_data(Y)
		tau_prop <- 0.5
		tau_target <- 0.30
		alpha_accept_n <- 0L
		alpha_accept_total <- 0L
		ALPHA <- if (K_ord >= 3L) {
			matrix(NA_real_, nrow = nscan/odens, ncol = K_ord - 2L,
			       dimnames = list(NULL, paste0("alpha_", seq.int(2L, K_ord - 1L))))
		} else NULL
	} else {
		alpha <- numeric(0)
		Y_int <- NULL
		tau_prop <- NA_real_
		K_ord <- 0L
		ALPHA <- NULL
	}

	# symmetric ordinal uses a dedicated sampler that draws z on the upper
	# triangle from a tn with the doubled precision (variance 1/2, mean =
	# (ez_ij + ez_ji)/2) and mirrors to the lower triangle. see r/rz_ord_sym_fc.r.
	# when ordinal_cutpoints == "explicit", the ordinal slot is left as the
	# default sampler; the mcmc loop overrides it inline via y_int + alpha.
	update_Z_fn <- switch(family,
		normal = rZ_nrm_fc,
		binary = function(Z, EZ, rho, Y) {
			# fused c++ kernel (same full conditional and sweep order as
			# rZ_bin_fc); fall back to the r sampler on error
			z_new <- tryCatch(rZ_bin_fused_cpp(Z, EZ, rho, Y),
				error = function(e) NULL)
			if (is.null(z_new)) return(rZ_bin_fc(Z, EZ, rho, Y))
			dimnames(z_new) <- dimnames(Z)
			z_new
		},
		ordinal = if (symmetric) function(Z, EZ, rho, Y) rZ_ord_sym_fc(Z, EZ, Y) else rZ_ord_fc,
		cbin = function(Z, EZ, rho, Y) rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs),
		frn = function(Z, EZ, rho, Y) rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs),
		poisson = rZ_pois_fc,
		stop(paste("Unknown family:", family))
	)
	
	simulate_Y_fn <- switch(family,
		normal = function(EZ, rho, s2) simY_nrm(EZ, rho, s2),
		binary = function(EZ, rho, s2) simY_bin(EZ, rho),
		ordinal = function(EZ, rho, s2) simY_ord(EZ, rho, Y),
		cbin = function(EZ, rho, s2) 1*(simY_frn(EZ, rho, odmax, YO=Y)>0),
		frn = function(EZ, rho, s2) simY_frn(EZ, rho, odmax, YO=Y),
		# the fitted model is y ~ Pois(exp(z)), z ~ N(EZ, s2*Sigma_rho);
		# posterior-predictive draws must include the latent z layer or
		# GOF rejects correctly-specified overdispersed models and YPM is
		# biased low by exp(-s2/2)
		poisson = function(EZ, rho, s2) simY_pois(simZ(EZ, rho, s2)),
		stop(paste("Unknown family:", family))
	)
	
	needs_s2 <- family %in% c("normal", "poisson")
	
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
					}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
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
	
	#### mcmc loop ####
	for (s in 1:(nscan + burn)) {
		
		EZ <- EZ_cache
		rho_eff <- rho
		Z <- if (use_explicit_cutpoints) {
			# explicit-cutpoint ordinal: rz_ord_(sym_)explicit_fc takes
			# (z, ez, y_int, alpha) and ignores rho (handled via the
			# precision-2 / variance-1/2 convention internally).
			if (symmetric) {
				rZ_ord_sym_explicit_fc(Z, EZ, Y_int, alpha)
			} else {
				rZ_ord_explicit_fc(Z, EZ, Y_int, alpha)
			}
		} else if(family %in% c("normal", "poisson")) {
			update_Z_fn(Z, EZ, rho_eff, s2, Y)
		} else {
			update_Z_fn(Z, EZ, rho_eff, Y)
		}

		# cowles (1996) mh update on alpha, z-marginalised. runs after the
		# z sweep so ez is current; robbins-monro on tau_prop only during
		# burn-in to preserve the stationary distribution.
		if (use_explicit_cutpoints) {
			alpha_step <- sample_alpha_cowles(
				alpha = alpha, Y_int = Y_int, EZ = EZ,
				tau_prop = tau_prop, symmetric = symmetric)
			alpha <- alpha_step$alpha
			alpha_accept_total <- alpha_accept_total + isTRUE(alpha_step$accept)
			alpha_accept_n <- alpha_accept_n + 1L
			if (s <= burn && alpha_accept_n >= 25L) {
				acc_rate <- alpha_accept_total / alpha_accept_n
				log_tau <- log(tau_prop) + (acc_rate - tau_target) / sqrt(s)
				tau_prop <- min(max(exp(log_tau), 0.01), 2.0)
				alpha_accept_n <- 0L
				alpha_accept_total <- 0L
			}
		}
		
		if (needs_s2) {
			if (symmetric) {
				# symmetric path: e_ij == e_ji by construction, so the
				# rs2_fc whitening (which decorrelates a length-2 dyad pair)
				# double-counts the degrees of freedom and gives s2 -> ve/2.
				# sample directly from the upper-triangle residual variance
				# with the same scale-free one-pseudo-observation prior as
				# rs2_fc (an absolute IG(2, 1) prior biases s2 toward 1
				# whenever Y is off unit scale).
				resid <- Z - EZ
				ut <- upper.tri(resid)
				e <- resid[ut]; e <- e[is.finite(e)]
				n_e <- length(e)
				if (n_e > 1) {
					ssr <- sum(e^2)
					nu0 <- 1
					s20 <- if (!is.null(prior$s20)) {
						prior$s20
					} else {
						ms <- ssr / n_e
						if (!is.finite(ms) || ms <= 0) 1 else ms
					}
					s2 <- 1 / rgamma(1, shape = (n_e + nu0) / 2,
					                 rate  = (ssr + nu0 * s20) / 2)
				}
			} else {
				s2 <- rs2_fc(Z, rho_eff, offset = EZ, s20 = prior$s20)
			}
		}

		if(dcor) {
			# rho is the slowest-mixing parameter (random-walk MH); applying
			# the kernel several times per scan is cheap relative to the Z
			# and beta updates and substantially improves its ESS
			for(k_rho in 1:4) {
				rho<-rrho_mh(Z, rho, s2, offset=EZ)
			}
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
				
				# rSuv_fc forms the prior scale internally as kappa0 * Suv0,
				# so Suv0 is passed through unscaled
				Suv0_use <- if(!is.null(prior$Suv0)) prior$Suv0 else diag(2*R)
				kappa0_use <- if(!is.null(prior$kappa0)) prior$kappa0 else (2*R + 2)
				Suv <- rSuv_fc(U, V, Suv0=Suv0_use, kappa0=kappa0_use)
				
				UV<-rUV_fc(E, U, V, Suv, rho, s2, offset=0) 
			}
			
			U<-UV$U ; V<-UV$V
			UV_cache <- U %*% t(V)
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
		}
		
		# normal/ordinal/poisson sample Sab via rSab_fc every scan.
		# binary/cbin/frn use the raSab_* joint (a,b,Sab) updates below
		# when rvar & cvar & !symmetric (binary refreshes Sab via rSab_fc
		# before its joint update); the other variance configurations fall
		# back to rSab_fc alone.
		if(is.element(family,c("normal","ordinal","poisson"))) {
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
				# same `-c(0, integer(0))` guard as in the colnames(beta) build
				excl_idx <- c(if (intercept) 1L else integer(0), rb, cb)
				bint <- if (intercept) beta[1L] else numeric(0)
				bdyad <- if (length(excl_idx) > 0L) beta[-excl_idx] else beta
				sbeta <- c(bint, bn, bdyad)
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

			# explicit-cutpoint alpha storage (alpha_2..alpha_{k-1}; alpha_1 = 0
			# is the identification anchor, not stored). k = 2 has zero free
			# cutpoints, in which case alpha is null and we skip the store.
			if (use_explicit_cutpoints && !is.null(ALPHA)) {
				ALPHA[sample_idx, ] <- alpha[-1L]
			}
			
			if(!is.null(posterior_opts)) {
				if(posterior_opts$save_UV && R > 0) {
					if(sample_idx %% posterior_opts$thin_UV == 0) {
						UV_sample_idx <- UV_sample_idx + 1
						if(UV_sample_idx <= dim(U_samples)[3]) {
							U_samples[,,UV_sample_idx] <- U
							if(!is.null(V_samples)) {
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
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
					} else if(is.list(custom_gof)) {
						for(i in seq_along(custom_gof)) {
							if(is.function(custom_gof[[i]])) {
								tryCatch({
									GOF[gof_idx, 5 + i] <- custom_gof[[i]](Ys)
								}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
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
		
	} # end mcmc   
	
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
		if (use_explicit_cutpoints && !is.null(ALPHA)) {
			ALPHA <- ALPHA[1:sample_idx, , drop=FALSE]
		}
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
		# stash the resolved scalar g and the prior list (with sampler-side
		# defaults filled in) so prior_summary() can show what was actually used
		.resolved_prior <- prior
		fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,L=NULL,UVPM=UVPM,
								YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
								mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R,
								n=n,p=dim(X)[3],rvar=rvar,cvar=cvar,dcor=dcor,
								g=g, prior=.resolved_prior,
								tryErrorChecks=tryErrorChecks,
								mh_counters=tryErrorChecks,
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
			# keep the R eigenvalues largest in magnitude and pair each with its
			# own eigenvector, so U L U' is the rank-R truncation of ULUPM
			eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
			U<-eULU$vec[,eR,drop=FALSE]
			L<-diag(eULU$val[eR], nrow=R)
			rownames(U)<-rownames(ULUPM)
		} else {
			U <- NULL
			L <- NULL
		}
		YPM<-.5*(YPM+t(YPM))
		start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=U, rho=rho, s2=s2, Sab=Sab)
		fit<-list(BETA=BETA,VC=VC,APM=APM,BPM=NULL,U=U,V=NULL,L=L,ULUPM=ULUPM,
							YPM=YPM,GOF=GOF,X=X,Y=Y,start_vals=start_vals_final,model.name=model.name,
							mode="unipartite",family=family,symmetric=symmetric,odmax=odmax,R=R,
							call=match.call(),n=n,p=dim(X)[3],rvar=rvar,cvar=cvar,dcor=dcor,
							g=g, prior=prior,                     # parity with the asymmetric path
							tryErrorChecks=tryErrorChecks,
							mh_counters=tryErrorChecks,
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

	# attach explicit-cutpoint posterior + bookkeeping. ordinal_cutpoints
	# is stored unconditionally for downstream methods (summary/predict)
	# so they can branch off it. alpha_post_mean prepends the anchor 0.
	fit$ordinal_cutpoints <- ordinal_cutpoints
	if (use_explicit_cutpoints) {
		fit$ordinal_levels <- ord_lvls
		if (!is.null(ALPHA)) {
			fit$ALPHA <- ALPHA
			fit$alpha_post_mean <- c(0, colMeans(ALPHA, na.rm = TRUE))
		} else {
			fit$alpha_post_mean <- 0
		}
	}

	# emit a single post-mcmc warning if the user-supplied custom_gof
	# function raised any errors during the run. without this the gof
	# matrix would carry na entries with no indication of why.
	if (exists(".cg_tracker", inherits = FALSE))
		.maybe_warn_custom_gof(.cg_tracker, where = "GOF")

	class(fit) <- "ame"

	if(use_sparse_matrices) {
		fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
	}

	fit
}
