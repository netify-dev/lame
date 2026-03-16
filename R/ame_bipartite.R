#' AME model fitting for bipartite (rectangular) networks
#' 
#' Internal function that handles AME model fitting for bipartite networks.
#' This function contains the core logic for rectangular adjacency matrices
#' with distinct row and column node sets.
#' 
#' @inheritParams ame
#' @param R_row integer: dimension of row node multiplicative effects
#' @param R_col integer: dimension of column node multiplicative effects
#' @return An ame object with posterior samples and estimates
#' @keywords internal
#' @noRd
ame_bipartite <- function(
	Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
	rvar = !(family=="rrl"), cvar = TRUE, R_row = 0, R_col = 0,
	family="normal", intercept=!is.element(family,c("rrl","ordinal")),
	odmax=NULL, prior=list(), g=NA,
	seed = 6886, nscan = 10000, burn = 500, odens = 25,
	verbose = TRUE, gof=TRUE, custom_gof=NULL,
	start_vals=NULL, periodic_save=FALSE, out_file=NULL,
	save_interval=0.25, model.name=NULL,
	posterior_opts = NULL, use_sparse_matrices = FALSE
	){

	#### parameter setup ####
	if(family == "rrl") {
		warning("rrl (row-ranked likelihood) is not well-defined for bipartite networks; results may be unreliable")
	}

	set.seed(seed)
	nA <- nrow(Y)
	nB <- ncol(Y)
	if(is.null(rownames(Y))) {
		rownames(Y) <- paste0("rowNode", 1:nA)
	}
	if(is.null(colnames(Y))) {
		colnames(Y) <- paste0("colNode", 1:nB)
	}
	
	n_eff <- floor((nscan - burn) / odens)
	
	Sab0 <- prior$Sab0
	if(is.null(Sab0)) Sab0 <- diag(2)
	eta0 <- prior$eta0
	if(is.null(eta0)) eta0 <- 4 + 3 * (nA + nB) / 200
	etaab <- prior$etaab
	if(is.null(etaab)) etaab <- eta0
	s20 <- prior$s20
	if(is.null(s20)) s20 <- 1
	s2u0 <- prior$s2u0
	if(is.null(s2u0)) s2u0 <- 1
	Suv0 <- prior$Suv0
	if(is.null(Suv0)) Suv0 <- diag(max(R_row, R_col))
	
	#### covariate setup ####
	if(!is.null(Xdyad)) {
		Xdyad <- as.array(Xdyad)
		if(length(dim(Xdyad)) == 2) {
			Xdyad <- array(Xdyad, c(nA, nB, 1))
		}
		p_dyad <- dim(Xdyad)[3]
	} else {
		p_dyad <- 0
	}
	
	if(!is.null(Xrow)) {
		Xrow <- as.matrix(Xrow)
		p_row <- ncol(Xrow)
	} else {
		p_row <- 0
	}
	
	if(!is.null(Xcol)) {
		Xcol <- as.matrix(Xcol)
		p_col <- ncol(Xcol)
	} else {
		p_col <- 0
	}
	
	p <- p_dyad + p_row + p_col + as.numeric(intercept)

	# covariate names for BETA columns
	beta_names <- character(0)
	if(intercept) beta_names <- c(beta_names, "intercept")
	if(p_dyad > 0) {
		dyad_nms <- dimnames(Xdyad)[[3]]
		if(is.null(dyad_nms)) dyad_nms <- paste0("Xdyad", 1:p_dyad)
		beta_names <- c(beta_names, paste0(dyad_nms, "_dyad"))
	}
	if(p_row > 0) {
		row_nms <- colnames(Xrow)
		if(is.null(row_nms)) row_nms <- paste0("Xrow", 1:p_row)
		beta_names <- c(beta_names, paste0(row_nms, "_row"))
	}
	if(p_col > 0) {
		col_nms <- colnames(Xcol)
		if(is.null(col_nms)) col_nms <- paste0("Xcol", 1:p_col)
		beta_names <- c(beta_names, paste0(col_nms, "_col"))
	}
	
	#### MCMC initialization ####
	if(is.null(start_vals)) {
		if(family == "normal") {
			Z <- Y
		} else if(family == "tobit") {
			Z <- Y
			Z[Y < 0] <- 0
			if(sum(Y > 0, na.rm = TRUE) > 0) {
				Z[Y == 0] <- -mean(Y[Y > 0], na.rm = TRUE)
			}
		} else if(family == "binary") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			Z[Y == 1] <- abs(Z[Y == 1])
			Z[Y == 0] <- -abs(Z[Y == 0])
		} else if(family == "ordinal") {
			Z <- matrix(qnorm((Y + 0.5) / (max(Y, na.rm = TRUE) + 1)), nA, nB)
		} else if(family == "cbin") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			Z[Y == 1] <- abs(Z[Y == 1])
			Z[Y == 0] <- -abs(Z[Y == 0])
		} else if(family == "frn") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			for(i in 1:nA) {
				if(any(Y[i,] > 0, na.rm = TRUE)) {
					ranked <- which(Y[i,] > 0)
					Z[i, ranked] <- sort(abs(rnorm(length(ranked))), decreasing = TRUE)
				}
			}
		} else if(family == "rrl") {
			Z <- matrix(0, nA, nB)
			for(i in 1:nA) {
				row_vals <- Y[i,]
				if(any(!is.na(row_vals))) {
					Z[i,] <- qnorm((rank(row_vals, na.last = "keep") - 0.5) / sum(!is.na(row_vals)))
				}
			}
		} else if(family == "poisson") {
			Z <- log(Y + 1)
		} else {
			Z <- Y
		}

		beta <- rep(0, p)
		a <- rep(0, nA)
		b <- rep(0, nB)
		Sab <- Sab0
		s2 <- 1
		
		if(R_row > 0) {
			U <- matrix(rnorm(nA * R_row, 0, 0.1), nA, R_row)
		} else {
			U <- NULL
		}
		
		if(R_col > 0) {
			V <- matrix(rnorm(nB * R_col, 0, 0.1), nB, R_col)
		} else {
			V <- NULL
		}
		
		if(R_row > 0 && R_col > 0) {
			G <- diag(min(R_row, R_col), R_row, R_col)
		} else {
			G <- NULL
		}
	} else {
		Z <- start_vals$Z
		beta <- start_vals$beta
		a <- start_vals$a
		b <- start_vals$b
		Sab <- start_vals$Sab
		s2 <- start_vals$s2
		U <- start_vals$U
		V <- start_vals$V
		G <- start_vals$G
	}
	
	if(family == "binary" && is.null(prior$eta0)) {
		ydist <- table(Y)
		ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
		YB <- 1 * (Y != ymode)
		ybar <- mean(YB, na.rm=TRUE)
		mu_bin <- qnorm(ybar)
		E <- (YB - ybar) / dnorm(qnorm(ybar))
		a_tmp <- rowMeans(E, na.rm=TRUE)
		b_tmp <- colMeans(E, na.rm=TRUE)
		a_tmp[is.na(a_tmp)] <- 0
		b_tmp[is.na(b_tmp)] <- 0
		vscale <- (var(a_tmp) + var(b_tmp)) / 2
		PHAT <- pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
		vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
		eta0 <- round(4 * vdfmlt)
		if(is.null(prior$etaab)) etaab <- eta0
	}
	
	BETA <- matrix(NA_real_, n_eff, p)
	if(length(beta_names) == p) colnames(BETA) <- beta_names
	VC <- matrix(NA_real_, n_eff, 5)
	colnames(VC) <- c("va", "cab", "vb", "ve", "rho")
	sample_idx <- 0

	U_samples <- NULL
	V_samples <- NULL
	a_samples <- NULL
	b_samples <- NULL
	
	if(!is.null(posterior_opts)) {
		if(posterior_opts$save_UV && R_row > 0 && R_col > 0) {
			n_UV_samples <- ceiling(n_eff / posterior_opts$thin_UV)
			U_samples <- array(NA_real_, dim = c(nA, R_row, n_UV_samples))
			V_samples <- array(NA_real_, dim = c(nB, R_col, n_UV_samples))
		}
		if(posterior_opts$save_ab && (rvar || cvar)) {
			n_ab_samples <- ceiling(n_eff / posterior_opts$thin_ab)
			if(rvar) a_samples <- matrix(NA_real_, nA, n_ab_samples)
			if(cvar) b_samples <- matrix(NA_real_, nB, n_ab_samples)
		}
	}
	
	UV_sample_idx <- 0
	ab_sample_idx <- 0
	APS <- rep(0, nA)
	BPS <- rep(0, nB)
	YPS <- array(0, dim = c(nA, nB))
	
	if(gof) {
		n_base_stats <- 3
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
		GOF <- matrix(NA, n_eff, n_total_stats)
		
		base_names <- c("sd.rowmean", "sd.colmean", "four.cycles")
		all_names <- c(base_names, custom_gof_names)
		colnames(GOF) <- all_names
	} else {
		GOF <- NULL
		custom_gof <- NULL
	}
	
	if(p > 0) {
		X <- array(0, c(nA, nB, p))
		idx <- 1
		
		if(intercept) {
			X[,,idx] <- 1
			idx <- idx + 1
		}
		
		if(p_dyad > 0) {
			for(k in 1:p_dyad) {
				X[,,idx] <- Xdyad[,,k]
				idx <- idx + 1
			}
		}
		
		if(p_row > 0) {
			for(k in 1:p_row) {
				X[,,idx] <- Xrow[,k]
				idx <- idx + 1
			}
		}
		
		if(p_col > 0) {
			for(k in 1:p_col) {
				for(j in 1:nB) {
					X[,j,idx] <- Xcol[j,k]
				}
				idx <- idx + 1
			}
		}
	} else {
		X <- NULL
	}
	
	if(verbose) {
		cli::cli_h2("Running MCMC for bipartite network")
		cli::cli_text("Dimensions: {.val {nA}} x {.val {nB}}")
		cli::cli_text("R_row = {.val {R_row}}, R_col = {.val {R_col}}")
		cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")
		
		if(gof) {
			n_gof_stats <- 3
			if(!is.null(custom_gof)) {
				if(is.function(custom_gof)) {
					n_gof_stats <- n_gof_stats + length(custom_gof(Y))
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
	
	update_Z_fn <- switch(family,
		normal = function(Z, EZ, s2, Y) Y,
		tobit = function(Z, EZ, s2, Y) {
			Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
			Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
			Z_new[Y > 0] <- Y[Y > 0]
			Z_new
		},
		binary = function(Z, EZ, s2, Y) {
			Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
			Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
			Z_new[Y == 1] <- pmax(Z_new[Y == 1], 0)
			Z_new
		},
		poisson = function(Z, EZ, s2, Y) {
			lambda <- exp(pmin(EZ, 10))
			not_na <- !is.na(Y)
			n_update <- sum(not_na)
			if(n_update > 0) {
				Z_new <- Z
				Z_new[not_na] <- log(Y[not_na] + rgamma(n_update, 1, lambda[not_na]))
				Z_new
			} else {
				Z
			}
		},
		ordinal = function(Z, EZ, s2, Y) {
			EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
		},
		cbin = function(Z, EZ, s2, Y) {
			Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
			Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
			Z_new[Y == 1] <- pmax(Z_new[Y == 1], 0)
			Z_new
		},
		frn = function(Z, EZ, s2, Y) {
			EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
		},
		rrl = function(Z, EZ, s2, Y) {
			EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
		},
		stop(paste("Unknown family:", family))
	)
	
	needs_s2 <- family %in% c("normal", "tobit", "poisson")
	
	UV_cache <- if(!is.null(U) && !is.null(V) && !is.null(G)) {
		U %*% G %*% t(V)
	} else {
		matrix(0, nA, nB)
	}

	# compute EZ_cache from current parameter values
	EZ_cache <- matrix(0, nA, nB)
	if(p > 0) for(k in 1:p) EZ_cache <- EZ_cache + beta[k] * X[,,k]
	if(rvar) EZ_cache <- EZ_cache + a
	if(cvar) for(j in 1:nB) EZ_cache[,j] <- EZ_cache[,j] + b[j]
	EZ_cache <- EZ_cache + UV_cache
	
	if(p > 0) {
		XtX <- matrix(0, p, p)
		for(k1 in 1:p) {
			for(k2 in k1:p) {
				XtX[k1,k2] <- sum(X[,,k1] * X[,,k2], na.rm = TRUE)
				if(k1 != k2) XtX[k2,k1] <- XtX[k1,k2]
			}
		}
	}
	
	iter_save <- 0

	#### MCMC loop ####
	for(iter in 1:(nscan + burn)) {
		EZ_curr <- EZ_cache

		Z <- update_Z_fn(Z, EZ_curr, s2, Y)
		
		if(p > 0) {
			resid <- Z
			if(rvar) resid <- resid - a
			if(cvar) {
				for(j in 1:nB) {
					resid[,j] <- resid[,j] - b[j]
				}
			}
			if(!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}
			
			XtX <- matrix(0, p, p)
			Xty <- rep(0, p)
			
			for(k1 in 1:p) {
				for(k2 in 1:p) {
					XtX[k1,k2] <- sum(X[,,k1] * X[,,k2], na.rm = TRUE)
				}
				Xty[k1] <- sum(X[,,k1] * resid, na.rm = TRUE)
			}
			
			if(is.na(g)) g <- nA * nB
			
			V_beta <- solve(XtX / s2 + diag(p) / (g * s20))
			m_beta <- V_beta %*% (Xty / s2)
			
			beta <- m_beta + t(chol(V_beta)) %*% rnorm(p)
		}
		
		if(rvar) {
			resid <- Z
			if(p > 0) {
				for(k in 1:p) {
					resid <- resid - beta[k] * X[,,k]
				}
			}
			if(cvar) {
				for(j in 1:nB) {
					resid[,j] <- resid[,j] - b[j]
				}
			}
			if(!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}
			
			for(i in 1:nA) {
				prec_a <- nB / s2 + 1 / Sab[1,1]
				mean_a <- sum(resid[i,], na.rm = TRUE) / s2 / prec_a
				a[i] <- rnorm(1, mean_a, sqrt(1/prec_a))
			}
		}
		
		if(cvar) {
			resid <- Z
			if(p > 0) {
				for(k in 1:p) {
					resid <- resid - beta[k] * X[,,k]
				}
			}
			if(rvar) resid <- resid - a
			if(!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}
			
			for(j in 1:nB) {
				prec_b <- nA / s2 + 1 / Sab[2,2]
				mean_b <- sum(resid[,j], na.rm = TRUE) / s2 / prec_b
				b[j] <- rnorm(1, mean_b, sqrt(1/prec_b))
			}
		}
		
		if(rvar || cvar) {
			# sample va and vb independently for bipartite (no cross-correlation)
			if(rvar) {
				Sab[1,1] <- 1/rgamma(1, (etaab + nA)/2, (etaab*Sab0[1,1] + sum(a^2))/2)
			}
			if(cvar) {
				Sab[2,2] <- 1/rgamma(1, (etaab + nB)/2, (etaab*Sab0[2,2] + sum(b^2))/2)
			}
			Sab[1,2] <- Sab[2,1] <- 0
		}
		
		if(R_row > 0 && R_col > 0) {
			resid <- Z
			if(p > 0) {
				for(k in 1:p) {
					resid <- resid - beta[k] * X[,,k]
				}
			}
			if(rvar) resid <- resid - a
			if(cvar) {
				for(j in 1:nB) {
					resid[,j] <- resid[,j] - b[j]
				}
			}
			
			for(i in 1:nA) {
				U_prop <- U
				U_prop[i,] <- U[i,] + rnorm(R_row, 0, 0.1)
				
				ll_curr <- -sum((resid[i,] - U[i,] %*% G %*% t(V))^2, na.rm = TRUE) / (2*s2)
				ll_prop <- -sum((resid[i,] - U_prop[i,] %*% G %*% t(V))^2, na.rm = TRUE) / (2*s2)
				
				if(log(runif(1)) < ll_prop - ll_curr) {
					U <- U_prop
				}
			}
			
			for(j in 1:nB) {
				V_prop <- V
				V_prop[j,] <- V[j,] + rnorm(R_col, 0, 0.1)
				
				ll_curr <- -sum((resid[,j] - U %*% G %*% V[j,])^2, na.rm = TRUE) / (2*s2)
				ll_prop <- -sum((resid[,j] - U %*% G %*% V_prop[j,])^2, na.rm = TRUE) / (2*s2)
				
				if(log(runif(1)) < ll_prop - ll_curr) {
					V <- V_prop
				}
			}
			
			UV_cache <- U %*% G %*% t(V)
		}
		
		if(family == "normal") {
			resid <- Z
			if(p > 0) {
				for(k in 1:p) {
					resid <- resid - beta[k] * X[,,k]
				}
			}
			if(rvar) resid <- resid - a
			if(cvar) {
				for(j in 1:nB) {
					resid[,j] <- resid[,j] - b[j]
				}
			}
			if(!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}
			
			shape <- (nA * nB) / 2
			rate <- sum(resid^2, na.rm = TRUE) / 2
			s2 <- 1 / rgamma(1, shape, rate)
		}
		
		# recompute EZ_cache from current parameter values
		EZ_cache <- matrix(0, nA, nB)
		if(p > 0) for(k in 1:p) EZ_cache <- EZ_cache + beta[k] * X[,,k]
		if(rvar) EZ_cache <- EZ_cache + a
		if(cvar) for(j in 1:nB) EZ_cache[,j] <- EZ_cache[,j] + b[j]
		EZ_cache <- EZ_cache + UV_cache
		
		if(iter > burn && (iter - burn) %% odens == 0) {
			iter_save <- iter_save + 1
			
			if(iter_save <= n_eff) {
				BETA[iter_save,] <- beta
				VC[iter_save,] <- c(Sab[1,1], Sab[1,2], Sab[2,2], s2, 0)  # no rho for bipartite
				
				if(rvar) APS <- APS + a
				if(cvar) BPS <- BPS + b
				if(!is.null(U) && !is.null(V) && !is.null(G)) {
					# uvps not accumulated for bipartite
				}
				
				# save posterior samples if requested
				if(!is.null(posterior_opts)) {
					# u/v samples (thinned)
					if(posterior_opts$save_UV && !is.null(U) && !is.null(V)) {
						if(iter_save %% posterior_opts$thin_UV == 0) {
							UV_sample_idx <- UV_sample_idx + 1
							if(UV_sample_idx <= dim(U_samples)[3]) {
								U_samples[,,UV_sample_idx] <- U
								V_samples[,,UV_sample_idx] <- V
							}
						}
					}
					
					# a/b samples (thinned)
					if(posterior_opts$save_ab) {
						if(iter_save %% posterior_opts$thin_ab == 0) {
							ab_sample_idx <- ab_sample_idx + 1
							if(rvar && !is.null(a_samples) && ab_sample_idx <= ncol(a_samples)) {
								a_samples[,ab_sample_idx] <- a
							}
							if(cvar && !is.null(b_samples) && ab_sample_idx <= ncol(b_samples)) {
								b_samples[,ab_sample_idx] <- b
							}
						}
					}
				}
			}
			
			EZ_curr <- matrix(0, nA, nB)
			if(p > 0) {
				for(k in 1:p) {
					EZ_curr <- EZ_curr + beta[k] * X[,,k]
				}
			}
			if(rvar) EZ_curr <- EZ_curr + a
			if(cvar) {
				for(j in 1:nB) {
					EZ_curr[,j] <- EZ_curr[,j] + b[j]
				}
			}
			if(!is.null(U) && !is.null(V) && !is.null(G)) {
				EZ_curr <- EZ_curr + U %*% G %*% t(V)
			}
			
			if(family == "normal") {
				Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
			} else if(family == "tobit") {
				Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
				Y_sim[Y_sim < 0] <- 0
			} else if(family == "binary") {
				Y_sim <- 1 * (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB) > 0)
			} else if(family == "ordinal") {
				Y_sim <- round(EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB))
				Y_sim[Y_sim < 0] <- 0
				if(any(!is.na(Y))) {
					max_cat <- max(Y, na.rm = TRUE)
					Y_sim[Y_sim > max_cat] <- max_cat
				}
			} else if(family == "cbin") {
				Y_sim <- 1 * (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB) > 0)
				if(!is.null(odmax)) {
					for(i in 1:nA) {
						if(sum(Y_sim[i,]) > odmax[i]) {
							# keep only top odmax[i] values
							thresh <- sort(EZ_curr[i,], decreasing = TRUE)[odmax[i]]
							Y_sim[i, EZ_curr[i,] < thresh] <- 0
						}
					}
				}
			} else if(family == "frn") {
				Y_sim <- matrix(0, nA, nB)
				if(!is.null(odmax)) {
					for(i in 1:nA) {
						ranks <- rank(-EZ_curr[i,], ties.method = "random")
						Y_sim[i, ranks <= odmax[i]] <- ranks[ranks <= odmax[i]]
					}
				}
			} else if(family == "rrl") {
				Y_sim <- matrix(0, nA, nB)
				for(i in 1:nA) {
					Y_sim[i,] <- rank(EZ_curr[i,] + rnorm(nB, 0, sqrt(s2)))
				}
			} else if(family == "poisson") {
				Y_sim <- matrix(rpois(nA * nB, exp(EZ_curr)), nA, nB)
			} else {
				Y_sim <- EZ_curr
			}
			
			if(iter_save <= n_eff) {
				YPS <- YPS + Y_sim
			}
			
			if(gof && iter_save <= n_eff) {
				gof_vals <- gof_stats_bipartite(Y_sim)
				GOF[iter_save, 1:3] <- gof_vals

				if(!is.null(custom_gof)) {
					if(is.function(custom_gof)) {
						tryCatch({
							custom_vals <- custom_gof(Y_sim)
							GOF[iter_save, (4):(3 + length(custom_vals))] <- custom_vals
						}, error = function(e) {})
					} else if(is.list(custom_gof)) {
						for(i in seq_along(custom_gof)) {
							if(is.function(custom_gof[[i]])) {
								tryCatch({
									GOF[iter_save, 3 + i] <- custom_gof[[i]](Y)
								}, error = function(e) {})
							}
						}
					}
				}
			}
		}
		
		if(verbose) {
			if(iter <= burn) {
				cli::cli_progress_update(id = NULL)
				if(iter == burn) {
					cli::cli_progress_done()
					cli::cli_progress_bar("Sampling", total = nscan, clear = TRUE)
				}
			} else if(iter > burn) {
				cli::cli_progress_update(id = NULL)
			}
		}
	}
	
	#### output assembly ####
	if(verbose) {
		cli::cli_progress_done()

		end_time <- Sys.time()
		elapsed <- difftime(end_time, start_time, units = "secs")
		cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
	}
	
	APM <- APS / n_eff
	names(APM) <- rownames(Y)
	BPM <- BPS / n_eff
	names(BPM) <- colnames(Y)
	YPM <- YPS / n_eff
	rownames(YPM) <- rownames(Y)
	colnames(YPM) <- colnames(Y)
	
	fit <- list(
		BETA = BETA,
		VC = VC,
		APM = APM,
		BPM = BPM,
		U = U,
		V = V,
		G = G,
		YPM = YPM,
		GOF = GOF,
		X = X,
		Y = Y,
		family = family,
		R_row = R_row,
		R_col = R_col,
		rvar = rvar,
		cvar = cvar,
		symmetric = FALSE,
		mode = "bipartite",
		nA = nA,
		nB = nB
	)
	
	if(!is.null(posterior_opts)) {
		if(!is.null(U_samples)) {
			actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
			if(actual_UV_samples > 0) {
				fit$U_samples <- U_samples[,,1:actual_UV_samples, drop=FALSE]
				fit$V_samples <- V_samples[,,1:actual_UV_samples, drop=FALSE]
			}
		}
		if(!is.null(a_samples)) {
			actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
			if(actual_ab_samples > 0 && rvar) {
				fit$a_samples <- a_samples[,1:actual_ab_samples, drop=FALSE]
			}
		}
		if(!is.null(b_samples)) {
			actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
			if(actual_ab_samples > 0 && cvar) {
				fit$b_samples <- b_samples[,1:actual_ab_samples, drop=FALSE]
			}
		}
	}
	
	class(fit) <- "ame"
	
	if(use_sparse_matrices) {
		fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
	}
	
	return(fit)
}