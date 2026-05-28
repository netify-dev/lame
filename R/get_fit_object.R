####
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
#' @param actorByYr List of actors by time point. In bipartite mode this is
#'   the per-year list of row actors.
#' @param colActorByYr Bipartite only. List of column actors by time point;
#'   defaults to \code{NULL} (unipartite).
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
#' @param Y_obs original observed network (stored for residuals computation)
#' @param G bipartite interaction matrix mapping row to column latent spaces
#' @param dynamic_beta logical or scalar; whether the BETA storage is 3-D
#'   (dynamic_beta path). Default \code{FALSE}.
#' @param beta_dynamic_mask logical vector marking which coefficients are dynamic.
#' @param beta_dynamic_groups character vector of per-coefficient block labels
#'   ("intercept", "dyad", "row", "col"); \code{""} for static coefficients.
#' @param rho_beta named numeric vector of per-block AR(1) rho values (one
#'   per dynamic block).
#' @param sigma_beta named numeric vector of per-block AR(1) innovation
#'   standard deviations.
#' @param RHO_BETA matrix of per-iteration rho_beta draws (rows = MCMC draw,
#'   cols = dynamic block).
#' @param SIGMA_BETA matrix of per-iteration sigma_beta draws.
#' @return Fitted AME object
#' @author Shahryar Minhas
#' @export get_fit_object
get_fit_object <- function(
	APS, BPS, UVPS, YPS,
	BETA, VC, GOF,
	Xlist, actorByYr, colActorByYr = NULL, start_vals,
	symmetric, tryErrorChecks,
	model.name=NULL,
	U=NULL, V=NULL, dynamic_uv=FALSE, dynamic_ab=FALSE,
	bip=FALSE,
	rho_ab=NULL, rho_uv=NULL,
	family=NULL, odmax=NULL, nA=NULL, nB=NULL, n_time=NULL,
	Y_obs=NULL, G=NULL,
	dynamic_beta=FALSE,
	beta_dynamic_mask=NULL, beta_dynamic_groups=NULL,
	rho_beta=NULL, sigma_beta=NULL,
	RHO_BETA=NULL, SIGMA_BETA=NULL
){
	# detect 3-D BETA (= dynamic_beta path). when 3-D, the per-period
	# posterior-mean beta is apply(BETA, c(2,3), mean) (p x T); the time-mean
	# beta is apply(BETA, 2, mean) (p,). these summaries feed the existing
	# downstream EZ/UVPM/YPM code.
	.beta_is_dynamic <- !is.null(BETA) && length(dim(BETA)) == 3L
	.beta_post_mean_per_t <- function(B) {
		if (length(dim(B)) == 3L) apply(B, c(2, 3), mean)
		else colMeans(B)
	}
	.beta_post_mean_overall <- function(B) {
		if (length(dim(B)) == 3L) apply(B, 2, mean)
		else colMeans(B)
	}

	####
	if(bip) {
		if(is.matrix(APS)) {
			row_actors <- rownames(APS)
			col_actors <- rownames(BPS)
		} else {
			row_actors <- names(APS)
			col_actors <- names(BPS)
		}
		# dynamic-ab bipartite stores APS/BPS as unnamed matrices; recover
		# the real row and column actor names from actorByYr/colActorByYr
		# (the MCMC's sorted actor order) so APM/BPM/U/V/YPM are not
		# labelled "Row1/Col1..." with the values silently in sorted order
		n_row_a <- if (is.matrix(APS)) nrow(APS) else length(APS)
		n_col_a <- if (is.matrix(BPS)) nrow(BPS) else length(BPS)
		if (is.null(row_actors)) {
			real_row <- sort(unique(unlist(actorByYr)))
			row_actors <- if (length(real_row) == n_row_a) real_row
			              else paste0("Row", seq_len(n_row_a))
		}
		if (is.null(col_actors)) {
			real_col <- sort(unique(unlist(colActorByYr)))
			col_actors <- if (length(real_col) == n_col_a) real_col
			              else paste0("Col", seq_len(n_col_a))
		}
		actors <- row_actors
	} else {
		if(is.matrix(APS)) {
			actors <- rownames(APS)
		} else {
			actors <- names(APS)
		}
		if(is.null(actors)) {
			n_a <- if(is.matrix(APS)) nrow(APS) else length(APS)
			# a dynamic-ab fit stores APS as an unnamed actor x time matrix;
			# recover the real actor names (the MCMC's sorted actor order)
			# from actorByYr so APM/EZ/YPM are not labelled "Actor1..."
			real_actors <- sort(unique(unlist(actorByYr)))
			actors <- if(length(real_actors) == n_a) {
				real_actors
			} else {
				paste0("Actor", 1:n_a)
			}
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

	####
	if(dynamic_ab && is.matrix(APS)) {
		APM <- APS
		BPM <- BPS
		if(bip) {
			rownames(APM) <- row_actors
			rownames(BPM) <- col_actors
			APM_avg <- rowMeans(APM); names(APM_avg) <- row_actors
			BPM_avg <- rowMeans(BPM); names(BPM_avg) <- col_actors
		} else {
			rownames(APM) <- rownames(BPM) <- actors
			APM_avg <- rowMeans(APM); names(APM_avg) <- actors
			BPM_avg <- rowMeans(BPM); names(BPM_avg) <- actors
		}
	} else {
		APM<-APS/nrow(VC)
		BPM<-BPS/nrow(VC)
		if(bip) {
			names(APM) <- row_actors
			names(BPM) <- col_actors
		} else {
			names(APM) <- names(BPM) <- actors
		}
		APM_avg <- APM
		BPM_avg <- BPM
	}
	YPM<-YPS/nrow(VC)
	if(!bip) {
		dimnames(YPM) <- list(actors, actors, pdLabs)
	} else {
		dimnames(YPM) <- list(row_actors, col_actors, pdLabs)
	}

	####
	if(!bip) {
		# per-period vs scalar beta-mean depending on dynamic_beta state
		beta_per_t <- .beta_post_mean_per_t(BETA)  # p x T (matrix) or p-vec
		beta_overall <- .beta_post_mean_overall(BETA)  # p-vec
		.beta_at_t <- function(t) {
			if (.beta_is_dynamic) beta_per_t[, t] else beta_overall
		}
		if(dynamic_uv && length(dim(UVPS)) == 3) {
			UVPM <- UVPS
			# build EZ per time slice so the time-varying U_t V_t' signal
			# flows into EZ; predict(type = "link") then returns one matrix
			# per period under dynamic_uv.
			N_t <- length(Xlist)
			n_actors <- dim(UVPS)[1]
			I_mat <- diag(n_actors)
			EZ <- array(NA_real_, dim = c(n_actors, n_actors, N_t))
			for (t in seq_len(N_t)) {
				a_t <- if (dynamic_ab && is.matrix(APM)) APM[, t] else APM_avg
				b_t <- if (dynamic_ab && is.matrix(BPM)) BPM[, t] else BPM_avg
				EZ_t <- get_EZ_cpp(list(Xlist[[t]]), .beta_at_t(t),
				                   outer(a_t, b_t, "+"),
				                   UVPS[, , t], I_mat)
				EZ[, , t] <- EZ_t[, , 1]
			}
		} else {
			UVPM <- UVPS
			if (length(dim(UVPM)) == 2) {
				# dynamic_ab with static UV: per-slice loop so each period
				# uses its own a_t / b_t.
				if (dynamic_ab && is.matrix(APM)) {
					N_t <- length(Xlist)
					n_actors <- nrow(UVPM)
					I_mat <- diag(n_actors)
					EZ <- array(NA_real_, dim = c(n_actors, n_actors, N_t))
					for (t in seq_len(N_t)) {
						EZ_t <- get_EZ_cpp(list(Xlist[[t]]), .beta_at_t(t),
						                   outer(APM[, t], BPM[, t], "+"),
						                   UVPM, I_mat)
						EZ[, , t] <- EZ_t[, , 1]
					}
				} else if (.beta_is_dynamic) {
					# dynamic_beta with static a/b: loop over t for per-period beta
					N_t <- length(Xlist)
					n_actors <- nrow(UVPM)
					I_mat <- diag(n_actors)
					EZ <- array(NA_real_, dim = c(n_actors, n_actors, N_t))
					for (t in seq_len(N_t)) {
						EZ_t <- get_EZ_cpp(list(Xlist[[t]]), .beta_at_t(t),
						                   outer(APM_avg, BPM_avg, "+"),
						                   UVPM, I_mat)
						EZ[, , t] <- EZ_t[, , 1]
					}
				} else {
					EZ <- get_EZ_cpp(Xlist,
					                 beta_overall,
					                 outer(APM_avg, BPM_avg, "+"),
					                 UVPM, diag(nrow(UVPM)))
				}
			}
		}

		dimnames(EZ) <- list(actors, actors, pdLabs)
	} else {
		UVPM <- UVPS
		if(!is.null(Xlist) && !is.null(G) && !is.null(nA) && !is.null(nB)) {
			N_t <- if(!is.null(n_time)) n_time else length(Xlist)
			# bipartite + dynamic_beta: use per-period beta in base_cube
			beta_per_t_bip <- .beta_post_mean_per_t(BETA)
			beta_overall_bip <- .beta_post_mean_overall(BETA)

			base_cube <- array(0, dim = c(nA, nB, N_t))
			for(t in 1:N_t) {
				beta_t_here <- if (.beta_is_dynamic) beta_per_t_bip[, t] else beta_overall_bip
				if(length(beta_t_here) > 0 && !is.null(Xlist[[t]])) {
					p_x <- min(length(beta_t_here), dim(Xlist[[t]])[3])
					for(k in 1:p_x) {
						base_cube[,,t] <- base_cube[,,t] + beta_t_here[k] * Xlist[[t]][,,k]
					}
				}
			}

			a_mat_ez <- matrix(APM_avg, nA, N_t)
			b_mat_ez <- matrix(BPM_avg, nB, N_t)

			RA <- if(!is.null(U) && length(dim(U)) == 3) dim(U)[2] else if(!is.null(U)) ncol(U) else 1
			RB <- if(!is.null(V) && length(dim(V)) == 3) dim(V)[2] else if(!is.null(V)) ncol(V) else 1
			if(!is.null(U) && length(dim(U)) == 2) {
				U_cube <- array(0, dim = c(nA, RA, N_t))
				for(t in 1:N_t) U_cube[,,t] <- U
			} else if(!is.null(U) && length(dim(U)) == 3) {
				U_cube <- U
			} else {
				U_cube <- array(0, dim = c(nA, 1, N_t))
			}
			if(!is.null(V) && length(dim(V)) == 2) {
				V_cube <- array(0, dim = c(nB, RB, N_t))
				for(t in 1:N_t) V_cube[,,t] <- V
			} else if(!is.null(V) && length(dim(V)) == 3) {
				V_cube <- V
			} else {
				V_cube <- array(0, dim = c(nB, 1, N_t))
			}

			EZ <- tryCatch(
				get_EZ_bip_cpp(base_cube, a_mat_ez, b_mat_ez, U_cube, V_cube, G),
				error = function(e) NULL
			)
		} else {
			EZ <- NULL
		}
	}

	####
	if(!bip) {
		if(dynamic_uv && length(dim(UVPM)) == 3) {
			dimnames(UVPM) <- list(actors, actors, pdLabs)
		} else if(length(dim(UVPM)) == 2) {
			rownames(UVPM)<-colnames(UVPM)<-actors
		}
	}
	# BETA is 2-D (matrix) or 3-D (array) -- the 3-D path is the dynamic_beta
	# storage and its first dim names should also be NULL'd if present
	if (length(dim(BETA)) == 2L) {
		rownames(BETA) <- NULL
	} else if (length(dim(BETA)) == 3L) {
		dn <- dimnames(BETA)
		if (!is.null(dn)) dimnames(BETA) <- list(NULL, dn[[2]], dn[[3]])
	}

	####
	if(!symmetric){
		if(dynamic_uv && !is.null(U) && !is.null(V)) {
			if(length(dim(U)) == 3) {
				if(bip) {
					dimnames(U) <- list(row_actors, NULL, pdLabs)
					dimnames(V) <- list(col_actors, NULL, pdLabs)
				} else {
					dimnames(U) <- list(actors, NULL, pdLabs)
					dimnames(V) <- list(actors, NULL, pdLabs)
				}
			} else {
				if(bip) {
					rownames(U) <- row_actors
					rownames(V) <- col_actors
				} else {
					rownames(U)<-rownames(V)<-actors
				}
			}
		} else {
			if(!bip) {
				UDV<-svd(UVPM)
				U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
				V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
				rownames(U)<-rownames(V)<-actors
			} else if(!is.null(U) && !is.null(V)) {
				if(length(dim(U)) == 3) {
					dimnames(U) <- list(row_actors, NULL, pdLabs)
					dimnames(V) <- list(col_actors, NULL, pdLabs)
				} else {
					rownames(U) <- row_actors
					rownames(V) <- col_actors
				}
			}
		}
	}
	
	if(symmetric){
		if(dynamic_uv && !is.null(U)) {
			if(length(dim(U)) == 3) {
				dimnames(U) <- list(actors, NULL, pdLabs)
				ULUPM <- UVPM
				L <- NULL
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

	####
	# transform EZ to count scale for poisson
	if(!is.null(family) && family == "poisson" && !is.null(EZ)) {
		if(is.array(EZ)) {
			EZ <- exp(EZ)
			EZ[EZ > 1e6] <- 1e6
		}
	}

	####
	# APM/BPM/U/V are indexed in the sorted actor order the MCMC runs in;
	# EZ/YPM must use the SAME per-slice order or a user combining APM with
	# YPM positionally gets silently misaligned actors. Sort each slice's
	# actor list so every per-actor object in the fit shares one ordering.
	if(!bip) {
		EZ <- array_to_list(EZ, lapply(actorByYr, sort), pdLabs)
		YPM <- array_to_list(YPM, lapply(actorByYr, sort), pdLabs)
	} else {
		if(!is.null(EZ) && is.array(EZ) && length(dim(EZ)) == 3) {
			N_t <- dim(EZ)[3]
			EZ_list <- vector("list", N_t)
			YPM_list <- vector("list", N_t)
			for(t in 1:N_t) {
				EZ_list[[t]] <- EZ[,,t]
				YPM_list[[t]] <- YPM[,,t]
			}
			names(EZ_list) <- names(YPM_list) <- pdLabs
			EZ <- EZ_list
			YPM <- YPM_list
		}
	}

	####
	if(is.array(GOF) && length(dim(GOF)) == 3) {
		n_gof <- dim(GOF)[1]
		gof_names <- dimnames(GOF)[[1]]
		if(is.null(gof_names)) {
			gof_names <- if(n_gof == 5) {
				c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
			} else {
				c("sd.rowmean", "sd.colmean", "four.cycles")
			}
		}
		GOF_list <- list()
		for(i in 1:n_gof) {
			GOF_list[[gof_names[i]]] <- GOF[i,,]
		}
		GOF <- GOF_list
	}

	####
	if(symmetric){
		fit <- list(
			BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
			YPM=YPM,GOF=GOF, start_vals=start_vals, tryErrorChecks=tryErrorChecks,
			model.name=model.name, family=family, symmetric=symmetric, odmax=odmax,
			mode=if(bip) "bipartite" else "unipartite")
		if(dynamic_ab && is.matrix(APM)) {
			fit$a_dynamic <- APM
			fit$APM <- APM_avg
		}
		if(!is.null(rho_ab)) fit$rho_ab <- rho_ab
		if(!is.null(rho_uv)) fit$rho_uv <- rho_uv
	}
	if(!symmetric){
		fit <- list(
			BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
			YPM=YPM,GOF=GOF, start_vals=start_vals, tryErrorChecks=tryErrorChecks,
			model.name=model.name, family=family, symmetric=symmetric, odmax=odmax,
			mode=if(bip) "bipartite" else "unipartite")
		if(dynamic_ab && is.matrix(APM)) {
			fit$a_dynamic <- APM
			fit$b_dynamic <- BPM
			fit$APM <- APM_avg
			fit$BPM <- BPM_avg
		}
		if(!is.null(rho_ab)) fit$rho_ab <- rho_ab
		if(!is.null(rho_uv)) fit$rho_uv <- rho_uv
	}
	if(!is.null(Y_obs)) fit$Y <- Y_obs
	if(bip) {
		fit$nA <- nA
		fit$nB <- nB
		if(!is.null(G)) fit$G <- G
	}
	if(!is.null(n_time)) fit$n_time <- n_time
	fit$dynamic_uv <- dynamic_uv
	fit$dynamic_ab <- dynamic_ab
	fit$dynamic_beta <- isTRUE(dynamic_beta) || .beta_is_dynamic
	# expose dynamic-beta metadata so downstream S3 methods can dispatch
	if (.beta_is_dynamic) {
		fit$beta_dynamic_mask   <- beta_dynamic_mask
		fit$beta_dynamic_groups <- beta_dynamic_groups
		# per-period posterior-mean beta convenience
		fit$beta_path <- apply(BETA, c(2, 3), mean)
		if (!is.null(RHO_BETA))   fit$RHO_BETA   <- RHO_BETA
		if (!is.null(SIGMA_BETA)) fit$SIGMA_BETA <- SIGMA_BETA
		if (!is.null(rho_beta))   fit$rho_beta   <- rho_beta
		if (!is.null(sigma_beta)) fit$sigma_beta <- sigma_beta
	}
	if(!is.null(Xlist)) fit$Xlist <- Xlist
	class(fit)<-"ame"
	return(fit)
}
####