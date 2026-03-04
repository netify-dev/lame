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
#' @param Y_obs original observed network (stored for residuals computation)
#' @param G bipartite interaction matrix mapping row to column latent spaces
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
	family=NULL, odmax=NULL, nA=NULL, nB=NULL, n_time=NULL,
	Y_obs=NULL, G=NULL
){

	####
	if(bip) {
		if(is.matrix(APS)) {
			row_actors <- rownames(APS)
			col_actors <- rownames(BPS)
			if(is.null(row_actors)) row_actors <- paste0("Row", 1:nrow(APS))
			if(is.null(col_actors)) col_actors <- paste0("Col", 1:nrow(BPS))
		} else {
			row_actors <- names(APS)
			col_actors <- names(BPS)
			if(is.null(row_actors)) row_actors <- paste0("Row", 1:length(APS))
			if(is.null(col_actors)) col_actors <- paste0("Col", 1:length(BPS))
		}
		actors <- row_actors
	} else {
		if(is.matrix(APS)) {
			actors <- rownames(APS)
		} else {
			actors <- names(APS)
		}
		if(is.null(actors)) {
			if(is.matrix(APS)) {
				actors <- paste0("Actor", 1:nrow(APS))
			} else {
				actors <- paste0("Actor", 1:length(APS))
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
		if(dynamic_uv && length(dim(UVPS)) == 3) {
			UVPM <- UVPS
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
		
		dimnames(EZ) <- list(actors, actors, pdLabs)
	} else {
		UVPM <- UVPS
		if(!is.null(Xlist) && !is.null(G) && !is.null(nA) && !is.null(nB)) {
			N_t <- if(!is.null(n_time)) n_time else length(Xlist)
			beta_mean <- colMeans(BETA)

			base_cube <- array(0, dim = c(nA, nB, N_t))
			for(t in 1:N_t) {
				if(length(beta_mean) > 0 && !is.null(Xlist[[t]])) {
					p_x <- min(length(beta_mean), dim(Xlist[[t]])[3])
					for(k in 1:p_x) {
						base_cube[,,t] <- base_cube[,,t] + beta_mean[k] * Xlist[[t]][,,k]
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
	rownames(BETA)<-NULL

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
	if(!bip) {
		EZ <- array_to_list(EZ, actorByYr, pdLabs)
		YPM <- array_to_list(YPM, actorByYr, pdLabs)
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
	if(!is.null(Xlist)) fit$Xlist <- Xlist
	class(fit)<-"ame"
	return(fit)
}
####