#' Bipartite network helper functions
#'
#' @name bipartite_helpers
#' @keywords internal

#' Initialize bipartite starting values
#' @return List of starting values for bipartite MCMC
#' @keywords internal
init_bipartite_startvals <- function(
	Y, family, nA, nB, RA, RB, T,
	Xlist = NULL, odmax = NULL
	){

	####
	# initialize Z
	Z <- array(dim = c(nA, nB, T))

	for(t in 1:T) {
		Yt <- Y[,,t]

		if(family == "normal") {
			Z[,,t] <- Yt
		} else if(family == "binary") {
			Z[,,t] <- matrix(qnorm((Yt + 0.5)/2), nA, nB)

			if(sum(Yt == 0, na.rm = TRUE) > 0 && sum(Yt == 1, na.rm = TRUE) > 0) {
				z0_max <- max(Z[,,t][Yt == 0], na.rm = TRUE)
				z1_min <- min(Z[,,t][Yt == 1], na.rm = TRUE)
				z_split <- (z0_max + z1_min) / 2
				Z[,,t] <- Z[,,t] - z_split
			}
		} else if(family == "poisson") {
			Z[,,t] <- log(Yt + 1)
		} else if(family %in% c("ordinal", "tobit")) {
			Z[,,t] <- Yt
		}
	}
	####

	####
	# initialize parameters
	a <- rep(0, nA)
	b <- rep(0, nB)

	if(RA > 0) {
		U <- array(rnorm(nA * RA * T, 0, 0.1), c(nA, RA, T))
	} else {
		U <- array(0, c(nA, 1, T))
	}

	if(RB > 0) {
		V <- array(rnorm(nB * RB * T, 0, 0.1), c(nB, RB, T))
	} else {
		V <- array(0, c(nB, 1, T))
	}

	if(RA > 0 && RB > 0) {
		min_rank <- min(RA, RB)
		G <- matrix(0, RA, RB)
		diag(G)[1:min_rank] <- 1
	} else {
		G <- matrix(0, max(1,RA), max(1,RB))
	}
	####

	####
	# beta and variance parameters
	if(!is.null(Xlist)) {
		p <- dim(Xlist[[1]])[3]
		beta <- rep(0, p)
	} else {
		beta <- numeric(0)
	}

	sigma2_a <- 1
	sigma2_b <- 1
	s2 <- 1
	rho <- NULL
	####

	list(
		Z = Z,
		a = a,
		b = b,
		U = U,
		V = V,
		G = G,
		beta = beta,
		sigma2_a = sigma2_a,
		sigma2_b = sigma2_b,
		s2 = s2,
		rho = rho
	)
}

#' Sample additive effects for bipartite networks
#' @return List with components \code{a} (row effects) and \code{b} (column effects)
#' @keywords internal
sample_ab_bipartite <- function(
	Z, EZ_without_ab, sigma2_a, sigma2_b, s2
	){
	dims <- dim(Z)
	nA <- dims[1]
	nB <- dims[2]
	T <- dims[3]

	R <- Z - EZ_without_ab

	####
	# sample row effects (a)
	a <- matrix(0, nA, T)
	prec_a <- T / s2 + 1 / sigma2_a

	for(i in 1:nA) {
		for(t in 1:T) {
			mean_a <- sum(R[i,,t], na.rm = TRUE) / (nB * s2)
			mean_a <- mean_a / prec_a
			a[i,t] <- rnorm(1, mean_a, sqrt(1/prec_a))
		}
	}

	# update residuals
	for(t in 1:T) {
		R[,,t] <- R[,,t] - a[,t]
	}
	####

	####
	# sample column effects (b)
	b <- matrix(0, nB, T)
	prec_b <- T / s2 + 1 / sigma2_b

	for(j in 1:nB) {
		for(t in 1:T) {
			mean_b <- sum(R[,j,t], na.rm = TRUE) / (nA * s2)
			mean_b <- mean_b / prec_b
			b[j,t] <- rnorm(1, mean_b, sqrt(1/prec_b))
		}
	}
	####

	list(a = a, b = b)
}

#' Update variance parameters for bipartite
#' @return List with components \code{sigma2_a} and \code{sigma2_b}
#' @keywords internal
update_variances_bipartite <- function(
	a, b, eta0_a = 2, eta0_b = 2,
	Sab0_aa = 1, Sab0_bb = 1
	){
	nA <- nrow(a)
	nB <- nrow(b)
	T <- ncol(a)

	# sample row variance
	shape_a <- (eta0_a + nA * T) / 2
	rate_a <- (eta0_a * Sab0_aa + sum(a^2)) / 2
	sigma2_a <- 1 / rgamma(1, shape_a, rate_a)

	# sample column variance
	shape_b <- (eta0_b + nB * T) / 2
	rate_b <- (eta0_b * Sab0_bb + sum(b^2)) / 2
	sigma2_b <- 1 / rgamma(1, shape_b, rate_b)

	list(sigma2_a = sigma2_a, sigma2_b = sigma2_b)
}

#' Compute GOF statistics for bipartite networks
#' @return List of observed and simulated GOF statistics
#' @keywords internal
compute_gof_bipartite <- function(Y_obs, Y_sim, family){
	T <- dim(Y_obs)[3]

	gof_stats <- list()

	####
	# density
	gof_stats$density_obs <- apply(Y_obs, 3, function(y) mean(y > 0, na.rm = TRUE))
	gof_stats$density_sim <- apply(Y_sim, 3, function(y) mean(y > 0, na.rm = TRUE))
	####

	####
	# degrees
	deg_bip <- compute_degrees_bip_cpp(Y_obs)
	gof_stats$row_degree_mean_obs <- colMeans(deg_bip$row_degrees)
	gof_stats$col_degree_mean_obs <- colMeans(deg_bip$col_degrees)

	deg_bip_sim <- compute_degrees_bip_cpp(Y_sim)
	gof_stats$row_degree_mean_sim <- colMeans(deg_bip_sim$row_degrees)
	gof_stats$col_degree_mean_sim <- colMeans(deg_bip_sim$col_degrees)
	####

	####
	# four-cycles and clustering
	gof_stats$four_cycles_obs <- as.numeric(count_four_cycles_bip_cpp(Y_obs))
	gof_stats$four_cycles_sim <- as.numeric(count_four_cycles_bip_cpp(Y_sim))

	gof_stats$clustering_obs <- as.numeric(clustering_coef_bip_cpp(Y_obs))
	gof_stats$clustering_sim <- as.numeric(clustering_coef_bip_cpp(Y_sim))
	####

	gof_stats
}

#' Convert bipartite list data to array format
#' @return List with components \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol} in array format
#' @keywords internal
list_to_array_bipartite <- function(
	rowActorSet, colActorSet, Y_list,
	Xdyad = NULL, Xrow = NULL, Xcol = NULL
	){
	nA <- length(rowActorSet)
	nB <- length(colActorSet)
	T <- length(Y_list)

	####
	# fill Y array
	Y_array <- array(NA, c(nA, nB, T),
		dimnames = list(rowActorSet, colActorSet, names(Y_list)))

	for(t in 1:T) {
		Yt <- Y_list[[t]]
		row_idx <- match(rownames(Yt), rowActorSet)
		col_idx <- match(colnames(Yt), colActorSet)

		for(i in seq_along(row_idx)) {
			for(j in seq_along(col_idx)) {
				if(!is.na(row_idx[i]) && !is.na(col_idx[j])) {
					Y_array[row_idx[i], col_idx[j], t] <- Yt[i, j]
				}
			}
		}
	}
	####

	####
	# coerce and reorder Xdyad to match sorted actor sets
	if(!is.null(Xdyad)) {
		Xdyad <- .coerce_Xdyad_3d(Xdyad)
		Xdyad_out <- vector("list", T)
		for(t in 1:T) {
			Xt <- Xdyad[[t]]
			if(is.null(Xt)) { Xdyad_out[[t]] <- NULL; next }
			pd_t <- dim(Xt)[3]
			Xt_new <- array(NA, dim = c(nA, nB, pd_t),
				dimnames = list(rowActorSet, colActorSet, dimnames(Xt)[[3]]))
			row_idx <- match(rownames(Xt), rowActorSet)
			col_idx <- match(colnames(Xt), colActorSet)
			for(i in seq_along(row_idx)) {
				for(j in seq_along(col_idx)) {
					if(!is.na(row_idx[i]) && !is.na(col_idx[j])) {
						Xt_new[row_idx[i], col_idx[j], ] <- Xt[i, j, ]
					}
				}
			}
			Xdyad_out[[t]] <- Xt_new
		}
		Xdyad <- Xdyad_out
	}
	####

	####
	# reorder Xrow to match sorted row actors and convert to 3D array
	if(!is.null(Xrow)) {
		pr <- ncol(Xrow[[1]])
		Xrow_array <- array(NA, dim = c(nA, pr, T),
			dimnames = list(rowActorSet, colnames(Xrow[[1]]), names(Y_list)))
		for(t in 1:T) {
			Xrt <- Xrow[[t]]
			if(is.null(Xrt)) next
			row_idx <- match(rownames(Xrt), rowActorSet)
			for(i in seq_along(row_idx)) {
				if(!is.na(row_idx[i])) Xrow_array[row_idx[i], , t] <- Xrt[i, ]
			}
		}
		Xrow <- Xrow_array
	}
	####

	####
	# reorder Xcol to match sorted column actors and convert to 3D array
	if(!is.null(Xcol)) {
		pc <- ncol(Xcol[[1]])
		Xcol_array <- array(NA, dim = c(nB, pc, T),
			dimnames = list(colActorSet, colnames(Xcol[[1]]), names(Y_list)))
		for(t in 1:T) {
			Xct <- Xcol[[t]]
			if(is.null(Xct)) next
			col_idx <- match(rownames(Xct), colActorSet)
			for(i in seq_along(col_idx)) {
				if(!is.na(col_idx[i])) Xcol_array[col_idx[i], , t] <- Xct[i, ]
			}
		}
		Xcol <- Xcol_array
	}
	####

	list(Y = Y_array, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)
}
