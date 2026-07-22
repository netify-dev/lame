#' Bipartite network helper functions
#'
#' @name bipartite_helpers
#' @keywords internal

#' Initialize bipartite starting values
#' @return List of starting values for bipartite MCMC
#' @keywords internal
init_bipartite_startvals <- function(
	Y, family, nA, nB, RA, RB, Tn,
	Xlist = NULL, odmax = NULL
	){

	####
	# initialize z. for the rank-family bipartite paths (cbin / frn)
	# we must seed z with a cone-respecting prior draw rather than leave it as
	# na; the row-cone samplers in r/rz_bipartite.r sample conditional on the
	# bounds, which are themselves a function of the surrounding z, and a fully-
	# na seed produces all-na bounds -> no movement. mirrors the cross-sectional
	# init at r/ame_bipartite.r:117-150.
	Z <- array(dim = c(nA, nB, Tn))

	for(t in 1:Tn) {
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
			# log(y + 1) for observed cells; the unipartite init in
			# r/get_start_vals.r imputes nan cells from the prior. mirror that
			# here so the sampler does not crash on the first sweep when y has
			# na cells (which is normal for partial-observation bipartite data).
			Z[,,t] <- log(Yt + 1)
			miss <- is.na(Z[,,t])
			if (any(miss)) {
				slice <- Z[,,t]
				slice[miss] <- stats::rnorm(sum(miss),
				                              mean = stats::median(slice, na.rm = TRUE),
				                              sd   = 1)
				Z[,,t] <- slice
			}
		} else if(family == "ordinal") {
			# rank-normal init on observed cells, ignoring na cells (which the
			# downstream `rz_ord_bip_fc` resamples from data-induced bounds)
			Z[,,t] <- matrix(qnorm((Yt + 0.5) / (max(Yt, na.rm = TRUE) + 1)),
			                  nA, nB)
			miss <- is.na(Z[,,t])
			if (any(miss)) Z[,,t][miss] <- stats::rnorm(sum(miss), 0, 1)
		} else if(family == "cbin") {
			# y=1 -> z > 0 (half-normal); y=0 -> z < 0; na -> unconstrained
			slice <- matrix(stats::rnorm(nA * nB), nA, nB)
			pos <- !is.na(Yt) & Yt == 1
			neg <- !is.na(Yt) & Yt == 0
			slice[pos] <-  abs(slice[pos])
			slice[neg] <- -abs(slice[neg])
			Z[,,t] <- slice
		} else if(family == "frn") {
			# rank-aware seed: top-ranked nominees positive and ordered;
			# non-nominees below all nominees. mirrors ame_bipartite.r:130-137.
			slice <- matrix(stats::rnorm(nA * nB), nA, nB)
			for (i in seq_len(nA)) {
				if (any(Yt[i, ] > 0, na.rm = TRUE)) {
					ranked <- which(Yt[i, ] > 0)
					slice[i, ranked] <- sort(abs(stats::rnorm(length(ranked))),
					                          decreasing = TRUE)
				}
			}
			Z[,,t] <- slice
		}
	}
	####

	####
	# initialize parameters
	a <- rep(0, nA)
	b <- rep(0, nB)

	if(RA > 0) {
		U <- array(rnorm(nA * RA * Tn, 0, 0.1), c(nA, RA, Tn))
	} else {
		U <- array(0, c(nA, 1, Tn))
	}

	if(RB > 0) {
		V <- array(rnorm(nB * RB * Tn, 0, 0.1), c(nB, RB, Tn))
	} else {
		V <- array(0, c(nB, 1, Tn))
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
	Tn <- dims[3]

	R <- Z - EZ_without_ab

	####
	# sample row effects (a)
	a <- matrix(0, nA, Tn)
	prec_a <- Tn / s2 + 1 / sigma2_a

	for(i in 1:nA) {
		for(t in 1:Tn) {
			mean_a <- sum(R[i,,t], na.rm = TRUE) / (nB * s2)
			mean_a <- mean_a / prec_a
			a[i,t] <- rnorm(1, mean_a, sqrt(1/prec_a))
		}
	}

	# update residuals
	for(t in 1:Tn) {
		R[,,t] <- R[,,t] - a[,t]
	}
	####

	####
	# sample column effects (b)
	b <- matrix(0, nB, Tn)
	prec_b <- Tn / s2 + 1 / sigma2_b

	for(j in 1:nB) {
		for(t in 1:Tn) {
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
	Tn <- ncol(a)

	# sample row variance
	shape_a <- (eta0_a + nA * Tn) / 2
	rate_a <- (eta0_a * Sab0_aa + sum(a^2)) / 2
	sigma2_a <- 1 / rgamma(1, shape_a, rate_a)

	# sample column variance
	shape_b <- (eta0_b + nB * Tn) / 2
	rate_b <- (eta0_b * Sab0_bb + sum(b^2)) / 2
	sigma2_b <- 1 / rgamma(1, shape_b, rate_b)

	list(sigma2_a = sigma2_a, sigma2_b = sigma2_b)
}

#' Compute GOF statistics for bipartite networks
#' @return List of observed and simulated GOF statistics
#' @keywords internal
compute_gof_bipartite <- function(Y_obs, Y_sim, family){
	Tn <- dim(Y_obs)[3]

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
	Tn <- length(Y_list)

	####
	# fill y array
	Y_array <- array(NA, c(nA, nB, Tn),
		dimnames = list(rowActorSet, colActorSet, names(Y_list)))

	for(t in 1:Tn) {
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

	# delegate to the package-internal `.lame_apply_suffix()` in r/globals.r
	# so the suffix rule is identical across the five naming sites (this
	# function, design_array, design_array_listwisedel, lame.r, ame_bipartite,
	# ame_als).
	.bip_apply_suffix <- .lame_apply_suffix

	####
	# coerce and reorder xdyad to match sorted actor sets
	if(!is.null(Xdyad)) {
		Xdyad <- .coerce_Xdyad_3d(Xdyad)
		Xdyad_out <- vector("list", Tn)
		for(t in 1:Tn) {
			Xt <- Xdyad[[t]]
			if(is.null(Xt)) { Xdyad_out[[t]] <- NULL; next }
			pd_t <- dim(Xt)[3]
			# append `_dyad` suffix to covariate-slice names so the resulting
			# beta carries the same `<name>_dyad` convention the unipartite
			# design produces (r/design_array.r)
			slice_nms <- dimnames(Xt)[[3]]
			if (is.null(slice_nms)) slice_nms <- paste0("Xdyad", seq_len(pd_t))
			slice_nms <- .bip_apply_suffix(slice_nms, "dyad")
			Xt_new <- array(NA, dim = c(nA, nB, pd_t),
				dimnames = list(rowActorSet, colActorSet, slice_nms))
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
	# reorder xrow to match sorted row actors and convert to 3d array
	if(!is.null(Xrow)) {
		pr <- ncol(Xrow[[1]])
		row_nms <- .bip_apply_suffix(
			if (is.null(colnames(Xrow[[1]]))) paste0("Xrow", seq_len(pr))
			else colnames(Xrow[[1]]), "row")
		Xrow_array <- array(NA, dim = c(nA, pr, Tn),
			dimnames = list(rowActorSet, row_nms, names(Y_list)))
		for(t in 1:Tn) {
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
	# reorder xcol to match sorted column actors and convert to 3d array
	if(!is.null(Xcol)) {
		pc <- ncol(Xcol[[1]])
		col_nms <- .bip_apply_suffix(
			if (is.null(colnames(Xcol[[1]]))) paste0("Xcol", seq_len(pc))
			else colnames(Xcol[[1]]), "col")
		Xcol_array <- array(NA, dim = c(nB, pc, Tn),
			dimnames = list(colActorSet, col_nms, names(Y_list)))
		for(t in 1:Tn) {
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
