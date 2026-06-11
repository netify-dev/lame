#' Helper functions for dynamic effects in LAME
#' 
#' Internal functions to support dynamic additive and multiplicative effects

#' Compute expected values with dynamic additive effects
#' @keywords internal
#' @noRd
get_EZ_unip_time <- function(Xlist, beta, ab, U, V, N = NULL) {
	U_is_dynamic <- length(dim(U)) == 3L
	V_is_dynamic <- length(dim(V)) == 3L
	ab_is_dynamic <- length(dim(ab)) == 3L
	if (!U_is_dynamic && !V_is_dynamic && !ab_is_dynamic) {
		return(get_EZ_cpp(Xlist, beta, ab, U, V))
	}
	if (is.null(N)) {
		N <- if (U_is_dynamic) dim(U)[3L] else if (V_is_dynamic) dim(V)[3L] else dim(ab)[3L]
	}
	nr <- if (U_is_dynamic) dim(U)[1L] else nrow(U)
	nc <- if (V_is_dynamic) dim(V)[1L] else nrow(V)
	out <- array(0, dim = c(nr, nc, N))
	for (t in seq_len(N)) {
		Ut <- if (U_is_dynamic) matrix(U[, , t], nrow = nr) else U
		Vt <- if (V_is_dynamic) matrix(V[, , t], nrow = nc) else V
		abt <- if (ab_is_dynamic) ab[, , t] else ab
		out[, , t] <- get_EZ_cpp(list(Xlist[[t]]), beta, abt, Ut, Vt)[, , 1L]
	}
	out
}

#' Compute expected values with dynamic additive effects
#' @keywords internal
#' @noRd
get_EZ_dynamic_ab <- function(Xlist, beta, a_mat, b_mat, U, V, N,
                              bip = FALSE, G = NULL, nA = NULL, nB = NULL) {
	if (bip) {
		# bipartite path: use get_ez_bip_cpp
		nA <- nA %||% nrow(a_mat)
		nB <- nB %||% nrow(b_mat)

		# build base cube from x * beta
		base_cube <- array(0, dim = c(nA, nB, N))
		for (t in 1:N) {
			if (length(beta) > 0 && !is.null(Xlist[[t]])) {
				n_covs <- min(length(beta), dim(Xlist[[t]])[3])
				for (k in 1:n_covs) {
					base_cube[,,t] <- base_cube[,,t] + beta[k] * Xlist[[t]][,,k]
				}
			}
		}

		if (is.null(dim(U)) || length(dim(U)) == 2) {
			RA <- if (is.null(dim(U))) 1L else ncol(U)
			U_cube <- array(0, dim = c(nA, RA, N))
			U_2d <- if (is.matrix(U)) U else matrix(U, nA, RA)
			for (t in 1:N) U_cube[,,t] <- U_2d
		} else {
			U_cube <- U
		}
		if (is.null(dim(V)) || length(dim(V)) == 2) {
			RB <- if (is.null(dim(V))) 1L else ncol(V)
			V_cube <- array(0, dim = c(nB, RB, N))
			V_2d <- if (is.matrix(V)) V else matrix(V, nB, RB)
			for (t in 1:N) V_cube[,,t] <- V_2d
		} else {
			V_cube <- V
		}

		if (is.null(G)) G <- matrix(0, dim(U_cube)[2], dim(V_cube)[2])

		EZ <- get_EZ_bip_cpp(base_cube, a_mat, b_mat, U_cube, V_cube, G)
	} else {
		# unipartite path
		EZ <- array(0, dim = c(nrow(a_mat), nrow(a_mat), N))
		for (t in 1:N) {
			ab_t <- outer(a_mat[,t], b_mat[,t], "+")
			Ut <- if (length(dim(U)) == 3L) matrix(U[, , t], nrow = dim(U)[1L]) else U
			Vt <- if (length(dim(V)) == 3L) matrix(V[, , t], nrow = dim(V)[1L]) else V
			EZ[,,t] <- get_EZ_cpp(list(Xlist[[t]]), beta, ab_t, Ut, Vt)[,,1]
		}
	}

	return(EZ)
}

#' Extract time-specific or averaged additive effects
#' @keywords internal
#' @noRd
get_ab_for_time <- function(a_mat, b_mat, t=NULL, average=FALSE) {
	if(is.null(a_mat) || is.null(b_mat)) {
		return(list(a=0, b=0))
	}
	
	if(average || is.null(t)) {
		# time-averaged
		return(list(a=rowMeans(a_mat), b=rowMeans(b_mat)))
	} else {
		# time-specific
		return(list(a=a_mat[,t], b=b_mat[,t]))
	}
}
