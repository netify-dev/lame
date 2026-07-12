# per-period bipartite interaction matrix g_t. when dynamic_g = true
# and the fit has ra, rb > 0, the linear predictor uses a separate
# g_t per period. the conditional posterior of vec(g_t) given u_t,
# v_t and the residual is gaussian:
#   vec(g_t) | . ~ n(m_t, s_t)
#   s_t^-1 = (v_t kron u_t)'(v_t kron u_t) / s^2 + prior precision.
# we canonicalise (u, g, v) per period by an svd of g_t.

#' @noRd
# sample g_t (ra x rb) conditional on u_t (na x ra), v_t (nb x rb), and
# the residual matrix e_t (na x nb). the prior on vec(g_t) is
# n(prior_mean, prior_cov * s^2) by default; we use prior_cov = i / lambda_g.
.sample_G_period <- function(E_t, U_t, V_t, s2,
                              lambda_G = 1.0,
                              prior_mean_vec = NULL) {
	RA <- ncol(U_t); RB <- ncol(V_t)
	if (RA == 0L || RB == 0L) return(matrix(0, RA, RB))
	# build design matrix x_des: (na*nb) x (ra*rb) so that
	# vec(e_t) = x_des %*% vec(g_t) (column-major vec)
	# each row corresponds to (i, j); the column for (r, s) of g is
	# u_t[i, r] * v_t[j, s].
	# equivalently: x_des = v_t kronecker u_t.
	X_des <- kronecker(V_t, U_t)
	y_vec <- as.numeric(E_t)
	ok <- is.finite(y_vec)
	X_obs <- X_des[ok, , drop = FALSE]
	y_obs <- y_vec[ok]

	prec <- crossprod(X_obs) / s2 + diag(1 / max(lambda_G, 1e-12),
	                                      RA * RB)
	rhs <- crossprod(X_obs, y_obs) / s2
	if (!is.null(prior_mean_vec) && length(prior_mean_vec) == RA * RB) {
		rhs <- rhs + prior_mean_vec / max(lambda_G, 1e-12)
	}
	L <- tryCatch(chol(prec), error = function(e) NULL)
	if (is.null(L)) {
		# tiny ridge and retry
		L <- chol(prec + diag(1e-6 * max(diag(prec), 1), nrow(prec)))
	}
	m_vec <- backsolve(L, forwardsolve(t(L), as.numeric(rhs)))
	z <- stats::rnorm(RA * RB)
	g_vec <- as.numeric(m_vec + backsolve(L, z))
	matrix(g_vec, RA, RB)
}

#' @noRd
# canonicalise (u, g, v) per period via svd of g_t so the parameterisation
# is identifiable. returns a list with canonical u/g/v cubes that
# preserve the linear predictor u %*% g %*% v'.
.canonicalize_UGV <- function(U_cube, G_cube, V_cube) {
	T_per <- if (length(dim(G_cube)) == 3L) dim(G_cube)[3L] else 1L
	U_can <- U_cube
	V_can <- V_cube
	G_can <- G_cube
	for (t in seq_len(T_per)) {
		G_t <- if (length(dim(G_cube)) == 3L) G_cube[, , t] else G_cube
		U_t <- if (length(dim(U_cube)) == 3L) U_cube[, , t] else U_cube
		V_t <- if (length(dim(V_cube)) == 3L) V_cube[, , t] else V_cube
		if (NCOL(U_t) == 0L || NCOL(V_t) == 0L) next
		sv <- tryCatch(svd(G_t), error = function(e) NULL)
		if (is.null(sv)) next
		# canonical: u_can_t = u_t %*% sv$u; g_can_t = diag(sv$d);
		#            v_can_t = v_t %*% sv$v
		# preserves u %*% g %*% v'
		if (length(dim(U_cube)) == 3L) {
			U_can[, , t] <- U_t %*% sv$u
			V_can[, , t] <- V_t %*% sv$v
			G_can[, , t] <- diag(sv$d, length(sv$d))
		} else if (t == 1L) {
			# static u/v: only canonicalise once
			U_can <- U_t %*% sv$u
			V_can <- V_t %*% sv$v
			G_can <- diag(sv$d, length(sv$d))
		}
	}
	list(U = U_can, G = G_can, V = V_can)
}
