# per-period bipartite interaction matrix g_t. when dynamic_g = true
# and the fit has ra, rb > 0, the linear predictor uses a separate
# g_t per period. the conditional posterior of vec(g_t) given u_t,
# v_t and the residual is gaussian:
#   vec(g_t) | . ~ n(m_t, s_t)
#   s_t^-1 = (v_t kron u_t)'(v_t kron u_t) / s^2 + prior precision.
# we canonicalise (u, g, v) per period by an svd of g_t.

#' @noRd
# Sample G_t (RA x RB) conditional on U_t (nA x RA), V_t (nB x RB), and
# the residual matrix E_t (nA x nB). The prior on vec(G_t) is
# N(prior_mean, prior_cov * s^2) by default; we use prior_cov = I / lambda_G.
.sample_G_period <- function(E_t, U_t, V_t, s2,
                              lambda_G = 1.0,
                              prior_mean_vec = NULL) {
	RA <- ncol(U_t); RB <- ncol(V_t)
	if (RA == 0L || RB == 0L) return(matrix(0, RA, RB))
	# Build design matrix X_des: (nA*nB) x (RA*RB) so that
	# vec(E_t) = X_des %*% vec(G_t) (column-major vec)
	# Each row corresponds to (i, j); the column for (r, s) of G is
	# U_t[i, r] * V_t[j, s].
	# Equivalently: X_des = V_t kronecker U_t.
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
# Canonicalise (U, G, V) per period via SVD of G_t so the parameterisation
# is identifiable. Returns a list with canonical U/G/V cubes that
# preserve the linear predictor U %*% G %*% V'.
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
		# canonical: U_can_t = U_t %*% sv$u; G_can_t = diag(sv$d);
		#            V_can_t = V_t %*% sv$v
		# preserves U %*% G %*% V'
		if (length(dim(U_cube)) == 3L) {
			U_can[, , t] <- U_t %*% sv$u
			V_can[, , t] <- V_t %*% sv$v
			G_can[, , t] <- diag(sv$d, length(sv$d))
		} else if (t == 1L) {
			# static U/V: only canonicalise once
			U_can <- U_t %*% sv$u
			V_can <- V_t %*% sv$v
			G_can <- diag(sv$d, length(sv$d))
		}
	}
	list(U = U_can, G = G_can, V = V_can)
}
