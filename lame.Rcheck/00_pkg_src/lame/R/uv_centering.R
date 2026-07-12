# reporting-only centred copies of u, v, and beta. canonical slots keep
# the sampled parameterisation so eta-reconstruction is unambiguous;
# the centred slots subtract column means and absorb the shift into
# the intercept column of beta.

#' @noRd
# centre u and v per period, absorbing the shift into the intercept of
# beta. returns a list(u_centered, v_centered, beta_centered).
#
# for static u/v (length(dim(u)) == 2), no per-period structure to centre;
# we still emit a copy as the "canonical" centred form for api symmetry.
.compute_uv_centered <- function(fit) {
	U <- fit$U
	V <- fit$V
	BETA <- fit$BETA
	if (is.null(U) || is.null(V) || NCOL(U) == 0L || NCOL(V) == 0L) {
		return(list(U_centered = U, V_centered = V, BETA_centered = BETA))
	}

	# bipartite products run through the g mixing matrix and the row /
	# column ranks may differ, so the square-uv' shift arithmetic below
	# does not apply; keep the sampled parameterisation for those fits.
	# exact [["G"]] indexing: $G would partial-match $GOF
	if (!is.null(fit[["G"]])) {
		return(list(U_centered = U, V_centered = V, BETA_centered = BETA))
	}

	beta_dyn <- length(dim(BETA)) == 3L
	beta_names <- if (beta_dyn) dimnames(BETA)[[2]] else colnames(BETA)
	has_int <- !is.null(beta_names) && beta_names[1L] == "intercept"

	if (length(dim(U)) == 3L) {
		# dynamic_uv: per-period u[, , t], v[, , t]
		T_per <- dim(U)[3L]
		U_c <- U
		V_c <- V
		# shift per period: mu_u_t = colmeans(u_t), mu_v_t = colmeans(v_t)
		# the u v' contribution is invariant under the joint shift
		#   u <- u - 1 mu_u', v unchanged, with eta correction = mu_u' (1 v') = mu_u' rowsums(v)
		# cleanest centring: subtract column means from each
		shift_eta <- matrix(0, nrow = nrow(U), ncol = nrow(V))
		shift_total_per_t <- numeric(T_per)
		for (t in seq_len(T_per)) {
			mu_u <- colMeans(U[, , t])
			mu_v <- colMeans(V[, , t])
			U_c[, , t] <- sweep(U[, , t], 2L, mu_u, "-")
			V_c[, , t] <- sweep(V[, , t], 2L, mu_v, "-")
			# residual shift absorbed into the intercept term: the per-cell
			# uv product changes by
			#   mu_u %*% t(v) + u %*% mu_v - mu_u %*% rep(1, n) * mu_v
			# averaged over the dyads, the per-period scalar shift is
			scalar_shift <- sum(mu_u * mu_v)
			shift_total_per_t[t] <- scalar_shift
		}
		# absorb the per-period scalar shift into the intercept column of beta
		BETA_c <- BETA
		if (has_int && beta_dyn) {
			for (t in seq_len(T_per)) {
				BETA_c[, 1L, t] <- BETA[, 1L, t] + shift_total_per_t[t]
			}
		}
		return(list(U_centered = U_c, V_centered = V_c, BETA_centered = BETA_c))
	}

	# static u/v: 2-d
	mu_u <- colMeans(U)
	mu_v <- colMeans(V)
	U_c <- sweep(U, 2L, mu_u, "-")
	V_c <- sweep(V, 2L, mu_v, "-")
	scalar_shift <- sum(mu_u * mu_v)
	BETA_c <- BETA
	if (has_int) {
		if (beta_dyn) {
			BETA_c[, 1L, ] <- BETA[, 1L, ] + scalar_shift
		} else {
			BETA_c[, 1L] <- BETA[, 1L] + scalar_shift
		}
	}
	list(U_centered = U_c, V_centered = V_c, BETA_centered = BETA_c)
}
