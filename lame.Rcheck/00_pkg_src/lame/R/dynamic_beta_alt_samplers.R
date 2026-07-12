# joint gaussian samplers for the dynamic_beta path when
# kind in {rw2, matern32}, plus a gap-aware ar(1) sampler used when
# time_index has unequal gaps. each one builds a t-by-t prior precision
# k and samples per coefficient from n(m, h^{-1}) where
#   h = blockdiag(x_t' x_t / s^2) + k / sigma_k^2
# rw2 uses the second-difference operator; matern 3/2 uses the
# state-space kernel; ar(1) gap uses the per-gap transition factors.
# triangular solves use forwardsolve(t(l), .) and backsolve(l, .).

#' @noRd
.rw2_precision <- function(T_per, jitter = 1e-6) {
	if (T_per < 3L) {
		# rw2 requires t >= 3; fall back to a diagonal ridge so the sampler
		# doesn't blow up.
		return(diag(jitter, T_per))
	}
	D <- matrix(0, T_per - 2L, T_per)
	for (i in seq_len(T_per - 2L)) {
		D[i, i]     <- 1
		D[i, i + 1L] <- -2
		D[i, i + 2L] <- 1
	}
	K <- crossprod(D)
	# rank deficient by 2; add a tiny ridge so the cholesky succeeds
	K + diag(jitter, T_per)
}

#' @noRd
.matern32_precision <- function(T_per, length_scale = NULL,
                                 time_positions = NULL, jitter = 1e-8) {
	if (is.null(length_scale) || !is.finite(length_scale) || length_scale <= 0) {
		length_scale <- max(2, T_per / 4)
	}
	tt <- if (is.null(time_positions)) seq_len(T_per) else as.numeric(time_positions)
	dst <- abs(outer(tt, tt, "-"))
	a <- sqrt(3) * dst / length_scale
	Sigma <- (1 + a) * exp(-a)
	# stable inverse via cholesky with a tiny ridge
	Sigma <- Sigma + diag(jitter, T_per)
	chol2inv(chol(Sigma))
}

#' @noRd
# gap-aware ar(1) precision. with strictly-increasing time positions
# and gaps delta_t = time_positions[t+1] - time_positions[t], the per-
# transition conditional variance is
#   q_t = sigma^2 * (1 - rho^{2 delta_t}) / (1 - rho^2)
# and the transition factor is f_t = rho^{delta_t}. the implied joint
# precision of (beta_1, ..., beta_t) (before dividing by sigma^2) is
# tridiagonal with the off-diagonal terms encoding -f_t / q_t per step
# and the diagonal pieces accumulating "from-prev" and "to-next" terms.
# the returned matrix is the kind-precision k such that the prior is
#   beta ~ n(0, sigma^2 k^{-1}).
.ar1_precision <- function(T_per, gaps = NULL, rho = 0, jitter = 1e-10) {
	if (T_per < 2L) return(matrix(1, 1L, 1L))
	if (is.null(gaps)) gaps <- rep(1, T_per - 1L)
	# clamp rho and per-step quantities
	rho <- max(min(rho, 0.999999), -0.999999)
	# stationary variance
	v_stat <- 1 / (1 - rho^2)
	K <- matrix(0, T_per, T_per)
	for (t in seq_len(T_per - 1L)) {
		dt <- as.numeric(gaps[t])
		Ft <- rho^dt
		qt <- (1 - rho^(2 * dt)) / (1 - rho^2)
		# contribution to (t, t), (t+1, t+1), (t, t+1), (t+1, t)
		# from the term -0.5 * (beta_{t+1} - ft * beta_t)^2 / qt
		K[t,     t]     <- K[t,     t]     + Ft^2 / qt
		K[t + 1, t + 1] <- K[t + 1, t + 1] + 1 / qt
		K[t,     t + 1] <- K[t,     t + 1] - Ft / qt
		K[t + 1, t]     <- K[t + 1, t]     - Ft / qt
	}
	# stationary initial condition: beta_1 ~ n(0, v_stat) -> precision 1/v_stat
	K[1, 1] <- K[1, 1] + 1 / v_stat
	K + diag(jitter, T_per)
}

#' @noRd
.dynamic_beta_kind_to_prec <- function(kind, T_per, length_scale = NULL,
                                        time_positions = NULL,
                                        gaps = NULL, rho = NULL) {
	# fall back to integer positions when not supplied
	if (is.null(time_positions)) time_positions <- seq_len(T_per)
	if (is.null(gaps)) gaps <- diff(time_positions)
	switch(kind,
		ar1      = .ar1_precision(T_per, gaps = gaps, rho = if (is.null(rho)) 0 else rho),
		rw1      = .ar1_precision(T_per, gaps = gaps, rho = 0.999999),
		rw2      = .rw2_precision(T_per),
		matern32 = .matern32_precision(T_per, length_scale,
		                                time_positions = time_positions),
		stop(sprintf("kind '%s' not handled by alt sampler", kind)))
}

#' @noRd
# build the per-period (x'x, x'r) sufficient statistics for the dynamic block.
# mirrors what `sample_beta_dynamic_cpp` does internally, in r, but only for
# the dynamic coefficients (the static block is handled separately).
#
# each (x_t, z_t, offset_t, beta_static, x^static_t) gives a contribution
#   r_t = z_t - offset_t - x^static_t %*% beta_static
#   xtx_t = x_t' x_t   (sum over observed dyads)
#   xty_t = x_t' r_t   (sum over observed dyads)
.sample_beta_dynamic_alt <- function(
		Xdyn_list, Xstat_list, Z_list, offset_list,
		beta_static, sigma_by_coef, Lambda,
		s2, kind, length_scale = NULL,
		time_positions = NULL, rho_by_coef = NULL,
		dyad_rho = 0, use_dyad_rho = FALSE,
		bipartite = FALSE, symmetric = FALSE) {

	T_per <- length(Z_list)
	p_dyn <- if (length(Xdyn_list) > 0L && !is.null(Xdyn_list[[1L]]))
		ncol(Xdyn_list[[1L]]) else length(sigma_by_coef)
	if (p_dyn == 0L) return(list(path = matrix(0, T_per, 0L), chol_fail = 0L))

	# canonicalise sigma_by_coef to per-coefficient
	sigma_by_coef <- as.numeric(sigma_by_coef)
	if (length(sigma_by_coef) == 1L)
		sigma_by_coef <- rep(sigma_by_coef, p_dyn)

	# per-period sufficient statistics
	XtX_by_t <- vector("list", T_per)
	Xtr_by_t <- vector("list", T_per)

	for (t in seq_len(T_per)) {
		Xd_t <- Xdyn_list[[t]]
		Zt   <- Z_list[[t]]
		off  <- offset_list[[t]]
		if (is.null(Xd_t) || !length(Xd_t)) {
			XtX_by_t[[t]] <- matrix(0, p_dyn, p_dyn)
			Xtr_by_t[[t]] <- numeric(p_dyn)
			next
		}
		r_t <- as.numeric(Zt) - as.numeric(off)
		if (length(beta_static) > 0L && length(Xstat_list) >= t &&
		    !is.null(Xstat_list[[t]]) && ncol(Xstat_list[[t]]) > 0L) {
			r_t <- r_t - as.numeric(Xstat_list[[t]] %*% beta_static)
		}
		# observed-dyad mask: na in z/offset -> drop
		ok <- is.finite(r_t)
		if (!any(ok)) {
			XtX_by_t[[t]] <- matrix(0, p_dyn, p_dyn)
			Xtr_by_t[[t]] <- numeric(p_dyn)
			next
		}
		X_obs <- Xd_t[ok, , drop = FALSE]
		r_obs <- r_t[ok]
			# data precision is 1/s2 per dyad; this sampler ignores dyad correlation
		XtX_by_t[[t]] <- crossprod(X_obs) / s2
		Xtr_by_t[[t]] <- as.numeric(crossprod(X_obs, r_obs)) / s2
	}

	# joint posterior precision per dynamic coefficient k. we treat
	# cross-coef terms as zero (diagonal lambda) -- matches the typical
	# default path. the cross-period precision comes from the kind prior.
	# when kind is ar1/rw1 the prior precision depends on rho; we build
	# a per-coefficient k when rho_by_coef is supplied.
	K_prior_default <- .dynamic_beta_kind_to_prec(
		kind, T_per, length_scale,
		time_positions = time_positions,
		rho = if (!is.null(rho_by_coef) && length(rho_by_coef) >= 1L)
			mean(rho_by_coef) else 0)

	path <- matrix(0, T_per, p_dyn)
	chol_fail <- 0L

	for (k in seq_len(p_dyn)) {
		# h_k: t x t precision
		# diagonal contribution from data per period: xtx_by_t[[t]][k,k]
		# (we ignore cross-coef terms because lambda is assumed diagonal)
		# cross-period contribution: k_prior / sigma_k^2
		# per-coefficient k when rho varies across coefs (ar1/rw1 only)
		K_prior <- if (!is.null(rho_by_coef) && kind %in% c("ar1", "rw1")) {
			.dynamic_beta_kind_to_prec(
				kind, T_per, length_scale,
				time_positions = time_positions,
				rho = if (kind == "rw1") 0.999999 else rho_by_coef[k])
		} else K_prior_default
		diag_data <- vapply(XtX_by_t, function(M) M[k, k], numeric(1))
		H_k <- diag(diag_data, T_per) + K_prior / max(sigma_by_coef[k]^2, 1e-12)
		# data-driven precision-weighted mean
		b_k <- vapply(Xtr_by_t, function(v) v[k], numeric(1))
		# subtract cross-coef contributions if lambda is non-diagonal:
		# data term adjusted by other-coef per-period draws if needed.
		# we approximate by holding the other-coef draws at zero on the
		# first sweep; iterate to convergence below if p_dyn > 1.
		# sample via cholesky: l l' = h_k
		L <- tryCatch(chol(H_k), error = function(e) NULL)
		if (is.null(L)) {
			# tiny ridge and retry
			H_k_ridge <- H_k + diag(1e-6 * max(diag(H_k), 1), T_per)
			L <- tryCatch(chol(H_k_ridge), error = function(e) NULL)
			if (is.null(L)) {
				chol_fail <- chol_fail + 1L
				next
			}
		}
		# m = h^{-1} b solved via the cholesky (r chol returns upper-triangular)
		m <- backsolve(L, forwardsolve(t(L), b_k))
		# sample noise eps ~ n(0, h^{-1}) via l^{-t} z, z ~ n(0, i)
		z <- stats::rnorm(T_per)
		eps <- backsolve(L, z)
		path[, k] <- m + eps
	}

	# refine cross-coefficient terms by gibbs sweep when p_dyn > 1, using the
	# correlated data contribution. two extra sweeps are enough in practice.
	if (p_dyn > 1L) {
		for (sweep in seq_len(2L)) {
			for (k in seq_len(p_dyn)) {
				K_prior <- if (!is.null(rho_by_coef) && kind %in% c("ar1", "rw1")) {
					.dynamic_beta_kind_to_prec(
						kind, T_per, length_scale,
						time_positions = time_positions,
						rho = if (kind == "rw1") 0.999999 else rho_by_coef[k])
				} else K_prior_default
				diag_data <- vapply(XtX_by_t, function(M) M[k, k], numeric(1))
				H_k <- diag(diag_data, T_per) + K_prior / max(sigma_by_coef[k]^2, 1e-12)
				# residual data term: subtract x'x[k, -k] * beta_{-k}
				adj_data <- numeric(T_per)
				for (t in seq_len(T_per)) {
					M <- XtX_by_t[[t]]
					adj_data[t] <- sum(M[k, -k] * path[t, -k])
				}
				b_k <- vapply(Xtr_by_t, function(v) v[k], numeric(1)) - adj_data
				L <- tryCatch(chol(H_k), error = function(e) NULL)
				if (is.null(L)) {
					chol_fail <- chol_fail + 1L
					next
				}
				m <- backsolve(L, forwardsolve(t(L), b_k))
				z <- stats::rnorm(T_per)
				eps <- backsolve(L, z)
				path[, k] <- m + eps
			}
		}
	}

	list(path = path, chol_fail = chol_fail)
}

# --------------------------------------------------------------------------
# hyperparameter update for non-ar(1) kinds: when kind is "rw2"/"matern32",
# rho is not meaningful. sigma^2 still has a conjugate inverse-gamma update
# conditional on the path. the contribution to the ig posterior is
#     shape += t/2 (or (t-2)/2 for rw2 because of rank-deficiency)
#     scale += 0.5 * t(beta_k) k beta_k
#
# we rely on the same prior_shape / prior_scale used elsewhere.
.sample_sigma_beta_alt <- function(beta_path, kind, length_scale = NULL,
                                    prior_shape, prior_scale) {
	T_per <- nrow(beta_path)
	p_dyn <- ncol(beta_path)
	K <- .dynamic_beta_kind_to_prec(kind, T_per, length_scale)
	# rank correction: rw2 has rank t - 2; matern32 has full rank
	df_per_coef <- switch(kind,
		rw2      = max(0L, T_per - 2L),
		matern32 = T_per,
		T_per)
	out <- numeric(p_dyn)
	for (k in seq_len(p_dyn)) {
		bk <- beta_path[, k]
		quad <- as.numeric(crossprod(bk, K %*% bk))
		shape_post <- prior_shape + df_per_coef / 2
		scale_post <- prior_scale + 0.5 * quad
		out[k] <- 1 / stats::rgamma(1L, shape = shape_post, rate = scale_post)
	}
	out
}
