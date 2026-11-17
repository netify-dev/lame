# the residual-variance prior is set before any update, since a rate from current residuals
# breaks the Gaussian mean conditional; normal uses response units, Poisson the log-rate scale
.lame_validate_s20 <- function(s20) {
	if (!is.numeric(s20) || length(s20) != 1L || !is.finite(s20) || s20 <= 0) {
		cli::cli_abort("{.code prior$s20} must be a finite positive scalar.")
	}
	as.numeric(s20)
}

.lame_residual_prior_scale <- function(Y, family, s20 = NULL) {
	if (!family %in% c("normal", "poisson")) return(NULL)
	if (!is.null(s20)) return(.lame_validate_s20(s20))
	# the response variance for normal, 1 on the Poisson log-rate scale (.lame_unit_var, normal_scale.R)
	.lame_unit_var(Y, family)
}

# prior scales g_j of the static coefficients on the dynamic-beta path, N(0, g_j * s2),
# indexed after scalar expansion
.lame_static_g_scale <- function(g, indices) {
	if (!length(indices)) return(numeric(0))
	out <- if (length(g) == 1L) rep(g, length(indices)) else g[indices]
	if (any(!is.finite(out) | out <= 0)) {
		cli::cli_abort("{.arg g} must supply a finite positive scale for every static coefficient.")
	}
	as.numeric(out)
}

# marginal precision of the coefficients `idx` under the zero-mean Gaussian prior with precision P
.lame_static_panel_precision <- function(P, idx) {
	rest <- setdiff(seq_len(nrow(P)), idx)
	out <- P[idx, idx, drop = FALSE]
	if (length(rest)) {
		out <- out - P[idx, rest, drop = FALSE] %*% solve(P[rest, rest, drop = FALSE], P[rest, idx, drop = FALSE])
	}
	(out + t(out)) / 2
}

# marginal prior variance of each coefficient of `idx` on its own under the zero-mean Gaussian prior with
# precision P: the variance it has as a static coefficient, which the long-run mean of an ar1 path and
# the first-period level of an rw2 or matern32 path take too
.lame_static_panel_var <- function(P, idx) {
	vapply(idx, function(k) 1 / .lame_static_panel_precision(P, k)[1L, 1L], numeric(1))
}

# prior precision of the static coefficients of a dynamic_beta fit: the panel precision P0, or
# N(0, g_j * s2) each when P0 is NULL
.lame_static_prior_prec <- function(P0, g, indices, s2) {
	if (!is.null(P0)) return(P0)
	diag(1 / (.lame_static_g_scale(g, indices) * s2), length(indices))
}

# terms the static coefficients of a dynamic_beta fit add to the s2 conditional (none under P0)
.lame_static_s2_terms <- function(P0, beta_static, g, indices) {
	if (!is.null(P0)) return(c(df = 0, ss = 0))
	.lame_s2_beta_terms(beta_static, .lame_static_g_scale(g, indices))
}

# prior of the static intercept, entry ks of the static block, given the other static coefficients
.lame_static_icpt_prior <- function(beta_static, ks, P0, g, indices, s2) {
	if (is.null(P0)) {
		return(list(kind = "static", value = beta_static[ks], prec = 1 / (.lame_static_g_scale(g, indices)[ks] * s2)))
	}
	q <- P0[ks, ks]
	list(kind = "static", value = beta_static[ks], prec = q, mean = -sum(P0[ks, -ks] * beta_static[-ks]) / q)
}

# terms that coefficients with prior variance proportional to s2 add to its conditional
.lame_s2_beta_terms <- function(beta, g) {
	if (!length(beta)) return(c(df = 0, ss = 0))
	if (length(g) == 1L) g <- rep(g, length(beta))
	if (length(g) != length(beta) || any(!is.finite(g) | g <= 0) ||
			any(!is.finite(beta))) cli::cli_abort("Invalid variance-scaled coefficient prior.")
	c(df = length(beta), ss = sum(beta^2 / g))
}

# linear positions of the diagonal cells of every period slice
.lame_diagonal_index <- function(Z) {
	n <- dim(Z)[1L]
	rep(seq.int(1L, by = n + 1L, length.out = n), dim(Z)[3L]) +
		rep(seq.int(0L, by = n * dim(Z)[2L], length.out = dim(Z)[3L]), each = n)
}

.lame_diagonal <- function(Z) {
	if (is.matrix(Z)) return(diag(Z))
	matrix(Z[.lame_diagonal_index(Z)], nrow = dim(Z)[1L], ncol = dim(Z)[3L])
}

# a collapsed update integrates the diagonal self-relations out; refresh them under the
# new covariance before a mean update reads them
.lame_refresh_diagonal <- function(Z, mean_diagonal, rho, s2) {
	if (is.matrix(Z)) {
		diag(Z) <- rnorm(nrow(Z), mean_diagonal, sqrt((1 + rho) * s2))
		return(Z)
	}
	n <- dim(Z)[1L]; nt <- dim(Z)[3L]
	mean_diagonal <- matrix(mean_diagonal, n, nt)
	draws <- vapply(
		seq_len(nt),
		function(t) rnorm(
			n,
			mean_diagonal[, t],
			sqrt((1 + .lame_rho_at(rho, t)) * s2)
		),
		numeric(n)
	)
	Z[.lame_diagonal_index(Z)] <- draws
	Z
}
