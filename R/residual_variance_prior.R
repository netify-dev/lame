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
