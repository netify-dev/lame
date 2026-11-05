# conditional updates for explicit cutpoints in directed correlated dyads
# flat prior on the ordered cutpoints, as in the Cowles update
.sample_alpha_given_z <- function(alpha, Y_int, Z) {
	if (length(alpha) < 2L)
		return(list(alpha = alpha, delta = numeric(0), accept = FALSE))
	obs <- !is.na(Y_int)
	if (any(!is.finite(Z[obs])) || any(!is.finite(alpha)) ||
			alpha[1L] != 0 || any(diff(alpha) <= 0))
		cli::cli_abort("Invalid explicit ordinal cutpoint state.")
	bounds <- c(-Inf, alpha, Inf)
	if (any(Y_int[obs] != as.integer(Y_int[obs])) ||
			any(Y_int[obs] < 1L | Y_int[obs] > length(alpha) + 1L) ||
			any(Z[obs] <= bounds[Y_int[obs]] | Z[obs] >= bounds[Y_int[obs] + 1L]))
		cli::cli_abort("Explicit ordinal latent state is outside its category support.")
	for (k in seq.int(2L, length(alpha))) {
		left <- Z[obs & Y_int == k]
		right <- Z[obs & Y_int == k + 1L]
		if (!length(left) || !length(right))
			cli::cli_abort("Explicit ordinal cutpoints require observed adjacent categories.")
		lo <- max(left); hi <- min(right)
		if (!is.finite(hi - lo) || hi <= lo)
			cli::cli_abort("Explicit ordinal cutpoint conditional is numerically unresolved.")
		alpha[k] <- lo + (hi - lo) * stats::runif(1L)
		if (!is.finite(alpha[k]) || alpha[k] <= lo || alpha[k] >= hi)
			cli::cli_abort("Explicit ordinal cutpoint draw is numerically unresolved.")
	}
	list(alpha = alpha, delta = .alpha_to_delta(alpha), accept = TRUE)
}

.rZ_ord_explicit_correlated_fc <- function(Z, EZ, Y_int, alpha, rho) {
	n <- nrow(Z)
	if (!is.matrix(Z) || ncol(Z) != n || !identical(dim(Z), dim(EZ)) ||
			!identical(dim(Z), dim(Y_int)) || length(rho) != 1L ||
			!is.finite(rho) || abs(rho) >= 1 || any(!is.finite(Z)) ||
			any(!is.finite(EZ)) || any(!is.finite(alpha)) ||
			alpha[1L] != 0 || any(diff(alpha) <= 0))
		cli::cli_abort("Invalid correlated explicit ordinal latent-state inputs.")
	obs <- !is.na(Y_int) & row(Y_int) != col(Y_int)
	if (any(Y_int[obs] != as.integer(Y_int[obs])) ||
			any(Y_int[obs] < 1L | Y_int[obs] > length(alpha) + 1L))
		cli::cli_abort("Invalid explicit ordinal category codes.")
	bounds <- c(-Inf, alpha, Inf)
	sd <- sqrt((1 - rho) * (1 + rho))
	# sweep the directed cells, each conditioning on its partner's latest value
	for (j in seq_len(n)) for (i in seq_len(n)) {
		if (i == j) next
		mean <- EZ[i, j] + rho * (Z[j, i] - EZ[j, i])
		if (!is.finite(mean)) cli::cli_abort("Correlated ordinal conditional mean is outside numerical range.")
		if (obs[i, j]) {
			lo <- bounds[Y_int[i, j]]; hi <- bounds[Y_int[i, j] + 1L]
			x <- sd * rtnorm_interval_logp(mean / sd, lo / sd, hi / sd)
			if (!is.finite(x) || x <= lo || x >= hi)
				cli::cli_abort("Correlated ordinal truncated-normal draw is numerically unresolved.")
		} else x <- stats::rnorm(1L, mean, sd)
		Z[i, j] <- x
	}
	# diagonal cells are nuisance draws, as in the GLS update of beta
	diag(Z) <- stats::rnorm(n, diag(EZ), sqrt(1 + rho))
	if (any(!is.finite(Z))) cli::cli_abort("Correlated ordinal latent draw is outside numerical range.")
	Z
}
