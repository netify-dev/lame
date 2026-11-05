# keep an AR coefficient proposal only inside the joint prior support
# (rejection preserves detailed balance; clipping to the boundary would not)
.lame_ar_rho_support <- function(
	candidate,
	current,
	sigma,
	stat_var_cap = Inf,
	bound = 0.99
) {
	if (!is.finite(candidate) || abs(candidate) > bound) return(current)
	if (is.finite(stat_var_cap) && sigma^2 > stat_var_cap * (1 - candidate^2))
		return(current)
	candidate
}

# reject residual variance draws below the bound implied by the stationary
# variance caps (a lazy Gibbs step); unbounded models return the draw unchanged
.lame_ar_s2_support <- function(
	candidate,
	current,
	rho_ab = NULL,
	sigma_ab = NULL,
	ab_cap = Inf,
	rho_uv = NULL,
	sigma_uv = NULL,
	uv_cap = Inf
) {
	lower <- 0
	if (is.finite(ab_cap) && !is.null(sigma_ab))
		lower <- max(lower, sigma_ab^2 / (ab_cap * (1 - rho_ab^2)))
	if (is.finite(uv_cap) && !is.null(sigma_uv))
		lower <- max(lower, sigma_uv^2 / (uv_cap * (1 - rho_uv^2)))
	if (is.finite(candidate) && candidate >= lower) candidate else current
}

.lame_ar_start_support <- function(
	rho,
	sigma,
	stat_var_cap = Inf,
	sigma_max = Inf,
	resume = FALSE
) {
	valid <- is.finite(rho) && abs(rho) <= .99 && is.finite(sigma) && sigma > 0 &&
		sigma <= sigma_max && (!is.finite(stat_var_cap) || sigma^2 <= stat_var_cap * (1 - rho^2))
	if (valid) return(list(rho = rho, sigma = sigma))
	if (resume) cli::cli_abort(
		"Checkpoint AR hyperparameters lie outside the configured joint prior support; start a new chain.", call = NULL
	)
	if (!is.finite(rho) || !is.finite(sigma) || sigma <= 0)
		cli::cli_abort(
			"Initial AR hyperparameters must have finite persistence and a positive innovation scale.", call = NULL
		)
	rho <- max(-.99, min(.99, rho))
	upper <- sigma_max
	if (is.finite(stat_var_cap)) upper <- min(upper, sqrt(stat_var_cap * (1 - rho^2)))
	sigma <- min(sigma, upper * (1 - 1e-8))
	list(rho = rho, sigma = sigma)
}
