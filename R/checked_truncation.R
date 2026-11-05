.lame_checked_truncated_gaussian <- function(mean, sd, lower, upper) {
	n <- length(mean)
	if (!is.numeric(mean) || any(!is.finite(mean)) || length(sd) != 1L ||
			!is.numeric(sd) || !is.finite(sd) || sd <= 0 ||
			!length(lower) || !length(upper) || anyNA(lower) || anyNA(upper)) {
		cli::cli_abort(
			"A truncated Gaussian conditional requires finite means, a positive finite scale and nonmissing bounds.",
			class = "lame_truncation_parameter_error"
		)
	}
	lower <- rep_len(lower, n)
	upper <- rep_len(upper, n)
	if (any(lower >= upper)) {
		cli::cli_abort(
			"The current latent state gives an infeasible truncation interval.",
			class = "lame_truncation_support_error"
		)
	}
	mu <- mean / sd
	lo <- lower / sd
	hi <- upper / sd
	if (any(!is.finite(mu)) || any(is.finite(lower) & !is.finite(lo)) ||
			any(is.finite(upper) & !is.finite(hi)) || any(lo >= hi)) {
		cli::cli_abort(
			"The truncated Gaussian conditional cannot be represented at its actual scale.",
			class = "lame_truncation_numerical_error"
		)
	}
	# draw on the standardized scale, with no floor on sd
	draw <- sd * rtnorm_interval_logp(mu, lo, hi)
	if (any(!is.finite(draw)) || any(draw <= lower | draw >= upper)) {
		cli::cli_abort(
			"A truncated Gaussian draw is nonfinite or outside its conditional support.",
			class = "lame_truncation_numerical_error"
		)
	}
	draw
}
