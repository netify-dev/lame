# Poisson log mass written as the log mass at the mode, log p(y | y), less
# y (expm1(delta) - delta) with delta = log(rate / y), so it stays exact near the mode and in the tails
.poisson_log_mass <- function(y, log_rate) {
	log_rate <- rep_len(log_rate, length(y))
	out <- rep(NA_real_, length(y))
	finite <- is.finite(log_rate) & is.finite(y)
	zero <- finite & y == 0
	out[zero] <- -exp(log_rate[zero])
	positive <- finite & y > 0
	yp <- y[positive]
	delta <- log_rate[positive] - log(yp)
	at_mode <- ifelse(
		yp == round(yp),
		stats::dpois(round(yp), round(yp), log = TRUE),
		yp * log(yp) - yp - lgamma(yp + 1)
	)
	out[positive] <- at_mode - yp * (expm1(delta) - delta)
	out[finite & y < 0] <- -Inf
	infinite <- !is.na(log_rate) & !is.finite(log_rate) & is.finite(y)
	out[infinite] <- -Inf
	out[infinite & log_rate < 0 & y == 0] <- 0
	out
}

# refuse negative or non-integer Poisson responses (no Poisson mass, so the latent Metropolis step
# could never move); Y is a matrix, array or list of matrices
.lame_check_poisson_counts <- function(Y) {
	num <- function(m) if (is.numeric(m) || is.logical(m)) as.vector(m) else NULL
	y <- if (is.list(Y) && !is.data.frame(Y)) unlist(lapply(Y, num), use.names = FALSE) else num(Y)
	y <- y[is.finite(y)]
	if (any(y < 0)) {
		cli::cli_abort(
			c(
				"{.arg family} = {.val poisson} requires non-negative counts.",
				"i" = "{.arg Y} has negative observed values."
			),
			class = "lame_poisson_counts"
		)
	}
	bad <- y != round(y)
	if (any(bad)) {
		ex <- utils::head(unique(y[bad]), 3L)
		cli::cli_abort(
			c(
				"{.arg family} = {.val poisson} requires integer counts.",
				"x" = "{.arg Y} has {sum(bad)} non-integer observed value{?s} (e.g. {.val {signif(ex, 6)}}).",
				"i" = "Round the counts if the fractions are numerical noise ({.code Y <- round(Y)}); for a continuous response use {.code family = \"normal\"}."
			),
			class = "lame_poisson_counts"
		)
	}
	invisible(NULL)
}

# invert a working log1p response; aborts only if the count mean overflows
.poisson_inverse_log1p <- function(eta) {
	# called for its overflow check only; the returned rate is not used
	.als_pois_rate(eta)
	pmax(expm1(eta), 0)
}
