#' Metropolis-Hastings update for the within-dyad correlation
#'
#' Draws a new value of the dyadic correlation rho from its full conditional
#' using a random-walk-style Metropolis-Hastings step. Each dyad's pair of
#' directed residuals is modelled as bivariate normal with unit variances
#' (after scaling by the dyadic variance) and correlation rho; the proposal is
#' a normal truncated to (-1, 1) and the acceptance ratio combines the
#' bivariate-normal log-likelihood, the prior of rho, and the
#' truncated-proposal normaliser.
#'
#' @param Z n x n normal relational matrix.
#' @param rho current value of rho.
#' @param s2 current value of the dyadic variance.
#' @param offset matrix matching \code{Z}; \code{Z - offset} is treated as
#'   dyadic noise (so it should absorb any additive/multiplicative effects).
#' @param asp \code{NULL} (the default), \code{TRUE} or \code{FALSE}: the prior
#'   of rho. \code{NULL} is the prior of \code{\link{ame}} and \code{\link{lame}},
#'   \eqn{2\,\mathrm{Beta}(2, 2) - 1}, with density proportional to
#'   \eqn{1 - \rho^2}; \code{TRUE} is the arc-sine prior, with density
#'   proportional to \eqn{(1 - \rho^2)^{-1/2}}; \code{FALSE} is the uniform
#'   prior.
#'
#' @details Exactly degenerate residual pairs (every pair exactly equal, or
#'   every pair exactly opposite) can make the correlation conditional improper
#'   on the full interval: with the default prior from four dyads, with the
#'   arc-sine prior from one and with the uniform prior from two. Such
#'   configurations raise an error before a proposal is drawn.
#'
#' @return a new (scalar) value of rho.
#' @keywords internal
#' @author lame authors
#' @export rrho_mh
rrho_mh <- function(Z, rho, s2 = 1, offset = 0, asp = NULL) {
	if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != ncol(Z)) {
		cli::cli_abort("{.arg Z} must be a square numeric matrix.")
	}
	if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
			length(s2) != 1L || !is.finite(s2) || s2 <= 0) {
		cli::cli_abort("{.arg rho} must be strictly between -1 and 1 and {.arg s2} must be finite and positive.")
	}
	if (!is.null(asp) && (!is.logical(asp) || length(asp) != 1L || is.na(asp))) {
		cli::cli_abort("{.arg asp} must be NULL or a single logical value.")
	}
	# exponent of (1 - rho^2) in the prior: package default, arc-sine or uniform
	power <- if (is.null(asp)) .lame_rho_prior_power() else if (asp) -0.5 else 0
	raw_resid <- Z - offset
	resid <- raw_resid / sqrt(s2)
	ut <- upper.tri(resid)
	dir_out <- resid[ut]
	dir_in <- t(resid)[ut]
	if (any(!is.finite(dir_out)) || any(!is.finite(dir_in))) {
		cli::cli_abort("{.arg Z} contains non-finite off-diagonal residuals.")
	}
	raw_out <- raw_resid[ut]
	raw_in <- t(raw_resid)[ut]
	zero_plus <- all(raw_out == -raw_in)
	zero_minus <- all(raw_out == raw_in)
	if (any((dir_out == 0 & raw_out != 0) | (dir_in == 0 & raw_in != 0))) {
		cli::cli_abort("Dyadic residuals cannot be represented at this variance.")
	}
	n_dyad <- length(dir_out)
	# a zero boundary contrast leaves a nonintegrable power at that endpoint
	if (n_dyad >= .lame_rho_improper_dyads(power)) {
		if (zero_plus || zero_minus) {
			cli::cli_abort("The dyadic correlation conditional is improper for these exactly degenerate residual pairs.")
		}
	}
	# norms of the orthogonal dyad contrasts, scaled to avoid underflow
	stable_norm <- function(x) {
		scale <- if (length(x)) max(abs(x)) else 0
		if (scale == 0) return(0)
		scale * sqrt(sum((x / scale)^2))
	}
	plus_norm <- stable_norm(dir_out / 2 + dir_in / 2)
	minus_norm <- stable_norm(dir_out / 2 - dir_in / 2)
	if (!is.finite(plus_norm) || !is.finite(minus_norm)) {
		cli::cli_abort("Dyadic residual sufficient statistics are not finite.")
	}
	if ((plus_norm == 0 && !zero_plus) ||
			(minus_norm == 0 && !zero_minus)) {
		cli::cli_abort("Dyadic residual contrasts cannot be represented at this scale.")
	}
	log_target <- function(r) {
		log_det <- log1p(r) + log1p(-r)
		-0.5 * n_dyad * log_det - (plus_norm / sqrt(1 + r))^2 -
			(minus_norm / sqrt(1 - r))^2 + power * log_det
	}
	old_target <- log_target(rho)
	if (!is.finite(old_target)) {
		cli::cli_abort("The current dyadic correlation has a non-finite log density.")
	}
	# proposal scale from the residuals alone
	prop_sd <- 1 / sqrt(max(1L, n_dyad))
	if (n_dyad >= 2L && sd(dir_out) > 0 && sd(dir_in) > 0) {
		scale <- 2 * (1 - cor(dir_out, dir_in)^2) / sqrt(n_dyad)
		if (is.finite(scale) && scale > 0) prop_sd <- scale
	}
	lo <- pnorm((-1 - rho) / prop_sd)
	hi <- pnorm((1 - rho) / prop_sd)
	rho_new <- rho + prop_sd * qnorm(runif(1, lo, hi))
	if (!is.finite(rho_new) || abs(rho_new) >= 1) return(rho)
	log_qnorm <- function(r) {
		log(pnorm((1 - r) / prop_sd) - pnorm((-1 - r) / prop_sd))
	}
	log_accept <- log_target(rho_new) - old_target + log_qnorm(rho) - log_qnorm(rho_new)
	if (log_accept > log(runif(1))) rho <- rho_new
	rho
}
