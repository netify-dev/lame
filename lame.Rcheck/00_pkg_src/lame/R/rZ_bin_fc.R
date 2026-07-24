#' Draw latent Z for the binary probit model
#'
#' Gibbs / Metropolis update of the latent normal matrix \code{Z} underlying a
#' binary relational matrix \code{Y}. Off-diagonal entries are drawn from the
#' dyadic full conditional \eqn{Z_{ij} \mid Z_{ji} \sim N(EZ_{ij} + \rho\,
#' (Z_{ji}-EZ_{ji}),\; 1-\rho^2)} truncated to the half-line implied by
#' \eqn{Y_{ij}} (positive when \eqn{Y_{ij}=1}, negative when \eqn{Y_{ij}=0},
#' unrestricted when \eqn{Y_{ij}} is missing). A correlated dyad-level
#' Metropolis proposal is then attempted to improve mixing, and the diagonal is
#' refreshed from its unconstrained conditional.
#'
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix, the new value of Z
#' @author lame authors
#' @keywords internal
#' @export rZ_bin_fc
rZ_bin_fc <-
function(Z, EZ, rho, Y) {

	n <- nrow(Z)
	cond_sd <- sqrt(1 - rho^2)

	# encode the truncation region for each cell by its observation:
	# a missing entry (coded -1) is unrestricted, a 0 forces Z < 0, a 1
	# forces Z > 0. Missing values are folded into the -1 code up front.
	obs <- Y
	obs[is.na(obs)] <- -1L

	above <- upper.tri(Z)
	below <- lower.tri(Z)

	# robust one-sided truncated-normal inverse-cdf sampler. Working on the
	# log probability scale keeps the draw well defined even when the
	# conditional mean sits many sd's outside the retained half-line, where a
	# plain qnorm(runif(pnorm(lo), pnorm(hi))) would underflow to a point mass.
	draw_half <- function(mu, lo, hi, keep) {
		m <- length(mu)
		u <- runif(m)
		if (is.infinite(lo) && is.infinite(hi)) {
			out <- mu + cond_sd * qnorm(u)
		} else if (is.infinite(lo)) {
			# retain (-Inf, hi]
			log_cdf <- pnorm((hi - mu) / cond_sd, log.p = TRUE)
			out <- mu + cond_sd * qnorm(log(u) + log_cdf, log.p = TRUE)
		} else {
			# retain [lo, Inf)
			log_ccdf <- pnorm((lo - mu) / cond_sd,
				lower.tail = FALSE, log.p = TRUE)
			out <- mu + cond_sd * qnorm(log(u) + log_ccdf,
				lower.tail = FALSE, log.p = TRUE)
		}
		# an endpoint u would produce +/-Inf; R's rng never returns it, but
		# guard anyway by holding the incoming value for such cells.
		spoil <- !is.finite(out)
		if (any(spoil)) out[spoil] <- keep[spoil]
		out
	}

	# single dyadic Gibbs sweep. The two triangles are visited in turn for
	# each observation code so that a cell always conditions on the current
	# value of its mirror entry.
	for (code in c(-1L, 0L, 1L)) {
		lo <- if (code == 1L) 0 else -Inf
		hi <- if (code == 0L) 0 else Inf
		for (half in c(1L, 2L)) {
			cells <- (if (half == 1L) above else below) & (obs == code)
			if (!any(cells)) next
			mu <- EZ[cells] + rho * (t(Z)[cells] - t(EZ)[cells])
			Z[cells] <- draw_half(mu, lo, hi, Z[cells])
		}
	}

	# dyad-level Metropolis proposal. c(1,1)*g_diag + t()*g_off applied to iid
	# standard normals yields a symmetric bivariate proposal with unit
	# variances and correlation rho; accept the pair only when both entries
	# land on the sign required by Y (missing dyads impose no sign).
	root_p <- sqrt(1 + rho)
	root_m <- sqrt(1 - rho)
	g_diag <- (root_p + root_m) / 2
	g_off <- (root_p - root_m) / 2
	noise <- matrix(rnorm(n * n), n, n)
	prop <- EZ + g_diag * noise + g_off * t(noise)
	compat <- (obs == -1L) | (sign(prop) == sign(obs - 0.5))
	diag(compat) <- TRUE
	compat <- compat & t(compat)
	Z[compat] <- prop[compat]

	# refresh the (unconstrained) diagonal from its conditional
	diag(Z) <- rnorm(n, diag(EZ), sqrt(1 + rho))

	Z
}
