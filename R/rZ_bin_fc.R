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

		# truncation code per cell: -1 (missing) is unrestricted, 0 forces Z < 0, 1 forces Z > 0
		obs <- Y
		obs[is.na(obs)] <- -1L

		above <- upper.tri(Z)
		below <- lower.tri(Z)

		# one-sided truncated normal by inverse cdf on the log scale, stable when the
		# conditional mean lies far outside the retained half-line
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
				log_ccdf <- pnorm(
					(lo - mu) / cond_sd,
					lower.tail = FALSE,
					log.p = TRUE
				)
				out <- mu + cond_sd * qnorm(
					log(u) + log_ccdf,
					lower.tail = FALSE,
					log.p = TRUE
				)
			}
			# hold the incoming value if an endpoint u gives +/-Inf
			spoil <- !is.finite(out)
			if (any(spoil)) out[spoil] <- keep[spoil]
			out
		}

		# Gibbs sweep over the two triangles in turn, so each cell conditions on the
		# current value of its mirror entry
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

		# dyad-level Metropolis move: a bivariate proposal with unit variances and
		# correlation rho, accepted when both entries match the sign of Y
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

		# diagonal from its unconstrained conditional
		diag(Z) <- rnorm(n, diag(EZ), sqrt(1 + rho))

		Z
	}

#' Sample Z under symmetric binary data (one latent per unordered pair)
#'
#' Symmetric analogue of \code{\link{rZ_bin_fc}}: each unordered actor pair
#' carries one latent value (\code{Z[i, j] = Z[j, i]}) drawn from a normal
#' with mean equal to the average linear predictor (EZ_ij + EZ_ji) / 2 and
#' variance 1, truncated to the half-line the observation requires (above
#' zero for a tie, below zero for a non-tie, unconstrained for a missing
#' pair). The lower triangle is mirrored after the upper-triangle draws so
#' the matrix is exactly symmetric at every sweep. The directed sampler run
#' at a near-unit dyadic correlation instead truncates the two cells of a
#' pair separately, which sharpens the likelihood of a symmetric tie and
#' inflates the fitted scale by a few percent.
#'
#' @param Z current symmetric latent matrix (n x n, Z = t(Z)).
#' @param EZ expected value of Z (regression, additive and multiplicative
#'   contributions).
#' @param Y observed symmetric binary matrix; the diagonal is ignored.
#' @param rho dyadic correlation carried by the fit; it only sets the
#'   variance of the diagonal draw, \code{1 + rho}, the convention of the
#'   coefficient kernel.
#' @return updated symmetric n x n latent matrix; the diagonal is refreshed
#'   with unconstrained normal draws around \code{diag(EZ)}.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_bin_sym_fc <- function(Z, EZ, Y, rho = 0) {
	n <- nrow(Z)
	ut <- upper.tri(Z)
	mu <- (EZ[ut] + t(EZ)[ut]) / 2
	y <- Y[ut]
	y_t <- t(Y)[ut]
	y[is.na(y)] <- y_t[is.na(y)]
	lb <- ifelse(!is.na(y) & y > 0, 0, -Inf)
	ub <- ifelse(!is.na(y) & y <= 0, 0, Inf)
	# tail-stable truncated draw (see rZ_ord_fc)
	z_new <- .lame_checked_truncated_gaussian(mu, 1, lb, ub)
	Z[ut] <- z_new
	Zt <- t(Z)
	Zt[ut] <- z_new
	Z <- t(Zt)
	diag(Z) <- stats::rnorm(n, diag(EZ), sqrt(1 + rho))
	Z
}
