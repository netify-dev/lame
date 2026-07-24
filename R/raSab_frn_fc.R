#' Simulate a and Sab from full conditional distributions under frn likelihood
#'
#' Simulate a and Sab from full conditional distributions under frn likelihood
#'
#'
#' @usage raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,
#' SS=round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
#' @param YL list of ranked individuals, from least to most preferred in each
#' row
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each individual
#' @param odobs observed outdegree
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @keywords internal
#' @author lame authors
#' @export raSab_frn_fc
raSab_frn_fc <-
function(Z, Y, YL, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL, SS=round(sqrt(nrow(Z)))) {

		n <- nrow(Z)
		if(is.null(Sab0)){ Sab0 <- diag(2) }
		if(is.null(eta0)){ eta0 <- 4 }

		## work on the row-effect-free residual R_ij = Z_ij - a_i, which does
		## not depend on a and stays fixed while a is refreshed below
		resid <- Z - outer(a, rep(1, n))
		rows <- seq_len(n)

		## Lower support point for a_i: the least-preferred nominee (YL[,1])
		## must keep a positive latent value, Z_ij = R_ij + a_i >= 0, hence
		## a_i >= -R[i, YL_i]. Rows with no nomination give NA -> unconstrained.
		lo <- -resid[cbind(rows, YL[, 1])]
		lo[is.na(lo)] <- -Inf

		## Upper support point for a_i: every observed non-nominee must stay at
		## a negative latent value, Z_ij = R_ij + a_i <= 0, so
		## a_i <= -max_{j non-nominated, observed} R_ij. Nominated cells and
		## missing cells are dropped from that maximum.
		drop <- is.na(Y) | (Y != 0)
		cand <- resid
		cand[drop] <- -Inf
		hi <- -apply(cand, 1, max, na.rm = TRUE)
		## a saturated out-degree (odobs == odmax) removes the ceiling
		hi[odobs == odmax] <- Inf

		for(iter in seq_len(SS)) {

				## full conditional of a given b under the bivariate N(0, Sab):
				## a | b ~ N( (S_ab / S_bb) b , S_aa - S_ab^2 / S_bb )
				slope <- Sab[1, 2] / Sab[2, 2]
				mu_a <- slope * b
				sd_a <- sqrt(Sab[1, 1] - slope * Sab[1, 2])

				## draw each a_i from that normal truncated to [lo_i, hi_i] via
				## the probability-integral-transform: map a uniform draw onto
				## the CDF interval, then invert
				p_lo <- pnorm((lo - mu_a) / sd_a)
				p_hi <- pnorm((hi - mu_a) / sd_a)
				u <- p_lo + runif(n) * (p_hi - p_lo)
				a <- mu_a + sd_a * qnorm(u)

				## conjugate inverse-Wishart refresh of the 2x2 covariance:
				## Sab ~ inv-Wishart(eta0 + n, eta0*Sab0 + [a b]'[a b]); sample it
				## as the inverse of a Wishart draw on the inverted scale
				post_scale <- eta0 * Sab0 + crossprod(cbind(a, b))
				Sab <- solve(rwish(solve(post_scale), eta0 + n))
		}

		list(Z = resid + outer(a, rep(1, n)), a = a, Sab = Sab)
}
