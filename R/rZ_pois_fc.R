#' Gibbs update for latent variable in a Poisson AME model
#'
#' Updates the latent variable Z in a Poisson AME model using a
#' Metropolis-Hastings step. The model assumes y_\{i,j\} ~ Poisson(exp(z_\{i,j\}))
#' where z_\{i,j\} is the latent variable representing the log mean.
#'
#' @usage rZ_pois_fc(Z, EZ, rho, s2, Y)
#' @param Z a square matrix, the current value of the latent variable
#' @param EZ expected value of Z (regression effects + random effects)
#' @param rho dyadic correlation
#' @param s2 dyadic variance (overdispersion parameter)
#' @param Y square relational matrix of observed counts
#' @return updated value of Z
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats dpois
#' @keywords internal
#' @export rZ_pois_fc
rZ_pois_fc <- function(Z, EZ, rho, s2, Y) {
	suppressWarnings(rZ_pois_fc_cpp(Z, EZ, rho, s2, Y))
}

#' Metropolis update for the latent variable of a symmetric Poisson AME model
#'
#' Symmetric analogue of \code{\link{rZ_pois_fc}}: each unordered pair
#' carries one latent value (\code{Z[i, j] = Z[j, i]}) with prior
#' \eqn{N(\mu_{ij}, s^2)}, where \eqn{\mu_{ij}} is the average of the two
#' directed linear predictors, and one Poisson observation. The upper
#' triangle is updated with an independent random-walk Metropolis step per
#' pair and mirrored to the lower triangle. Under the near-unit dyadic
#' correlation the directed kernel carries on symmetric fits, independent
#' per-cell proposals are almost always rejected because they break the tie
#' between the two cells; proposing the pair as one value removes that
#' penalty.
#'
#' @param Z current symmetric latent matrix (n x n, Z = t(Z)).
#' @param EZ expected value of Z (regression, additive and multiplicative
#'   contributions).
#' @param s2 dyadic variance (overdispersion parameter).
#' @param Y observed symmetric count matrix; the diagonal is ignored.
#' @param log_exposure optional log exposure offset added to the log rate.
#' @param rho dyadic correlation carried by the fit; it only sets the
#'   variance of the diagonal draw, \code{(1 + rho) s2}, the convention of
#'   the coefficient kernel.
#' @return updated symmetric n x n latent matrix; the diagonal is refreshed
#'   with unconstrained normal draws around \code{diag(EZ)}.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_pois_sym_fc <- function(Z, EZ, s2, Y, log_exposure = 0, rho = 0) {
	n <- nrow(Z)
	ut <- upper.tri(Z)
	mu <- (EZ[ut] + t(EZ)[ut]) / 2
	z <- Z[ut]
	y <- Y[ut]
	y_t <- t(Y)[ut]
	y[is.na(y)] <- y_t[is.na(y)]
	obs <- is.finite(y)
	sd_z <- sqrt(s2)
	z_prop <- z + stats::rnorm(length(z), 0, sd_z)
	loglik <- function(zz) {
		out <- numeric(length(zz))
		out[obs] <- .poisson_log_mass(y[obs], zz[obs] + log_exposure)
		out
	}
	lr <- loglik(z_prop) - loglik(z) +
		stats::dnorm(z_prop, mu, sd_z, log = TRUE) -
		stats::dnorm(z, mu, sd_z, log = TRUE)
	accept <- !is.na(lr) & (log(stats::runif(length(z))) < lr)
	z[accept] <- z_prop[accept]
	Z[ut] <- z
	Zt <- t(Z)
	Zt[ut] <- z
	Z <- t(Zt)
	# diagonal: unconstrained draw around EZ, as in the other latent samplers, so residuals
	# never see missing values
	diag(Z) <- stats::rnorm(n, diag(EZ), sqrt((1 + rho) * s2))
	Z
}
