#' Get fitted object from MCMC results
#'
#' @param start_vals List object that is null or contains
#' starting values
#' @param Y dependent variable in array format
#' @param family character vector (e.g. 'bin', 'nrm') specifying
#' family type
#' @param xP number of exogenous covariates
#' @param rvar logical indicating whether to include sender random
#' effects
#' @param cvar logical indicating whether to include receiver random
#' effects
#' @param R Number of dimensions for multiplicative effects
#' @param odmax vector of maximum ranks for cbin/frn families (optional)
#' @details In addition to the slots listed below, any of the dynamic
#'   hyperparameter scalars \code{rho_uv}, \code{sigma_uv}, \code{rho_ab},
#'   and \code{sigma_ab} present in \code{start_vals} (for example from
#'   \code{\link{als_start_vals}}) are passed through so a dynamic MCMC run
#'   can warm-start them, as are the adapted step sizes
#'   \code{level_move_sd} and \code{dev_move_sd} of the binary scale move, so
#'   that a run continued with \code{\link{lame_resume}} keeps the step it
#'   had tuned rather than re-adapting from the default.
#' @param bip logical; \code{TRUE} for a bipartite (rectangular) panel, in
#'   which case the Poisson start values keep the diagonal of a square matrix
#'   as a real cell
#' @param unit_var the unit of the response, the variance a normal response
#'   is measured in (1 where the link fixes the scale of the latent response):
#'   the starting dyadic variance, the spread of the starting coefficients and
#'   the floor of the starting covariance of the effects are written in it
#' @return List of starting values for MCMC
#' @author Shahryar Minhas
#' @keywords internal
#' @export get_start_vals

get_start_vals <- function(
	start_vals,
	Y,
	family,
	xP,
	rvar,
	cvar,
	R,
	odmax = NULL,
	bip = FALSE,
	unit_var = 1
) {

	N <- dim(Y)[3]

	# partial start_vals are merged over the defaults computed below, so a caller
	# can warm-start just `beta`
	user_start <- start_vals
	if (!is.null(user_start) && !is.list(user_start)) {
		cli::cli_abort("{.arg start_vals} must be a list (or NULL).")
	}
	start_vals <- NULL

	if (is.null(start_vals)) {

		# initialize Z
		Z <- array(dim = dim(Y))
		for (t in 1:N) {
			if (family == "normal") {Z[, , t] <- Y[, , t] }
			if (family == "ordinal") {Z[, , t] <- matrix(zscores(Y[, , t]), nrow(Y[, , t]), ncol(Y[, , t]))}
			if (family == "binary") {
				Z[, , t] <- matrix(zscores(Y[, , t]), nrow(Y[, , t]), nrow(Y[, , t]))
				zy_max <- if (sum(Y[, , t] == 0, na.rm = TRUE) > 0) {
					max(Z[, , t][Y[, , t] == 0], na.rm = TRUE)
				} else {
					min(Z[, , t][Y[, , t] == 1], na.rm = TRUE) - 1e-6
				}
				zy_min <- if (sum(Y[, , t] == 1, na.rm = TRUE) > 0) {
					min(Z[, , t][Y[, , t] == 1], na.rm = TRUE)
				} else {
					max(Z[, , t][Y[, , t] == 0], na.rm = TRUE) + 1e-6
				}
				z01 <- .5 * (zy_max + zy_min)
				Z[, , t] <- Z[, , t] - z01
			}
			if (family == "poisson") {
				# zero cells start at the period's mean rate (see .ame_pois_init_z in ame_unipartite.R)
				Z[, , t] <- .ame_pois_init_z(Y[, , t], unipartite = !bip)
			}

			if (is.element(family, c("cbin", "frn"))) {
				if (is.null(odmax)) {
					cli::cli_abort("odmax must be provided for family=cbin/frn")
				}
				Z[, , t] <- Y[, , t]
				for (i in 1:nrow(Y[, , t])) {
					yi <- Y[i, , t]
					zi <- zscores(yi)
					rnkd <- which(!is.na(yi) & yi > 0)
					if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
						zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 1e-3
					}

					if (length(rnkd) < odmax[i]) {
						urnkd <- which(!is.na(yi) & yi == 0)
						if (length(urnkd) > 0 && max(zi[urnkd]) > 0) { zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 1e-3 }
					}

					Z[i, , t] <- zi
				}
			}
		}

		# fill missing entries by row and column mean imputation
		ZA <- Z
		for (t in 1:N) {
			mu <- mean(Z[, , t], na.rm = TRUE)
			a <- rowMeans(Z[, , t], na.rm = TRUE) ; b <- colMeans(Z[, , t], na.rm = TRUE)
			a[is.na(a)] <- 0 ; b[is.na(b)] <- 0
			ZA[, , t] <- mu + outer(a, b, "+")
		}
		Z[is.na(Z)] <- ZA[is.na(Z)]

		# initialize the other parameters
		if (xP > 0) {
			beta <- rnorm(xP, mean = 0, sd = 0.1 * sqrt(unit_var))
			z_mean <- mean(Z, na.rm = TRUE)
			if (!is.finite(z_mean)) {
				z_mean <- 0
			}
			beta[1] <- z_mean
			# a binary panel starts at the intercept that reproduces its tie rate given
			# the starting effects, so a sparse panel does not start several probit units away
			if (identical(family, "binary")) {
				p_obs <- mean(Y > 0, na.rm = TRUE)
				if (is.finite(p_obs) && p_obs > 0 && p_obs < 1) {
					ab0 <- outer(a, b, "+")
					lev <- tryCatch(
						stats::uniroot(
							function(m) mean(stats::pnorm(m + ab0), na.rm = TRUE) - p_obs,
							c(-20, 20)
						)$root,
						error = function(e) stats::qnorm(p_obs)
					)
					if (is.finite(lev)) beta[1] <- lev
				}
			}
		} else {
			beta <- numeric(0)
		}
		s2 <- unit_var
		rho <- 0
		Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar))
		Sab <- .lame_start_sab_covariance(Sab, c(rvar, cvar), floor = 1e-6 * unit_var)
		U <- V <- matrix(0, nrow(Y[, , 1]), R)

	}

	# merge any user-supplied slots over the defaults computed above
	if (!is.null(user_start)) {
		# exact `[[` lookups: `$b` would partial-match `beta` and `$rho` would match
		# `rho_uv` / `rho_ab`
		if (!is.null(user_start[["Z"]])) Z <- user_start[["Z"]]
		if (!is.null(user_start[["beta"]])) beta <- user_start[["beta"]]
		if (!is.null(user_start[["a"]])) a <- user_start[["a"]]
		if (!is.null(user_start[["b"]])) b <- user_start[["b"]]
		if (!is.null(user_start[["U"]])) U <- user_start[["U"]]
		if (!is.null(user_start[["V"]])) V <- user_start[["V"]]
		if (!is.null(user_start[["rho"]])) rho <- user_start[["rho"]]
		if (!is.null(user_start[["s2"]])) s2 <- user_start[["s2"]]
		if (!is.null(user_start[["Sab"]])) Sab <- user_start[["Sab"]]
	}

	out <- list(
		Z = Z,
		beta = beta,
		a = a,
		b = b,
		U = U,
		V = V,
		rho = rho,
		s2 = s2,
		Sab = Sab
	)

	# dynamic hyperparameter scalars (e.g. from als_start_vals) pass through;
	# lame() reads them when the matching dynamic option is on
	for (nm in c(
		"rho_uv", "sigma_uv", "rho_ab", "sigma_ab",
		"level_move_sd", "dev_move_sd"
	)) {
		val <- user_start[[nm]]
		if (!is.null(val)) {
			if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
				cli::cli_abort("{.arg start_vals${nm}} must be a single finite number.")
			}
			out[[nm]] <- val
		}
	}

	return(out)

}

# signed starting factors for a static symmetric fit: the leading R
# eigenvectors (by absolute eigenvalue) of the period-averaged residual
# without its additive part, scaled as ame()'s SVD start. V = U L carries each
# eigenvalue's sign, which rUV_sym_fc_cpp reads as the starting weight; NULL
# when the residual has no signal
.lame_signed_sym_start <- function(Z, R, unit_sd = 1, canonical_sign = FALSE) {
	Zbar <- if (length(dim(Z)) == 3L) apply(Z, c(1L, 2L), mean, na.rm = TRUE) else Z
	if (!is.matrix(Zbar) || nrow(Zbar) != ncol(Zbar) || nrow(Zbar) < R) return(NULL)
	Zbar[!is.finite(Zbar)] <- NA
	diag(Zbar) <- NA
	E0 <- Zbar - mean(Zbar, na.rm = TRUE)
	r0 <- rowMeans(E0, na.rm = TRUE)
	r0[!is.finite(r0)] <- 0
	E0 <- E0 - outer(r0, r0, "+")
	E0[is.na(E0)] <- 0
	E0 <- (E0 + t(E0)) / 2
	if (!any(E0 != 0)) return(NULL)
	eig <- eigen(E0, symmetric = TRUE)
	k <- order(-abs(eig$values))[seq_len(R)]
	if (!all(is.finite(eig$values[k])) || any(eig$values[k] == 0)) return(NULL)
	# eigenvector signs are arbitrary and can change with the scale of the matrix, so a normal
	# response takes the sign that makes the largest coordinate positive and y and c y start alike
	if (canonical_sign) {
		for (j in k) {
			v <- eig$vectors[, j]
			if (v[which.max(abs(v))] < 0) eig$vectors[, j] <- -v
		}
	}
	# unit-free coordinates; the weights carry the unit of the response (V = U L, |L| = unit_sd)
	U <- eig$vectors[, k, drop = FALSE] %*% diag(sqrt(0.5 * abs(eig$values[k])) / sqrt(unit_sd), nrow = R)
	list(U = U, V = U %*% diag(sign(eig$values[k]) * unit_sd, nrow = R))
}

# start values with the names dropped from the coefficient and effect slots,
# since the samplers would copy a label into every scalar computed from them;
# the positions and the interaction matrix keep their dimnames
.start_vals_plain <- function(start_vals) {
	if (!is.list(start_vals)) return(start_vals)
	for (nm in intersect(c("beta", "a", "b"), names(start_vals))) {
		x <- start_vals[[nm]]
		if (is.numeric(x)) start_vals[[nm]] <- unname(x)
	}
	start_vals
}
