.rUV_dynamic_validate <- function(U, V, ET, rho_uv, sigma_uv, s2,
                                  symmetric = FALSE, extra = list()) {
	dU <- dim(U)
	dV <- dim(V)
	dE <- dim(ET)
	if (length(dU) != 3L || length(dV) != 3L) {
		cli::cli_abort("{.arg U} and {.arg V} must be 3D arrays.")
	}
	if (dU[2L] != dV[2L]) {
		cli::cli_abort(c(
			"{.arg U} and {.arg V} must have the same latent rank.",
			"i" = "Got dim(U)[2] = {dU[2L]} and dim(V)[2] = {dV[2L]}."))
	}
	if (dU[3L] != dV[3L]) {
		cli::cli_abort(c(
			"{.arg U} and {.arg V} must have the same number of time periods.",
			"i" = "Got dim(U)[3] = {dU[3L]} and dim(V)[3] = {dV[3L]}."))
	}
	if (length(dE) != 3L ||
	    !identical(as.integer(dE), as.integer(c(dU[1L], dV[1L], dU[3L])))) {
		cli::cli_abort(c(
			"{.arg ET} must have dimensions nrow(U) by nrow(V) by n_time.",
			"i" = "Expected [{paste(c(dU[1L], dV[1L], dU[3L]), collapse = ', ')}], got [{paste(dE, collapse = ', ')}]."))
	}
	if (isTRUE(symmetric) && dU[1L] != dV[1L]) {
		cli::cli_abort("{.arg symmetric} = TRUE requires U and V to have the same number of actors.")
	}
	scalars <- c(rho_uv = rho_uv, sigma_uv = sigma_uv, s2 = s2, extra)
	for (nm in names(scalars)) {
		x <- scalars[[nm]]
		if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
			cli::cli_abort("{.arg {nm}} must be a finite numeric scalar.")
		}
	}
	if (sigma_uv <= 0 || s2 <= 0) {
		cli::cli_abort("{.arg sigma_uv} and {.arg s2} must be positive.")
	}
	invisible(list(dU = dU, dV = dV, dE = dE))
}

#' Gibbs sampling of dynamic U and V with AR(1) evolution
#'
#' Updates latent factor positions U and V that evolve over time according to
#' an AR(1) process: u_\{i,t\} = rho * u_\{i,t-1\} + epsilon_\{i,t\}
#'
#' @usage rUV_dynamic_fc(U, V, ET, rho_uv, sigma_uv, s2, shrink=TRUE, symmetric=FALSE)
#' @param U 3D array of current U positions (n x R x T)
#' @param V 3D array of current V positions (n x R x T)
#' @param ET 3D array of residuals (n x n x T)
#' @param rho_uv AR(1) autoregressive parameter for latent positions
#' @param sigma_uv Innovation standard deviation for latent positions
#' @param s2 dyadic variance
#' @param shrink whether to apply shrinkage (default TRUE)
#' @param symmetric whether the network is symmetric (default FALSE)
#' @return list with updated U and V arrays
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rUV_dynamic_fc
rUV_dynamic_fc <- function(U, V, ET, rho_uv, sigma_uv, s2, shrink=TRUE, symmetric=FALSE) {
	.rUV_dynamic_validate(U, V, ET, rho_uv, sigma_uv, s2,
	                      symmetric = symmetric)
	rUV_dynamic_fc_cpp(U, V, ET, rho_uv, sigma_uv, s2, shrink, symmetric)
}

#' Gibbs sampling of dynamic U and V with snap-shift dynamics
#'
#' Like \code{rUV_dynamic_fc} but lets each actor's latent position either drift
#' under the AR(1) prior or jump under a diffuse N(0, kappa^2) snap prior, chosen
#' per actor-period by a Gaussian log-marginal model selection.
#'
#' @param U 3D array of current U positions (n x R x T)
#' @param V 3D array of current V positions (n x R x T)
#' @param ET 3D array of residuals (n x n x T)
#' @param rho_uv AR(1) autoregressive parameter for the drift prior
#' @param sigma_uv Innovation standard deviation for the drift prior
#' @param s2 dyadic variance
#' @param kappa diffuse snap-prior standard deviation
#' @param pi_snap prior snap probability
#' @param delta_u_current current sender-side snap indicators. If \code{NULL},
#'   a zero matrix is used.
#' @param delta_v_current current receiver-side snap indicators. If \code{NULL},
#'   a zero matrix is used.
#' @param shrink whether to apply shrinkage (default TRUE)
#' @param symmetric whether the network is symmetric (default FALSE)
#' @return list with updated U, V arrays and delta_u, delta_v snap indicators
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rUV_dynamic_snap_fc
rUV_dynamic_snap_fc <- function(U, V, ET, rho_uv, sigma_uv, s2, kappa, pi_snap,
		delta_u_current = NULL, delta_v_current = NULL,
		shrink = TRUE, symmetric = FALSE) {
	dU <- dim(U)
	dV <- dim(V)
	.rUV_dynamic_validate(U, V, ET, rho_uv, sigma_uv, s2,
	                      symmetric = symmetric,
	                      extra = list(kappa = kappa, pi_snap = pi_snap))
	if (kappa <= 0 || pi_snap <= 0 || pi_snap >= 1) {
		cli::cli_abort("{.arg kappa} must be positive and {.arg pi_snap} must be in (0, 1).")
	}
	if (is.null(delta_u_current)) {
		delta_u_current <- matrix(0, dU[1L], dU[3L])
	}
	if (is.null(delta_v_current)) {
		delta_v_current <- matrix(0, dV[1L], dV[3L])
	}
	if (!identical(dim(delta_u_current), c(dU[1L], dU[3L]))) {
		cli::cli_abort("{.arg delta_u_current} must have dimensions nrow(U) by n_time.")
	}
	if (!identical(dim(delta_v_current), c(dV[1L], dV[3L]))) {
		cli::cli_abort("{.arg delta_v_current} must have dimensions nrow(V) by n_time.")
	}
	delta_u_current[, 1L] <- 0
	delta_v_current[, 1L] <- 0
	rUV_dynamic_snap_fc_cpp(U, V, ET, rho_uv, sigma_uv, s2, kappa, pi_snap,
	                        delta_u_current, delta_v_current, shrink, symmetric)
}

#' Gibbs sampling of dynamic U and V with heavy-tailed (Student-t) innovations
#'
#' Like \code{rUV_dynamic_fc} but the AR(1) innovations are Student-t via a
#' scale-mixture of normals, a continuous heavy-tailed alternative to snap-shift.
#'
#' @param U 3D array of current U positions (n x R x T)
#' @param V 3D array of current V positions (n x R x T)
#' @param ET 3D array of residuals (n x n x T)
#' @param rho_uv AR(1) autoregressive parameter
#' @param sigma_uv Innovation scale
#' @param s2 dyadic variance
#' @param nu Student-t degrees of freedom
#' @param lambda_u current local scales for U (n x T)
#' @param lambda_v current local scales for V (n x T)
#' @param shrink whether to apply shrinkage (default TRUE)
#' @param symmetric whether the network is symmetric (default FALSE)
#' @return list with updated U, V arrays and lambda_u, lambda_v local scales
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rUV_dynamic_t_fc
rUV_dynamic_t_fc <- function(U, V, ET, rho_uv, sigma_uv, s2, nu,
		lambda_u = NULL, lambda_v = NULL, shrink = TRUE, symmetric = FALSE) {
	dU <- dim(U)
	dV <- dim(V)
	.rUV_dynamic_validate(U, V, ET, rho_uv, sigma_uv, s2,
	                      symmetric = symmetric,
	                      extra = list(nu = nu))
	if (nu <= 2) {
		cli::cli_abort("{.arg nu} must be greater than 2.")
	}
	if (is.null(lambda_u)) {
		lambda_u <- matrix(1, dU[1L], dU[3L])
	}
	if (is.null(lambda_v)) {
		lambda_v <- matrix(1, dV[1L], dV[3L])
	}
	if (!identical(dim(lambda_u), c(dU[1L], dU[3L]))) {
		cli::cli_abort("{.arg lambda_u} must have dimensions nrow(U) by n_time.")
	}
	if (!identical(dim(lambda_v), c(dV[1L], dV[3L]))) {
		cli::cli_abort("{.arg lambda_v} must have dimensions nrow(V) by n_time.")
	}
	rUV_dynamic_t_fc_cpp(U, V, ET, rho_uv, sigma_uv, s2, nu, lambda_u, lambda_v, shrink, symmetric)
}
