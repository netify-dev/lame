# h-step-ahead forecasting. propagates (beta_t, a_t, b_t, u_t, v_t)
# forward via the ar(1) / rw1 recursion, sampling hyperparameters
# (rho_beta, sigma_beta, ...) from their posterior per draw. emits a
# warning when the posterior 95% interval for rho_beta touches 0.99.

#' @noRd
.propagate_ar1_vec <- function(state_T, rho, sigma, h) {
	# state_T: length-p vector of current beta at period T (one draw)
	# rho, sigma: scalar AR(1) parameters
	# h: number of periods to forecast
	# returns: h x p matrix of future states
	p <- length(state_T)
	out <- matrix(NA_real_, nrow = h, ncol = p)
	cur <- as.numeric(state_T)
	for (k in seq_len(h)) {
		cur <- rho * cur + stats::rnorm(p, mean = 0, sd = sigma)
		out[k, ] <- cur
	}
	out
}

#' @noRd
.propagate_per_coef <- function(state_T, rho_by_coef, sigma_by_coef, h) {
	# state_T:        length-p vector
	# rho_by_coef:    length-p vector of AR(1) rho values
	# sigma_by_coef:  length-p vector of AR(1) sigma values
	# returns: h x p matrix
	p <- length(state_T)
	out <- matrix(NA_real_, nrow = h, ncol = p)
	cur <- as.numeric(state_T)
	for (k in seq_len(h)) {
		eps <- stats::rnorm(p, mean = 0, sd = sigma_by_coef)
		cur <- rho_by_coef * cur + eps
		out[k, ] <- cur
	}
	out
}

#' Forecast h periods ahead from a fitted lame object
#'
#' Propagates the AR(1) / RW1 state-space model forward by drawing from
#' the joint posterior of \eqn{(\beta_t, a_t, b_t, U_t, V_t)} and the
#' implied AR(1) hyperparameters. Returned as a list of length \code{h},
#' each element an \code{n x n} (or \code{nA x nB}) matrix on the
#' response or link scale, averaged over posterior draws (or returned as
#' an array of all draws when \code{by_draw = TRUE}).
#'
#' Forecast variance grows linearly in \code{h} when \eqn{\rho_\beta} is
#' close to 1. A warning fires if the posterior 95\% interval for
#' \eqn{\rho_\beta} contains 0.99.
#'
#' @param object A fitted \code{lame} object with at least one dynamic
#'   component active.
#' @param h Number of periods to forecast (must be a positive integer).
#' @param newdata Optional list of \code{h} dyadic covariate arrays for the
#'   future periods. If \code{NULL}, the covariates at the last fitted
#'   period are reused (which assumes covariates are constant or that the
#'   user just wants a "baseline" forecast).
#' @param type One of \code{"link"} or \code{"response"}. \code{"response"}
#'   applies the family inverse-link.
#' @param n_draws Number of posterior draws to use. Default uses all stored
#'   draws.
#' @param by_draw If \code{TRUE}, return a 4-D array \code{[n, n, h, n_draws]}
#'   instead of summarising to per-period means.
#' @param seed Optional RNG seed.
#'
#' @return List of \code{h} matrices (when \code{by_draw = FALSE}), or a
#'   4-D array (when \code{by_draw = TRUE}).
#'
#' @noRd
.forecast_lame <- function(object,
                           h           = 1L,
                           newdata     = NULL,
                           newexposure = NULL,
                           type        = c("link", "response"),
                           n_draws     = NULL,
                           by_draw     = FALSE,
                           seed        = NULL) {
	type <- match.arg(type)
	h <- as.integer(h)
	if (!is.null(seed)) set.seed(seed)
	if (h < 1L) cli::cli_abort("`h` must be at least 1.")

	# at least one dynamic flag must be active for forecasting to make sense
	has_dyn <- isTRUE(object$dynamic_beta) || isTRUE(object$dynamic_ab) ||
	           isTRUE(object$dynamic_uv)
	if (!has_dyn) {
		cli::cli_warn(c(
			"Forecasting from a fit with NO dynamic components produces a constant path.",
			"i" = "Refit with one of {.code dynamic_beta}, {.code dynamic_ab}, {.code dynamic_uv} for genuine forecasts."))
	}

	# dimensions
	n_a <- if (!is.null(object$nA)) object$nA else nrow(object$YPM[[1]] %||% object$YPM)
	n_b <- if (!is.null(object$nB)) object$nB else ncol(object$YPM[[1]] %||% object$YPM)
	bip <- !is.null(object$mode) && identical(object$mode, "bipartite")

	# pick draws
	n_iter <- dim(object$BETA)[1L]
	if (is.null(n_draws)) n_draws <- n_iter
	n_draws <- min(n_draws, n_iter)
	idx <- sample.int(n_iter, n_draws, replace = (n_draws > n_iter))

	# warn once per fit when rho_beta posterior is near the unit root
	if (!is.null(object$RHO_BETA)) {
		q <- stats::quantile(object$RHO_BETA, probs = c(0.025, 0.975), na.rm = TRUE)
		if (q[2L] >= 0.99) {
			warned <- getOption("lame.forecast_unit_root_warned", character())
			fit_key <- digest::digest(object$RHO_BETA)
			if (!fit_key %in% warned) {
				kind <- object$dynamic_beta_kind %||% "ar1"
				var_msg <- switch(kind,
					rw1 = "Under RW1 the {.arg h}-step forecast variance grows linearly in {.arg h}; treat horizons beyond {.val 3} cautiously.",
					rw2 = "Under RW2 the {.arg h}-step forecast propagation uses the AR(1) recursion at {.code rho = 1} so the variance grows linearly in {.arg h}; treat horizons beyond {.val 3} cautiously.",
					matern32 = "Under Matern 3/2 the {.arg h}-step forecast propagation falls back to the AR(1) recursion at {.code rho = 1} so the variance grows linearly in {.arg h}; treat horizons beyond {.val 3} cautiously.",
					"Under AR(1) the {.arg h}-step forecast variance saturates at the stationary level {.code sigma_beta^2 / (1 - rho_beta^2)}, which itself explodes as {.code rho_beta -> 1}; horizons beyond {.val 3} become uninformative in this regime."
				)
				msgs <- c("!" = "Posterior 95% interval for {.code rho_beta} reaches {.val {round(q[2], 3)}}.",
				          "i" = var_msg)
				if (identical(kind, "ar1")) {
					msgs <- c(msgs,
					          "i" = "Consider refitting with {.code dynamic_beta_kind = \"rw1\"} if the underlying process is genuinely unit-root.")
				}
				cli::cli_warn(msgs)
				options(lame.forecast_unit_root_warned = c(warned, fit_key))
			}
		}
	}

	# get last-period state for each draw
	beta_is_dyn <- length(dim(object$BETA)) == 3L
	# returns: an n x n x h x n_draws array
	out <- array(NA_real_, c(n_a, n_b, h, n_draws))

	# pre-cache the last-period covariates as the default extrapolation
	last_X <- if (!is.null(object$Xlist) && length(object$Xlist) > 0L)
		object$Xlist[[length(object$Xlist)]] else NULL
	# user-supplied future covariates take precedence
	X_future <- if (!is.null(newdata)) newdata else replicate(h, last_X, simplify = FALSE)
	if (length(X_future) < h) {
		# cycle
		X_future <- X_future[((seq_len(h) - 1L) %% length(X_future)) + 1L]
	}
	# Align user-supplied newdata to the model's internal design. The fitted
	# Xlist (and therefore beta_path / BETA) carries slices in the canonical
	# layout [intercept?, *_row, *_col, *_dyad]; a user passing `newdata`
	# supplies only the *substantive* dyadic covariates. We rebuild the full
	# design per horizon via the shared .build_full_design() helper -- dyad
	# slices come from newdata (by name when labelled), the intercept + nodal
	# (row/col) slices are held at the last observed period's fitted values --
	# so the positional eta loop below lines up coefficient-for-slice. This is
	# the same name-safe convention predict.ame() / predict.lame() use; the
	# previous code only prepended an intercept and so dropped the nodal
	# fixed effects when a forecast model also had Xrow / Xcol.
	if (!is.null(newdata)) {
		beta_names <- if (length(dim(object$BETA)) == 3L)
			dimnames(object$BETA)[[2L]] else colnames(object$BETA)
		# the held-fixed source for intercept + nodal slices is the last
		# observed period's design
		fitted_design_last <- last_X
		X_future <- lapply(X_future, function(Xk) {
			if (is.null(Xk)) return(Xk)
			.build_full_design(Xk, beta_names, fitted_design_last, n_a, n_b)
		})
	}

	# rho_beta / sigma_beta per draw (with the right shape). Use posterior
	# mean if iteration-level matrices aren't available.
	get_rho_sigma_beta <- function(i) {
		if (!is.null(object$RHO_BETA) && !is.null(object$SIGMA_BETA)) {
			rho_b <- object$RHO_BETA[i, , drop = TRUE]
			sig_b <- object$SIGMA_BETA[i, , drop = TRUE]
		} else if (!is.null(object$rho_beta) && !is.null(object$sigma_beta)) {
			rho_b <- object$rho_beta
			sig_b <- object$sigma_beta
		} else {
			rho_b <- 0.8
			sig_b <- 0.25
		}
		list(rho = rho_b, sigma = sig_b)
	}

	# rho_ab / sigma_ab
	get_rho_sigma_ab <- function(i) {
		rho_a <- if (!is.null(object$rho_ab)) {
			if (is.numeric(object$rho_ab) && length(object$rho_ab) >= i) object$rho_ab[i] else mean(object$rho_ab)
		} else 0
		sig_a <- if (!is.null(object$sigma_ab)) object$sigma_ab else 0.1
		list(rho = rho_a, sigma = sig_a)
	}

	# rho_uv / sigma_uv
	get_rho_sigma_uv <- function(i) {
		rho_u <- if (!is.null(object$rho_uv)) {
			if (is.numeric(object$rho_uv) && length(object$rho_uv) >= i) object$rho_uv[i] else mean(object$rho_uv)
		} else 0
		sig_u <- if (!is.null(object$sigma_uv)) object$sigma_uv else 0.1
		list(rho = rho_u, sigma = sig_u)
	}

	# iterate over draws
	for (s_idx in seq_along(idx)) {
		i <- idx[s_idx]

		# beta path forward
		if (beta_is_dyn) {
			# beta_T: last-period beta for this draw (length p)
			beta_T <- object$BETA[i, , dim(object$BETA)[3L]]
			rs <- get_rho_sigma_beta(i)
			# expand block-level (rho, sigma) to per-coef
			groups <- object$beta_dynamic_groups
			if (is.null(groups)) groups <- rep("dyad", length(beta_T))
			# fast-path: scalar rho/sigma
			if (is.null(names(rs$rho))) {
				rho_by_coef <- rep(rs$rho, length(beta_T))
				sig_by_coef <- rep(rs$sigma, length(beta_T))
			} else {
				# named per-block; broadcast to per-coef via group label
				rho_by_coef <- ifelse(groups %in% names(rs$rho),
				                      rs$rho[groups],
				                      0)  # static coefs: rho effectively 0; we'll handle below
				sig_by_coef <- ifelse(groups %in% names(rs$sigma),
				                      rs$sigma[groups],
				                      0)
			}
			# mask static coefs: they don't evolve. Set rho = 1, sigma = 0.
			static_mask <- if (!is.null(object$beta_dynamic_mask))
				!object$beta_dynamic_mask else rep(FALSE, length(beta_T))
			rho_by_coef[static_mask] <- 1
			sig_by_coef[static_mask] <- 0

			beta_path <- .propagate_per_coef(beta_T, rho_by_coef, sig_by_coef, h)
		} else {
			# static beta: replicate the single beta vector across h periods
			beta_T <- object$BETA[i, ]
			beta_path <- matrix(rep(beta_T, h), nrow = h, byrow = TRUE)
		}

		# a, b forward
		if (isTRUE(object$dynamic_ab) && !is.null(object$a_dynamic) && !is.null(object$b_dynamic)) {
			a_T <- object$a_dynamic[, dim(object$a_dynamic)[2L]]
			b_T <- object$b_dynamic[, dim(object$b_dynamic)[2L]]
			rs_ab <- get_rho_sigma_ab(i)
			a_path <- .propagate_ar1_vec(a_T, rs_ab$rho, rs_ab$sigma, h)
			b_path <- .propagate_ar1_vec(b_T, rs_ab$rho, rs_ab$sigma, h)
		} else {
			a_T <- object$APM
			b_T <- object$BPM
			a_path <- matrix(rep(a_T, h), nrow = h, byrow = TRUE)
			b_path <- matrix(rep(b_T, h), nrow = h, byrow = TRUE)
		}

		# U, V forward. Skip UV propagation entirely when R = 0 (no latent
		# factors): the matrices have zero columns OR the bipartite fit stores
		# placeholder cubes with second dim 0 / 1 that the predict path cannot
		# usefully propagate.
		R_fit <- max(c(object$R %||% 0L,
		               object$R_row %||% 0L, object$R_col %||% 0L))
		has_uv <- !is.null(object$U) && !is.null(object$V) &&
		          (length(dim(object$U)) >= 2L) &&
		          (NCOL(object$U) > 0L) && (NCOL(object$V) > 0L) &&
		          (R_fit > 0L)
		if (has_uv && isTRUE(object$dynamic_uv) && length(dim(object$U)) == 3L) {
			U_T <- object$U[, , dim(object$U)[3L]]
			V_T <- object$V[, , dim(object$V)[3L]]
			rs_uv <- get_rho_sigma_uv(i)
			# propagate every element independently (the AR(1) is element-wise)
			U_path <- array(NA_real_, c(nrow(U_T), ncol(U_T), h))
			V_path <- array(NA_real_, c(nrow(V_T), ncol(V_T), h))
			cur_U <- U_T; cur_V <- V_T
			for (k in seq_len(h)) {
				cur_U <- rs_uv$rho * cur_U + matrix(stats::rnorm(prod(dim(U_T)), 0, rs_uv$sigma),
				                                    nrow(U_T), ncol(U_T))
				cur_V <- rs_uv$rho * cur_V + matrix(stats::rnorm(prod(dim(V_T)), 0, rs_uv$sigma),
				                                    nrow(V_T), ncol(V_T))
				U_path[, , k] <- cur_U
				V_path[, , k] <- cur_V
			}
		} else if (has_uv) {
			U_T <- if (length(dim(object$U)) == 3L) object$U[, , dim(object$U)[3L]] else object$U
			V_T <- if (length(dim(object$V)) == 3L) object$V[, , dim(object$V)[3L]] else object$V
			U_path <- replicate(h, U_T)
			V_path <- replicate(h, V_T)
		} else {
			U_path <- NULL
			V_path <- NULL
		}

		# compute eta at each future period
		for (k in seq_len(h)) {
			eta <- matrix(0, n_a, n_b)
			# X beta
			Xk <- X_future[[k]]
			if (!is.null(Xk)) {
				p_k <- if (length(dim(Xk)) == 3L) dim(Xk)[3L] else 1L
				for (j in seq_len(min(p_k, ncol(beta_path)))) {
					if (length(dim(Xk)) == 3L) {
						eta <- eta + beta_path[k, j] * Xk[, , j]
					} else {
						eta <- eta + beta_path[k, j] * Xk
					}
				}
			}
			# a + b
			eta <- eta + outer(a_path[k, ], b_path[k, ], "+")
			# UV
			if (!is.null(U_path) && !is.null(V_path)) {
				Uk <- if (length(dim(U_path)) == 3L) U_path[, , k] else U_path
				Vk <- if (length(dim(V_path)) == 3L) V_path[, , k] else V_path
				if (bip && !is.null(object$G)) {
					eta <- eta + Uk %*% object$G %*% t(Vk)
				} else {
					eta <- eta + Uk %*% t(Vk)
				}
			}
			out[, , k, s_idx] <- eta
		}
	}

	# response-scale transform per family. for poisson, multiply by the
	# future-period exposure when supplied; default to the last observed
	# exposure when the fit has period_exposure but the user passed no
	# newexposure, and to 1 when neither is present.
	family <- object$family %||% "normal"
	if (identical(family, "poisson") && type == "response") {
		# determine future-period exposure
		if (is.null(newexposure)) {
			if (!is.null(object$period_exposure)) {
				newexposure <- rep(utils::tail(object$period_exposure, 1L), h)
			} else {
				newexposure <- rep(1, h)
			}
		}
		if (length(newexposure) != h) {
			cli::cli_abort("{.arg newexposure} must have length {.val {h}}.")
		}
		if (any(!is.finite(newexposure)) || any(newexposure < 0)) {
			cli::cli_abort("{.arg newexposure} must be non-negative and finite.")
		}
		# apply per-period exposure: out[, , k, draw] *= exposure[k] on the
		# response scale (after exp(eta))
		eta_capped <- pmin(out, log(.Machine$double.xmax))
		out <- exp(eta_capped)
		for (k in seq_len(h)) out[, , k, ] <- out[, , k, ] * newexposure[k]
	} else {
		transform_eta <- function(eta) {
			switch(family,
			       normal = eta,
			       nrm = eta,
			       binary = stats::pnorm(eta),
			       bin = stats::pnorm(eta),
			       cbin = stats::pnorm(eta),
			       tobit = pmax(0, eta),
			       poisson = exp(pmin(eta, log(.Machine$double.xmax))),
			       eta)
		}
		if (type == "response") {
			out <- array(transform_eta(as.numeric(out)),
			             dim = dim(out), dimnames = dimnames(out))
		}
	}

	if (isTRUE(by_draw)) return(out)
	# per-period mean across draws
	lapply(seq_len(h), function(k) apply(out[, , k, , drop = FALSE], c(1, 2), mean))
}
