# h-step-ahead forecasting: propagates (beta_t, a_t, b_t, u_t, v_t) forward by the
# AR(1) / random-walk recursion, with hyperparameters taken from each posterior draw

#' @noRd
.propagate_ar1_vec <- function(state_T, rho, sigma, h, center = FALSE) {
	# h x p matrix of future states from the current state, with scalar AR(1) rho and sigma
	p <- length(state_T)
	out <- matrix(NA_real_, nrow = h, ncol = p)
	cur <- as.numeric(state_T)
	for (k in seq_len(h)) {
		eps <- stats::rnorm(p, mean = 0, sd = sigma)
		if (isTRUE(center)) eps <- eps - mean(eps)
		cur <- rho * cur + eps
		out[k, ] <- cur
	}
	out
}

#' @noRd
.propagate_per_coef <- function(state_T, rho_by_coef, sigma_by_coef, h) {
	# h x p matrix of future states with per-coefficient AR(1) rho and sigma
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
#' Each forecast uses a matched posterior draw of the coefficients,
#' effects, actor levels, positions and bipartite coupling. New fits retain
#' this terminal state automatically; full trajectories need not be saved.
#' Older fits can use aligned saved effect and position draws. A fit lacking
#' either source must be refitted before forecasting. AR(1), snap and
#' heavy-tailed position forecasts use their fitted transition laws.
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
#' @param predictive If \code{TRUE}, draw the outcome itself rather than its
#'   conditional mean: each posterior draw's linear predictor gains that draw's
#'   observation noise (\code{normal}: \eqn{N(0, ve)} per cell; \code{poisson}:
#'   a count from \eqn{Poisson(e^{z})} with \eqn{z \sim N(\eta, ve)}). The
#'   explicit ordinal models draw categories using unit Gaussian observation
#'   noise and the cutpoints of the matched posterior draw. Other families
#'   retain the conditional-mean forecast convention.
#' @param seed Optional RNG seed.
#'
#' @return List of \code{h} matrices (when \code{by_draw = FALSE}), or a
#'   4-D array (when \code{by_draw = TRUE}).
#'
#' @noRd
.forecast_lame <- function(
	object,
	h = 1L,
	newdata = NULL,
	newexposure = NULL,
	type = c("link", "response"),
	n_draws = NULL,
	by_draw = FALSE,
	predictive = FALSE,
	seed = NULL,
	newdesign = NULL
) {
	# a withheld product (.ae_withhold_product) is restored to its variational means
	object <- .ae_product_restore(object)
	type <- match.arg(type)
	h <- as.integer(h)
	if (!is.null(seed)) set.seed(seed)
	if (h < 1L) cli::cli_abort("`h` must be at least 1.")

	if (!is.null(object$dynamic_beta_per_actor) ||
			(isTRUE(object$dynamic_beta) && object$dynamic_beta_kind %in% c("rw2", "matern32"))) {
		cli::cli_abort(
			"Forecasting is not implemented for this coefficient transition model.",
			class = "lame_forecast_model"
		)
	}
	# forecasting needs at least one dynamic component
	has_dyn <- isTRUE(object$dynamic_beta) || isTRUE(object$dynamic_ab) ||
		isTRUE(object$dynamic_uv)
	if (!has_dyn) {
		cli::cli_warn(c(
			"Forecasting from a fit with no dynamic components produces a constant path.",
			"i" = "Refit with one of {.code dynamic_beta}, {.code dynamic_ab}, {.code dynamic_uv} for time-varying forecasts."
		))
	}

	# actor-union sizes (APM/BPM), since with changing composition period 1 can hold fewer actors
	n_a <- if (!is.null(object$nA)) object$nA
	else if (!is.null(object$APM)) length(object$APM)
	else nrow(object$YPM[[1]] %||% object$YPM)
	# a symmetric fit stores one effect per actor and no BPM, so columns mirror rows
	n_b <- if (!is.null(object$nB)) object$nB
	else if (!is.null(object$BPM)) length(object$BPM)
	else if (isTRUE(object$symmetric) && !is.null(object$APM)) length(object$APM)
	else ncol(object$YPM[[1]] %||% object$YPM)
	bip <- !is.null(object$mode) && identical(object$mode, "bipartite")
	row_actors <- names(object$APM) %||% rownames(object$YPM[[1]] %||% object$YPM)
	col_actors <- names(object$BPM) %||%
		(if (isTRUE(object$symmetric)) names(object$APM) else NULL) %||%
		colnames(object$YPM[[1]] %||% object$YPM)
	if (length(row_actors) != n_a) row_actors <- NULL
	if (length(col_actors) != n_b) col_actors <- NULL

	# every parameter block uses the same retained posterior row
	joint <- .lame_forecast_state(object)
	available <- joint$draw_index
	if (is.null(n_draws)) n_draws <- length(available)
	if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
			n_draws < 1 || n_draws != floor(n_draws)) {
		cli::cli_abort("{.arg n_draws} must be a positive integer.")
	}
	n_draws <- min(n_draws, length(available))
	state_idx <- sample.int(length(available), n_draws)
	idx <- available[state_idx]

	# warn once per fit when rho_beta is near the unit root
	if (!is.null(object$RHO_BETA)) {
		q <- stats::quantile(object$RHO_BETA, probs = c(0.025, 0.975), na.rm = TRUE)
		if (q[2L] >= 0.99) {
			warned <- .lame_state$forecast_unit_root_warned %||% character()
			fit_key <- if (requireNamespace("digest", quietly = TRUE)) {
				digest::digest(object$RHO_BETA)
			} else {
				paste(
					length(object$RHO_BETA),
					signif(sum(object$RHO_BETA, na.rm = TRUE), 12),
					sep = "_"
				)
			}
			if (!fit_key %in% warned) {
				kind <- object$dynamic_beta_kind %||% "ar1"
				var_msg <- switch(
					kind,
					rw1 = "Under RW1 the {.arg h}-step forecast variance grows linearly in {.arg h}, so intervals widen at every step.",
					rw2 = "Under RW2 the {.arg h}-step forecast propagation uses the AR(1) recursion at {.code rho = 1} so the variance grows linearly in {.arg h} and intervals widen at every step.",
					matern32 = "Under Matern 3/2 the {.arg h}-step forecast propagation falls back to the AR(1) recursion at {.code rho = 1} so the variance grows linearly in {.arg h} and intervals widen at every step.",
					"Under AR(1) the {.arg h}-step forecast variance saturates at the stationary level {.code sigma_beta^2 / (1 - rho_beta^2)}, which itself explodes as {.code rho_beta -> 1}; horizons beyond {.val 3} become uninformative in this regime."
				)
				msgs <- c(
					"!" = "Posterior 95% interval for {.code rho_beta} reaches {.val {round(q[2], 3)}}.",
					"i" = var_msg
				)
				if (identical(kind, "ar1")) {
					msgs <- c(
						msgs,
						"i" = "Consider refitting with {.code dynamic_beta_kind = \"rw1\"} if the underlying process is genuinely unit-root."
					)
				}
				cli::cli_warn(msgs)
				.lame_state$forecast_unit_root_warned <- c(warned, fit_key)
			}
		}
	}

	# last-period state for each draw
	beta_is_dyn <- length(dim(object$BETA)) == 3L
	beta_transition <- if (beta_is_dyn) .lame_beta_transition_scale(object) else NULL
	# n x n x h x n_draws array, actor-named on the first two dimensions
	out <- array(
		NA_real_,
		c(n_a, n_b, h, n_draws),
		dimnames = list(row_actors, col_actors, NULL, NULL)
	)

	# canonical full designs keep nodal roles, symmetric merging and actor order
	last_X <- if (length(object$Xlist)) utils::tail(object$Xlist, 1L)[[1L]] else NULL
	beta_names <- dimnames(object$BETA)[[2L]]
	X_future <- if (!is.null(newdesign)) newdesign else if (!is.null(newdata)) newdata else
		replicate(h, last_X, simplify = FALSE)
	if (length(X_future) < h) X_future <- X_future[((seq_len(h) - 1L) %% length(X_future)) + 1L]
	X_future <- lapply(X_future, function(Xk) {
		if (is.null(Xk)) return(Xk)
		if (!is.null(newdata) && is.null(newdesign)) {
			Xk <- .lame_forecast_design(Xk, NULL, row_actors, col_actors, n_a, n_b)
			Xk <- .build_full_design(Xk, beta_names, last_X, n_a, n_b)
		}
		Xk <- .lame_forecast_design(Xk, beta_names, row_actors, col_actors, n_a, n_b)
		.lame_tempered_prediction_design(object, Xk)
	})

	# aligned transition rows at the same retained state index
	get_rho_sigma_beta <- function(i) {
		list(
			rho = object$RHO_BETA[i, , drop = TRUE],
			sigma = object$SIGMA_BETA[i, , drop = TRUE]
		)
	}

	get_rho_sigma_ab <- function(i) {
		rho_a <- if (!is.null(object$rho_ab)) {
			if (is.numeric(object$rho_ab) && length(object$rho_ab) >= i) object$rho_ab[i] else mean(object$rho_ab)
		} else 0
		sig_a <- if (!is.null(object$sigma_ab)) {
			if (is.numeric(object$sigma_ab) && length(object$sigma_ab) >= i) object$sigma_ab[i] else mean(object$sigma_ab)
		} else 0.1
		list(rho = rho_a, sigma = sig_a)
	}

	get_rho_sigma_uv <- function(i) {
		rho_u <- if (!is.null(object$rho_uv)) {
			if (is.numeric(object$rho_uv) && length(object$rho_uv) >= i) object$rho_uv[i] else mean(object$rho_uv)
		} else 0
		sig_u <- if (!is.null(object$sigma_uv)) {
			if (is.numeric(object$sigma_uv) && length(object$sigma_uv) >= i) object$sigma_uv[i] else mean(object$sigma_uv)
		} else 0.1
		list(rho = rho_u, sigma = sig_u)
	}

	# multiplicative term, skipped when there are no latent factors (R = 0)
	R_fit <- max(c(
		object[["R", exact = TRUE]] %||% 0L,
		object$R_row %||% 0L,
		object$R_col %||% 0L
	))
	U_fit <- object$U
	V_fit <- object[["V"]]
	# symmetric fits store only U (the term is u l u'), so V mirrors U
	if (is.null(V_fit) && isTRUE(object$symmetric)) V_fit <- U_fit
	has_uv <- !is.null(U_fit) && !is.null(V_fit) &&
		(length(dim(U_fit)) >= 2L) &&
		(NCOL(U_fit) > 0L) && (NCOL(V_fit) > 0L) &&
		(R_fit > 0L)
	uv_is_dynamic <- has_uv && isTRUE(object$dynamic_uv) &&
		length(dim(U_fit)) == 3L
	# V mirrors U for symmetric fits and for dynamic fits whose V is not a cube
	sym_uv <- has_uv && (isTRUE(object$symmetric) ||
			(uv_is_dynamic && length(dim(V_fit)) != 3L))
	# keep the n x R matrix shape when R = 1
	last_uv_slice <- function(A) {
		if (length(dim(A)) == 3L)
			matrix(A[, , dim(A)[3L]], nrow = dim(A)[1L], ncol = dim(A)[2L])
		else A
	}
	U_T <- if (has_uv) last_uv_slice(U_fit) else NULL
	V_T <- if (!has_uv) NULL else if (sym_uv) U_T else last_uv_slice(V_fit)
	# bipartite fits carry the interaction matrix
	G_fc <- if (has_uv && bip) object[["G"]] else NULL
	# a symmetric fit's signed diagonal L plays the role of the interaction matrix
	if (has_uv && !bip && sym_uv && !is.null(object$L) &&
			length(object$L) == ncol(U_T)) {
		G_fc <- diag(as.numeric(object$L), nrow = ncol(U_T))
	}
	uv_product <- function(Um, Vm) {
		if (!is.null(G_fc)) Um %*% G_fc %*% t(Vm) else Um %*% t(Vm)
	}

	for (s_idx in seq_along(idx)) {
		i <- idx[s_idx]
		state_i <- state_idx[s_idx]
		if (has_uv) {
			U_T <- matrix(joint$U[, , state_i], n_a)
			V_T <- matrix(joint$V[, , state_i], n_b)
			G_fc <- if (bip) matrix(joint$G[, , state_i], ncol(U_T)) else NULL
			if (sym_uv) {
				den <- colSums(U_T^2)
				L_i <- ifelse(den > 0, colSums(U_T * V_T) / den, 0)
				G_fc <- diag(L_i, ncol(U_T))
				V_T <- U_T
			}
		}

		# beta path forward
		if (beta_is_dyn) {
			beta_T <- object$BETA[i, , dim(object$BETA)[3L]]
			rs <- get_rho_sigma_beta(i)
			# expand block-level (rho, sigma) to per-coefficient values
			groups <- object$beta_dynamic_groups
			if (is.null(groups)) groups <- rep("dyad", length(beta_T))
			if (is.null(names(rs$rho))) {
				rho_by_coef <- rep(rs$rho, length(beta_T))
				sig_by_coef <- rep(rs$sigma, length(beta_T))
			} else {
				# every dynamic coefficient needs its block's transition draws
				dyn_coef <- if (!is.null(object$beta_dynamic_mask)) {
					object$beta_dynamic_mask
				} else {
					rep(TRUE, length(beta_T))
				}
				has_draws <- groups %in% names(rs$rho) & groups %in% names(rs$sigma)
				missing_groups <- unique(groups[dyn_coef & !has_draws])
				if (length(missing_groups)) {
					cli::cli_abort(
						c(
							"Dynamic coefficient block{?s} {.val {missing_groups}} {?has/have} no transition draws in {.code fit$RHO_BETA} and {.code fit$SIGMA_BETA}.",
							"i" = "The forecast needs the rho and sigma draws of every dynamic coefficient block."
						),
						class = "lame_forecast_block_missing"
					)
				}
				# per-block values broadcast to coefficients by group label
				rho_by_coef <- ifelse(
					groups %in% names(rs$rho),
					rs$rho[groups],
					0
				)
				sig_by_coef <- ifelse(
					groups %in% names(rs$sigma),
					rs$sigma[groups],
					0
				)
			}
			# static coefficients do not evolve: rho = 1, sigma = 0
			static_mask <- if (!is.null(object$beta_dynamic_mask))
				!object$beta_dynamic_mask else rep(FALSE, length(beta_T))
			rho_by_coef[static_mask] <- 1
			sig_by_coef[static_mask] <- 0
			sig_by_coef <- sig_by_coef * beta_transition$sd

			beta_path <- .propagate_per_coef(beta_T, rho_by_coef, sig_by_coef, h)
		} else {
			# static beta repeated across the h periods
			beta_T <- object$BETA[i, ]
			beta_path <- matrix(rep(beta_T, h), nrow = h, byrow = TRUE)
		}

		# keep the posterior dependence between terminal states and parameters
		a_T <- joint$a[, state_i]
		b_T <- joint$b[, state_i]
		if (isTRUE(object$dynamic_ab)) {
			a_lev <- joint$a_level[, state_i]
			b_lev <- joint$b_level[, state_i]
			rs_ab <- get_rho_sigma_ab(i)
			center_ab <- object$intercept %||%
				object$start_vals$checkpoint_state$layout$intercept %||%
				("intercept" %in% dimnames(object$BETA)[[2L]])
			a_path <- if (isFALSE(object$rvar)) matrix(0, h, length(a_T)) else
				.propagate_ar1_vec(a_T - a_lev, rs_ab$rho, rs_ab$sigma, h, center_ab) +
				matrix(a_lev, h, length(a_lev), byrow = TRUE)
			b_path <- if (isTRUE(object$symmetric)) a_path else if (isFALSE(object$cvar)) {
				matrix(0, h, length(b_T))
			} else .propagate_ar1_vec(b_T - b_lev, rs_ab$rho, rs_ab$sigma, h, center_ab) +
				matrix(b_lev, h, length(b_lev), byrow = TRUE)
		} else {
			a_path <- matrix(a_T, h, length(a_T), byrow = TRUE)
			b_path <- matrix(b_T, h, length(b_T), byrow = TRUE)
		}

		if (uv_is_dynamic) {
			rs_uv <- get_rho_sigma_uv(i)
			kind <- object$dynamic_uv_kind %||% "ar1"
			pi_i <- if (!is.null(object$pi_snap)) object$pi_snap[i] else NULL
			U_path <- array(NA_real_, c(dim(U_T), h))
			V_path <- array(NA_real_, c(dim(V_T), h))
			cur_U <- U_T; cur_V <- V_T
			for (k in seq_len(h)) {
				cur_U <- .lame_forecast_uv_step(
					cur_U,
					rs_uv$rho,
					rs_uv$sigma,
					kind,
					pi_i,
					object$prior$snap_kappa,
					object$prior$t_nu %||% 4
				)
				cur_V <- if (sym_uv) cur_U else .lame_forecast_uv_step(
					cur_V,
					rs_uv$rho,
					rs_uv$sigma,
					kind,
					pi_i,
					object$prior$snap_kappa,
					object$prior$t_nu %||% 4
				)
				U_path[, , k] <- cur_U; V_path[, , k] <- cur_V
			}
		} else if (has_uv) {
			U_path <- array(rep(U_T, h), c(dim(U_T), h))
			V_path <- array(rep(V_T, h), c(dim(V_T), h))
		} else {
			U_path <- V_path <- NULL
		}

		# linear predictor at each future period
		for (k in seq_len(h)) {
			eta <- matrix(0, n_a, n_b)
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
			eta <- eta + outer(a_path[k, ], b_path[k, ], "+")
			if (!is.null(U_path) && !is.null(V_path)) {
				Uk <- if (length(dim(U_path)) == 3L)
					matrix(U_path[, , k], dim(U_path)[1L], dim(U_path)[2L])
				else U_path
				Vk <- if (length(dim(V_path)) == 3L)
					matrix(V_path[, , k], dim(V_path)[1L], dim(V_path)[2L])
				else V_path
				prod_k <- uv_product(Uk, Vk)
				eta <- eta + prod_k
			}
			if (isTRUE(object$symmetric) && !bip) eta <- (eta + t(eta)) / 2
			out[, , k, s_idx] <- eta
		}
	}

	# response scale by family; Poisson multiplies by the future exposure, defaulting
	# to the last observed exposure, else 1
	family <- object$family %||% "normal"
	# only families with observation noise on the response scale draw an outcome
	explicit_ordinal <- identical(family, "ordinal") && identical(object$ordinal_cutpoints, "explicit")
	predictive <- isTRUE(predictive) && (family %in% c("normal", "nrm", "poisson") || explicit_ordinal)
	if (predictive && isTRUE(object$dynamic_rho)) {
		cli::cli_abort(
			c(
				"Outcome forecasts are not implemented for period-specific reciprocity without a future reciprocity model.",
				"i" = "Use conditional-mean forecasts, or simulate fitted periods with type = 'posterior'."
			),
			class = "lame_forecast_model"
		)
	}
	if (identical(family, "poisson") && type == "response") {
		newexposure <- newexposure %||% rep(utils::tail(object$period_exposure %||% 1, 1L), h)
		if (length(newexposure) != h)
			cli::cli_abort("{.arg newexposure} must have length {.val {h}}.")
		if (any(!is.finite(newexposure) | newexposure < 0))
			cli::cli_abort("{.arg newexposure} must be non-negative and finite.")
	}
	if (predictive) {
		# observation variance and reciprocity at the selected posterior iteration;
		# symmetric outcomes share noise, bipartite cells do not
		observation_fit <- object
		if (family == "nrm") observation_fit$family <- "normal"
		if (family == "poisson") observation_fit$period_exposure <- newexposure
		for (d in seq_len(n_draws)) for (k in seq_len(h)) {
			eta <- matrix(out[, , k, d], n_a, n_b, dimnames = list(row_actors, col_actors))
			out[, , k, d] <- .predictive_observation(eta, observation_fit, idx[d], k)$Y
		}
	} else if (explicit_ordinal && type == "response") {
		for (d in seq_len(n_draws)) for (k in seq_len(h))
			out[, , k, d] <- .ordinal_response_mean(
				out[, , k, d],
				.ordinal_saved_cutpoints(object, idx[d]),
				object$ordinal_levels
			)
	} else if (identical(family, "poisson") && type == "response") {
		if (!is.matrix(object$VC) || !"ve" %in% colnames(object$VC))
			cli::cli_abort("Poisson forecasts require saved observation variance draws.")
		s2_draw <- object$VC[idx, "ve"]
		if (any(!is.finite(s2_draw) | s2_draw < 0))
			cli::cli_abort("Poisson forecasts require finite nonnegative observation variance draws.")
		# point forecast is the mean count at the mean linear predictor and variance
		# (as .poisson_fitted_mean); by_draw keeps each draw's own mean
		if (!isTRUE(by_draw)) {
			s2_bar <- mean(s2_draw)
			return(lapply(seq_len(h), function(k) {
				eta_bar <- apply(out[, , k, , drop = FALSE], c(1, 2), mean)
				mu <- if (newexposure[k] == 0) eta_bar * 0 else
					exp(eta_bar + s2_bar / 2 + log(newexposure[k]))
				if (any(!is.finite(mu[!is.na(eta_bar)])))
					cli::cli_abort("The Poisson forecast mean overflowed; no clipped mean was substituted.")
				mu
			}))
		}
		for (d in seq_len(n_draws)) for (k in seq_len(h)) {
			eta <- matrix(out[, , k, d], n_a, n_b)
			mu <- if (newexposure[k] == 0) matrix(0, n_a, n_b) else
				exp(eta + s2_draw[d] / 2 + log(newexposure[k]))
			if (any(!is.finite(mu[!is.na(eta)])))
				cli::cli_abort("The Poisson forecast mean overflowed; no clipped mean was substituted.")
			mu[is.na(eta)] <- NA_real_
			out[, , k, d] <- mu
		}
	} else if (type == "response" && family %in% c("binary", "bin", "cbin")) {
		out <- array(stats::pnorm(out), dim = dim(out), dimnames = dimnames(out))
	}

	if (isTRUE(by_draw)) {
		attr(out, "draw_index") <- idx
		return(out)
	}
	lapply(seq_len(h), function(k) apply(out[, , k, , drop = FALSE], c(1, 2), mean))
}
