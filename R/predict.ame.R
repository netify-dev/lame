#' Predict method for AME models
#'
#' @description
#' Generate predictions from fitted AME models, including point estimates
#' and predictive distributions.
#'
#' @param object Fitted AME model object
#' @param newdata Optional dyadic covariates for prediction: a single
#'   3-D array (\code{n x n x p_dyad}) whose slices match the dyadic
#'   covariates the model was fit with (same order; do not include the
#'   intercept or nodal-covariate slices). Nodal covariates cannot be
#'   changed at prediction time -- their fitted contribution is reused.
#'   When omitted, the fitted design is used. Explicit ordinal response predictions
#'   average numeric category expectations over matched predictor and cutpoint
#'   draws. Newdata predictions require saved draws for enabled additive and
#'   latent effects; disabled effects and zero latent rank require no such draws.
#' @param type Character; type of prediction:
#'   - "response": predicted values on response scale (default)
#'   - "link": predicted values on link scale
#'   - "distribution": full posterior predictive distribution
#' @param n_samples For type="distribution", number of posterior samples
#' @param include_uncertainty Logical; include joint parameter uncertainty (default TRUE).
#'   Requires aligned additive and latent draws saved through
#'   \code{posterior_options(save_ab = TRUE, save_UV = TRUE)}. With FALSE,
#'   only observation noise is simulated around the fitted mean predictor.
#' @param ... Additional arguments (not used)
#'
#' @return Depending on type:
#'   - "response"/"link": Matrix of predictions
#'   - "distribution": Array of posterior predictive samples
#'
#' @details For MCMC Poisson models the response is the mean count of a cell,
#'   \eqn{\exp(\eta + \sigma^2 / 2)}, at the posterior mean \eqn{\eta} of its
#'   linear predictor and the posterior mean \eqn{\sigma^2} of the residual
#'   variance (the model is \eqn{Y \sim \mathrm{Poisson}(\exp(z))},
#'   \eqn{z \sim N(\eta, \sigma^2)}). The stored \code{YPM} holds it, and
#'   \code{newdata} predictions use the same form at the new linear predictor.
#'   \code{type = "distribution"} draws counts from the posterior predictive
#'   distribution instead.
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @examples
#' \donttest{
#' # Fit model
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Point predictions
#' Y_pred <- predict(fit)
#'
#' # Predictions on link scale
#' Y_link <- predict(fit, type = "link")
#' }
#'
#' @export
predict.ame <- function(
	object,
	newdata = NULL,
	type = c("response", "link", "distribution"),
	n_samples = 100,
	include_uncertainty = TRUE,
	...
) {
	.lame_require_draws(object, "predict")

	type <- match.arg(type)

	# network dimensions
	n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
	m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))

	# newdata must carry one slice per dyadic coefficient
	if (!is.null(newdata)) {
		if (!is.array(newdata) || length(dim(newdata)) != 3) {
			cli::cli_abort("newdata must be a 3-dimensional array")
		}
		if (dim(newdata)[1] != n || dim(newdata)[2] != m) {
			cli::cli_abort("newdata dimensions don't match model ({.val {n}} x {.val {m}})")
		}
		# expected number of dyadic covariate slices
		beta_names <- colnames(object$BETA)
		has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
		# dyadic covariates are the beta columns ending in "_dyad" or ".dyad"
		n_dyad_expected <- if (is.null(beta_names)) ncol(object$BETA) - has_intercept
		else sum(grepl("_dyad$|\\.dyad$", beta_names))
		if (n_dyad_expected > 0L && dim(newdata)[3] != n_dyad_expected) {
			cli::cli_abort(c(
				"newdata has {.val {dim(newdata)[3]}} covariate slice{?s} but the model expects {.val {n_dyad_expected}}.",
				"i" = "Supply the same dyadic covariates the model was fit with."
			))
		}
		# realign named newdata to the fit's sorted actor order, which APM, BPM and UVPM follow
		actor_order <- .fit_actor_order(object)
		X <- .reorder_newdata_actors(newdata, actor_order$rows, actor_order$cols)
	} else {
		X <- NULL
		actor_order <- .fit_actor_order(object)
	}

	if (type == "distribution") {
		return(predict_distribution(object, X, n_samples, include_uncertainty))
	} else {
		if (!is.null(X)) {
			# linear predictor: dyad slices from newdata, intercept and nodal slices held at their
			# fitted values (.xbeta_newdata(), shared with predict.lame() and the forecaster)
			beta_mean <- if (length(dim(object$BETA)) == 3L)
				apply(object$BETA, 2, mean) else colMeans(object$BETA)
			beta_names <- if (length(dim(object$BETA)) == 3L)
				dimnames(object$BETA)[[2]] else colnames(object$BETA)
			names(beta_mean) <- beta_names
			XB <- .xbeta_newdata(beta_mean, X, object[["X"]], n, m)

			# additive effects; symmetric fits store one effect (BPM is NULL) that enters both margins
			if (!is.null(object$APM)) {
				for (i in 1:n) XB[i, ] <- XB[i, ] + object$APM[i]
			}
			if (!is.null(object$BPM)) {
				for (j in 1:m) XB[, j] <- XB[, j] + object$BPM[j]
			} else if (isTRUE(object$symmetric) && !is.null(object$APM)) {
				for (j in 1:m) XB[, j] <- XB[, j] + object$APM[j]
			}

			# multiplicative effects
			UVPM <- reconstruct_UVPM(object)
			if (!is.null(UVPM)) {
				XB <- XB + UVPM
			}

			EZ <- XB
			# label with the fit's sorted actor order, matching fitted() and YPM
			row_nm <- rownames(object$YPM) %||% actor_order$rows
			col_nm <- colnames(object$YPM) %||% actor_order$cols
			if (!is.null(row_nm) && length(row_nm) == nrow(EZ))
				rownames(EZ) <- row_nm
			if (!is.null(col_nm) && length(col_nm) == ncol(EZ))
				colnames(EZ) <- col_nm
			# inverse link for type = "response", else EZ on the link scale
			if (type == "response") {
				out <- if (identical(object$family, "ordinal") &&
						identical(object$ordinal_cutpoints, "explicit")) {
					design <- .build_full_design(X, beta_names, object[["X"]], n, m)
					.ordinal_joint_response_mean(object, design)
				} else .fit_to_response(object, EZ)
				dimnames(out) <- dimnames(EZ)
				return(out)
			}
			return(EZ)
		} else {
			if (type == "link") {
				EZ <- reconstruct_EZ(object)
				return(EZ)
			} else {
				# YPM is the stored response-scale posterior summary
				return(object$YPM)
			}
		}
	}
}

#' Generate posterior predictive distribution
#' @noRd
predict_distribution <- function(object, X, n_samples, include_uncertainty) {
	if (length(n_samples) != 1L || !is.finite(n_samples) || n_samples < 1L || n_samples != as.integer(n_samples))
		cli::cli_abort("n_samples must be a positive integer.")
	if (include_uncertainty) {
		state <- .joint_predictive_state(object)
		design <- if (is.null(X)) NULL else .build_full_design(
			X,
			colnames(object$BETA),
			object[["X"]],
			nrow(object$YPM),
			ncol(object$YPM)
		)
		out <- .simulate_joint(object, state, n_samples, designs = if (is.null(design)) NULL else list(design))
		y <- simplify2array(out$Y)
		attr(y, "draw_index") <- out$draw_index
		return(y)
	}
	eta <- predict(object, newdata = X, type = "link")
	fixed <- object
	fixed$BETA <- matrix(colMeans(object$BETA), nrow = 1L, dimnames = list(NULL, colnames(object$BETA)))
	fixed$VC <- matrix(colMeans(object$VC), nrow = 1L, dimnames = list(NULL, colnames(object$VC)))
	if (!is.null(object$ALPHA)) fixed$ALPHA <- matrix(colMeans(object$ALPHA), nrow = 1L)
	simplify2array(lapply(seq_len(n_samples), function(i) .predictive_observation(eta, fixed, 1L)$Y))
}

#' Transform linear predictor to response scale
#' @noRd
transform_to_response <- function(Z, family, s2 = 0, exposure = 1) {
	switch(
		family,
		normal = Z,
		binary = pnorm(Z),
		# mean count of a Poisson-lognormal cell: Z is the unexposed log rate, s2 the
		# variance of its lognormal layer, exposure the offset
		poisson = if (exposure == 0) Z * 0 else exp(Z + s2 / 2 + log(exposure)),
		ordinal = Z, # link scale: no cutpoints here
		cbin = pnorm(Z),
		frn = Z,
		Z # default passthrough
	)
}

#' Extract fitted values from AME model
#'
#' Returns the stored response-scale summary (YPM).
#' For binary models, these are predicted probabilities between 0 and 1.
#' For normal models, these are predicted continuous values.
#' For MCMC Poisson models, these are mean counts,
#' \eqn{\exp(\eta + \sigma^2 / 2)}, at the posterior mean \eqn{\eta} of the
#' linear predictor and the posterior mean \eqn{\sigma^2} of the residual
#' variance, the form the fast estimators report. Predictive quantiles and
#' event probabilities come from \code{predict(type = "distribution")}.
#'
#' @param object Fitted AME model object (class "ame").
#' @param ... Additional arguments (not used).
#'
#' @return An n x n matrix (unipartite) or nA x nB matrix (bipartite) of
#'   fitted values on the response scale. Diagonal entries are \code{NA} for
#'   unipartite networks.
#'
#' @seealso \code{\link{predict.ame}} for predictions with type control,
#'   \code{\link{residuals.ame}} for residuals
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
fitted.ame <- function(object, ...) {
	.lame_require_draws(object, "fitted")
	return(object$YPM)
}

#' Extract residuals from AME model
#'
#' Computes residuals as the difference between observed values and fitted
#' values. For \code{type = "response"}, returns \code{Y - fitted(object)}.
#' For \code{type = "pearson"}, returns response residuals scaled by the
#' standard deviation the family assigns to a response with that fitted value
#' \eqn{\mu} (\code{fitted(object)}, which for a Poisson fit already carries
#' the exposure):
#' \tabular{ll}{
#'   \strong{family} \tab \strong{variance} \cr
#'   normal \tab \eqn{\sigma^2_e}, the posterior mean residual variance \cr
#'   binary, cbin \tab \eqn{\mu (1 - \mu)} (a cbin cell is a nomination with
#'     that marginal probability) \cr
#'   poisson \tab \eqn{\mu + \mu^2 (e^{\sigma^2_e} - 1)}, the variance of a
#'     count when \eqn{y \sim \mathrm{Poisson}(e^z)} and
#'     \eqn{z \sim N(\eta, \sigma^2_e)} \cr
#'   ordinal \tab the variance of the category distribution at the
#'     posterior-mean linear predictor; with explicit cutpoints the
#'     posterior-mean cutpoints, with data-induced cutpoints those that
#'     reproduce the observed category shares (the rule the posterior
#'     predictive draws use) \cr
#'   frn \tab not defined: a rank within its row has a variance that depends on
#'     the other cells of the row, and the call is an error of class
#'     \code{lame_pearson_residuals}
#' }
#' A cell whose variance is zero has no standardised residual and is
#' \code{NA}.
#'
#' @param object Fitted AME model object (class "ame").
#' @param type Character; \code{"response"} (default) for raw residuals or
#'   \code{"pearson"} for standardized residuals.
#' @param ... Additional arguments (not used).
#'
#' @return An n x n matrix (unipartite) or nA x nB matrix (bipartite) of
#'   residuals. Entries where the original data was \code{NA} remain \code{NA}.
#'
#' @seealso \code{\link{fitted.ame}}, \code{\link{predict.ame}}
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
residuals.ame <- function(object, type = c("response", "pearson"), ...) {
	.lame_require_draws(object, "residuals")
	type <- match.arg(type)

	# exact lookup: `$Y` would partial-match `$YPM` on fits without Y (amen fits) and zero every residual
	Y <- object[["Y", exact = TRUE]] %||% object[["Y_obs", exact = TRUE]]
	if (is.null(Y)) {
		cli::cli_abort(c(
			"This fit does not store the observed {.var Y}, so residuals cannot be computed.",
			"i" = "Refit with {.fn lame::ame}, which retains {.var Y} on the fit object."
		))
	}

	if (type == "response") {
		return(Y - object$YPM)
	} else {
		# Pearson residuals scaled by the family's variance at the fitted value
		eta <- if (identical(object[["family", exact = TRUE]], "ordinal")) reconstruct_EZ(object)
		return(.pearson_residual(object, Y, object$YPM, eta))
	}
}

#' Predict method for LAME models
#'
#' @description
#' Generate predictions from a fitted longitudinal AME model. Returns a list of
#' matrices (one per time point) on the requested scale.
#'
#' @param object Fitted LAME model object.
#' @param newdata Optional dyadic covariate arrays. With \code{h = 0}, a list of
#'   \code{T} arrays (\code{[n_row, n_col, p]} each, same actors as the fit) for
#'   counterfactual in-sample predictions; \code{NULL} returns the training-data
#'   predictions. With \code{h > 0}, a list of \code{h} arrays for the future
#'   periods; \code{NULL} reuses the last fitted period's covariates, which
#'   scores the forecast against stale covariates whenever they vary over time.
#' @param type Character; \code{"response"} (default) or \code{"link"}.
#' @param ... Additional arguments (not used).
#'
#' @return With \code{interval = "none"} and \code{by_draw = FALSE}, a list of
#'   prediction matrices (one per period, or per horizon when \code{h > 0}).
#'   With \code{by_draw = TRUE}, an \code{[n_row, n_col, h, n_draws]} array.
#'   With \code{interval = "credible"} or \code{"prediction"}, a list of
#'   \code{h} lists holding \code{$lower}, \code{$median} and \code{$upper}
#'   matrices.
#'
#' @details For \code{h > 0} only \eqn{U_t V_t'} is identified, so each
#'   posterior draw of the coordinates carries its own rotation and sign and
#'   the product of the posterior-mean coordinates falls well short of the
#'   posterior mean of the product (\code{UVPM}, or \code{ULUPM} on symmetric
#'   fits). Under the element-wise AR(1),
#'   \eqn{E[u_{T+h} v_{T+h}'] = \rho_{uv}^{2h} E[u_T v_T']} off the diagonal,
#'   so the forecast carries \code{UVPM} at period \eqn{T} forward by
#'   \eqn{\rho_{uv}^{2h}} and each draw's propagated coordinates enter only as
#'   the deviation around that mean. Static positions are carried forward
#'   unchanged. Dynamic additive effects are propagated as
#'   their actor level plus an AR(1) deviation.
#'
#'   For MCMC Poisson models the response is the mean count of a cell in a
#'   period, \eqn{e_t \exp(\eta + \sigma^2 / 2)}, with \eqn{e_t} the exposure
#'   of the period (\code{period_exposure}; one without), \eqn{\eta} the
#'   posterior mean of the linear predictor and \eqn{\sigma^2} the posterior
#'   mean of the residual variance. The stored \code{YPM} holds it in sample,
#'   \code{newdata} predictions use the same form at the new linear predictor,
#'   and a forecast point (\code{h > 0}) uses the means of the forecast linear
#'   predictor and of the variance over the forecast draws, and the exposure
#'   \code{newexposure}. With \code{by_draw = TRUE} the draws' own conditional
#'   means, \eqn{e \exp(\eta_s + \sigma_s^2 / 2)}, are returned and
#'   \code{interval = "credible"} summarizes them; with
#'   \code{interval = "prediction"} the summaries are of simulated counts.
#'
#' @param h Integer >= 0: forecast horizon. When \code{h = 0} (default),
#'   returns in-sample predictions as before. When \code{h > 0}, propagates
#'   the AR(1) (or RW1) state-space model forward by \code{h} periods and
#'   returns a list of \code{h} matrices (one per future period). Requires
#'   at least one dynamic component on the fit. Warns when posterior
#'   \eqn{\rho_\beta} is near 1.
#'   Coefficient innovations retain the fitted design-based scale. Dynamic
#'   coefficient forecasts require equally spaced AR(1) or RW1 periods;
#'   irregular-period, RW2, Matern and actor-specific coefficient laws are
#'   refused until their future transition is implemented.
#' @param by_draw When \code{TRUE} and \code{h > 0}, returns an
#'   \code{n x n x h x n_draws} array of per-draw forecasts instead of
#'   per-period means.
#' @param interval One of \code{"none"} (default), \code{"credible"} or
#'   \code{"prediction"}. When \code{h > 0}, the per-period output is a list
#'   of length-3 lists with \code{$lower}, \code{$median}, \code{$upper}
#'   matrices computed at the \code{probs} quantiles across posterior draws.
#'   \code{"credible"} is the interval of the conditional mean
#'   \eqn{E[Y_{T+h} \mid \theta_{T+h}]} given each propagated parameter state:
#'   it carries the posterior uncertainty in that state
#'   and adds no observation noise, so it is narrower than the spread of
#'   realised outcomes and under-covers them (on simulated 30-actor normal
#'   panels its nominal 95 percent interval covered 51 to 61 percent of
#'   held-out cells). \code{"prediction"} adds the outcome noise on top of
#'   that mean: a draw from \eqn{N(\eta, ve)} per posterior draw and cell for
#'   \code{"normal"}, a count from \eqn{Poisson(e^{z})} with
#'   \eqn{z \sim N(\eta, ve)} for \code{"poisson"}. For explicit ordinal models,
#'   credible intervals summarize numeric category expectations using matched
#'   saved cutpoints; prediction intervals summarize category draws including
#'   unit Gaussian observation noise. For \code{"binary"}, \code{"cbin"},
#'   \code{"frn"} and data-induced ordinal models, this option retains the
#'   existing credible interval convention and reports that choice. Ignored for in-sample
#'   (\code{h = 0}) predictions.
#'   Outcome forecasts preserve symmetric and within-dyad observation dependence.
#'   Their unipartite diagonals are missing; bipartite diagonal cells are retained.
#'   They are refused for \code{dynamic_rho = TRUE}, which supplies no law
#'   for future reciprocity; conditional-mean forecasts remain available.
#' @param probs Length-2 vector of lower / upper quantiles for the
#'   interval when \code{interval} is \code{"credible"} or
#'   \code{"prediction"}. Default \code{c(0.025, 0.975)}.
#' @param newexposure Optional length-\code{h} non-negative numeric
#'   vector of future-period exposures (Poisson only). When omitted
#'   and the fit has \code{period_exposure} stored, defaults to the
#'   last observed exposure; when both are absent, defaults to 1.
#' @param n_draws Number of posterior draws to use when \code{h > 0}.
#'   Default \code{NULL} uses all stored draws. Ignored for in-sample
#'   (\code{h = 0}) predictions.
#' @param seed Optional RNG seed for the \code{h > 0} forecast draws,
#'   making forecasts reproducible. Ignored when \code{h = 0}.
#' @method predict lame
#' @export
predict.lame <- function(
	object,
	newdata = NULL,
	type = c("response", "link"),
	h = 0L,
	by_draw = FALSE,
	interval = c("none", "credible", "prediction"),
	probs = c(0.025, 0.975),
	newexposure = NULL,
	n_draws = NULL,
	seed = NULL,
	...
) {
	.lame_require_draws(object, "predict")
	type <- match.arg(type)
	interval <- match.arg(interval)
	family <- object$family %||% "normal"

	# h > 0 forecasts the AR(1)/rw1 state forward through .forecast_lame(); negative h is refused
	if (!is.null(h) && length(h) == 1L && is.finite(h) && h < 0L) {
		cli::cli_abort("{.arg h} must be at least 0; got {.val {h}}.")
	}
	explicit_ordinal <- identical(family, "ordinal") && identical(object$ordinal_cutpoints, "explicit")
	want_interval <- interval %in% c("credible", "prediction")
	if (!is.null(h) && h > 0L) {
		# "credible" is the interval of the conditional mean; "prediction" draws the outcome,
		# which only families with response-scale noise need
		predictive <- identical(interval, "prediction")
		if (predictive && !family %in% c("normal", "nrm", "poisson") && !explicit_ordinal) {
			cli::cli_inform(c(
				"i" = "{.code interval = \"prediction\"} on family {.val {family}} reports the interval on the tie probability.",
				"i" = "A draw of the outcome itself is 0 or 1 (or a rank), so its own prediction interval carries no information; this is the {.code interval = \"credible\"} interval of the probability."
			))
			predictive <- FALSE
		}
		if (predictive && identical(type, "link") &&
				!family %in% c("normal", "nrm")) {
			cli::cli_abort(c(
				"{.code interval = \"prediction\"} needs {.code type = \"response\"} on family {.val {family}}.",
				"i" = "The outcome noise is on the response scale; the link scale has no outcome to draw."
			))
		}
		# draws are kept (by_draw) so the interval quantiles can be taken
		raw <- .forecast_lame(
			object,
			h = as.integer(h),
			newdata = newdata,
			newexposure = newexposure,
			type = if (type == "response") "response" else "link",
			by_draw = isTRUE(by_draw) || want_interval,
			predictive = predictive,
			n_draws = n_draws,
			seed = seed
		)
		if (!want_interval) return(raw)
		# raw is [n_a, n_b, h, n_draws]; summarise lower/median/upper per step
		if (length(dim(raw)) != 4L) {
			cli::cli_abort("Internal: expected a 4-D forecast array for interval summarisation.")
		}
		if (length(probs) != 2L || any(!is.finite(probs)) ||
				any(probs <= 0) || any(probs >= 1) || probs[1L] >= probs[2L]) {
			cli::cli_abort("{.arg probs} must be a strictly-ordered length-2 vector in (0, 1).")
		}
		h_int <- as.integer(h)
		return(lapply(seq_len(h_int), function(k) {
			slab <- raw[, , k, , drop = FALSE]
			# drop the singleton third dim, keeping the actor dimnames as in the point path
			slab <- array(
				slab,
				c(dim(slab)[1L], dim(slab)[2L], dim(slab)[4L]),
				dimnames = list(dimnames(raw)[[1L]], dimnames(raw)[[2L]], NULL)
			)
			list(
				lower = apply(
					slab,
					c(1L, 2L),
					stats::quantile,
					probs = probs[1L],
					na.rm = TRUE
				),
				median = apply(
					slab,
					c(1L, 2L),
					stats::quantile,
					probs = 0.5,
					na.rm = TRUE
				),
				upper = apply(
					slab,
					c(1L, 2L),
					stats::quantile,
					probs = probs[2L],
					na.rm = TRUE
				)
			)
		}))
	} else if (want_interval) {
		cli::cli_warn(c(
			"{.code interval = \"{interval}\"} only has an effect for {.code h > 0} forecasts.",
			"i" = "Ignoring for in-sample predictions."
		))
	}

	# training-data predictions
	if (is.null(newdata)) {
		if (type == "response") {
			# return the stored response-scale summary, as in fitted.lame
			if (!is.null(object$YPM)) return(object$YPM)
		}
		EZ <- object$EZ
		if (is.null(EZ)) {
			cli::cli_warn("EZ not available in model object")
			return(object$YPM)
		}
		if (type == "link") return(EZ)
		out <- lapply(seq_along(EZ), function(t) .fit_to_response(object, EZ[[t]], t))
		names(out) <- names(EZ)
		return(out)
	}

	# newdata counterfactual predictions
	actor_slopes <- !is.null(object$dynamic_beta_per_actor)
	if (actor_slopes && (!identical(object$per_actor_model_version, 2L) ||
			!is.matrix(object$theta_actor_mean) || is.null(object$per_actor_design_idx))) {
		cli::cli_abort(c(
			"Counterfactual predictions require corrected actor-specific coefficient summaries.",
			"i" = "Refit with the current package and {.code keep_per_actor = \"summary\"} or {.code \"draws\"}."
		))
	}
	# rebuild EZ per period from the new dyadic covariates, other terms held fixed
	if (!is.list(newdata)) {
		cli::cli_abort("{.arg newdata} must be a list of dyadic covariate arrays, one per time slice.")
	}
	n_t <- length(object$EZ)
	if (length(newdata) != n_t) {
		cli::cli_abort(c(
			"{.arg newdata} has {length(newdata)} slice{?s} but the fit has {.val {n_t}}.",
			"i" = "Supply one dyadic covariate array per time period."
		))
	}
	# posterior-mean coefficients per period: p x T for dynamic_beta, the static mean repeated otherwise
	if (length(dim(object$BETA)) == 3L) {
		beta_per_t <- apply(object$BETA, c(2, 3), mean)
		beta_names <- dimnames(object$BETA)[[2]]
	} else {
		beta_per_t <- matrix(
			colMeans(object$BETA),
			ncol(object$BETA),
			length(object$EZ %||% list(1))
		)
		beta_names <- colnames(object$BETA)
	}
	is_dyad <- if (is.null(beta_names)) {
		c(FALSE, rep(TRUE, nrow(beta_per_t) - 1L))[seq_len(nrow(beta_per_t))]
	} else {
		grepl("_dyad$|\\.dyad$", beta_names)
	}
	# per-period multiplicative term from the stored posterior-mean product, the identified
	# quantity that fit$EZ, fitted() and the forecaster use
	uv_pm <- object[["UVPM", exact = TRUE]] %||% object[["ULUPM", exact = TRUE]]
	UVPM <- if (is.null(uv_pm)) {
		NULL
	} else if (length(dim(uv_pm)) == 3L) {
		lapply(seq_len(dim(uv_pm)[3L]), function(k) uv_pm[, , k])
	} else {
		rep(list(uv_pm), n_t)
	}
	# fitted per-period design (intercept, nodal and dyad slices); NULL for a dyad-only fit,
	# where .build_full_design() adds the intercept slice
	fitted_design_list <- object[["Xlist", exact = TRUE]] %||% object[["X", exact = TRUE]]
	# realign each newdata slice by name to the fit's sorted actor order
	actor_order <- .fit_actor_order(object)
	ez_new <- vector("list", n_t)
	for (t in seq_len(n_t)) {
		Xn <- newdata[[t]]
		if (length(dim(Xn)) == 2L) Xn <- array(Xn, c(dim(Xn), 1L))
		Xn <- .reorder_newdata_actors(Xn, actor_order$rows, actor_order$cols)
		Xn <- .lame_tempered_prediction_design(object, Xn)
		n_row <- dim(Xn)[1L]; n_col <- dim(Xn)[2L]
		fd_t <- if (is.list(fitted_design_list)) fitted_design_list[[t]] else fitted_design_list
		# a slice may cover a subset of the actors (changing composition): effects, product and
		# nodal design are selected by name; an unnamed slice of another size cannot be placed
		sel <- .newdata_actor_select(Xn, actor_order, n_row, n_col, t)
		if (!is.null(fd_t) && length(dim(fd_t)) == 3L &&
				(dim(fd_t)[1L] != n_row || dim(fd_t)[2L] != n_col)) {
			fd_t <- fd_t[sel$rows, sel$cols, , drop = FALSE]
		}
		# full design: dyad slices from newdata, intercept and nodal slices at their fitted values
		Xfull_t <- .build_full_design(Xn, beta_names, fd_t, n_row, n_col)
		if (explicit_ordinal && type == "response") {
			xb <- .ordinal_joint_response_mean(object, Xfull_t, t, sel$rows, sel$cols)
			dimnames(xb) <- list(
				sel$row_names %||% actor_order$rows[sel$rows],
				sel$col_names %||% actor_order$cols[sel$cols]
			)
			ez_new[[t]] <- xb
			next
		}
		xb <- matrix(0, n_row, n_col)
		for (k in seq_len(dim(Xfull_t)[3L])) xb <- xb + beta_per_t[k, t] * Xfull_t[, , k]
		if (actor_slopes) {
			margin <- if (identical(object$dynamic_beta_per_actor, "row")) 1L else 2L
			actors <- if (margin == 1L) sel$rows else sel$cols
			x_actor <- matrix(Xfull_t[, , object$per_actor_design_idx], n_row, n_col)
			xb <- xb + sweep(x_actor, margin, object$theta_actor_mean[actors, t], "*")
		}
		# additive effects, per period when dynamic_ab
		a_t <- if (!is.null(object$a_dynamic)) object$a_dynamic[, t] else object$APM
		b_t <- if (!is.null(object$b_dynamic)) object$b_dynamic[, t] else object$BPM
		# symmetric fits store no BPM: the receiver margin is the sender margin
		if (is.null(b_t) && isTRUE(object$symmetric)) b_t <- a_t
		if (!is.null(a_t) && length(a_t) != n_row) a_t <- a_t[sel$rows]
		if (!is.null(b_t) && length(b_t) != n_col) b_t <- b_t[sel$cols]
		if (!is.null(a_t)) for (i in seq_along(a_t)) xb[i, ] <- xb[i, ] + a_t[i]
		if (!is.null(b_t)) for (j in seq_along(b_t)) xb[, j] <- xb[, j] + b_t[j]
		if (!is.null(UVPM) && length(UVPM) >= t) {
			uv_t <- UVPM[[t]]
			if (!identical(dim(uv_t), c(n_row, n_col))) uv_t <- uv_t[sel$rows, sel$cols, drop = FALSE]
			xb <- xb + uv_t
		}
		# the stored predictor of a symmetric fit is symmetrised
		if (isTRUE(object$symmetric) && nrow(xb) == ncol(xb)) {
			xb <- (xb + t(xb)) / 2
		}
		# label with the sorted actor order of the prediction
		if (!is.null(actor_order$rows) && length(actor_order$rows) == n_row)
			rownames(xb) <- actor_order$rows
		if (!is.null(actor_order$cols) && length(actor_order$cols) == n_col)
			colnames(xb) <- actor_order$cols
		if (!is.null(sel$row_names)) rownames(xb) <- sel$row_names
		if (!is.null(sel$col_names)) colnames(xb) <- sel$col_names
		ez_new[[t]] <- xb
	}
	if (type == "link" || explicit_ordinal) return(ez_new)
	lapply(seq_along(ez_new), function(t) .fit_to_response(object, ez_new[[t]], t))
}

# positions of a slice's actors in the fit's actor order: all of them, a named subset
# (changing composition), or an error for an unnamed slice of another size
.newdata_actor_select <- function(Xn, actor_order, n_row, n_col, t) {
	dn <- dimnames(Xn)
	rn <- if (length(dn) >= 1L) dn[[1L]] else NULL
	cn <- if (length(dn) >= 2L) dn[[2L]] else NULL
	pick <- function(nm, order, n, what) {
		if (is.null(order) || length(order) == n) return(list(idx = seq_len(n), names = NULL))
		if (!is.null(nm) && all(nm %in% order)) return(list(idx = match(nm, order), names = nm))
		cli::cli_abort(c(
			"{.arg newdata} period {t} has {n} {what} actor{?s} but the fit has {length(order)}.",
			"i" = "Name the {what}s of every {.arg newdata} slice with the fitted actor names so a subset of actors can be placed."
		))
	}
	r <- pick(rn, actor_order$rows, n_row, "row")
	c <- pick(cn, actor_order$cols, n_col, "column")
	list(rows = r$idx, cols = c$idx, row_names = r$names, col_names = c$names)
}

#' Extract fitted values from LAME model
#'
#' Returns the stored response-scale summaries. For MCMC Poisson models these
#' are mean counts, \eqn{e_t \exp(\eta + \sigma^2 / 2)}, at the posterior mean
#' \eqn{\eta} of the linear predictor and the posterior mean \eqn{\sigma^2} of
#' the residual variance, with \eqn{e_t} the exposure of the period
#' (\code{period_exposure}; one without); the form the fast estimators report.
#' Predictive quantiles and event probabilities come from \code{simulate()} and
#' from \code{predict(fit, h = 1, interval = "prediction")}.
#'
#' @param object Fitted LAME model
#' @param ... Additional arguments
#'
#' @return List of fitted value matrices (one per time point)
#'
#' @method fitted lame
#' @export
fitted.lame <- function(object, ...) {
	.lame_require_draws(object, "fitted")
	return(object$YPM)
}

#' Extract residuals from LAME model
#'
#' Response residuals are \code{Y - fitted(object)} in each period. Pearson
#' residuals scale them by the standard deviation the family assigns to a
#' response with that fitted value, as in \code{\link{residuals.ame}}: the
#' residual variance for the normal family, \eqn{\mu (1 - \mu)} for binary and
#' cbin, \eqn{\mu + \mu^2 (e^{\sigma^2_e} - 1)} for Poisson (the fitted mean
#' \eqn{\mu} carries the exposure of the period), and the variance of the
#' category distribution at the posterior-mean linear predictor for ordinal
#' fits. The frn family has no cell-wise variance and is an error of class
#' \code{lame_pearson_residuals}. A cell whose variance is zero is \code{NA}.
#'
#' @param object Fitted LAME model
#' @param type Type of residuals ("response" or "pearson")
#' @param ... Additional arguments
#'
#' @return List of residual matrices (one per time point)
#'
#' @method residuals lame
#' @export
residuals.lame <- function(object, type = c("response", "pearson"), ...) {
	.lame_require_draws(object, "residuals")
	type <- match.arg(type)

	# exact lookup: `$Y` would partial-match `$YPM` on fits without Y
	Y <- object[["Y", exact = TRUE]]
	if (is.null(Y)) {
		cli::cli_abort(c(
			"This fit does not store the observed {.var Y}, so residuals cannot be computed.",
			"i" = "Refit with {.fn lame::lame}, which retains {.var Y} on the fit object."
		))
	}

	# 3-D array to a list of matrices
	if (is.array(Y) && length(dim(Y)) == 3) {
		n_time <- dim(Y)[3]
		Y <- lapply(seq_len(n_time), function(t) Y[, , t])
	}

	# the ordinal variance is that of the category distribution at the
	# posterior-mean linear predictor of the period
	eta <- if (type == "pearson" && identical(object[["family", exact = TRUE]], "ordinal"))
		reconstruct_EZ(object)
	if (is.null(eta)) eta <- vector("list", length(object$YPM))

	# lame stores Y and YPM as lists
	mapply(function(y, ypm, eta_t) {
		# under changing composition Y covers every actor and YPM only the period's own:
		# take Y's cells for those actors
		if (!identical(dim(y), dim(ypm)) && !is.null(rownames(y)) &&
				!is.null(rownames(ypm)) && !is.null(colnames(ypm)) &&
				all(rownames(ypm) %in% rownames(y)) && all(colnames(ypm) %in% colnames(y))) {
			y <- y[rownames(ypm), colnames(ypm), drop = FALSE]
		}
		if (type == "response") {
			y - ypm
		} else {
			if (!is.null(eta_t) && !identical(dim(eta_t), dim(ypm)) && !is.null(rownames(eta_t)) &&
					!is.null(rownames(ypm)) && !is.null(colnames(eta_t)) && !is.null(colnames(ypm))) {
				eta_t <- eta_t[rownames(ypm), colnames(ypm), drop = FALSE]
			}
			.pearson_residual(object, y, ypm, eta_t)
		}
	}, Y, object$YPM, eta, SIMPLIFY = FALSE)
}
