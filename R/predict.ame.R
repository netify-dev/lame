####
#' Predict method for AME models
#' 
#' @description
#' Generate predictions from fitted AME models, including point estimates
#' and predictive distributions.
#' 
#' @param object Fitted AME model object
#' @param newdata Optional new data (covariates) for prediction
#' @param type Character; type of prediction:
#'   - "response": predicted values on response scale (default)
#'   - "link": predicted values on link scale
#'   - "distribution": full posterior predictive distribution
#' @param n_samples For type="distribution", number of posterior samples
#' @param include_uncertainty Logical; include parameter uncertainty (default TRUE)
#' @param ... Additional arguments (not used)
#' 
#' @return Depending on type:
#'   - "response"/"link": Matrix of predictions
#'   - "distribution": Array of posterior predictive samples
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @examples
#' \donttest{
#' # Fit model
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
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
	include_uncertainty = TRUE, ...
	){
	
	type <- match.arg(type)
	
	# network dimensions
	n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
	m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))
	
	# validate new covariates including the covariate-count (3rd dim), which
	# must match the dyadic coefficients of `object$beta` so a wrong-shape
	# `newdata` is rejected rather than producing wrong predictions.
	if(!is.null(newdata)) {
		if(!is.array(newdata) || length(dim(newdata)) != 3) {
			cli::cli_abort("newdata must be a 3-dimensional array")
		}
		if(dim(newdata)[1] != n || dim(newdata)[2] != m) {
			cli::cli_abort("newdata dimensions don't match model ({.val {n}} x {.val {m}})")
		}
		# number of dyadic covariate slices expected
		beta_names <- colnames(object$BETA)
		has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
		# dyadic covariates are the beta columns ending in "_dyad" or ".dyad"
		n_dyad_expected <- if (is.null(beta_names)) ncol(object$BETA) - has_intercept
		                   else sum(grepl("_dyad$|\\.dyad$", beta_names))
		if (n_dyad_expected > 0L && dim(newdata)[3] != n_dyad_expected) {
			cli::cli_abort(c(
				"newdata has {.val {dim(newdata)[3]}} covariate slice{?s} but the model expects {.val {n_dyad_expected}}.",
				"i" = "Supply the same dyadic covariates the model was fit with."))
		}
		X <- newdata
	} else {
		X <- NULL
	}
	
	if(type == "distribution") {
		return(predict_distribution(object, X, n_samples, include_uncertainty))
	} else {
		if(!is.null(X)) {
			# linear predictor from new covariates. `x` (= newdata) carries
			# only the dyadic covariate slices; the intercept and any nodal
			# (xrow / xcol) contributions are held at their fitted values.
			# .xbeta_newdata() rebuilds the full design in the model's
			# canonical slice order (dyad slices from newdata, the rest from
			# fit$x) and contracts it against the named coefficient vector --
			# the single name-safe convention shared with predict.lame() and
			# the forecaster. (the old positional `beta[k+1]` fallback
			# silently multiplied dyadic covariates by nodal coefficients.)
			beta_mean <- if (length(dim(object$BETA)) == 3L)
				apply(object$BETA, 2, mean) else colMeans(object$BETA)
			beta_names <- if (length(dim(object$BETA)) == 3L)
				dimnames(object$BETA)[[2]] else colnames(object$BETA)
			names(beta_mean) <- beta_names
			XB <- .xbeta_newdata(beta_mean, X, object$X, n, m)

			# add sender/receiver random effects
			if(!is.null(object$APM)) {
				for(i in 1:n) XB[i,] <- XB[i,] + object$APM[i]
			}
			if(!is.null(object$BPM)) {
				for(j in 1:m) XB[,j] <- XB[,j] + object$BPM[j]
			}
			
			# multiplicative effects
			UVPM <- reconstruct_UVPM(object)
			if(!is.null(UVPM)) {
				XB <- XB + UVPM
			}

			EZ <- XB
			# apply the inverse link when type = "response"; otherwise
			# return ez on the link scale.
			if (type == "response") {
				return(transform_to_response(EZ, object$family))
			}
			return(EZ)
		} else {
			if(type == "link") {
				EZ <- reconstruct_EZ(object)
				return(EZ)
			} else {
				# ypm is the posterior mean on the response scale, computed
				# during mcmc by averaging simulated outcomes. returning it
				# directly avoids jensen's inequality bias that would arise
				# from transforming the posterior mean of the linear predictor.
				return(object$YPM)
			}
		}
	}
}

####
#' Generate posterior predictive distribution
#' @noRd
predict_distribution <- function(object, X, n_samples, include_uncertainty) {

	n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
	m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))

	beta_names <- colnames(object$BETA)
	has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"

	Y_pred <- array(NA, dim = c(n, m, n_samples))

	if(include_uncertainty && !is.null(object$BETA)) {
		# n_mcmc = number of post-burnin iterations (first dim) for both 2-d and 3-d beta
		n_mcmc <- dim(object$BETA)[1]
		if(n_mcmc < n_samples) {
			cli::cli_warn("Only {.val {n_mcmc}} MCMC samples available, using all")
			n_samples <- n_mcmc
		}
		beta_is_dyn <- length(dim(object$BETA)) == 3L

		sample_idx <- sample(n_mcmc, n_samples, replace = (n_samples > n_mcmc))

		for(i in 1:n_samples) {
			if (beta_is_dyn) {
				# for predict.ame() (cross-sectional), 3-d beta shouldn't exist;
				# defensive: use period 1 of the path
				beta <- object$BETA[sample_idx[i], , 1]
			} else {
				beta <- object$BETA[sample_idx[i], ]
			}
			if(!is.null(X)) {
				# name-safe per-draw newdata reconstruction (same helper as
				# predict.ame): dyad slices from newdata, intercept + nodal
				# held at fitted design values.
				names(beta) <- beta_names
				XB <- .xbeta_newdata(beta, X, object$X, n, m)
			} else {
				if(has_intercept) {
					XB <- matrix(beta[1], n, m)
				} else {
					XB <- matrix(0, n, m)
				}
			}
			
			if(!is.null(object$APM)) {
				XB <- XB + object$APM
			}
			if(!is.null(object$BPM)) {
				XB <- sweep(XB, 2, object$BPM, "+")
			}
			
			if(!is.null(object$U) && !is.null(object$V)) {
				if(!is.null(object$U_samples) && !is.null(object$V_samples)) {
					idx <- min(i, dim(object$U_samples)[3])
					UV <- object$U_samples[,,idx] %*% t(object$V_samples[,,idx])
				} else {
					# no uv samples stored, use posterior means
					UV <- object$U %*% t(object$V)
				}
				XB <- XB + UV
			}
			
			# observation noise
			s2 <- object$VC[sample_idx[i], ncol(object$VC)]
			Z <- XB + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
			
			# apply inverse link
			Y_pred[,,i] <- transform_to_response(Z, object$family)
		}
	} else {
		EZ <- if(!is.null(X)) {
			beta_mean <- if (length(dim(object$BETA)) == 3L)
				apply(object$BETA, 2, mean) else colMeans(object$BETA)
			names(beta_mean) <- beta_names
			# name-safe newdata reconstruction (dyad from newdata; intercept +
			# nodal held at fitted values) -- same helper as everywhere else.
			.xbeta_newdata(beta_mean, X, object$X, n, m)
		} else {
			reconstruct_EZ(object)
		}
		
		s2 <- mean(object$VC[, ncol(object$VC)])
		
		for(i in 1:n_samples) {
			Z <- EZ + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
			Y_pred[,,i] <- transform_to_response(Z, object$family)
		}
	}
	
	return(Y_pred)
}

####
#' Transform linear predictor to response scale
#' @noRd
transform_to_response <- function(Z, family) {
	switch(family,
		normal = Z,
		binary = pnorm(Z),
		poisson = exp(Z),
		tobit = pmax(0, Z),
		ordinal = Z,  # proper ordinal transform needs cutpoints we don't have here
		cbin = pnorm(Z),
		frn = Z,
		rrl = Z,
		Z  # default passthrough
	)
}

####
#' Extract fitted values from AME model
#'
#' Returns the posterior mean of the network on the response scale (YPM).
#' For binary models, these are predicted probabilities between 0 and 1.
#' For normal models, these are predicted continuous values.
#' For Poisson models, these are predicted counts.
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
	return(object$YPM)
}

####
#' Extract residuals from AME model
#'
#' Computes residuals as the difference between observed values and fitted
#' values. For \code{type = "response"}, returns \code{Y - fitted(object)}.
#' For \code{type = "pearson"}, returns response residuals scaled by the
#' standard deviation implied by the family (e.g., \code{sqrt(p*(1-p))} for
#' binary).
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
	type <- match.arg(type)

	if(!is.null(object$Y)) {
		Y <- object$Y
	} else if(!is.null(object$Y_obs)) {
		Y <- object$Y_obs
	} else {
		cli::cli_warn("original Y not stored in model object, residuals may be incorrect")
		Y <- object$YPM
	}

	if(type == "response") {
		return(Y - object$YPM)
	} else {
		# pearson residuals scaled by variance
		s2 <- mean(object$VC[, ncol(object$VC)])
		return((Y - object$YPM) / sqrt(s2))
	}
}

####
#' Predict method for LAME models
#'
#' @description
#' Generate predictions from a fitted longitudinal AME model. Returns a list of
#' matrices (one per time point) on the requested scale.
#'
#' @param object Fitted LAME model object.
#' @param newdata Optional list of \code{T} dyadic covariate arrays
#'   (\code{[n_row, n_col, p]} each, with the same actors as the fit) to
#'   compute counterfactual predictions. When \code{NULL}, the training-data
#'   predictions are returned.
#' @param type Character; \code{"response"} (default) or \code{"link"}.
#' @param ... Additional arguments (not used).
#'
#' @return List of prediction matrices (one per time point).
#'
#' @param h Integer >= 0: forecast horizon. When \code{h = 0} (default),
#'   returns in-sample predictions as before. When \code{h > 0}, propagates
#'   the AR(1) (or RW1) state-space model forward by \code{h} periods and
#'   returns a list of \code{h} matrices (one per future period). Requires
#'   at least one dynamic component on the fit. Warns when posterior
#'   \eqn{\rho_\beta} is near 1.
#' @param by_draw When \code{TRUE} and \code{h > 0}, returns an
#'   \code{n x n x h x n_draws} array of per-draw forecasts instead of
#'   per-period means.
#' @param interval One of \code{"none"} (default) or \code{"credible"}.
#'   When \code{h > 0} and \code{"credible"}, the per-period output is
#'   a list of length-3 lists with \code{$lower}, \code{$median},
#'   \code{$upper} matrices computed at the \code{probs} quantiles
#'   across posterior draws. Ignored for in-sample (\code{h = 0})
#'   predictions.
#' @param probs Length-2 vector of lower / upper quantiles for the
#'   credible interval when \code{interval = "credible"}. Default
#'   \code{c(0.025, 0.975)}.
#' @param newexposure Optional length-\code{h} non-negative numeric
#'   vector of future-period exposures (Poisson only). When omitted
#'   and the fit has \code{period_exposure} stored, defaults to the
#'   last observed exposure; when both are absent, defaults to 1.
#' @method predict lame
#' @export
predict.lame <- function(object, newdata = NULL,
                         type = c("response", "link"),
                         h = 0L, by_draw = FALSE,
                         interval = c("none", "credible"),
                         probs = c(0.025, 0.975),
                         newexposure = NULL, ...) {
	type <- match.arg(type)
	interval <- match.arg(interval)
	family <- object$family %||% "normal"

	# forecast horizon. h > 0 delegates to .forecast_lame() which
	# propagates the ar(1)/rw1 state-space model forward. abort on
	# negative h so a typo does not silently return in-sample predictions.
	if (!is.null(h) && length(h) == 1L && is.finite(h) && h < 0L) {
		cli::cli_abort("{.arg h} must be at least 0; got {.val {h}}.")
	}
	if (!is.null(h) && h > 0L) {
		# force by_draw = true internally when an interval is requested so
		# quantiles can be summarised across draws.
		want_interval <- identical(interval, "credible")
		raw <- .forecast_lame(object,
		                      h           = as.integer(h),
		                      newdata     = newdata,
		                      newexposure = newexposure,
		                      type        = if (type == "response") "response" else "link",
		                      by_draw     = isTRUE(by_draw) || want_interval)
		if (!want_interval) return(raw)
		# raw is [n_a, n_b, h, n_draws]; collapse to lower/median/upper per k
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
			# collapse the singleton 3rd dim
			slab <- array(slab, c(dim(slab)[1L], dim(slab)[2L], dim(slab)[4L]))
			list(
				lower  = apply(slab, c(1L, 2L), stats::quantile,
				               probs = probs[1L], na.rm = TRUE),
				median = apply(slab, c(1L, 2L), stats::quantile,
				               probs = 0.5, na.rm = TRUE),
				upper  = apply(slab, c(1L, 2L), stats::quantile,
				               probs = probs[2L], na.rm = TRUE)
			)
		}))
	} else if (identical(interval, "credible")) {
		cli::cli_warn(c(
			"{.code interval = \"credible\"} only has an effect for {.code h > 0} forecasts.",
			"i" = "Ignoring for in-sample predictions."))
	}

	# ---- training-data predictions ----
	if (is.null(newdata)) {
		if (type == "response") {
			# match fitted.lame: ypm is the posterior mean on the response
			# scale, avoiding jensen bias from transforming the posterior
			# mean of the linear predictor
			if (!is.null(object$YPM)) return(object$YPM)
		}
		EZ <- object$EZ
		if (is.null(EZ)) {
			cli::cli_warn("EZ not available in model object")
			return(object$YPM)
		}
		if (type == "link") return(EZ)
		return(lapply(EZ, function(ez) transform_to_response(ez, family)))
	}

	# ---- newdata counterfactual predictions ----
	# rebuild ez per slice using new dyadic covariates while holding intercept,
	# additive effects and the multiplicative term fixed
	if (!is.list(newdata)) {
		cli::cli_abort("{.arg newdata} must be a list of dyadic covariate arrays, one per time slice.")
	}
	n_t <- length(object$EZ)
	if (length(newdata) != n_t) {
		cli::cli_abort(c(
			"{.arg newdata} has {length(newdata)} slice{?s} but the fit has {.val {n_t}}.",
			"i" = "Supply one dyadic covariate array per time period."))
	}
	# beta_mean_per_t handles 3-d dynamic_beta beta: a p x t matrix where
	# row k of column t is the posterior-mean coefficient k at period t. for
	# 2-d beta we replicate the single mean across t periods.
	if (length(dim(object$BETA)) == 3L) {
		beta_per_t <- apply(object$BETA, c(2, 3), mean)
		beta_names <- dimnames(object$BETA)[[2]]
	} else {
		beta_per_t <- matrix(colMeans(object$BETA),
		                     ncol(object$BETA),
		                     length(object$EZ %||% list(1)))
		beta_names <- colnames(object$BETA)
	}
	has_int <- !is.null(beta_names) && beta_names[1] == "intercept"
	is_dyad <- if (is.null(beta_names)) {
		c(FALSE, rep(TRUE, nrow(beta_per_t) - 1L))[seq_len(nrow(beta_per_t))]
	} else {
		grepl("_dyad$|\\.dyad$", beta_names)
	}
	n_dyad <- sum(is_dyad)
	# build a per-time-slice u v' product; tolerate static (2d), dynamic (3d),
	# absent (null), and zero-column r=0 cases
	UVPM <- if (!is.null(object$U) && !is.null(object$V) &&
	            ncol(object$U) > 0 && ncol(object$V) > 0) {
		if (length(dim(object$U)) == 3L) {
			Map(function(U, V) U %*% t(V),
			    asplit(object$U, 3), asplit(object$V, 3))
		} else {
			uv_static <- object$U %*% t(object$V)
			rep(list(uv_static), n_t)
		}
	} else NULL
	# the fitted per-period design (intercept + nodal + dyad slices, canonical
	# order) used to hold the intercept and any xrow / xcol contributions
	# fixed while only the dyadic covariates change. may be a per-period list
	# (f$x), or null for a dyad-only fit where .build_full_design synthesises
	# the intercept slice itself.
	fitted_design_list <- object$X
	# lame sorts actors alphabetically on ingest, so the stored additive /
	# multiplicative effects are in sorted order. realign each newdata array to
	# that order by name (no-op when names already match or are absent) before
	# building the design, otherwise a user passing newdata in the original
	# actor order is silently misaligned against a_t / b_t / uvpm.
	actor_order <- .fit_actor_order(object)
	ez_new <- vector("list", n_t)
	for (t in seq_len(n_t)) {
		Xn <- newdata[[t]]
		if (length(dim(Xn)) == 2L) Xn <- array(Xn, c(dim(Xn), 1L))
		Xn <- .reorder_newdata_actors(Xn, actor_order$rows, actor_order$cols)
		n_row <- dim(Xn)[1L]; n_col <- dim(Xn)[2L]
		fd_t <- if (is.list(fitted_design_list)) fitted_design_list[[t]] else fitted_design_list
		# name-safe full-design reconstruction: dyad slices from newdata, the
		# intercept + nodal (row/col) slices held at their fitted values. uses
		# the same helper as predict.ame() and the forecaster so the nodal
		# beta'x is included for fits with xrow / xcol.
		Xfull_t <- .build_full_design(Xn, beta_names, fd_t, n_row, n_col)
		xb <- matrix(0, n_row, n_col)
		for (k in seq_len(dim(Xfull_t)[3L])) xb <- xb + beta_per_t[k, t] * Xfull_t[, , k]
		# additive effects (per-time when dynamic_ab; else broadcast static)
		a_t <- if (!is.null(object$a_dynamic)) object$a_dynamic[, t] else object$APM
		b_t <- if (!is.null(object$b_dynamic)) object$b_dynamic[, t] else object$BPM
		if (!is.null(a_t)) for (i in seq_along(a_t)) xb[i, ] <- xb[i, ] + a_t[i]
		if (!is.null(b_t)) for (j in seq_along(b_t)) xb[, j] <- xb[, j] + b_t[j]
		if (!is.null(UVPM) && length(UVPM) >= t) xb <- xb + UVPM[[t]]
		# label with the canonical (sorted) actor order the prediction is in
		if (!is.null(actor_order$rows) && length(actor_order$rows) == n_row)
			rownames(xb) <- actor_order$rows
		if (!is.null(actor_order$cols) && length(actor_order$cols) == n_col)
			colnames(xb) <- actor_order$cols
		ez_new[[t]] <- xb
	}
	if (type == "link") return(ez_new)
	lapply(ez_new, function(ez) transform_to_response(ez, family))
}

####
#' Extract fitted values from LAME model
#'
#' @param object Fitted LAME model
#' @param ... Additional arguments
#'
#' @return List of fitted value matrices (one per time point)
#'
#' @method fitted lame
#' @export
fitted.lame <- function(object, ...) {
	return(object$YPM)
}

####
#' Extract residuals from LAME model
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
	type <- match.arg(type)

	if(!is.null(object$Y)) {
		Y <- object$Y
	} else {
		cli::cli_warn("original Y not stored in model object, residuals may be incorrect")
		Y <- object$YPM
	}

	# convert 3d array to list of matrices if needed
	if(is.array(Y) && length(dim(Y)) == 3) {
		n_time <- dim(Y)[3]
		Y <- lapply(seq_len(n_time), function(t) Y[,,t])
	}

	# lame stores y and ypm as lists
	s2 <- mean(object$VC[, ncol(object$VC)])

	mapply(function(y, ypm) {
		if(type == "response") {
			y - ypm
		} else {
			(y - ypm) / sqrt(s2)
		}
	}, Y, object$YPM, SIMPLIFY = FALSE)
}
####
