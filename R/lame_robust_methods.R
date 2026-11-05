#' Summarize a robust longitudinal network model
#'
#' Summaries use retained posterior draws. Coefficient summaries include the
#' static covariate coefficients and every time intercept; `summary()` also
#' reports estimated variance components and AR coefficients. Unused factor
#' hyperparameters are excluded when the rank is zero.
#'
#' @param object,x A `lame_robust` fit, or its summary for the summary print method.
#' @param parm Coefficient names or positive column indices; omit for all.
#' @param level Central credible-interval probability, strictly between zero and one.
#' @param ... Additional arguments are accepted only by `fitted()`, which passes
#'   them to [lame::predict.lame_robust()].
#' @return `coef()` returns posterior coefficient means. `vcov()` returns their
#'   posterior covariance. `confint()` returns central posterior quantiles.
#'   `fitted()` returns fitted mean networks; `nobs()` counts observed
#'   off-diagonal responses. `summary()` returns a `summary.lame_robust` object
#'   containing a posterior summary table and resolved priors. Print methods
#'   return their input invisibly.
#' @seealso [lame::lame_robust()], [lame::compute_mcmc_diagnostics()]
#' @name summary.lame_robust
NULL

.lame_robust_unused <- function(dots) {
	if (length(dots)) cli::cli_abort("Unsupported argument{?s}: {.arg {names(dots)}}.", call = NULL)
}

.lame_robust_design <- function(object, h, newdata) {
	D <- object$data
	if (length(h) != 1L || !is.numeric(h) || is.complex(h) || !is.finite(h) || h < 0 || h != floor(h) || h > .Machine$integer.max)
		cli::cli_abort("h must be a nonnegative integer.", call = NULL)
	if (h == 0L) {
		if (!is.null(newdata)) cli::cli_abort("newdata is only used for positive forecast horizons.", call = NULL)
		return(list(X = D$X, periods = D$tt, names = D$years))
	}
	if (D$p == 0L) {
		if (!is.null(newdata)) cli::cli_abort("This model has no dyadic covariates; omit newdata.", call = NULL)
		return(list(X = matrix(numeric(), D$n * D$n * h, 0L), periods = as.integer(h), names = paste0("h", seq_len(h))))
	}
	if (is.matrix(newdata) || length(dim(newdata)) == 3L) newdata <- list(newdata)
	if (!is.list(newdata) || length(newdata) != h)
		cli::cli_abort("newdata must supply one dyadic-covariate array per forecast period.", call = NULL)
	future_names <- names(newdata)
	if (is.null(future_names)) future_names <- paste0("h", seq_len(h))
	if (anyNA(future_names) || any(!nzchar(future_names)) || anyDuplicated(future_names))
		cli::cli_abort("Forecast period names must be nonmissing and unique.", call = NULL)
	# align the design dimensions through the fitting boundary; a duplicated
	# period allows a one-step horizon
	base <- matrix(0, D$n, D$n, dimnames = list(D$actor_names, D$actor_names))
	newdata <- lapply(newdata, function(x) {
		if (is.matrix(x)) x <- array(x, c(dim(x), 1L), dimnames = c(dimnames(x), list(NULL)))
		if (length(dim(x)) == 3L && dim(x)[3L] == D$p && is.null(dimnames(x)[[3L]])) {
			dn <- dimnames(x)
			if (is.null(dn)) dn <- vector("list", 3L)
			dn[[3L]] <- colnames(D$X)
			dimnames(x) <- dn
		}
		x
	})
	x <- if (h == 1L) c(newdata, newdata) else newdata
	names(x) <- NULL
	y <- rep(list(base), length(x))
	valid <- .lame_robust_input(y, x, object$R)$Xdyad
	if (dim(valid[[1L]])[3L] != D$p || !setequal(dimnames(valid[[1L]])[[3L]], colnames(D$X)))
		cli::cli_abort("Forecast covariates must match the fitted covariate names and count.", call = NULL)
	valid <- lapply(valid[seq_len(h)], function(x) x[, , colnames(D$X), drop = FALSE])
	X <- do.call(rbind, lapply(valid, function(x) matrix(x, D$n * D$n, D$p)))
	list(X = X, periods = as.integer(h), names = future_names)
}

.lame_robust_indices <- function(object, nsim = NULL) {
	n <- nrow(object$draws$hyper)
	if (is.null(nsim)) return(seq_len(n))
	if (length(nsim) != 1L || !is.numeric(nsim) || is.complex(nsim) || !is.finite(nsim) || nsim < 1 ||
			nsim != floor(nsim) || nsim > .Machine$integer.max)
		cli::cli_abort("nsim must be a positive integer or NULL for every retained draw.", call = NULL)
	if (nsim == n) return(seq_len(n))
	sort(sample.int(n, as.integer(nsim), replace = nsim > n))
}

.lame_robust_prediction <- function(object, design, index, h, innovations, noise) {
	D <- object$data
	z <- object$draws
	n <- D$n
	rank <- D$rank
	nt <- design$periods
	result <- array(
		NA_real_,
		c(n, n, nt, length(index)),
		dimnames = list(D$actor_names, D$actor_names, design$names, NULL)
	)
	for (j in seq_along(index)) {
		s <- index[j]
		v <- z$hyper[s, ]
		if (h > 0L) {
			mu <- z$mu[s, D$tt]
			a <- z$A[, D$tt, s]
			b <- z$B[, D$tt, s]
			U <- matrix(z$U[, , D$tt, s], n, rank)
			V <- matrix(z$V[, , D$tt, s], n, rank)
		}
		for (t in seq_len(nt)) {
			if (h == 0L) {
				mu <- z$mu[s, t]
				a <- z$A[, t, s]
				b <- z$B[, t, s]
				U <- matrix(z$U[, , t, s], n, rank)
				V <- matrix(z$V[, , t, s], n, rank)
			} else {
				mu <- v["rho_mu"] * mu
				a <- z$a0[, s] + v["rho_ab"] * (a - z$a0[, s])
				b <- z$b0[, s] + v["rho_ab"] * (b - z$b0[, s])
				U <- v["rho_uv"] * U
				V <- v["rho_uv"] * V
				if (innovations) {
					mu <- mu + stats::rnorm(1L, sd = sqrt(v["var_mu"]))
					df_actor <- object$actor_df %||% Inf
					tau <- if (is.finite(df_actor)) stats::rgamma(2L, df_actor / 2, df_actor / 2) else c(1, 1)
					ea <- stats::rnorm(n, sd = sqrt(v["var_ab"] / tau[1L]))
					eb <- stats::rnorm(n, sd = sqrt(v["var_ab"] / tau[2L]))
					a <- a + ea - mean(ea)
					b <- b + eb - mean(eb)
					if (rank) {
						U <- U + matrix(stats::rnorm(n * rank, sd = sqrt(v["var_uv"])), n, rank)
						V <- V + matrix(stats::rnorm(n * rank, sd = sqrt(v["var_uv"])), n, rank)
					}
				}
			}
			xbeta <- matrix(0, n, n)
			if (D$p) {
				rows <- seq_len(n * n) + (t - 1L) * n * n
				xbeta <- matrix(design$X[rows, , drop = FALSE] %*% z$beta[s, ], n, n)
			}
			eta <- xbeta + mu + outer(a, b, "+")
			if (rank) eta <- eta + tcrossprod(U, V)
			if (noise) {
				df_observation <- object$observation_df %||% Inf
				omega <- if (is.finite(df_observation)) matrix(
					stats::rgamma(
						n * n,
						df_observation / 2,
						df_observation / 2
					),
					n,
					n
				) else 1
				eta <- eta + matrix(
					stats::rnorm(
						n * n,
						sd = sqrt(v["ve"] / (z$lambda[, , s] * omega))
					),
					n,
					n
				)
			}
			diag(eta) <- NA_real_
			result[, , t, j] <- eta
		}
	}
	attr(result, "draw_index") <- index
	result
}

#' Predict and simulate robust longitudinal networks
#'
#' Posterior means and joint draws from a [lame::lame_robust()] fit. Prediction uses
#' the conditional mean; simulation also includes observation noise with the
#' persistent directed-pair precision from the same posterior iteration.
#'
#' @param object A `lame_robust` fit.
#' @param h Nonnegative integer forecast horizon. Zero returns observed periods.
#' @param newdata For a positive horizon with covariates, a list of one dyadic
#'   design array per future period. For one period, a matrix or array is also
#'   accepted. Actor and covariate names are aligned to the fit. All future
#'   covariates must be supplied and finite off the diagonal. Omit otherwise.
#' @param type Either `"response"` or `"link"`; both denote the continuous mean.
#' @param by_draw Logical; return joint posterior mean-network draws. Forecast
#'   draws include future latent innovations. The default posterior mean
#'   integrates those innovations analytically before averaging over draws.
#' @param nsim Number of posterior draws. `NULL` uses every retained iteration.
#'   Simulation defaults to one. Smaller requests sample without replacement;
#'   larger requests sample with replacement. Selected indices are sorted.
#' @param seed Optional seed preserving the caller's RNG state.
#' @param ... No additional arguments are supported.
#'
#' @details
#' Future intercepts, actor deviations, and factors follow their sampled AR(1)
#' laws. Actor innovations respect the sum-to-zero constraints. Finite
#' `actor_df` draws a new precision for each actor role and forecast period. Persistent
#' actor levels and pair precisions remain the same within each simulated
#' trajectory. Simulation at `h = 0` replicates observations conditional on
#' the saved mean paths and persistent pair scales. Finite `observation_df`
#' draws fresh temporary residual precisions for every replicated observation,
#' including at `h = 0`; fitted observation precisions are not reused.
#' Simulation at a positive horizon draws both future
#' states and observation noise. Future latent uncertainty is included when
#' `by_draw = TRUE`; its Monte Carlo average need not equal the analytically
#' integrated default point forecast exactly.
#'
#' @return Prediction returns a sender-by-receiver-by-period array, or a
#'   sender-by-receiver-by-period-by-draw array when `by_draw = TRUE`.
#'   Simulation always returns the four-dimensional form. Diagonals are `NA`.
#'   Draw arrays carry a `draw_index` attribute identifying the joint posterior
#'   iteration. Singleton dimensions are preserved.
#' @name predict.lame_robust
#' @export
predict.lame_robust <- function(
	object,
	h = 0L,
	newdata = NULL,
	type = c("response", "link"),
	by_draw = FALSE,
	nsim = NULL,
	seed = NULL,
	...
) {
	.lame_robust_unused(list(...))
	type <- match.arg(type)
	if (length(by_draw) != 1L || !is.logical(by_draw) || is.na(by_draw))
		cli::cli_abort("by_draw must be TRUE or FALSE.", call = NULL)
	design <- .lame_robust_design(object, h, newdata)
	.lame_robust_seed(seed, function() {
		index <- .lame_robust_indices(object, nsim)
		z <- .lame_robust_prediction(
			object,
			design,
			index,
			h,
			innovations = isTRUE(by_draw) && h > 0L,
			noise = FALSE
		)
		if (by_draw) return(z)
		answer <- apply(z, c(1L, 2L, 3L), mean)
		array(answer, dim(z)[1:3], dimnames = dimnames(z)[1:3])
	})
}

#' @rdname predict.lame_robust
#' @export
simulate.lame_robust <- function(
	object,
	nsim = 1L,
	seed = NULL,
	h = 0L,
	newdata = NULL,
	...
) {
	.lame_robust_unused(list(...))
	design <- .lame_robust_design(object, h, newdata)
	.lame_robust_seed(seed, function() {
		index <- .lame_robust_indices(object, nsim)
		.lame_robust_prediction(
			object,
			design,
			index,
			h,
			innovations = h > 0L,
			noise = TRUE
		)
	})
}

.lame_robust_coefficient_draws <- function(object) {
	B <- object$draws$beta
	mu <- object$draws$mu
	colnames(mu) <- paste0("intercept[", object$data$years, "]")
	cbind(B, mu)
}

#' @rdname summary.lame_robust
#' @export
coef.lame_robust <- function(object, ...) {
	.lame_robust_unused(list(...))
	colMeans(.lame_robust_coefficient_draws(object))
}

#' @rdname summary.lame_robust
#' @export
vcov.lame_robust <- function(object, ...) {
	.lame_robust_unused(list(...))
	stats::cov(.lame_robust_coefficient_draws(object))
}

#' @rdname summary.lame_robust
#' @export
confint.lame_robust <- function(object, parm, level = 0.95, ...) {
	.lame_robust_unused(list(...))
	if (length(level) != 1L || !is.numeric(level) || is.complex(level) || !is.finite(level) || level <= 0 || level >= 1)
		cli::cli_abort("level must be between zero and one.", call = NULL)
	B <- .lame_robust_coefficient_draws(object)
	if (!missing(parm)) {
		valid_names <- is.character(parm) && length(parm) > 0L && !anyNA(parm) && all(parm %in% colnames(B))
		valid_indices <- is.numeric(parm) && !is.complex(parm) && length(parm) > 0L &&
			all(is.finite(parm)) && all(parm == floor(parm) & parm >= 1L & parm <= ncol(B))
		if (!valid_names && !valid_indices)
			cli::cli_abort("parm must contain existing coefficient names or positive column indices.", call = NULL)
		B <- B[, parm, drop = FALSE]
	}
	probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
	answer <- t(apply(B, 2L, stats::quantile, probs = probs, names = FALSE))
	colnames(answer) <- paste0(format(100 * probs, trim = TRUE), "%")
	answer
}

#' @rdname summary.lame_robust
#' @export
fitted.lame_robust <- function(object, ...) predict.lame_robust(object, ...)
#' @rdname summary.lame_robust
#' @export
nobs.lame_robust <- function(object, ...) {
	.lame_robust_unused(list(...))
	sum(object$data$observed)
}

#' @rdname summary.lame_robust
#' @method summary lame_robust
#' @export
summary.lame_robust <- function(object, ...) {
	.lame_robust_unused(list(...))
	B <- .lame_robust_coefficient_draws(object)
	H <- object$draws$hyper
	if (!object$R) H <- H[, !colnames(H) %in% c("rho_uv", "var_uv"), drop = FALSE]
	x <- cbind(B, H)
	quantiles <- t(apply(x, 2L, stats::quantile, c(.025, .5, .975), names = FALSE))
	tab <- cbind(
		mean = colMeans(x),
		sd = apply(x, 2L, stats::sd),
		q025 = quantiles[, 1L],
		median = quantiles[, 2L],
		q975 = quantiles[, 3L]
	)
	structure(
		list(
			call = object$call,
			family = object$family,
			df = object$df,
			actor_df = object$actor_df %||% Inf,
			observation_df = object$observation_df %||% Inf,
			n_actors = object$data$n,
			n_periods = object$data$tt,
			R = object$R,
			n_draws = nrow(x),
			coefficients = tab,
			prior = object$prior
		),
		class = "summary.lame_robust"
	)
}

#' @rdname summary.lame_robust
#' @method print lame_robust
#' @export
print.lame_robust <- function(x, ...) {
	.lame_robust_unused(list(...))
	.lame_cli_stdout({
		n_draws <- nrow(x$draws$hyper)
		cli::cli_h1("Robust longitudinal network model")
		info <- c(
			"*" = "{.val {x$data$n}} actors, {.val {x$data$tt}} periods, rank {.val {x$R}}, {.val {n_draws}} retained draws",
			"*" = if (is.finite(x$df)) {
				"Persistent directed-pair variance, df = {.val {x$df}}"
			} else if (is.finite(x$observation_df %||% Inf)) {
				"Common residual scale across directed pairs"
			} else {
				"Homogeneous Gaussian residual variance"
			}
		)
		if (is.finite(x$observation_df %||% Inf)) {
			info <- c(info, "*" = "Additional observation precision, df = {.val {x$observation_df}}")
		}
		if (is.finite(x$actor_df %||% Inf)) {
			info <- c(info, "*" = "Actor innovation precision varies by role and period, df = {.val {x$actor_df}}")
		}
		cli::cli_bullets(info)
	})
	invisible(x)
}

#' @rdname summary.lame_robust
#' @method print summary.lame_robust
#' @export
print.summary.lame_robust <- function(x, ...) {
	.lame_robust_unused(list(...))
	.lame_cli_stdout({
		cli::cli_h1("Robust longitudinal network model: posterior summary")
		cli::cli_bullets(c(
			"*" = "Family: {.field {x$family}}",
			"*" = "Retained draws: {.val {x$n_draws}}"
		))
		cli::cli_h3("Parameters")
		cli::cli_verbatim(utils::capture.output(print(x$coefficients)))
	})
	invisible(x)
}
