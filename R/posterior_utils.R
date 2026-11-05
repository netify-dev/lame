#' Options for saving posterior samples during MCMC
#'
#' @param save_UV Logical; whether to save samples of U and V matrices (default FALSE)
#' @param save_ab Logical; whether to save samples of additive effects a and b (default FALSE)
#' @param thin_UV Integer; thinning interval for U/V samples (default 10)
#' @param thin_ab Integer; thinning interval for a/b samples (default 10)
#'
#' @return List of posterior saving options, to be passed to the
#'   \code{posterior_opts} argument of \code{\link{ame}} or \code{\link{lame}}:
#'   \code{ame(Y, ..., posterior_opts = posterior_options(save_UV = TRUE))} or
#'   \code{lame(Y, ..., posterior_opts = posterior_options(save_UV = TRUE))}.
#'   With \code{save_UV = TRUE}, static fits store \code{U_samples} /
#'   \code{V_samples} (plus \code{G_samples} for bipartite fits) and
#'   \code{dynamic_uv} fits store the per-period trajectory in
#'   \code{U_draws} / \code{V_draws}; see \code{\link{latent_positions}} and
#'   \code{\link{simulate_posterior}} for summarising them. With
#'   \code{save_ab = TRUE}, static effects use \code{a_samples} /
#'   \code{b_samples} (the same draws in both for a symmetric fit, which has
#'   one effect per actor); longitudinal dynamic effects use
#'   actor-by-period-by-draw arrays \code{a_draws} / \code{b_draws}. \code{thin_ab} is measured
#'   in retained coefficient draws; \code{thin_ab = 1} keeps them aligned.
#'
#' @examples
#' opts <- posterior_options(save_UV = TRUE, thin_UV = 5)
#' opts
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
posterior_options <- function(
	save_UV = FALSE,
	save_ab = FALSE,
	thin_UV = 10,
	thin_ab = 10
) {
	list(
		save_UV = save_UV,
		save_ab = save_ab,
		thin_UV = thin_UV,
		thin_ab = thin_ab
	)
}

#' Simulate posterior distributions from fitted AME model
#'
#' @param fit Fitted ame model object
#' @param component Character; which component to simulate: "UV", "ab", "beta", "Y"
#' @param n_samples Number of posterior samples to return. Defaults to
#'   \code{NULL}, which uses every saved draw when the fit stores them
#'   (see \code{\link{posterior_options}}) and 100 otherwise. A smaller
#'   value draws a random subset, which widens Monte Carlo error in the
#'   tails of any interval computed from the result.
#' @param seed Random seed for reproducibility
#'
#' @return Array or matrix of posterior samples
#'
#' @details
#' This function can simulate posterior distributions even when they weren't
#' saved during MCMC, by using the posterior means and variance components.
#'
#' For more accurate posteriors, use posterior_options() during model fitting
#' to save the actual MCMC samples. Longitudinal additive trajectories require
#' saved draws; their joint posterior is not reconstructed from summaries.
#'
#' @examples
#' \donttest{
#' # Fit a model with multiplicative effects
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Get posterior samples of regression coefficients
#' beta_post <- simulate_posterior(fit, "beta", n_samples = 50)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
simulate_posterior <- function(
	fit,
	component = c("UV", "ab", "beta", "Y"),
	n_samples = NULL,
	seed = NULL
) {

	component <- match.arg(component)

	if (!is.null(seed)) set.seed(seed)

	# longitudinal posterior predictive networks use their own simulator
	if (inherits(fit, "lame") && identical(component, "Y")) {
		cli::cli_abort(c(
			"{.code component = \"{component}\"} is not available for longitudinal fits.",
			"i" = "Use {.fn simulate} on the fit for posterior predictive networks."
		))
	}
	has_samples <- switch(
		component,
		UV = (!is.null(fit$U_samples) && !is.null(fit$V_samples)) ||
			(!is.null(fit$U_draws) && !is.null(fit$V_draws)),
		ab = (!is.null(fit$a_samples) && !is.null(fit$b_samples)) ||
			(!is.null(fit$a_draws) && !is.null(fit$b_draws)),
		beta = !is.null(fit$BETA),
		Y = FALSE # always simulated
	)

	if (inherits(fit, "lame") && identical(component, "ab") && !has_samples) {
		cli::cli_abort(c(
			"Additive posterior draws were not saved in this longitudinal fit.",
			"i" = "Refit with {.code posterior_options(save_ab = TRUE)} to retain joint actor trajectories."
		))
	}

	if (has_samples && component != "Y") {
		cli::cli_text("Using saved MCMC samples for {.field {component}}")

		if (component == "UV") {
			# dynamic_uv fits store 4-D trajectories in U_draws/V_draws
			if (!is.null(fit$U_samples) && !is.null(fit$V_samples)) {
				return(sample_from_saved_UV(fit, n_samples))
			}
			return(sample_from_saved_UV_dynamic(fit, n_samples))
		} else if (component == "ab") {
			return(sample_from_saved_ab(fit, n_samples))
		} else if (component == "beta") {
			return(sample_from_saved_beta(fit, n_samples))
		}

	} else {
		# no stored draws to size against: default to 100
		if (is.null(n_samples)) n_samples <- 100
		if (component == "UV") {
			simulate_UV_posterior(fit, n_samples)
		} else if (component == "ab") {
			simulate_ab_posterior(fit, n_samples)
		} else if (component == "beta") {
			return(simulate_beta_posterior(fit, n_samples))
		} else if (component == "Y") {
			return(simulate_Y_posterior(fit, n_samples))
		}
	}
}

#' Sample from saved U/V posteriors
#' @noRd
sample_from_saved_UV <- function(fit, n_samples) {
	n_saved <- dim(fit$U_samples)[3]
	idx <- .saved_draw_index(n_samples, n_saved)
	n_samples <- length(idx)

	n_row <- nrow(fit$U_samples)
	n_col <- nrow(fit$V_samples)

	UV_samples <- array(NA, dim = c(n_row, n_col, n_samples))

	# bipartite fits have the multiplicative term U G V' with a per-draw G
	has_G <- !is.null(fit$G_samples)
	for (i in 1:n_samples) {
		Ui <- matrix(fit$U_samples[, , idx[i]], nrow = n_row)
		Vi <- matrix(fit$V_samples[, , idx[i]], nrow = n_col)
		UV_samples[, , i] <- if (has_G) {
			Ui %*% fit$G_samples[, , idx[i]] %*% t(Vi)
		} else {
			Ui %*% t(Vi)
		}
	}

	return(UV_samples)
}

#' Sample per-period latent similarity from a dynamic_uv fit
#'
#' dynamic_uv fits store the full latent trajectory as 4-D draws
#' (\code{U_draws}/\code{V_draws}: actor x rank x time x draw) instead of the
#' 3-D \code{U_samples}/\code{V_samples} of a static fit. Returns an
#' actor x actor x time x draw array of \eqn{U_t V_t'} (rows x cols x time x
#' draw for bipartite; symmetric fits store \eqn{V = U L} per draw, so the
#' product is the \eqn{U L U'} similarity).
#' @noRd
sample_from_saved_UV_dynamic <- function(fit, n_samples) {
	dU <- dim(fit$U_draws) # n_row x R x T x draws
	idx <- .saved_draw_index(n_samples, dU[4])
	n_samples <- length(idx)
	n_row <- dU[1]
	n_col <- dim(fit$V_draws)[1]
	TT <- dU[3]

	# bipartite dynamic fits: each period is U_t G V_t' with one G per draw
	has_G <- !is.null(fit$G_samples)
	UV <- array(NA_real_, dim = c(n_row, n_col, TT, n_samples))
	for (i in seq_len(n_samples)) {
		Gi <- if (has_G) fit$G_samples[, , idx[i]] else NULL
		for (t in seq_len(TT)) {
			# matrix() keeps the rank dimension when R == 1
			Ut <- matrix(fit$U_draws[, , t, idx[i]], nrow = n_row)
			Vt <- matrix(fit$V_draws[, , t, idx[i]], nrow = n_col)
			UV[, , t, i] <- if (has_G) Ut %*% Gi %*% t(Vt) else Ut %*% t(Vt)
		}
	}
	UV
}

#' Simulate U/V from posterior -- requires saved MCMC samples
#' @noRd
simulate_UV_posterior <- function(fit, n_samples) {

	cli::cli_abort(c(
		"Cannot simulate UV posterior without saved MCMC samples.",
		"i" = "Re-fit the model with {.code posterior_opts = posterior_options(save_UV = TRUE)}",
		"i" = "Then call {.fn simulate_posterior} again."
	))
}

#' Sample from saved a/b posteriors
#' @noRd
sample_from_saved_ab <- function(fit, n_samples) {
	if (!is.null(fit$a_draws) && !is.null(fit$b_draws)) {
		n_saved <- dim(fit$a_draws)[3L]
		idx <- .saved_draw_index(n_samples, n_saved)
		return(list(
			a = fit$a_draws[, , idx, drop = FALSE],
			b = fit$b_draws[, , idx, drop = FALSE]
		))
	}
	n_saved <- ncol(fit$a_samples)
	idx <- .saved_draw_index(n_samples, n_saved)

	ab_samples <- list(
		a = fit$a_samples[, idx, drop = FALSE],
		b = fit$b_samples[, idx, drop = FALSE]
	)

	return(ab_samples)
}

#' Simulate a/b from posterior -- requires saved MCMC samples
#' @noRd
simulate_ab_posterior <- function(fit, n_samples) {

	cli::cli_abort(c(
		"Cannot simulate ab posterior without saved MCMC samples.",
		"i" = "Re-fit the model with {.code posterior_opts = posterior_options(save_ab = TRUE)}",
		"i" = "Then call {.fn simulate_posterior} again."
	))
}

#' Sample from saved beta posteriors
#' @noRd
sample_from_saved_beta <- function(fit, n_samples) {
	d <- dim(fit$BETA)
	idx <- .saved_draw_index(n_samples, d[1])
	if (length(d) == 3L) {
		# dynamic_beta: BETA is draws x coefs x time
		return(fit$BETA[idx, , , drop = FALSE])
	}
	return(fit$BETA[idx, , drop = FALSE])
}

#' Index into stored posterior draws
#'
#' Returns every draw in order when \code{n_samples} is \code{NULL} or at
#' least the number stored, so intervals use the full posterior sample; a
#' smaller request draws a random subset without replacement.
#' @noRd
.saved_draw_index <- function(n_samples, n_saved) {
	if (is.null(n_samples) || n_samples >= n_saved) {
		if (!is.null(n_samples) && n_samples > n_saved) {
			cli::cli_warn(c(
				"Requested {.val {n_samples}} samples but only {.val {n_saved}} are stored.",
				"i" = "Returning all {.val {n_saved}} stored draws."
			))
		}
		return(seq_len(n_saved))
	}
	sort(sample(n_saved, n_samples))
}

#' Simulate beta from approximate posterior
#' @noRd
simulate_beta_posterior <- function(fit, n_samples) {

	if (is.null(fit$BETA)) {
		return(matrix(0, n_samples, 0))
	}

	beta_mean <- colMeans(fit$BETA)
	beta_cov <- cov(fit$BETA)

	if (requireNamespace("MASS", quietly = TRUE)) {
		beta_samples <- MASS::mvrnorm(n_samples, beta_mean, beta_cov)
	} else {
		p <- length(beta_mean)
		beta_samples <- matrix(NA, n_samples, p)
		for (i in 1:p) {
			beta_samples[, i] <- rnorm(n_samples, beta_mean[i], sqrt(beta_cov[i, i]))
		}
	}

	return(beta_samples)
}

#' Simulate Y from posterior predictive distribution
#' @noRd
simulate_Y_posterior <- function(fit, n_samples) {
	out <- simulate(fit, nsim = n_samples, type = "posterior")
	y <- simplify2array(out$Y)
	attr(y, "draw_index") <- out$draw_index
	y
}

#' Extract posterior quantiles for model components
#'
#' @param fit Fitted ame model object
#' @param component Character; which component: "beta", "UV", "ab". The
#'   additive effects "ab" come back as a list of the row (\code{a}) and
#'   column (\code{b}) quantiles; a symmetric fit has one effect per actor, so
#'   \code{b} equals \code{a}.
#' @param probs Numeric vector of probabilities for quantiles
#'
#' @return Matrix or array of posterior quantiles
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' X <- array(rnorm(n * n), dim = c(n, n, 1))
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, Xdyad = X, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' posterior_quantiles(fit, component = "beta")
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
posterior_quantiles <- function(
	fit,
	component = c("beta", "UV", "ab"),
	probs = c(0.025, 0.5, 0.975)
) {

	component <- match.arg(component)

	if (component == "beta") {
		if (is.null(fit$BETA)) {
			cli::cli_abort("No regression coefficients in model")
		}

		quantiles <- apply(fit$BETA, 2, quantile, probs = probs, na.rm = TRUE)
		rownames(quantiles) <- paste0("q", probs * 100)

		return(quantiles)

	} else if (component == "UV") {
		if (is.null(fit$U) || NCOL(fit$U) == 0L) {
			cli::cli_abort("No multiplicative effects in this fit ({.code R = 0}).")
		}
		if (!((!is.null(fit$U_samples) && !is.null(fit$V_samples)) || (!is.null(fit$U_draws) && !is.null(fit$V_draws)))) {
			cli::cli_abort(c(
				"Quantiles of the multiplicative effects need the saved position draws.",
				"i" = "Refit with {.code posterior_opts = posterior_options(save_UV = TRUE)}."
			))
		}

		UV_samples <- simulate_posterior(fit, "UV", n_samples = 1000)
		n <- dim(UV_samples)[1]
		m <- dim(UV_samples)[2]

		UV_quantiles <- array(NA, dim = c(n, m, length(probs)))

		for (i in 1:n) {
			for (j in 1:m) {
				UV_quantiles[i, j, ] <- quantile(
					UV_samples[i, j, ],
					probs = probs,
					na.rm = TRUE
				)
			}
		}

		dimnames(UV_quantiles)[[3]] <- paste0("q", probs * 100)

		return(UV_quantiles)

	} else if (component == "ab") {
		# a symmetric fit has one effect per actor and no separate column effects
		if (is.null(fit$APM) || (is.null(fit$BPM) && !isTRUE(fit$symmetric))) {
			cli::cli_abort("No additive effects in model")
		}

		ab_samples <- simulate_posterior(fit, "ab")

		a_quantiles <- apply(
			ab_samples$a,
			seq_len(length(dim(ab_samples$a)) - 1L),
			quantile,
			probs = probs,
			na.rm = TRUE
		)
		b_quantiles <- apply(
			ab_samples$b,
			seq_len(length(dim(ab_samples$b)) - 1L),
			quantile,
			probs = probs,
			na.rm = TRUE
		)

		rownames(a_quantiles) <- paste0("q", probs * 100)
		rownames(b_quantiles) <- paste0("q", probs * 100)

		return(list(a = a_quantiles, b = b_quantiles))
	}
}
