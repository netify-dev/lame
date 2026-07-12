# exact rolling-origin leave-future-out cv. for each leave-out period
# t, refits on y[, , 1:(t-1)] and scores y[, , t]. defaults to the
# last 3 periods because each leave-out needs its own refit.

#' Exact rolling-origin leave-future-out cross-validation
#'
#' For a fitted \code{lame} object with \eqn{T} periods, refits the
#' model on the first \eqn{t - 1} periods (for each \eqn{t} in
#' \code{periods}) and computes the expected log predictive density
#' (elpd) of period \eqn{t} under the refit. Returns the summed elpd
#' across all leave-out periods, along with per-period contributions.
#'
#' @param fit A fitted \code{lame} object.
#' @param periods Integer vector of leave-out periods to evaluate.
#'   Default is the last 3 periods (\code{tail(seq_len(T), 3L)}). Each
#'   period \eqn{t} must satisfy \eqn{t \ge 2}.
#' @param refit Logical; if \code{TRUE} (default), refits on the
#'   training window. If \code{FALSE}, uses the original posterior
#'   means (a much rougher approximation).
#' @param ... Passed to the refit \code{lame()} call (typically
#'   \code{nscan}, \code{burn}, \code{odens}, \code{verbose}).
#'
#' @return A list with \code{elpd_lfo} (total summed elpd), \code{pointwise}
#'   (per-dyad log-density at each leave-out period; a list of numeric
#'   vectors, one per period -- \code{unlist(pointwise)} gives a flat vector
#'   suitable for \code{loo::loo_compare()}-style stacking),
#'   \code{p_lfo} (effective number of parameters), \code{per_period}
#'   (data frame with \code{period}, \code{elpd}, \code{n_obs}), and
#'   \code{periods} (the periods evaluated).
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE)
#' lfo_res <- lfo(fit, periods = 4L, refit = TRUE,
#'                nscan = 50, burn = 10, odens = 5, verbose = FALSE)
#' print(lfo_res)
#' }
#'
#' @importFrom stats predict
#' @export
lfo <- function(fit, periods = NULL, refit = TRUE, ...) {
	if (!inherits(fit, "lame")) {
		cli::cli_abort("{.arg fit} must be a fitted {.cls lame} object.")
	}

	# extract per-period y and x -- fit$y can be a list or a 3-d array
	Y_src <- fit$Y_obs %||% fit$Y
	Y_list <- if (is.list(Y_src)) {
		Y_src
	} else if (is.array(Y_src) && length(dim(Y_src)) == 3L) {
		lapply(seq_len(dim(Y_src)[3L]), function(t) Y_src[, , t])
	} else {
		cli::cli_abort("Could not extract per-period Y from the fit.")
	}
	T_per <- length(Y_list)
	X_list <- fit$Xlist
	if (is.null(X_list)) X_list <- rep(list(NULL), T_per)

	if (is.null(periods)) {
		periods <- utils::tail(seq_len(T_per), 3L)
	}
	periods <- sort(unique(as.integer(periods)))
	if (any(periods < 2L) || any(periods > T_per)) {
		cli::cli_abort("Every leave-out period must be in {.val 2:T}.")
	}

	# extract fit-time arguments for refit
	family <- fit$family %||% "normal"
	R_rank <- if (!is.null(fit$U)) NCOL(fit$U) else 0L
	dyn_beta <- fit$dynamic_beta %||% FALSE
	dyn_beta_kind <- fit$dynamic_beta_kind %||% "ar1"

	per_period <- data.frame(period = periods, elpd = NA_real_, n_obs = NA_integer_)
	# per-dyad pointwise log-density at each leave-out period (one element of
	# the list per period). concatenated, these are what loo::loo_compare()
	# style workflows would consume as the "pointwise elpd" vector.
	pointwise_list <- vector("list", length(periods))
	names(pointwise_list) <- as.character(periods)

	for (rr in seq_along(periods)) {
		t_out <- periods[rr]
		Y_train <- Y_list[seq_len(t_out - 1L)]
		X_train <- X_list[seq_len(t_out - 1L)]
		Y_test  <- Y_list[[t_out]]
		X_test  <- X_list[[t_out]]

		refit_args <- list(
			Y = Y_train, Xdyad = X_train,
			family = family, R = R_rank,
			dynamic_beta = dyn_beta,
			dynamic_beta_kind = dyn_beta_kind,
			verbose = FALSE, plot = FALSE
		)
		# allow user to override / pass nscan, burn, odens, etc.
		usr <- list(...)
		refit_args[names(usr)] <- usr

		if (isTRUE(refit)) {
			refit_obj <- do.call(lame, refit_args)
		} else {
			refit_obj <- fit
		}

		# 1-step forecast under the refit
		fc_mean <- predict(refit_obj, h = 1L, type = "link")[[1L]]
		# pointwise elpd on y_test using the family-implied density
		obs_idx <- which(!is.na(Y_test))
		y_obs <- Y_test[obs_idx]
		ez_obs <- fc_mean[obs_idx]

		# variance scale: use the posterior mean s2 from the refit
		s2 <- if (!is.null(refit_obj$VC) && ncol(refit_obj$VC) > 0L)
			mean(refit_obj$VC[, ncol(refit_obj$VC)], na.rm = TRUE) else 1
		s2 <- max(s2, 1e-8)

		# compute per-dyad log-density first; sum over dyads for the
		# per-period elpd. exposing the pointwise vector is what lets
		# downstream workflows feed lfo() output into `loo::loo_compare()`
		# style stacking or pareto-k diagnostics.
		ll_dyad <- switch(family,
			normal  = stats::dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE),
			binary  = {
				p_hat <- stats::pnorm(ez_obs)
				p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
				y_obs * log(p_hat) + (1 - y_obs) * log(1 - p_hat)
			},
			cbin    = {
				p_hat <- stats::pnorm(ez_obs)
				p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
				y_obs * log(p_hat) + (1 - y_obs) * log(1 - p_hat)
			},
			poisson = {
				lam <- pmin(exp(ez_obs), 1e8)
				stats::dpois(round(y_obs), lambda = lam, log = TRUE)
			},
			# fall back to gaussian on link for any unsupported family
			stats::dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
		)
		pointwise_list[[rr]] <- ll_dyad
		per_period$elpd[rr]  <- sum(ll_dyad)
		per_period$n_obs[rr] <- length(obs_idx)
	}

	elpd_lfo <- sum(per_period$elpd)
	# rough effective param count via the spread of per-period elpd
	p_lfo <- max(0, stats::var(per_period$elpd) * (nrow(per_period) - 1L) / 2)

	res <- list(
		elpd_lfo    = elpd_lfo,
		p_lfo       = p_lfo,
		per_period  = per_period,
		periods     = periods,
		refit       = refit,
		# pointwise (per-dyad) log-density at each leave-out period; a list
		# with one numeric vector per period. concatenate into one long
		# vector with unlist() to feed loo::loo_compare()-style stacking.
		pointwise   = pointwise_list
	)
	class(res) <- c("lfo_lame", "list")
	res
}

#' Print method for lfo() results
#' @param x A \code{lfo_lame} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   summary of the leave-future-out cross-validation results.
#' @export
#' @method print lfo_lame
print.lfo_lame <- function(x, ...) {
	cli::cli_h2("Exact rolling-origin LFO")
	cli::cli_text("Periods evaluated: {.val {x$periods}}")
	cli::cli_text("Refit per leave-out: {.val {x$refit}}")
	cli::cli_text("Total elpd_lfo: {.val {signif(x$elpd_lfo, 5)}}")
	print(x$per_period, row.names = FALSE)
	invisible(x)
}
