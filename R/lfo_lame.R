# exact rolling-origin leave-future-out CV: each leave-out period t is scored by a refit on periods 1:(t-1)

# row-wise log(mean(exp(x))) on the log scale, so one underflowing draw cannot sink the average
#' @noRd
.lfo_log_mean_exp_rows <- function(m) {
	if (nrow(m) == 0L) return(numeric(0))
	mx <- apply(m, 1L, max)
	mx[!is.finite(mx)] <- 0
	mx + log(rowMeans(exp(m - mx)))
}

# align a held-out period to the refit's forecast grid; with changing actors, cells the two do not
# share are dropped
#' @noRd
.lfo_align_heldout <- function(Y_t, row_nm, col_nm, n_r, n_c, t_out) {
	if (identical(dim(Y_t), c(n_r, n_c)) &&
			(is.null(rownames(Y_t)) || is.null(colnames(Y_t)) || is.null(row_nm) || is.null(col_nm) ||
					(identical(rownames(Y_t), row_nm) && identical(colnames(Y_t), col_nm)))) return(Y_t)
	if (is.null(row_nm) || is.null(col_nm) ||
			is.null(rownames(Y_t)) || is.null(colnames(Y_t))) {
		cli::cli_abort(c(
			"Period {.val {t_out}} is {.val {nrow(Y_t)}} x {.val {ncol(Y_t)}} but the refit forecasts a {.val {n_r}} x {.val {n_c}} grid.",
			"i" = "Actor names on {.arg Y} are needed to align periods with different actor sets."
		))
	}
	out <- matrix(NA_real_, n_r, n_c, dimnames = list(row_nm, col_nm))
	ri <- intersect(row_nm, rownames(Y_t))
	ci <- intersect(col_nm, colnames(Y_t))
	if (length(ri) > 0L && length(ci) > 0L) out[ri, ci] <- Y_t[ri, ci]
	out
}

# subset only the time margin; static designs remain static
.lfo_subset_covariate <- function(x, periods, nodal = FALSE) {
	if (is.null(x)) return(NULL)
	if (is.list(x)) return(x[periods])
	d <- dim(x)
	if (length(d) == if (nodal) 3L else 4L) {
		idx <- lapply(d, seq_len); idx[[length(d)]] <- periods
		return(do.call(`[`, c(list(x), idx, list(drop = FALSE))))
	}
	x
}

#' Exact rolling-origin leave-future-out cross-validation
#'
#' For a fitted \code{lame} object with \eqn{T} periods, refits the
#' model on the first \eqn{t - 1} periods (for each \eqn{t} in
#' \code{periods}) and computes the expected log predictive density
#' (elpd) of period \eqn{t} under the refit. Returns the summed elpd
#' across all leave-out periods, along with per-period contributions.
#'
#' The one-step forecast is taken with the leave-out period's own dyadic
#' covariates supplied as \code{newdata}. This matters because
#' \code{predict(fit, h = 1)} called without \code{newdata} reuses the last
#' fitted period's covariate array, which scores the held-out period against
#' the wrong covariates whenever they vary over time.
#'
#' Each dyad is scored by its posterior predictive log density,
#' \eqn{\log \frac{1}{S} \sum_s p(y \mid \eta_s)}, averaging the family
#' density over the posterior draws of the forecast linear predictor rather
#' than evaluating it once at the posterior-mean predictor. For
#' normal and Poisson models, each predictor uses its matching residual
#' variance draw. Poisson scoring integrates the lognormal observation
#' effect and uses held-out exposure. Supported families are normal, binary
#' and Poisson; symmetric networks score each unordered dyad once.
#' Refits retain the original covariate roles, selected dynamic coefficients
#' and user priors, recomputing empirical defaults on training data only.
#'
#' @param fit A fitted \code{lame} object.
#' @param periods Integer vector of leave-out periods to evaluate.
#'   Default is the last 3 periods (\code{tail(seq_len(T), 3L)}). Each
#'   period \eqn{t} must satisfy \eqn{t \ge 2}.
#' @param refit Logical; if \code{TRUE} (default), refits on the
#'   training window. If \code{FALSE}, scores each period under the
#'   original fit, which was estimated on every period including the one
#'   being scored, so the elpd is optimistic.
#' @param ... Passed to the refit \code{lame()} call (typically
#'   \code{nscan}, \code{burn}, \code{odens}, \code{verbose}).
#'
#' @return A list with \code{elpd_lfo} (total summed elpd), \code{pointwise}
#'   (per-dyad log-density at each leave-out period; a list of numeric
#'   vectors, one per period -- \code{unlist(pointwise)} gives a flat vector
#'   suitable for \code{loo::loo_compare()}-style stacking),
#'   \code{p_lfo} (\code{NA}; an effective parameter count is not
#'   estimated by this rolling-origin calculation), \code{per_period}
#'   (data frame with \code{period}, \code{elpd}, \code{n_obs}), and
#'   \code{periods} (the periods evaluated).
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE, chains = 1)
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
	# short training windows often put rho_beta at the unit root; count those notes and report them once
	unit_root_hits <- 0L
	muffle_unit_root <- function(expr) {
		withCallingHandlers(expr, warning = function(w) {
			if (grepl("rho_beta", conditionMessage(w), fixed = TRUE)) {
				unit_root_hits <<- unit_root_hits + 1L
				invokeRestart("muffleWarning")
			}
		})
	}

	# per-period Y; fit$Y can be a list or a 3-d array
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

	family <- fit$family %||% "normal"
	if (!family %in% c("normal", "binary", "poisson")) {
		cli::cli_abort(
			"LFO scoring supports normal, binary and Poisson outcomes.",
			class = "lame_lfo_family"
		)
	}
	if (isTRUE(refit) && (is.null(fit$refit_options) || is.null(fit$refit_data))) {
		cli::cli_abort(
			c(
				"This fit lacks the original model specification for faithful LFO refits.",
				"i" = "Refit with the current version of {.fn lame} before calling {.fn lfo}."
			),
			class = "lame_lfo_refit"
		)
	}

	per_period <- data.frame(period = periods, elpd = NA_real_, n_obs = NA_integer_)
	# per-dyad pointwise log-density at each leave-out period
	pointwise_list <- vector("list", length(periods))
	names(pointwise_list) <- as.character(periods)

	for (rr in seq_along(periods)) {
		t_out <- periods[rr]
		Y_train <- Y_list[seq_len(t_out - 1L)]
		Y_test <- Y_list[[t_out]]
		X_test <- X_list[[t_out]]

		refit_args <- fit$refit_options %||% list()
		refit_args$Y <- Y_train
		for (nm in c("Xdyad", "Xrow", "Xcol")) {
			refit_args[nm] <- list(.lfo_subset_covariate(
				fit$refit_data[[nm]],
				seq_len(t_out - 1L),
				nodal = nm != "Xdyad"
			))
		}
		for (nm in c("time_index", "period_exposure")) {
			if (!is.null(refit_args[[nm]])) refit_args[[nm]] <- refit_args[[nm]][seq_len(t_out - 1L)]
		}
		refit_args$chains <- 1L
		refit_args$verbose <- FALSE
		refit_args$plot <- FALSE

		# user overrides such as nscan, burn, odens
		usr <- list(...)
		refit_args[names(usr)] <- usr

		if (isTRUE(refit)) {
			refit_obj <- muffle_unit_root(do.call(lame, refit_args))
		} else {
			refit_obj <- fit
		}

		# the held-out design supplies nodal and dyadic covariates; the forecaster maps named actors
		eta_draws <- muffle_unit_root(.forecast_lame(
			refit_obj,
			h = 1L,
			type = "link",
			newdesign = list(X_test),
			by_draw = TRUE
		))

		draw_index <- attr(eta_draws, "draw_index")
		dm <- dim(eta_draws)
		fc_rows <- dimnames(eta_draws)[[1L]]
		fc_cols <- dimnames(eta_draws)[[2L]]
		Y_score <- .lfo_align_heldout(
			Y_test,
			fc_rows,
			fc_cols,
			dm[1L],
			dm[2L],
			t_out
		)
		# [n_row, n_col, h = 1, draw] to one column of dyads per draw
		dim(eta_draws) <- c(dm[1L] * dm[2L], dm[4L])
		if (!identical(refit_obj$mode, "bipartite")) diag(Y_score) <- NA_real_
		if (isTRUE(refit_obj$symmetric)) Y_score[lower.tri(Y_score, diag = TRUE)] <- NA_real_
		obs_idx <- which(is.finite(Y_score))
		y_obs <- Y_score[obs_idx]
		if (family == "binary" && any(!y_obs %in% c(0, 1)))
			cli::cli_abort("Binary held-out outcomes must be zero or one.")
		if (family == "poisson" && any(y_obs < 0 | y_obs != floor(y_obs)))
			cli::cli_abort("Poisson held-out outcomes must be nonnegative integer counts.")
		eta_obs <- eta_draws[obs_idx, , drop = FALSE]
		if (any(!is.finite(eta_obs))) cli::cli_abort("Held-out cells require finite predictive draws.")
		n_obs_t <- length(y_obs)

		# residual variance from the same posterior draws as the forecasts
		s2 <- if (family == "binary") rep(1, ncol(eta_obs)) else {
			if (is.null(draw_index) || is.null(refit_obj$VC) || !"ve" %in% colnames(refit_obj$VC))
				cli::cli_abort("LFO scoring requires matched residual variance draws.")
			variance <- refit_obj$VC[draw_index, "ve"]
			if (any(!is.finite(variance) | variance <= 0))
				cli::cli_abort("LFO scoring requires finite positive residual variance draws.")
			variance
		}
		exposure <- (fit$period_exposure %||% rep(1, T_per))[t_out]
		if (family == "poisson" && (length(exposure) != 1L || !is.finite(exposure) || exposure < 0))
			cli::cli_abort("Held-out Poisson exposure must be finite and nonnegative.")
		lp <- matrix(NA_real_, n_obs_t, ncol(eta_obs))
		for (ss in seq_len(ncol(eta_obs))) {
			eta <- eta_obs[, ss]
			lp[, ss] <- switch(
				family,
				normal = stats::dnorm(y_obs, mean = eta, sd = sqrt(s2[ss]), log = TRUE),
				binary = ifelse(
					y_obs == 1,
					stats::pnorm(eta, log.p = TRUE),
					stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
				),
				poisson = if (exposure == 0) ifelse(y_obs == 0, 0, -Inf) else
					.lame_pois_marginal_loglik(y_obs, eta + log(exposure), s2[ss])
			)
		}

		lp <- matrix(lp, nrow = n_obs_t)
		ll_dyad <- .lfo_log_mean_exp_rows(lp)
		pointwise_list[[rr]] <- ll_dyad
		per_period$elpd[rr] <- sum(ll_dyad)
		per_period$n_obs[rr] <- n_obs_t
	}

	elpd_lfo <- sum(per_period$elpd)
	# between-period score variation is not an effective parameter count, so none is reported
	p_lfo <- NA_real_

	if (unit_root_hits > 0L) {
		cli::cli_inform(
			c(
				"i" = "{unit_root_hits} of the {nrow(per_period)} window refits had a {.code rho_beta} posterior at the unit root, as short training windows often do; their one-step forecasts carry the unit-root variance."
			)
		)
	}
	res <- list(
		elpd_lfo = elpd_lfo,
		p_lfo = p_lfo,
		per_period = per_period,
		periods = periods,
		refit = refit,
		# per-dyad log-density for each leave-out period
		pointwise = pointwise_list
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
	.lame_cli_stdout({
		cli::cli_h2("Exact rolling-origin LFO")
		cli::cli_text("Periods evaluated: {.val {x$periods}}")
		cli::cli_text("Refit per leave-out: {.val {x$refit}}")
		cli::cli_text("Total elpd_lfo: {.val {signif(x$elpd_lfo, 5)}}")
		cli::cli_verbatim(utils::capture.output(print(x$per_period, row.names = FALSE)))
	})
	invisible(x)
}
