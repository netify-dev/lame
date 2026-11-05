# broom-style tidy() and glance() methods and a long-format prediction-draws helper

#' Tidy method for fitted `ame` / `lame` objects
#'
#' Returns a data frame with one row per estimated coefficient,
#' compatible with the \pkg{broom} idiom. For \code{dynamic_beta} fits
#' (3-D \code{BETA}), returns one row per coefficient \emph{per period}
#' with a \code{period} column. Standard errors are posterior standard
#' deviations; \code{statistic} is \code{estimate / std.error}.
#'
#' \strong{Note on \code{p.value}.} This column is included for
#' \pkg{broom} compatibility; it is a posterior summary, not a classical test:
#' the two-sided Normal tail probability of the posterior mean over the marginal
#' posterior standard deviation, matching \code{summary(fit)}. The
#' \code{conf.low} / \code{conf.high} columns are the credible interval. The
#' posterior probability of the sign is available directly from \code{x$BETA},
#' for example \code{mean(sign(BETA) == sign(mean(BETA)))}.
#'
#' Loaded as an S3 method against \code{generics::tidy} when the
#' \pkg{generics} package is available; works as
#' \code{tidy(fit)} either way once \pkg{broom} is loaded.
#'
#' @param x A fitted \code{ame} / \code{lame} object.
#' @param conf.int Logical; include 95\% credible interval columns
#'   (\code{conf.low}, \code{conf.high}). Default \code{TRUE}.
#' @param conf.level Confidence level for the interval. Default
#'   \code{0.95}.
#' @param ... Ignored.
#'
#' @return Data frame with columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic}, \code{p.value},
#'   \code{conf.low}, \code{conf.high}, and (for dynamic_beta fits)
#'   \code{period}.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'             nscan = 100, burn = 20, odens = 5, verbose = FALSE, chains = 1)
#' tidy(fit)
#' }
#'
#' @export
tidy.ame <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	.lame_require_draws(x, "tidy")
	alpha <- (1 - conf.level) / 2
	probs <- c(alpha, 1 - alpha)

	B <- x$BETA
	if (is.null(B)) {
		cli::cli_abort("Fit has no {.code $BETA} to tidy.")
	}
	# no coefficients: an empty table
	if (dim(B)[2L] == 0L) {
		out <- data.frame(
			term = character(0),
			estimate = numeric(0),
			std.error = numeric(0),
			statistic = numeric(0),
			p.value = numeric(0),
			conf.low = numeric(0),
			conf.high = numeric(0),
			stringsAsFactors = FALSE
		)
		if (length(dim(B)) == 3L) out$period <- character(0)
		return(out)
	}

	if (length(dim(B)) == 3L) {
		# dynamic_beta: [n_iter, p, t]
		coef_nms <- dimnames(B)[[2]] %||% paste0("v", seq_len(dim(B)[2]))
		per_nms <- dimnames(B)[[3]] %||% paste0("t", seq_len(dim(B)[3]))
		out_list <- list()
		for (t_idx in seq_along(per_nms)) {
			slice <- B[, , t_idx, drop = FALSE]
			slice <- matrix(slice, nrow = dim(slice)[1L], ncol = dim(slice)[2L])
			est <- colMeans(slice)
			se <- apply(slice, 2L, stats::sd)
			ci <- apply(slice, 2L, stats::quantile, probs = probs, na.rm = TRUE)
			out_list[[t_idx]] <- data.frame(
				term = coef_nms,
				period = per_nms[t_idx],
				estimate = est,
				std.error = se,
				statistic = est / pmax(se, .Machine$double.eps),
				p.value = 2 * (1 - stats::pnorm(abs(est / pmax(se, .Machine$double.eps)))),
				conf.low = ci[1L, ],
				conf.high = ci[2L, ],
				stringsAsFactors = FALSE
			)
		}
		out <- do.call(rbind, out_list)
		rownames(out) <- NULL
	} else {
		coef_nms <- colnames(B) %||% paste0("v", seq_len(ncol(B)))
		est <- colMeans(B)
		se <- apply(B, 2L, stats::sd)
		ci <- apply(B, 2L, stats::quantile, probs = probs, na.rm = TRUE)
		stat <- est / pmax(se, .Machine$double.eps)
		out <- data.frame(
			term = coef_nms,
			estimate = est,
			std.error = se,
			statistic = stat,
			p.value = 2 * (1 - stats::pnorm(abs(stat))),
			conf.low = ci[1L, ],
			conf.high = ci[2L, ],
			stringsAsFactors = FALSE
		)
		rownames(out) <- NULL
	}
	# a saturated nomination cap leaves row levels unfixed: withhold the
	# intercept and row-covariate coefficients
	coef_hide <- .lame_ab_coef_hidden_rows(.lame_ab_coef_reason(x), unique(out$term))
	if (length(coef_hide)) {
		num_cols <- intersect(
			c(
				"estimate", "std.error", "statistic", "p.value",
				"conf.low", "conf.high"
			),
			names(out)
		)
		out[out$term %in% coef_hide, num_cols] <- NA_real_
	}

	if (!isTRUE(conf.int)) {
		out$conf.low <- NULL; out$conf.high <- NULL
	}
	# a tibble when tibble is installed
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

#' @rdname tidy.ame
#' @export
tidy.lame <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	.lame_require_draws(x, "tidy")
	tidy.ame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

# fallback generic so `tidy(fit)` works without broom or generics; .onLoad
# registers the methods with generics::tidy when it is available

#' S3 generic for `tidy`
#'
#' Light-weight fallback generic so that calls of the form
#' \code{tidy(fit)} dispatch through R's S3 system even when the
#' \pkg{broom} or \pkg{generics} packages are not loaded. When either
#' is loaded, its generic resolves first; this generic only fires for
#' bare-namespace use.
#'
#' @param x An object to tidy.
#' @param ... Passed to the relevant method.
#' @return A data frame; method-specific schema.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' tidy(fit)
#' @export
tidy <- function(x, ...) UseMethod("tidy")

# register the tidy and glance methods with generics and broom when installed,
# without a hard dependency
#' @noRd
.register_tidy_methods <- function() {
	if (requireNamespace("generics", quietly = TRUE)) {
		try(
			registerS3method(
				"tidy",
				"ame",
				tidy.ame,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"lame",
				tidy.lame,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"ame_als",
				tidy.ame_als,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"lame_als",
				tidy.lame_als,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"boot_ame",
				tidy.boot_ame,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"ame",
				glance.ame,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"lame",
				glance.lame,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"ame_als",
				glance.ame_als,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"lame_als",
				glance.lame_als,
				envir = asNamespace("generics")
			),
			silent = TRUE
		)
	}
	if (requireNamespace("broom", quietly = TRUE)) {
		try(
			registerS3method(
				"tidy",
				"ame",
				tidy.ame,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"lame",
				tidy.lame,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"ame_als",
				tidy.ame_als,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"lame_als",
				tidy.lame_als,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"tidy",
				"boot_ame",
				tidy.boot_ame,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"ame",
				glance.ame,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"lame",
				glance.lame,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"ame_als",
				glance.ame_als,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
		try(
			registerS3method(
				"glance",
				"lame_als",
				glance.lame_als,
				envir = asNamespace("broom")
			),
			silent = TRUE
		)
	}
}

#' Long-format draws of the linear predictor for marginaleffects-style use
#'
#' Returns a long-format data frame with one row per
#' \code{(draw, i, j, period)} combination, giving the per-draw linear
#' predictor (or response-scale prediction) at each dyad and period.
#' Only the regression coefficients vary across \code{.draw}: the
#' additive (\code{a}, \code{b}) and multiplicative (\code{U}, \code{V})
#' effects are held at their posterior means, so the spread across draws
#' reflects coefficient uncertainty only, not the full posterior of the
#' linear predictor.
#' Intended for \pkg{marginaleffects}- / \pkg{tidybayes}-style
#' downstream summarisation: column names follow the
#' \code{.draw} / \code{.chain} / \code{.iteration} / \code{.value}
#' convention so that \code{tidybayes::spread_draws()} and
#' \code{marginaleffects::posterior_draws()} auto-dispatch on the
#' returned data frame. Actor and period names from \code{fit$Y}'s
#' dimnames are carried forward into the \code{actor_i},
#' \code{actor_j}, \code{period_label} columns.
#'
#' @param object A fitted \code{lame} object.
#' @param newdata Optional list of \code{T} dyadic covariate arrays for
#'   counterfactual evaluation. Default \code{NULL} uses in-sample
#'   covariates.
#' @param type One of \code{"link"} (default) or \code{"response"}. For a
#'   Poisson fit the response of a draw is its mean count,
#'   \eqn{e_t \exp(\eta + \sigma^2 / 2)}, with the draw's own residual variance
#'   and the exposure of the period.
#' @param n_draws Number of posterior draws to use. Default 100.
#' @param seed Optional RNG seed.
#'
#' @return A long-format data frame with columns \code{.chain},
#'   \code{.iteration}, \code{.draw}, \code{period}, \code{period_label},
#'   \code{i}, \code{j}, \code{actor_i}, \code{actor_j}, \code{.value}.
#'   Returned as a tibble when the \pkg{tibble} package is available;
#'   otherwise a plain data frame. The \code{.iteration} column identifies
#'   the retained coefficient iteration; \code{.draw} numbers the requested
#'   draws, including repeated selections.
#'
#' @examples
#' set.seed(1)
#' n <- 12; Tt <- 3
#' Y <- lapply(seq_len(Tt), function(t) {
#'   m <- matrix(rnorm(n * n), n, n); diag(m) <- NA; m
#' })
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' draws <- prediction_draws_long(fit, n_draws = 10, seed = 1)
#' head(draws)
#'
#' @export
prediction_draws_long <- function(
	object,
	newdata = NULL,
	type = c("link", "response"),
	n_draws = 100L,
	seed = NULL
) {
	type <- match.arg(type)
	if (!is.null(seed)) set.seed(seed)
	if (!inherits(object, "lame")) {
		cli::cli_abort("Only {.cls lame} fits are supported.")
	}
	if (is.null(object$EZ) || length(object$EZ) < 1L) {
		cli::cli_abort("Fit has no per-period EZ; cannot draw predictions.")
	}
	T_per <- length(object$EZ)
	n_a <- nrow(object$EZ[[1L]]); n_b <- ncol(object$EZ[[1L]])

	n_iter <- dim(object$BETA)[1L]
	idx <- sample.int(n_iter, n_draws, replace = (n_draws > n_iter))

	beta_is_dyn <- length(dim(object$BETA)) == 3L

	# coefficient names, used to rebuild a newdata design
	beta_names <- if (beta_is_dyn) dimnames(object$BETA)[[2]] else colnames(object$BETA)

	# covariates: newdata, else the stored design
	X_src <- if (!is.null(newdata)) newdata else object$Xlist
	if (is.null(X_src) || length(X_src) < T_per) {
		# repeat the last period's covariates
		X_src <- replicate(T_per, X_src[[length(X_src)]], simplify = FALSE)
	}
	# newdata carries dyadic slices only: rebuild the full design per period,
	# with intercept and nodal terms from the fit
	X_full <- lapply(seq_len(T_per), function(t) {
		Xt <- X_src[[t]]
		if (is.null(dim(Xt))) Xt <- array(Xt, c(n_a, n_b, 1L))
		if (length(dim(Xt)) == 2L) Xt <- array(Xt, c(dim(Xt), 1L))
		if (is.null(newdata)) return(Xt)
		actor_order <- .fit_actor_order(object)
		Xt <- .reorder_newdata_actors(Xt, actor_order$rows, actor_order$cols)
		fd_t <- if (is.list(object$Xlist)) object$Xlist[[t]] else object$Xlist
		.build_full_design(Xt, beta_names, fd_t, n_a, n_b)
	})
	# stored posterior-mean product (ULUPM on symmetric fits)
	uv_pm <- object[["UVPM", exact = TRUE]] %||% object[["ULUPM", exact = TRUE]]

	# actor and period labels from the dimnames of Y, else indices
	actor_i_names <- NULL
	actor_j_names <- NULL
	period_names <- NULL
	if (!is.null(object$Y)) {
		if (is.array(object$Y) && length(dim(object$Y)) >= 2L) {
			dn <- dimnames(object$Y)
			if (!is.null(dn)) {
				actor_i_names <- dn[[1L]]
				actor_j_names <- dn[[2L]]
				if (length(dn) >= 3L) period_names <- dn[[3L]]
			}
		} else if (is.list(object$Y) && length(object$Y) > 0L) {
			dn1 <- dimnames(object$Y[[1L]])
			if (!is.null(dn1)) {
				actor_i_names <- dn1[[1L]]
				actor_j_names <- dn1[[2L]]
			}
			period_names <- names(object$Y)
		}
	}
	if (is.null(actor_i_names)) actor_i_names <- paste0("i", seq_len(n_a))
	if (is.null(actor_j_names)) actor_j_names <- paste0("j", seq_len(n_b))
	if (is.null(period_names)) period_names <- paste0("t", seq_len(T_per))

	# a Poisson response draw is the mean count exp(eta + s2 / 2) with the
	# draw's own residual variance
	ve_draws <- if (type == "response" && identical(object$family, "poisson")) {
		.poisson_residual_variance(object)
		object$VC[, "ve"]
	}

	# linear predictor per draw and period
	rows <- vector("list", n_draws * T_per)
	rix <- 1L

	for (s_idx in seq_along(idx)) {
		i_draw <- idx[s_idx]
		for (t_idx in seq_len(T_per)) {
			if (beta_is_dyn) {
				bt <- object$BETA[i_draw, , t_idx]
			} else {
				bt <- object$BETA[i_draw, ]
			}
			Xt <- X_full[[t_idx]]
			# design slices follow the coefficient order, so the product is positional
			eta <- matrix(0, n_a, n_b)
			for (k in seq_len(min(length(bt), dim(Xt)[3L]))) {
				eta <- eta + bt[k] * Xt[, , k]
			}
			# additive effects, per period under dynamic_ab
			a_t <- if (!is.null(object$a_dynamic)) object$a_dynamic[, t_idx] else object$APM
			b_t <- if (!is.null(object$b_dynamic)) object$b_dynamic[, t_idx] else object$BPM
			# symmetric fits store one additive effect per actor and no BPM
			if (is.null(b_t) && isTRUE(object$symmetric)) b_t <- a_t
			if (!is.null(a_t)) eta <- eta + outer(a_t, rep(0, n_b), "+")
			if (!is.null(b_t)) eta <- sweep(eta, 2L, b_t, "+")
			# multiplicative term from the stored posterior-mean product, as in fit$EZ
			if (!is.null(uv_pm)) {
				eta <- eta + (if (length(dim(uv_pm)) == 3L) uv_pm[, , t_idx] else uv_pm)
			}
			if (isTRUE(object$symmetric) && n_a == n_b) eta <- (eta + t(eta)) / 2
			if (type == "response") {
				eta <- .fit_to_response(object, eta, t_idx, s2 = ve_draws[i_draw])
			}
			# long rows with tidybayes column names (.chain, .iteration, .draw, .value)
			grid <- expand.grid(
				i = seq_len(n_a),
				j = seq_len(n_b),
				KEEP.OUT.ATTRS = FALSE
			)
			rows[[rix]] <- data.frame(
				.chain = 1L,
				.iteration = i_draw,
				.draw = s_idx,
				period = t_idx,
				period_label = period_names[t_idx],
				i = grid$i,
				j = grid$j,
				actor_i = actor_i_names[grid$i],
				actor_j = actor_j_names[grid$j],
				.value = as.numeric(eta),
				stringsAsFactors = FALSE
			)
			rix <- rix + 1L
		}
	}

	out <- do.call(rbind, rows)
	# a tibble when tibble is installed
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

# glance: one-row model summary for broom and modelsummary

#' Glance method for fitted `ame` / `lame` objects
#'
#' One-row data frame summarising model-level statistics, in the
#' \pkg{broom} idiom. Used by \code{modelsummary::modelsummary()} and
#' similar tabling tools to populate the lower goodness-of-fit panel of
#' a regression table.
#'
#' @param x A fitted \code{ame} / \code{lame} object.
#' @param ... Ignored.
#'
#' @return A one-row data frame with columns:
#' \itemize{
#'   \item \code{nobs} -- number of observed dyads (NA cells excluded).
#'   \item \code{n_actors} -- number of distinct actors; for bipartite
#'     fits this is the row-actor count (a 12 x 9 fit reports 12).
#'   \item \code{n_row_actors}, \code{n_col_actors} -- row and column
#'     actor counts for bipartite fits (\code{NA} for unipartite).
#'   \item \code{n_periods} -- number of time periods (1 for \code{ame}).
#'   \item \code{n_stored} -- number of stored MCMC draws.
#'   \item \code{family} -- outcome family, such as \code{"normal"} or
#'     \code{"binary"}.
#'   \item \code{mode} -- \code{"unipartite"} or \code{"bipartite"}.
#'   \item \code{R} -- latent-space dimension (or max of \code{R_row},
#'     \code{R_col} for bipartite).
#'   \item \code{dynamic_uv}, \code{dynamic_ab}, \code{dynamic_beta} --
#'     logicals; whether each component is time-varying.
#'   \item \code{elpd_loo} -- leave-one-out expected log predictive
#'     density. Populated only when a \code{loo} object has been cached
#'     on the fit: fit with \code{save_log_lik = TRUE}, then attach it
#'     via \code{fit$loo <- loo(fit)}. NA otherwise (calling
#'     \code{loo(fit)} alone does not modify the fit object).
#' }
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'             nscan = 100, burn = 20, odens = 5, verbose = FALSE, chains = 1)
#' glance(fit)
#' }
#' @export
glance.ame <- function(x, ...) {
	.lame_require_draws(x, "glance")
	# exact lookup, since `$Y` partial-matches `$YPM`
	Y <- x[["Y", exact = TRUE]] %||% x$YPM
	nobs <- if (is.list(Y)) {
		sum(vapply(Y, function(yt) sum(is.finite(yt)), integer(1L)))
	} else if (is.array(Y)) {
		sum(is.finite(Y))
	} else NA_integer_

	is_bip <- identical(x$mode, "bipartite")
	nA_bip <- if (is_bip) length(x$APM %||% x$row_names %||% integer(0)) else NA_integer_
	nB_bip <- if (is_bip) length(x$BPM %||% x$col_names %||% integer(0)) else NA_integer_
	# for bipartite fits n_actors is the row count; n_col_actors holds the columns
	n_act <- if (is_bip) {
		as.integer(nA_bip)
	} else {
		length(x$APM %||% x$BPM %||% integer(0))
	}

	n_periods <- if (is.list(Y)) length(Y) else
	if (length(dim(Y)) == 3L) dim(Y)[3] else 1L
	n_stored <- if (is.matrix(x$BETA)) nrow(x$BETA) else
	if (length(dim(x$BETA)) == 3L) dim(x$BETA)[1L] else NA_integer_

	# elpd_loo only when a loo result is cached on the fit, else NA
	elpd <- NA_real_
	if (!is.null(x$loo) && inherits(x$loo, "loo")) {
		est <- x$loo$estimates
		if (is.matrix(est) && "elpd_loo" %in% rownames(est)) {
			elpd <- as.numeric(est["elpd_loo", "Estimate"])
		}
	}

	# an amen fit carries no $family or $R: take the rank from U, family stays NA
	is_legacy <- is.null(x[["family", exact = TRUE]]) &&
		is.null(x[["R", exact = TRUE]])
	R_used <- if (is_bip) {
		max(c(x$R_row %||% 0L, x$R_col %||% 0L, x$R %||% 0L))
	} else if (is_legacy) {
		if (is.matrix(x$U)) ncol(x$U) else 0L
	} else {
		x$R %||% 0L
	}

	# dynamic flags from the fit, else inferred from the shapes of U, the
	# additive effects and BETA
	dyn_uv_flag <- if (!is.null(x$dynamic_uv)) isTRUE(x$dynamic_uv) else
		(!is.null(x$U) && length(dim(x$U)) == 3L)
	dyn_ab_flag <- if (!is.null(x$dynamic_ab)) isTRUE(x$dynamic_ab) else
		(!is.null(x$a_dynamic) || !is.null(x$b_dynamic) ||
				(!is.null(x$APS) && length(dim(x$APS)) == 3L) ||
				(!is.null(x$BPS) && length(dim(x$BPS)) == 3L))
	dyn_beta_flag <- if (!is.null(x$dynamic_beta)) isTRUE(x$dynamic_beta) else
		(!is.null(x$BETA) && length(dim(x$BETA)) == 3L)

	out <- data.frame(
		nobs = as.integer(nobs),
		n_actors = as.integer(n_act),
		n_row_actors = as.integer(nA_bip),
		n_col_actors = as.integer(nB_bip),
		n_periods = as.integer(n_periods),
		n_stored = as.integer(n_stored),
		family = x$family %||% NA_character_,
		mode = if (is_bip) "bipartite" else "unipartite",
		R = as.integer(R_used),
		dynamic_uv = dyn_uv_flag,
		dynamic_ab = dyn_ab_flag,
		dynamic_beta = dyn_beta_flag,
		elpd_loo = elpd,
		stringsAsFactors = FALSE
	)
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

#' @rdname glance.ame
#' @export
glance.lame <- function(x, ...) glance.ame(x, ...)

# fallback generic so `glance(fit)` works without broom or generics
#' S3 generic for `glance`
#'
#' Light-weight fallback so \code{glance(fit)} dispatches through S3
#' even when \pkg{broom} or \pkg{generics} is not loaded.
#'
#' @param x An object to glance at.
#' @param ... Passed to the method.
#' @return A one-row data frame.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' glance(fit)
#' @export
glance <- function(x, ...) UseMethod("glance")


#' Tidy method for fitted `ame_als` / `lame_als` objects
#'
#' Returns a data frame with one row per regression coefficient,
#' compatible with the \pkg{broom} idiom, so that ALS fits compose with
#' \pkg{modelsummary} / \pkg{kableExtra} pipelines next to MCMC fits.
#' Standard errors come from the bootstrap object attached to
#' \code{x$bootstrap} when present (preferred, fully propagated), and
#' otherwise from \code{\link{vcov.ame_als}}: the Laplace covariance of an
#' empirical-Bayes fit (the default for normal, binary and Poisson fits) or
#' the conditional sandwich covariance of a fixed-effects fit.
#' \code{statistic} is \code{estimate / std.error};
#' \code{p.value} is the Normal-approximation two-sided tail
#' \eqn{2(1 - \Phi(|z|))} from the bootstrap or sandwich standard error. It is
#' a Wald-style summary for the point estimator, not a posterior probability.
#'
#' One row is returned for every regression coefficient of \code{coef(fit)}:
#' the intercept, the dyadic covariates and the node covariates. The
#' node-covariate coefficients of a fixed-effects fit have no sandwich
#' covariance, so they keep their estimate with \code{NA} uncertainty. Additive
#' (\code{a}, \code{b}) and multiplicative (\code{U}, \code{V}) parameters are
#' not included; use \code{\link{ame_als_bootstrap}} and inspect the bootstrap
#' object directly if you need them.
#'
#' @param x A fitted \code{ame_als} / \code{lame_als} object.
#' @param conf.int Logical; include \code{conf.low} / \code{conf.high}
#'   columns. Default \code{TRUE}.
#' @param conf.level Confidence level. Default \code{0.95}.
#' @param ... Passed to \code{\link{vcov.ame_als}}
#'   (e.g. \code{cluster = "dyad"}).
#' @return Data frame with columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic}, \code{p.value},
#'   \code{conf.low}, \code{conf.high}, plus a \code{se_source} column
#'   recording \code{"bootstrap"}, \code{"laplace"} (an empirical-Bayes
#'   normal, binary or Poisson fit, whose covariance covers every coefficient,
#'   node-covariate coefficients included) or \code{"sandwich"}.
#' @examples
#' \donttest{
#' data(YX_nrm)
#' fit <- ame_als(Y = YX_nrm$Y, Xdyad = YX_nrm$X[, , 2:3],
#'                family = "normal", R = 1, verbose = FALSE)
#' tidy(fit)
#' }
#' @export
tidy.ame_als <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	if (length(conf.level) != 1L || !is.finite(conf.level) ||
			conf.level <= 0 || conf.level >= 1) {
		cli::cli_abort("{.arg conf.level} must be a single number in (0, 1).")
	}
	# regression coefficients only, with bootstrap summaries when available,
	# else the sandwich or Laplace covariance
	if (!is.null(x$bootstrap)) {
		coef_nms <- names(x$coefficients)
		boot_coefs <- x$bootstrap$coefs
		# the coefficient columns of the bootstrap replicates
		nms <- intersect(coef_nms, colnames(boot_coefs))
		if (length(nms) == 0L) {
			nms <- coef_nms
		}
		# point estimates are the fit's coefficients, not the bootstrap mean
		est <- x$coefficients[nms]
		if (!is.null(boot_coefs) && all(nms %in% colnames(boot_coefs))) {
			se <- apply(
				boot_coefs[, nms, drop = FALSE],
				2L,
				stats::sd,
				na.rm = TRUE
			)
			pr <- .ae_boot_probs(x$bootstrap, conf.level)
			ci_lo <- apply(
				boot_coefs[, nms, drop = FALSE],
				2L,
				stats::quantile,
				probs = pr[1],
				na.rm = TRUE
			)
			ci_hi <- apply(
				boot_coefs[, nms, drop = FALSE],
				2L,
				stats::quantile,
				probs = pr[2],
				na.rm = TRUE
			)
			ci <- cbind(ci_lo, ci_hi)
			rownames(ci) <- nms
			# parameters the bootstrap withholds
			w <- intersect(.ae_boot_withheld_names(x$bootstrap, "coefficient"), nms)
			se[w] <- NA_real_
			ci[w, ] <- NA_real_
		} else {
			# standard errors from the confint() half-width
			ci_all <- suppressMessages(stats::confint(x, level = conf.level, ...))
			ci <- ci_all[intersect(rownames(ci_all), nms), , drop = FALSE]
			nms <- rownames(ci)
			est <- x$coefficients[nms]
			a <- (1 - conf.level) / 2
			z <- stats::qnorm(1 - a)
			se <- (ci[, 2L] - ci[, 1L]) / (2 * z)
		}
		se_source <- "bootstrap"
	} else {
		V <- suppressMessages(stats::vcov(x, ...))
		# terms outside the covariance (node covariates) get NA uncertainty
		nms <- names(x$coefficients)
		est <- x$coefficients
		se <- stats::setNames(rep(NA_real_, length(nms)), nms)
		se[rownames(V)] <- sqrt(diag(V))
		a <- (1 - conf.level) / 2
		zq <- stats::qnorm(c(a, 1 - a))
		ci <- cbind(est + zq[1] * se, est + zq[2] * se)
		rownames(ci) <- nms
		se_source <- if (!is.null(x$eb$vcov)) "laplace" else "sandwich"
	}
	stat <- est / pmax(se, .Machine$double.eps)
	out <- data.frame(
		term = nms,
		estimate = unname(est),
		std.error = unname(se),
		statistic = unname(stat),
		p.value = 2 * (1 - stats::pnorm(abs(unname(stat)))),
		conf.low = unname(ci[, 1L]),
		conf.high = unname(ci[, 2L]),
		se_source = se_source,
		stringsAsFactors = FALSE
	)
	if (!isTRUE(conf.int)) {
		out$conf.low <- NULL; out$conf.high <- NULL
	}
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

#' @rdname tidy.ame_als
#' @export
tidy.lame_als <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	tidy.ame_als(x, conf.int = conf.int, conf.level = conf.level, ...)
}

#' Glance method for fitted `ame_als` / `lame_als` objects
#'
#' One-row data frame summarising an ALS fit, compatible with
#' \pkg{broom} / \pkg{modelsummary}. Reports observation count, actor
#' count, latent dimension, family, mode, ALS convergence flag, and
#' iteration count.
#'
#' @param x A fitted \code{ame_als} / \code{lame_als} object.
#' @param ... Ignored.
#' @return One-row data frame with columns \code{nobs},
#'   \code{n_actors}, \code{n_periods}, \code{family}, \code{mode},
#'   \code{R}, \code{converged}, \code{iterations}, \code{se_source}.
#' @examples
#' \donttest{
#' data(YX_nrm)
#' fit <- ame_als(Y = YX_nrm$Y, Xdyad = YX_nrm$X[, , 2:3],
#'                family = "normal", R = 1, verbose = FALSE)
#' glance(fit)
#' }
#' @export
glance.ame_als <- function(x, ...) {
	Y <- x$Y
	nobs <- if (is.list(Y)) {
		sum(vapply(Y, function(yt) sum(is.finite(yt)), integer(1L)))
	} else if (is.array(Y)) {
		sum(is.finite(Y))
	} else NA_integer_

	is_bip <- identical(x$mode, "bipartite")
	a_len <- length(x$a %||% x$a_row %||% integer(0))
	b_len <- length(x$b %||% x$b_col %||% integer(0))
	# for bipartite fits n_actors is the row count, as in glance.ame
	n_act <- if (is_bip) as.integer(a_len) else as.integer(max(a_len, b_len))
	nA_bip <- if (is_bip) as.integer(a_len) else NA_integer_
	nB_bip <- if (is_bip) as.integer(b_len) else NA_integer_

	n_periods <- if (isTRUE(x$longitudinal)) {
		if (is.list(Y)) length(Y) else
		if (length(dim(Y)) == 3L) dim(Y)[3] else 1L
	} else 1L

	se_src <- if (!is.null(x$bootstrap)) {
		"bootstrap"
	} else if (!is.null(x$eb$vcov)) {
		"laplace"
	} else {
		"sandwich"
	}
	nnm <- x$non_normal_method %||% NA_character_
	out <- data.frame(
		nobs = as.integer(nobs),
		n_actors = n_act,
		n_row_actors = nA_bip,
		n_col_actors = nB_bip,
		n_periods = as.integer(n_periods),
		family = x$family %||% NA_character_,
		mode = if (is_bip) "bipartite" else "unipartite",
		R = as.integer(x$R %||% 0L),
		non_normal_method = nnm,
		converged = isTRUE(x$converged),
		iterations = as.integer(x$iterations %||% NA_integer_),
		se_source = se_src,
		stringsAsFactors = FALSE
	)
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

#' @rdname glance.ame_als
#' @export
glance.lame_als <- function(x, ...) glance.ame_als(x, ...)


#' Tidy method for a standalone bootstrap object (`boot_ame`)
#'
#' `ame_als_bootstrap()` returns an object of class \code{"boot_ame"}
#' (not \code{"ame_als"}); this tidy method exposes the bootstrap
#' estimates as a broom-style data frame so the standalone object
#' composes with \pkg{modelsummary} the same way an embedded
#' \code{ame_als(..., bootstrap = N)} fit does.
#'
#' @param x A \code{boot_ame} object.
#' @param conf.level Confidence level. Default \code{0.95}.
#' @param ... Ignored.
#' @return Data frame with columns \code{term, estimate, std.error,
#'   statistic, p.value, conf.low, conf.high, se_source}.
#' @export
tidy.boot_ame <- function(x, conf.level = 0.95, ...) {
	if (length(conf.level) != 1L || !is.finite(conf.level) ||
			conf.level <= 0 || conf.level >= 1) {
		cli::cli_abort("{.arg conf.level} must be a single number in (0, 1).")
	}
	nms <- x$param_names %||% colnames(x$coefs)
	est <- x$point_est
	se <- x$se
	# stored intervals are 95%; other levels come from the replicates
	if (abs(conf.level - 0.95) < 1e-8 && !is.null(x$ci_lo)) {
		ci_lo <- x$ci_lo; ci_hi <- x$ci_hi
	} else {
		pr <- .ae_boot_probs(x, conf.level)
		ci_lo <- apply(x$coefs, 2L, stats::quantile, probs = pr[1], na.rm = TRUE)
		ci_hi <- apply(x$coefs, 2L, stats::quantile, probs = pr[2], na.rm = TRUE)
	}
	# parameters the bootstrap withholds
	w <- .ae_boot_withheld_names(x, "coefficient")
	se[nms %in% w] <- NA_real_
	ci_lo[nms %in% w] <- NA_real_
	ci_hi[nms %in% w] <- NA_real_
	stat <- est / pmax(se, .Machine$double.eps)
	out <- data.frame(
		term = nms,
		estimate = unname(est),
		std.error = unname(se),
		statistic = unname(stat),
		p.value = 2 * (1 - stats::pnorm(abs(unname(stat)))),
		conf.low = unname(ci_lo),
		conf.high = unname(ci_hi),
		se_source = "bootstrap",
		stringsAsFactors = FALSE
	)
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}
