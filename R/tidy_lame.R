# broom-style tidy() method plus a long-format prediction-draws helper.

#' Tidy method for fitted `ame` / `lame` objects
#'
#' Returns a data frame with one row per estimated coefficient,
#' compatible with the \pkg{broom} idiom. For \code{dynamic_beta} fits
#' (3-D \code{BETA}), returns one row per coefficient \emph{per period}
#' with a \code{period} column. Standard errors are posterior standard
#' deviations; \code{statistic} is \code{estimate / std.error}.
#'
#' \strong{Note on \code{p.value}.} This column is included for
#' \pkg{broom} compatibility but is \emph{not} a classical test. It is a
#' two-sided Normal approximation based on the posterior mean and marginal
#' posterior standard deviation, matching the calculation in
#' \code{summary(fit)}. Use it as a compact signal that the marginal posterior
#' is far from zero, and report it alongside the \code{conf.low} /
#' \code{conf.high} credible interval. When sign certainty matters, compute it
#' directly from \code{x$BETA}, for example
#' \code{mean(sign(BETA) == sign(mean(BETA)))}.
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
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             nscan = 100, burn = 20, odens = 5, verbose = FALSE)
#' tidy(fit)
#' }
#'
#' @export
tidy.ame <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	alpha <- (1 - conf.level) / 2
	probs <- c(alpha, 1 - alpha)

	B <- x$BETA
	if (is.null(B)) {
		cli::cli_abort("Fit has no {.code $BETA} to tidy.")
	}

	if (length(dim(B)) == 3L) {
		# dynamic_beta: [n_iter, p, t]
		coef_nms <- dimnames(B)[[2]] %||% paste0("v", seq_len(dim(B)[2]))
		per_nms  <- dimnames(B)[[3]] %||% paste0("t", seq_len(dim(B)[3]))
		out_list <- list()
		for (t_idx in seq_along(per_nms)) {
			slice <- B[, , t_idx, drop = FALSE]
			slice <- matrix(slice, nrow = dim(slice)[1L], ncol = dim(slice)[2L])
			est <- colMeans(slice)
			se  <- apply(slice, 2L, stats::sd)
			ci  <- apply(slice, 2L, stats::quantile, probs = probs, na.rm = TRUE)
			out_list[[t_idx]] <- data.frame(
				term      = coef_nms,
				period    = per_nms[t_idx],
				estimate  = est,
				std.error = se,
				statistic = est / pmax(se, .Machine$double.eps),
				p.value   = 2 * (1 - stats::pnorm(abs(est / pmax(se, .Machine$double.eps)))),
				conf.low  = ci[1L, ],
				conf.high = ci[2L, ],
				stringsAsFactors = FALSE
			)
		}
		out <- do.call(rbind, out_list)
		rownames(out) <- NULL
	} else {
		coef_nms <- colnames(B) %||% paste0("v", seq_len(ncol(B)))
		est <- colMeans(B)
		se  <- apply(B, 2L, stats::sd)
		ci  <- apply(B, 2L, stats::quantile, probs = probs, na.rm = TRUE)
		stat <- est / pmax(se, .Machine$double.eps)
		out <- data.frame(
			term      = coef_nms,
			estimate  = est,
			std.error = se,
			statistic = stat,
			p.value   = 2 * (1 - stats::pnorm(abs(stat))),
			conf.low  = ci[1L, ],
			conf.high = ci[2L, ],
			stringsAsFactors = FALSE
		)
		rownames(out) <- NULL
	}

	if (!isTRUE(conf.int)) {
		out$conf.low <- NULL; out$conf.high <- NULL
	}
	# return as a tibble when tibble is available so the result composes
	# naturally with dplyr / ggplot2 verbs in a tidyverse pipeline
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

#' @rdname tidy.ame
#' @export
tidy.lame <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	tidy.ame(x, conf.int = conf.int, conf.level = conf.level, ...)
}

# generics::tidy is the canonical generic for broom-style methods.
# provide a fallback generic so `tidy(fit)` works even if neither
# `broom` nor `generics` is loaded. when `generics` (or `broom`,
# which imports the generic) is loaded, .onload registers tidy.ame
# / tidy.lame against generics::tidy so that `broom::tidy(fit)`
# also dispatches correctly.

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
#' @export
tidy <- function(x, ...) UseMethod("tidy")

# register tidy.ame / tidy.lame against generics::tidy (and broom::tidy
# when broom is loaded) so that broom::tidy(fit) — which imports
# generics::tidy — finds our methods. this is the standard pattern for
# package-provided broom methods without making generics a hard dep.
#' @noRd
.register_tidy_methods <- function() {
	if (requireNamespace("generics", quietly = TRUE)) {
		try(registerS3method("tidy", "ame", tidy.ame,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("tidy", "lame", tidy.lame,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("tidy", "ame_als", tidy.ame_als,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("tidy", "lame_als", tidy.lame_als,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("tidy", "boot_ame", tidy.boot_ame,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("glance", "ame", glance.ame,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("glance", "lame", glance.lame,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("glance", "ame_als", glance.ame_als,
		                     envir = asNamespace("generics")), silent = TRUE)
		try(registerS3method("glance", "lame_als", glance.lame_als,
		                     envir = asNamespace("generics")), silent = TRUE)
	}
	if (requireNamespace("broom", quietly = TRUE)) {
		try(registerS3method("tidy", "ame", tidy.ame,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("tidy", "lame", tidy.lame,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("tidy", "ame_als", tidy.ame_als,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("tidy", "lame_als", tidy.lame_als,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("tidy", "boot_ame", tidy.boot_ame,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("glance", "ame", glance.ame,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("glance", "lame", glance.lame,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("glance", "ame_als", glance.ame_als,
		                     envir = asNamespace("broom")), silent = TRUE)
		try(registerS3method("glance", "lame_als", glance.lame_als,
		                     envir = asNamespace("broom")), silent = TRUE)
	}
}

#' Long-format draws of the linear predictor for marginaleffects-style use
#'
#' Returns a long-format data frame with one row per
#' \code{(draw, i, j, period)} combination, giving the per-draw linear
#' predictor (or response-scale prediction) at each dyad and period.
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
#' @param type One of \code{"link"} (default) or \code{"response"}.
#' @param n_draws Number of posterior draws to use. Default 100.
#' @param seed Optional RNG seed.
#'
#' @return A long-format data frame with columns \code{.chain},
#'   \code{.iteration}, \code{.draw}, \code{period}, \code{period_label},
#'   \code{i}, \code{j}, \code{actor_i}, \code{actor_j}, \code{.value}.
#'   Returned as a tibble when the \pkg{tibble} package is available;
#'   otherwise a plain data frame.
#'
#' @export
prediction_draws_long <- function(object,
                                  newdata = NULL,
                                  type    = c("link", "response"),
                                  n_draws = 100L,
                                  seed    = NULL) {
	type <- match.arg(type)
	if (!is.null(seed)) set.seed(seed)
	if (!inherits(object, "lame")) {
		cli::cli_abort("Only {.cls lame} fits are supported.")
	}
	# extract dims
	if (is.null(object$EZ) || length(object$EZ) < 1L) {
		cli::cli_abort("Fit has no per-period EZ; cannot draw predictions.")
	}
	T_per <- length(object$EZ)
	n_a <- nrow(object$EZ[[1L]]); n_b <- ncol(object$EZ[[1L]])

	# pick draw indices
	n_iter <- dim(object$BETA)[1L]
	idx <- sample.int(n_iter, n_draws, replace = (n_draws > n_iter))

	beta_is_dyn <- length(dim(object$BETA)) == 3L

	# beta names + dyadic mask
	beta_names <- if (beta_is_dyn) dimnames(object$BETA)[[2]] else colnames(object$BETA)
	has_int <- !is.null(beta_names) && beta_names[1L] == "intercept"
	is_dyad <- if (is.null(beta_names))
		c(FALSE, rep(TRUE, ncol(object$BETA) - 1L))
		else grepl("_dyad$|\\.dyad$", beta_names)
	n_dyad <- sum(is_dyad)

	# x source
	X_src <- if (!is.null(newdata)) newdata else object$Xlist
	if (is.null(X_src) || length(X_src) < T_per) {
		# fall back to repeating the last period's covariates
		X_src <- replicate(T_per, X_src[[length(X_src)]], simplify = FALSE)
	}

	# resolve actor and period labels from the fit's y dimnames so the long
	# table carries the same ids the user supplied at fit time. falls back
	# to integer indices when dimnames are missing.
	actor_i_names <- NULL
	actor_j_names <- NULL
	period_names  <- NULL
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
	if (is.null(period_names))  period_names  <- paste0("t", seq_len(T_per))

	# per-time per-draw eta
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
			Xt <- X_src[[t_idx]]
			if (is.null(dim(Xt))) Xt <- array(Xt, c(n_a, n_b, 1L))
			if (length(dim(Xt)) == 2L) Xt <- array(Xt, c(dim(Xt), 1L))
			eta <- matrix(0, n_a, n_b)
			dyad_beta_t <- bt[is_dyad]
			for (k in seq_len(min(n_dyad, dim(Xt)[3L]))) {
				eta <- eta + dyad_beta_t[k] * Xt[, , k]
			}
			if (has_int) eta <- eta + bt[1L]
			# additive effects (per-time when dynamic_ab; else broadcast)
			a_t <- if (!is.null(object$a_dynamic)) object$a_dynamic[, t_idx] else object$APM
			b_t <- if (!is.null(object$b_dynamic)) object$b_dynamic[, t_idx] else object$BPM
			if (!is.null(a_t)) eta <- eta + outer(a_t, rep(0, n_b), "+")
			if (!is.null(b_t)) eta <- sweep(eta, 2L, b_t, "+")
			# uv
			if (!is.null(object$U) && !is.null(object$V) &&
			    NCOL(object$U) > 0L && NCOL(object$V) > 0L) {
				if (length(dim(object$U)) == 3L) {
					Ut <- object$U[, , t_idx]
					Vt <- object$V[, , t_idx]
				} else {
					Ut <- object$U; Vt <- object$V
				}
				eta <- eta + Ut %*% t(Vt)
			}
			if (type == "response") {
				eta <- transform_to_response(eta, object$family %||% "normal")
			}
			# stack into long. column names follow tidybayes / marginaleffects
			# convention (.chain / .iteration / .draw / .value) so downstream
			# packages auto-dispatch. actor_i / actor_j / period_label
			# preserve the ids the user supplied at fit time.
			grid <- expand.grid(i = seq_len(n_a), j = seq_len(n_b),
			                    KEEP.OUT.ATTRS = FALSE)
			rows[[rix]] <- data.frame(
				.chain        = 1L,
				.iteration    = s_idx,
				.draw         = s_idx,
				period        = t_idx,
				period_label  = period_names[t_idx],
				i             = grid$i,
				j             = grid$j,
				actor_i       = actor_i_names[grid$i],
				actor_j       = actor_j_names[grid$j],
				.value        = as.numeric(eta),
				stringsAsFactors = FALSE
			)
			rix <- rix + 1L
		}
	}

	out <- do.call(rbind, rows)
	# return as a tibble when tibble is available so dplyr / tidyverse
	# verbs feel natural; otherwise fall back to a plain data.frame
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}

# --------------------------------------------------------------------------
# glance: model-level one-row summary for broom / modelsummary integration
# --------------------------------------------------------------------------

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
#'   \item \code{nobs} — number of observed dyads (NA cells excluded).
#'   \item \code{n_actors} — number of distinct actors (for bipartite,
#'     row + column actors).
#'   \item \code{n_periods} — number of time periods (1 for \code{ame}).
#'   \item \code{n_stored} — number of stored MCMC draws.
#'   \item \code{family} — outcome family ({"normal","binary",...}).
#'   \item \code{mode} — {"unipartite","bipartite"}.
#'   \item \code{R} — latent-space dimension (or max of \code{R_row},
#'     \code{R_col} for bipartite).
#'   \item \code{dynamic_uv}, \code{dynamic_ab}, \code{dynamic_beta} —
#'     logicals; whether each component is time-varying.
#'   \item \code{elpd_loo} — leave-one-out expected log predictive
#'     density, if a cached \code{loo} object is on the fit (NA if not
#'     fit with \code{save_log_lik = TRUE}).
#' }
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             nscan = 100, burn = 20, odens = 5, verbose = FALSE)
#' glance(fit)
#' }
#' @export
glance.ame <- function(x, ...) {
	Y <- x$Y %||% x$YPM
	nobs <- if (is.list(Y)) {
		sum(vapply(Y, function(yt) sum(is.finite(yt)), integer(1L)))
	} else if (is.array(Y)) {
		sum(is.finite(Y))
	} else NA_integer_

	is_bip <- identical(x$mode, "bipartite")
	nA_bip <- if (is_bip) length(x$APM %||% x$row_names %||% integer(0)) else NA_integer_
	nB_bip <- if (is_bip) length(x$BPM %||% x$col_names %||% integer(0)) else NA_integer_
	# n_actors for bipartite is the row count (the "primary" actor side);
	# the column actor count is surfaced separately as n_col_actors. a user
	# reading n_actors on a 12x9 bipartite fit gets 12 (not 21 from
	# `na + nb`), which matches the dyadic-network convention and what
	# `nrow(fit$y[[1]])` would tell them.
	n_act <- if (is_bip) {
		as.integer(nA_bip)
	} else {
		length(x$APM %||% x$BPM %||% integer(0))
	}

	n_periods <- if (is.list(Y)) length(Y) else
		if (length(dim(Y)) == 3L) dim(Y)[3] else 1L
	n_stored <- if (is.matrix(x$BETA)) nrow(x$BETA) else
		if (length(dim(x$BETA)) == 3L) dim(x$BETA)[1L] else NA_integer_

	# the elpd_loo is only present when the user explicitly cached it on the
	# fit object via loo::loo(); leave na otherwise so modelsummary does not
	# imply it was computed.
	elpd <- NA_real_
	if (!is.null(x$loo) && inherits(x$loo, "loo")) {
		est <- x$loo$estimates
		if (is.matrix(est) && "elpd_loo" %in% rownames(est)) {
			elpd <- as.numeric(est["elpd_loo", "Estimate"])
		}
	}

	R_used <- if (is_bip) {
		max(c(x$R_row %||% 0L, x$R_col %||% 0L, x$R %||% 0L))
	} else {
		x$R %||% 0L
	}

	# reconstruct the dynamic flags. lame() fits preserve $dynamic_uv /
	# $dynamic_ab / $dynamic_beta as logicals on the fit object, so prefer
	# those when present (most reliable). otherwise (e.g. ame() cross-
	# sectional fits, or a fit produced before these fields were stored)
	# fall back to shape inference: 3-d u means time-varying latent
	# positions; a 2-d $a_dynamic or $b_dynamic matrix (or 3-d $aps /
	# $bps, in older fit objects) means time-varying additive effects;
	# 3-d beta means dynamic_beta is active.
	dyn_uv_flag <- if (!is.null(x$dynamic_uv)) isTRUE(x$dynamic_uv) else
		(!is.null(x$U) && length(dim(x$U)) == 3L)
	dyn_ab_flag <- if (!is.null(x$dynamic_ab)) isTRUE(x$dynamic_ab) else
		(!is.null(x$a_dynamic) || !is.null(x$b_dynamic) ||
		 (!is.null(x$APS) && length(dim(x$APS)) == 3L) ||
		 (!is.null(x$BPS) && length(dim(x$BPS)) == 3L))
	dyn_beta_flag <- if (!is.null(x$dynamic_beta)) isTRUE(x$dynamic_beta) else
		(!is.null(x$BETA) && length(dim(x$BETA)) == 3L)

	out <- data.frame(
		nobs           = as.integer(nobs),
		n_actors       = as.integer(n_act),
		n_row_actors   = as.integer(nA_bip),
		n_col_actors   = as.integer(nB_bip),
		n_periods      = as.integer(n_periods),
		n_stored       = as.integer(n_stored),
		family         = x$family %||% NA_character_,
		mode           = if (is_bip) "bipartite" else "unipartite",
		R              = as.integer(R_used),
		dynamic_uv     = dyn_uv_flag,
		dynamic_ab     = dyn_ab_flag,
		dynamic_beta   = dyn_beta_flag,
		elpd_loo       = elpd,
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

# fallback generic so `glance(fit)` works without broom/generics loaded
#' S3 generic for `glance`
#'
#' Light-weight fallback so \code{glance(fit)} dispatches through S3
#' even when \pkg{broom} or \pkg{generics} is not loaded.
#'
#' @param x An object to glance at.
#' @param ... Passed to the method.
#' @return A one-row data frame.
#' @export
glance <- function(x, ...) UseMethod("glance")


#' Tidy method for fitted `ame_als` / `lame_als` objects
#'
#' Returns a data frame with one row per regression coefficient,
#' compatible with the \pkg{broom} idiom, so that ALS fits compose with
#' \pkg{modelsummary} / \pkg{kableExtra} pipelines next to MCMC fits.
#' Standard errors come from the sandwich covariance
#' (\code{\link{vcov.ame_als}}) by default, or from the bootstrap object
#' attached to \code{x$bootstrap} when present (preferred, fully
#' propagated). \code{statistic} is \code{estimate / std.error};
#' \code{p.value} is the Normal-approximation two-sided tail
#' \eqn{2(1 - \Phi(|z|))} from the bootstrap or sandwich standard error. It is
#' a Wald-style summary for the point estimator, not a posterior probability.
#'
#' Only the intercept and dyadic-covariate coefficients are returned,
#' matching \code{coef(fit)} on the sandwich-covered subset. Additive
#' (\code{a}, \code{b}), multiplicative (\code{U}, \code{V}), and
#' node-covariate parameters are not included; use
#' \code{\link{ame_als_bootstrap}} and inspect the bootstrap object
#' directly if you need them.
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
#'   recording \code{"bootstrap"} or \code{"sandwich"}.
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
#' fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
#'                family = "binary", R = 1, verbose = FALSE)
#' tidy(fit)
#' }
#' @export
tidy.ame_als <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	if (length(conf.level) != 1L || !is.finite(conf.level) ||
	    conf.level <= 0 || conf.level >= 1) {
		cli::cli_abort("{.arg conf.level} must be a single number in (0, 1).")
	}
	# prefer bootstrap-based summaries when available; fall back to the
	# sandwich vcov for a fast point-and-ci table when not. the tidy
	# contract is the intercept and dyadic-covariate coefficients only;
	# the a / b / u / v / vc bootstrap intervals are reachable via
	# `confint(fit$bootstrap)` directly for users who want them.
	if (!is.null(x$bootstrap)) {
		# restrict to the regression-coefficient subset: same set
		# `vcov.ame_als()` and `coef(x)` use. the bootstrap object stores
		# replicate coefficients in `$coefs`.
		coef_nms <- names(x$coefficients)
		boot_coefs <- x$bootstrap$coefs
		# constrain to the rows that actually have a regression coefficient,
		# excluding any a / b / u / v / vc entries that the bootstrap object
		# also stores
		nms <- intersect(coef_nms, colnames(boot_coefs))
		if (length(nms) == 0L) {
			# fall back to all regression coefficient names if column
			# alignment is wrong
			nms <- coef_nms
		}
		est <- x$coefficients[nms]
		# point estimates can shift slightly between fits; use the original
		# fit's coefficients, not the bootstrap mean
		if (!is.null(boot_coefs) && all(nms %in% colnames(boot_coefs))) {
			se <- apply(boot_coefs[, nms, drop = FALSE], 2L,
			            stats::sd, na.rm = TRUE)
			a <- (1 - conf.level) / 2
			ci_lo <- apply(boot_coefs[, nms, drop = FALSE], 2L,
			               stats::quantile, probs = a, na.rm = TRUE)
			ci_hi <- apply(boot_coefs[, nms, drop = FALSE], 2L,
			               stats::quantile, probs = 1 - a, na.rm = TRUE)
			ci <- cbind(ci_lo, ci_hi)
			rownames(ci) <- nms
		} else {
			# fall back path: derive se from confint() half-width
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
		nms <- rownames(V)
		est <- x$coefficients[nms]
		se  <- sqrt(diag(V))
		a   <- (1 - conf.level) / 2
		zq  <- stats::qnorm(c(a, 1 - a))
		ci  <- cbind(est + zq[1] * se, est + zq[2] * se)
		rownames(ci) <- nms
		se_source <- "sandwich"
	}
	stat <- est / pmax(se, .Machine$double.eps)
	out <- data.frame(
		term       = nms,
		estimate   = unname(est),
		std.error  = unname(se),
		statistic  = unname(stat),
		p.value    = 2 * (1 - stats::pnorm(abs(unname(stat)))),
		conf.low   = unname(ci[, 1L]),
		conf.high  = unname(ci[, 2L]),
		se_source  = se_source,
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
#' data(YX_bin_list)
#' Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
#' fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
#'                family = "binary", R = 1, verbose = FALSE)
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
	# ame_als stores per-actor effects under different field names than ame()
	a_len <- length(x$a %||% x$a_row %||% integer(0))
	b_len <- length(x$b %||% x$b_col %||% integer(0))
	# n_actors: row-side actors for bipartite (matches glance.ame), separate
	# n_row_actors / n_col_actors columns below for the column side.
	n_act <- if (is_bip) as.integer(a_len) else as.integer(max(a_len, b_len))
	nA_bip <- if (is_bip) as.integer(a_len) else NA_integer_
	nB_bip <- if (is_bip) as.integer(b_len) else NA_integer_

	n_periods <- if (isTRUE(x$longitudinal)) {
		if (is.list(Y)) length(Y) else
		if (length(dim(Y)) == 3L) dim(Y)[3] else 1L
	} else 1L

	se_src <- if (!is.null(x$bootstrap)) "bootstrap" else "sandwich"
	# non_normal_method distinguishes the irls path from the rank-normal
	# transform. als supports normal / binary / poisson; tobit / ordinal
	# inference goes through the mcmc path.
	nnm <- x$non_normal_method %||% NA_character_
	out <- data.frame(
		nobs            = as.integer(nobs),
		n_actors        = n_act,
		n_row_actors    = nA_bip,
		n_col_actors    = nB_bip,
		n_periods       = as.integer(n_periods),
		family          = x$family %||% NA_character_,
		mode            = if (is_bip) "bipartite" else "unipartite",
		R               = as.integer(x$R %||% 0L),
		non_normal_method = nnm,
		converged       = isTRUE(x$converged),
		iterations      = as.integer(x$iterations %||% NA_integer_),
		se_source       = se_src,
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
	se  <- x$se
	a   <- (1 - conf.level) / 2
	# default ci_lo / ci_hi on the boot_ame object are at 95%; recompute
	# from the replicate matrix when the user asks for a different level
	if (abs(conf.level - 0.95) < 1e-8 && !is.null(x$ci_lo)) {
		ci_lo <- x$ci_lo; ci_hi <- x$ci_hi
	} else {
		ci_lo <- apply(x$coefs, 2L, stats::quantile, probs = a, na.rm = TRUE)
		ci_hi <- apply(x$coefs, 2L, stats::quantile, probs = 1 - a, na.rm = TRUE)
	}
	stat <- est / pmax(se, .Machine$double.eps)
	out <- data.frame(
		term       = nms,
		estimate   = unname(est),
		std.error  = unname(se),
		statistic  = unname(stat),
		p.value    = 2 * (1 - stats::pnorm(abs(unname(stat)))),
		conf.low   = unname(ci_lo),
		conf.high  = unname(ci_hi),
		se_source  = "bootstrap",
		stringsAsFactors = FALSE
	)
	if (requireNamespace("tibble", quietly = TRUE)) {
		out <- tibble::as_tibble(out)
	}
	out
}
