# ame_als_methods_extra.R
#
# additional S3 methods for ame_als so simulate / gof_plot / ab_plot /
# prior_summary / nobs work uniformly on MCMC and ALS fits.

####
#' Number of observed dyads in an ame_als fit
#'
#' Mirrors \code{\link{nobs.ame}} so the MCMC and ALS fits expose the same
#' accessor.
#'
#' @param object an \code{ame_als} fit.
#' @param ... ignored.
#' @return Integer: number of observed dyads.
#' @method nobs ame_als
#' @export
nobs.ame_als <- function(object, ...) {
	Y <- object$Y
	if (is.null(Y)) cli::cli_abort("Cannot compute {.fn nobs}: {.field Y} not stored on the fit.")
	if (is.list(Y)) sum(vapply(Y, function(y) sum(is.finite(y)), integer(1)))
	else sum(is.finite(Y))
}

####
#' Simulate networks from a fitted ame_als model
#'
#' Posterior-predictive equivalent for the fast ALS estimator: draws \code{nsim}
#' replicates of \code{Y} from the parametric model implied by the fitted
#' \code{mu}, \code{beta}, \code{a}, \code{b}, \code{U}, \code{V}, and the
#' family-appropriate noise distribution. Returned object has class
#' \code{"ame.sim"} so it is compatible with \code{plot_ppc_*}-style consumers
#' (where applicable).
#'
#' Caveat: ALS is a point estimator. \code{simulate.ame_als} therefore holds
#' \code{mu, beta, a, b, U, V} FIXED at the point estimate and only the
#' noise is resampled. For uncertainty over the parameters themselves use
#' \code{\link{ame_als_bootstrap}} (whose replicates each carry their own
#' resampled \code{Y}) and combine those.
#'
#' @param object an \code{ame_als} fit.
#' @param nsim integer; number of replicates to simulate (default 1).
#' @param seed optional RNG seed.
#' @param ... ignored.
#' @return An object of class \code{"ame.sim"} with element \code{Y} -- a list
#'   of \code{nsim} simulated outcome arrays (one matrix for a cross-section,
#'   one list-of-matrices per slice for a longitudinal fit).
#' @method simulate ame_als
#' @export
simulate.ame_als <- function(object, nsim = 1, seed = NULL, ...) {
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			.old_seed <- get(".Random.seed", envir = globalenv())
			on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
			        add = TRUE)
		}
		set.seed(seed)
	}
	nr <- object$dims$n_row; nc <- object$dims$n_col
	Tt <- object$n_time
	family <- object$family
	# linear-predictor base (mu + Xbeta + outer(a, b) + UV') per time slice
	EZ_list <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		EZ_t <- object$EZ
		if (is.list(EZ_t)) EZ_t <- EZ_t[[t]] else EZ_t <- EZ_t
		EZ_list[[t]] <- EZ_t
	}
	# residual SD for normal / tobit
	sigma_e <- sqrt(max(object$VC["ve"], 0, na.rm = TRUE))
	sim_list <- vector("list", nsim)
	for (s in seq_len(nsim)) {
		Yt_list <- vector("list", Tt)
		for (t in seq_len(Tt)) {
			EZ_t <- EZ_list[[t]]
			Yt <- switch(family,
				normal  = EZ_t + matrix(stats::rnorm(nr * nc, 0, sigma_e), nr, nc),
				binary  = 1 * (EZ_t + matrix(stats::rnorm(nr * nc), nr, nc) > 0),
				ordinal = EZ_t + matrix(stats::rnorm(nr * nc), nr, nc),
				poisson = matrix(stats::rpois(nr * nc, lambda = pmax(exp(EZ_t), 0)),
				                  nr, nc),
				EZ_t)
			if (!identical(object$mode, "bipartite") && nr == nc) diag(Yt) <- NA_real_
			# symmetric fit -> draws must also be symmetric (noise broke it)
			if (isTRUE(object$symmetric) && nr == nc) {
				Yt[lower.tri(Yt)] <- t(Yt)[lower.tri(Yt)]
			}
			dimnames(Yt) <- list(object$row_names, object$col_names)
			Yt_list[[t]] <- Yt
		}
		# cross-sectional fit (Tt == 1): return a single matrix per replicate
		sim_list[[s]] <- if (Tt == 1L) Yt_list[[1]] else Yt_list
	}
	out <- list(Y = sim_list, n = nsim, family = family,
	            mode = object$mode %||% "unipartite",
	            symmetric = isTRUE(object$symmetric),
	            n_time = Tt,
	            note = "Conditional on point-estimate parameters; resampling only the noise. For full uncertainty use ame_als_bootstrap().")
	class(out) <- c("ame_als.sim", "ame.sim")
	out
}

####
#' Goodness-of-fit posterior predictive check for an ame_als fit
#'
#' Bootstrap-based analogue of the MCMC \code{\link{gof_plot}}: draws
#' \code{nsim} simulated networks from the fitted ALS model, computes the
#' standard network statistics (\code{sd.rowmean}, \code{sd.colmean},
#' \code{dyad.dep}, \code{cycle.dep}, \code{trans.dep} for unipartite;
#' \code{sd.rowmean}, \code{sd.colmean}, \code{four.cycles} for bipartite) on
#' each replicate, and overlays the observed value on a histogram of replicate
#' values.
#'
#' Uses \code{\link{simulate.ame_als}} for the replicates and is therefore
#' subject to the same caveat: the noise is resampled, but the point estimates
#' of \code{mu, beta, a, b, U, V} are held fixed. The MCMC \code{gof_plot.ame}
#' integrates over the posterior of those parameters too, so the ALS GOF is a
#' weaker check.
#'
#' @param fit an \code{ame_als} fit.
#' @param nsim integer; number of replicates (default 100).
#' @param seed optional RNG seed.
#' @param ... reserved.
#' @return A \code{ggplot} (or \code{patchwork}) object.
#' @export
gof_plot.ame_als <- function(fit, nsim = 100, seed = NULL, ...) {
	if (!inherits(fit, "ame_als")) {
		cli::cli_abort("{.arg fit} must be an {.cls ame_als} object.")
	}
	mode_arg <- fit$mode %||% "unipartite"
	sims <- simulate(fit, nsim = nsim, seed = seed)
	# observed-statistic vector
	Y_obs <- if (is.list(fit$Y)) fit$Y[[1]] else fit$Y
	if (length(dim(Y_obs)) == 3L) Y_obs <- Y_obs[, , 1]
	obs <- gof_stats(Y_obs, mode = mode_arg)
	# per-replicate stats
	rep_stats <- vapply(sims$Y, function(Yt) {
		Y1 <- if (is.list(Yt)) Yt[[1]] else Yt
		if (length(dim(Y1)) == 3L) Y1 <- Y1[, , 1]
		gof_stats(Y1, mode = mode_arg)
	}, numeric(length(obs)))
	if (is.null(dim(rep_stats))) rep_stats <- matrix(rep_stats, nrow = 1)
	rownames(rep_stats) <- names(obs)
	df <- do.call(rbind, lapply(seq_along(obs), function(j) {
		data.frame(
			statistic = names(obs)[j],
			observed  = unname(obs[j]),
			value     = unname(rep_stats[j, ]),
			stringsAsFactors = FALSE)
	}))
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg ggplot2} is required for {.fn gof_plot.ame_als}.")
	}
	# dual-encode observed cue (Okabe-Ito #D55E00 + dashed linetype) for
	# consistency with the MCMC gof_plot; the colour-blind-safe vermillion
	# also reads as "orange" so legend/prose match.
	ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
		ggplot2::geom_histogram(ggplot2::aes(fill = "Bootstrap replicates"),
		                        bins = 30) +
		ggplot2::geom_vline(ggplot2::aes(xintercept = .data$observed,
		                                 colour = "Observed",
		                                 linetype = "Observed"),
		                    linewidth = 1) +
		ggplot2::scale_fill_manual(NULL,
		                           values = c("Bootstrap replicates" = "grey60")) +
		ggplot2::scale_colour_manual(NULL,
		                             values = c("Observed" = "#D55E00")) +
		ggplot2::scale_linetype_manual(NULL,
		                               values = c("Observed" = "dashed")) +
		ggplot2::facet_wrap(~ .data$statistic, scales = "free") +
		ggplot2::labs(title = "ame_als Goodness-of-Fit (Bootstrap-Based Posterior Predictive)",
		              subtitle = paste0("R = ", nsim,
		                                " parametric draws; dashed orange line = observed"),
		              x = NULL, y = "Count") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border     = ggplot2::element_blank(),
			axis.ticks       = ggplot2::element_blank(),
			legend.position  = "top",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x     = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y     = ggplot2::element_text(color = "white", hjust = 0)
		)
}

####
#' Additive-effects plot for an ame_als fit
#'
#' Sender (a) and receiver (b) additive-effect point estimates as a sorted
#' lollipop chart. Mirrors \code{\link{ab_plot}} but without posterior
#' intervals: this is a point estimator. If the fit was produced via
#' \code{ame_als(..., bootstrap = N)}, bootstrap 95\% intervals are drawn as
#' error bars.
#'
#' @param fit an \code{ame_als} fit.
#' @param effect \code{"sender"} (default), \code{"receiver"}, or \code{"both"}.
#' @param top_n integer; show only the top / bottom \code{top_n} actors by
#'   absolute effect (default \code{Inf}, show all).
#' @param ... reserved.
#' @return A \code{ggplot} object.
#' @export
ab_plot.ame_als <- function(fit, effect = c("sender", "receiver", "both"),
                             top_n = Inf, ...) {
	if (!inherits(fit, "ame_als")) {
		cli::cli_abort("{.arg fit} must be an {.cls ame_als} object.")
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg ggplot2} is required for {.fn ab_plot.ame_als}.")
	}
	effect <- match.arg(effect)
	a <- fit$a; b <- fit$b
	bt <- fit$bootstrap
	mk <- function(vec, lbl) {
		if (is.null(vec)) return(NULL)
		nm <- names(vec) %||% paste0("a", seq_along(vec))
		df <- data.frame(actor = nm, value = unname(vec),
		                 effect = lbl, stringsAsFactors = FALSE)
		# bootstrap intervals when available
		if (!is.null(bt) && !is.null(bt[[paste0(lbl, "_q")]])) {
			q <- bt[[paste0(lbl, "_q")]]
			df$lo <- q[, 1]; df$hi <- q[, 2]
		} else {
			df$lo <- NA_real_; df$hi <- NA_real_
		}
		# top-n filtering
		if (is.finite(top_n) && nrow(df) > top_n) {
			ord <- order(abs(df$value), decreasing = TRUE)
			df <- df[ord[seq_len(top_n)], , drop = FALSE]
		}
		df$actor <- factor(df$actor, levels = df$actor[order(df$value)])
		df
	}
	dfs <- list(); if (effect %in% c("sender", "both")) dfs$a <- mk(a, "sender")
	if (effect %in% c("receiver", "both")) dfs$b <- mk(b, "receiver")
	df_all <- do.call(rbind, dfs)
	if (is.null(df_all) || nrow(df_all) == 0L) {
		cli::cli_abort("No additive effects to plot for {.arg effect} = {.val {effect}}.")
	}
	p <- ggplot2::ggplot(df_all,
	                    ggplot2::aes(x = .data$value, y = .data$actor)) +
		ggplot2::geom_point()
	if (!all(is.na(df_all$lo))) {
		p <- p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo,
		                                              xmax = .data$hi),
		                                  height = 0.2, alpha = 0.6)
	}
	p + ggplot2::geom_vline(xintercept = 0, linetype = 3) +
		ggplot2::facet_wrap(~ effect, scales = "free_y") +
		ggplot2::labs(title = "ame_als Additive Effects",
		              subtitle = if (!is.null(bt)) "With Bootstrap 95% CIs" else
		                          "Point Estimate (no CI; refit with bootstrap = N for intervals)",
		              x = "Effect (Linear-Predictor Scale)", y = NULL) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border     = ggplot2::element_blank(),
			axis.ticks       = ggplot2::element_blank(),
			legend.position  = "top",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x     = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y     = ggplot2::element_text(color = "white", hjust = 0)
		)
}

####
#' Print the priors used by an AME / LAME / ame_als fit
#'
#' Bayesian-hygiene helper modelled on \code{rstanarm::prior_summary()}.
#' Prints the priors that were actually in effect for a fit, with defaults
#' filled in. For \code{ame_als} fits it just states that the estimator is a
#' frequentist point estimator and there are no priors.
#'
#' @param object a fitted \code{ame}, \code{lame}, or \code{ame_als} object.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
prior_summary <- function(object, ...) UseMethod("prior_summary")

#' @rdname prior_summary
#' @method prior_summary default
#' @export
prior_summary.default <- function(object, ...) {
	cli::cli_abort("No {.fn prior_summary} method for an object of class {.cls {class(object)[1]}}.")
}

#' @rdname prior_summary
#' @method prior_summary ame
#' @export
prior_summary.ame <- function(object, ...) {
	fit_kind <- if (inherits(object, "lame")) "lame" else "ame"
	cli::cli_h2("Priors in effect ({fit_kind} fit)")
	cli::cli_text("Regression coefficients: {.code beta ~ N(0, g * sigma^2 * (X'X)^-1)} (g-prior).")
	g_used <- object$g
	if (is.null(g_used) || all(is.na(g_used))) {
		cli::cli_text("  {.code g} (top-level argument) = {.emph not exposed on this fit object} (sampler-filled default; refit with a newer build to expose the numeric value)")
	} else if (length(g_used) == 1L) {
		cli::cli_text("  {.code g} (top-level argument) = {.val {round(g_used, 4)}}")
	} else {
		cli::cli_text("  {.code g} (vector, top-level argument) = {.val {paste(round(g_used, 4), collapse=', ')}}")
	}
	cli::cli_text("Note: {.code g} is set at top level on {.fn {fit_kind}} -- not inside {.arg prior = list(...)}.")
	pr <- object$prior %||% list()
	# print every documented prior key whose value is on the fit
	for (nm in c("Sab0", "eta0", "etaab", "Suv0", "kappa0", "s20", "s2u0")) {
		v <- pr[[nm]]
		if (!is.null(v)) {
			val_str <- if (is.matrix(v)) {
				paste0("matrix(", paste(round(as.vector(v), 4), collapse = ", "), ")")
			} else {
				paste(round(unlist(v), 4), collapse = ", ")
			}
			cli::cli_text("  {.code {nm}} = {.val {val_str}}")
		}
	}
	if (length(pr) == 0L) {
		cli::cli_text(cli::col_grey("Other prior hyperparameters: sampler-filled defaults (see ?ame Prior Distributions)."))
	}
	invisible(object)
}

#' @rdname prior_summary
#' @method prior_summary lame
#' @export
prior_summary.lame <- function(object, ...) {
	prior_summary.ame(object, ...)
	pr <- object$prior %||% list()
	if (isTRUE(object$dynamic_uv) || isTRUE(object$dynamic_ab) ||
	    isTRUE(object$dynamic_beta)) {
		cli::cli_h3("Dynamic AR(1) priors")
		for (nm in c("rho_uv_mean", "rho_uv_sd",
		             "rho_ab_mean", "rho_ab_sd",
		             "sigma_uv_shape", "sigma_uv_scale",
		             "sigma_ab_shape", "sigma_ab_scale")) {
			v <- pr[[nm]]
			if (!is.null(v)) cli::cli_text("  {.code {nm}} = {.val {round(v, 4)}}")
		}
		if (isTRUE(object$dynamic_beta)) {
			for (nm in c("rho_beta_mean", "rho_beta_sd",
			             "sigma_beta_shape", "sigma_beta_scale",
			             "sigma_beta_init",
			             "rho_beta_lower", "rho_beta_upper",
			             "beta0_mean", "beta0_var")) {
				v <- pr[[nm]]
				if (!is.null(v)) cli::cli_text("  {.code {nm}} = {.val {round(v, 4)}}")
			}
		}
	}
	invisible(object)
}

#' @rdname prior_summary
#' @method prior_summary ame_als
#' @export
prior_summary.ame_als <- function(object, ...) {
	cli::cli_h2("Priors (ame_als fit)")
	cli::cli_text("The ALS estimator is a frequentist iterative block coordinate descent point estimator. There are no priors.")
	cli::cli_text("For uncertainty, refit with {.code bootstrap = N} (parametric / block bootstrap).")
	cli::cli_text("For Bayesian inference (priors, credible intervals, posterior summaries), use {.fn ame} / {.fn lame} (MCMC).")
	invisible(object)
}
