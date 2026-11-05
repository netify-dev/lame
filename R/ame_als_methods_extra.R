# extra S3 methods for ame_als so simulate, gof_plot, ab_plot, prior_summary and
# nobs work the same way on MCMC and ALS fits

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

#' Simulate networks from a fitted ame_als model
#'
#' Draws \code{nsim} replicates of \code{Y} from the fitted ALS mean structure:
#' \code{mu}, \code{beta}, \code{a}, \code{b}, \code{U}, \code{V}, and the
#' family-appropriate noise distribution. For directed unipartite normal,
#' binary and Poisson fits, reciprocal dyads use the fitted \code{rho}. Binary IRLS fits
#' use their stated marginal link and a Gaussian copula with the fitted
#' dyadic correlation. Square bipartite networks retain their diagonal cells.
#' The returned
#' object has class
#' \code{"ame.sim"} so it is compatible with \code{plot_ppc_*}-style consumers
#' (where applicable).
#'
#' ALS is a point estimator. \code{simulate.ame_als} holds
#' \code{mu, beta, a, b, U, V} at the point estimate and resamples the noise.
#' For uncertainty over the parameters themselves, use
#' \code{\link{ame_als_bootstrap}} (whose replicates each carry their own
#' resampled \code{Y}) and combine those simulations.
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
		.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
		.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
		on.exit({
			if (.had_seed) {
				assign(".Random.seed", .old_seed, envir = globalenv())
			} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
				rm(".Random.seed", envir = globalenv())
			}
		}, add = TRUE)
		set.seed(seed)
	}
	nr <- object$dims$n_row; nc <- object$dims$n_col
	Tt <- object$n_time
	family <- object$family
	rho <- as.numeric(object$VC["rho"] %||% 0)
	if (!is.finite(rho) || isTRUE(object$symmetric) ||
			identical(object$mode, "bipartite") || nr != nc) {
		rho <- 0
	}
	if (.ae_rho_required(object)) {
		cli::cli_abort(
			c(
				"Dyadic dependence is not estimated from the available reciprocal pairs.",
				"i" = "Use {.fn ame} or {.fn lame} (MCMC) to integrate uncertainty in dyadic dependence."
			),
			class = "lame_als_prediction"
		)
	}
	rho <- max(min(rho, 1), -1)
	draw_error <- function(sd = 1) {
		E <- matrix(stats::rnorm(nr * nc, 0, sd), nr, nc)
		if (abs(rho) > 1e-12) {
			for (i in seq_len(nr - 1L)) {
				for (j in (i + 1L):nr) {
					z1 <- stats::rnorm(1)
					z2 <- stats::rnorm(1)
					E[i, j] <- sd * z1
					E[j, i] <- sd * (rho * z1 + sqrt(1 - rho^2) * z2)
				}
			}
		}
		if (nr == nc && !identical(object$mode, "bipartite")) diag(E) <- NA_real_
		E
	}
	# linear-predictor base per time slice
	EZ_list <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		EZ_t <- object$EZ
		if (is.list(EZ_t)) EZ_t <- EZ_t[[t]] else EZ_t <- EZ_t
		EZ_list[[t]] <- EZ_t
	}
	# residual sd for the Gaussian response
	sigma_e <- sqrt(max(object$VC["ve"], 0, na.rm = TRUE))
	# an IRLS Poisson fit reports EZ on the latent scale of y ~ Pois(exp(z)),
	# z ~ N(EZ, ve); the transform path has no latent layer
	pois_latent_sd <- if (identical(family, "poisson") &&
			identical(object$link, "log")) sigma_e else 0
	sim_list <- vector("list", nsim)
	for (s in seq_len(nsim)) {
		Yt_list <- vector("list", Tt)
		for (t in seq_len(Tt)) {
			EZ_t <- EZ_list[[t]]
			Yt <- switch(
				family,
				normal = EZ_t + draw_error(sigma_e),
				binary = {
					threshold <- if (identical(object$link, "logit")) {
						-stats::qnorm(stats::plogis(EZ_t, log.p = TRUE), log.p = TRUE)
					} else -EZ_t
					1 * (draw_error(1) > threshold)
				},
				ordinal = EZ_t + draw_error(1),
				# the latent layer is drawn so the replicates are not equidispersed
				poisson = matrix(
					stats::rpois(
						nr * nc,
						lambda = .als_pois_rate(if (pois_latent_sd > 0) {
							EZ_t + draw_error(pois_latent_sd)
						} else EZ_t)
					),
					nr,
					nc
				),
				EZ_t
			)
			if (!identical(object$mode, "bipartite") && nr == nc) diag(Yt) <- NA_real_
			# a symmetric fit gives symmetric draws
			if (isTRUE(object$symmetric) && nr == nc) {
				Yt[lower.tri(Yt)] <- t(Yt)[lower.tri(Yt)]
			}
			dimnames(Yt) <- list(object$row_names, object$col_names)
			Yt_list[[t]] <- Yt
		}
		# a cross-sectional fit returns one matrix per replicate
		sim_list[[s]] <- if (Tt == 1L) Yt_list[[1]] else Yt_list
	}
	out <- list(
		Y = sim_list,
		n = nsim,
		family = family,
		mode = object$mode %||% "unipartite",
		symmetric = isTRUE(object$symmetric),
		n_time = Tt,
		note = "Conditional on point-estimate parameters; resampling only the noise. For full uncertainty use ame_als_bootstrap()."
	)
	class(out) <- c("ame_als.sim", "ame.sim")
	out
}

#' Goodness-of-fit check for an ame_als fit
#'
#' Bootstrap-style analogue of the MCMC \code{\link{gof_plot}}: draws
#' \code{nsim} simulated networks from the fitted ALS model, computes the
#' standard network statistics (\code{sd.rowmean}, \code{sd.colmean},
#' \code{dyad.dep}, \code{cycle.dep}, \code{trans.dep} for unipartite;
#' \code{sd.rowmean}, \code{sd.colmean}, \code{four.cycles} for bipartite,
#' replacing the count with weighted \code{four.cycle.dep} for Gaussian
#' outcomes) on
#' each replicate, and compares the observed value with the replicates.
#'
#' A cross-sectional fit gives one histogram per statistic with the observed
#' value overlaid. A longitudinal fit gives, as \code{\link{gof_plot}} does for
#' an MCMC \code{lame} fit, the statistics of \emph{every} time slice: for each
#' statistic, the observed value in each period against the median and the
#' central 95 percent of the simulated values in that period.
#'
#' Uses \code{\link{simulate.ame_als}} for the replicates and is therefore
#' conditional on the fitted point estimates: the noise is resampled, but
#' \code{mu, beta, a, b, U, V} are held fixed. The MCMC \code{gof_plot.ame}
#' also averages over posterior movement in those parameters.
#'
#' @param fit an \code{ame_als} fit.
#' @param nsim integer; number of replicates (default 100).
#' @param seed optional RNG seed.
#' @param ... reserved.
#' @return A \code{ggplot} (or \code{patchwork}) object.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' gof_plot(fit, nsim = 20, seed = 1)
#' @export
gof_plot.ame_als <- function(fit, nsim = 100, seed = NULL, ...) {
	if (!inherits(fit, "ame_als")) {
		cli::cli_abort("{.arg fit} must be an {.cls ame_als} object.")
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg ggplot2} is required for {.fn gof_plot.ame_als}.")
	}
	mode_arg <- fit$mode %||% "unipartite"
	Tt <- fit$n_time %||% 1L
	sims <- simulate(fit, nsim = nsim, seed = seed)
	# slice t of an outcome stored as a matrix, a 3-D array or a list of slices
	slice_of <- function(Y, t) {
		if (is.list(Y)) Y[[t]] else if (length(dim(Y)) == 3L) Y[, , t] else Y
	}
	# statistics of every slice for the observed network and each replicate
	obs <- lapply(seq_len(Tt), function(t) {
		gof_stats(slice_of(fit$Y, t), mode = mode_arg, family = fit$family)
	})
	rep_stats <- lapply(seq_len(Tt), function(t) {
		Y_obs <- slice_of(fit$Y, t)
		m <- vapply(sims$Y, function(Yr) {
			Y1 <- slice_of(Yr, t)
			Y1[!is.finite(Y_obs)] <- NA_real_
			gof_stats(Y1, mode = mode_arg, family = fit$family)
		}, numeric(length(obs[[t]])))
		if (is.null(dim(m))) m <- matrix(m, nrow = length(obs[[t]]))
		rownames(m) <- names(obs[[t]])
		m
	})
	if (Tt > 1L) {
		return(.ae_gof_plot_slices(obs, rep_stats, nsim))
	}
	obs <- obs[[1L]]
	rep_stats <- rep_stats[[1L]]
	df <- do.call(rbind, lapply(seq_along(obs), function(j) {
		data.frame(
			statistic = names(obs)[j],
			observed = unname(obs[j]),
			value = unname(rep_stats[j, ]),
			stringsAsFactors = FALSE
		)
	}))
	# observed cue drawn with both colour and linetype
	ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
		ggplot2::geom_histogram(
			ggplot2::aes(fill = "Bootstrap replicates"),
			bins = 30
		) +
		ggplot2::geom_vline(
			ggplot2::aes(
				xintercept = .data$observed,
				colour = "Observed",
				linetype = "Observed"
			),
			linewidth = 1
		) +
		ggplot2::scale_fill_manual(
			NULL,
			values = c("Bootstrap replicates" = "grey60")
		) +
		ggplot2::scale_colour_manual(
			NULL,
			values = c("Observed" = "black")
		) +
		ggplot2::scale_linetype_manual(
			NULL,
			values = c("Observed" = "dashed")
		) +
		ggplot2::facet_wrap(
			~ .data$statistic,
			scales = "free",
			labeller = ggplot2::as_labeller(gof_label)
		) +
		ggplot2::labs(
			title = "Goodness of fit (bootstrap-based)",
			subtitle = paste0(nsim, " parametric draws; dashed line = observed"),
			x = NULL,
			y = "Count"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			plot.caption = ggplot2::element_text(hjust = 0),
			axis.ticks = ggplot2::element_blank(),
			legend.position = "top",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y = ggplot2::element_text(color = "white", hjust = 0)
		)
}

# per-period GOF of a longitudinal ame_als fit, drawn as gof_plot() draws a lame
# fit: observed value against the median and central 95 percent of the replicates
.ae_gof_plot_slices <- function(obs, rep_stats, nsim) {
	stat_names <- names(obs[[1L]])
	plot_data <- do.call(rbind, lapply(seq_along(obs), function(t) {
		do.call(rbind, lapply(stat_names, function(nm) {
			v <- rep_stats[[t]][nm, ]
			v <- v[is.finite(v)]
			o <- unname(obs[[t]][nm])
			if (!length(v) || !is.finite(o)) return(NULL)
			data.frame(
				time = t,
				statistic = nm,
				observed = o,
				median = stats::median(v),
				lower = stats::quantile(v, 0.025, names = FALSE),
				upper = stats::quantile(v, 0.975, names = FALSE),
				stringsAsFactors = FALSE
			)
		}))
	}))
	if (is.null(plot_data) || nrow(plot_data) == 0L) {
		cli::cli_abort("No finite GOF values are available to plot.")
	}
	plot_data$statistic <- gof_label(plot_data$statistic)
	ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$time)) +
		ggplot2::geom_ribbon(
			ggplot2::aes(
				ymin = .data$lower,
				ymax = .data$upper,
				fill = "Simulated 95%"
			),
			alpha = 0.3
		) +
		ggplot2::geom_line(
			ggplot2::aes(
				y = .data$median,
				colour = "Simulated median",
				linetype = "Simulated median"
			),
			linewidth = 1
		) +
		ggplot2::geom_line(
			ggplot2::aes(
				y = .data$observed,
				colour = "Observed",
				linetype = "Observed"
			),
			linewidth = 1
		) +
		ggplot2::geom_point(ggplot2::aes(y = .data$observed, colour = "Observed"), size = 2) +
		ggplot2::scale_fill_manual(NULL, values = c("Simulated 95%" = "grey60")) +
		ggplot2::scale_colour_manual(
			NULL,
			values = c(
				"Observed" = "black",
				"Simulated median" = "grey25"
			)
		) +
		ggplot2::scale_linetype_manual(
			NULL,
			values = c(
				"Observed" = "solid",
				"Simulated median" = "dashed"
			)
		) +
		ggplot2::facet_wrap(~ .data$statistic, scales = "free_y", ncol = 2) +
		ggplot2::labs(
			title = "Goodness of fit over time (bootstrap-based)",
			subtitle = paste0(nsim, " parametric draws per period; band = central 95%"),
			x = "Period",
			y = "Value"
		) +
		gof_theme() +
		ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines")) +
		ggplot2::guides(
			fill = ggplot2::guide_legend(order = 1),
			colour = ggplot2::guide_legend(order = 2),
			linetype = ggplot2::guide_legend(order = 2)
		)
}

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
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' ab_plot(fit, effect = "sender")
#' @export
ab_plot.ame_als <- function(
	fit,
	effect = c("sender", "receiver", "both"),
	top_n = Inf,
	...
) {
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
		df <- data.frame(
			actor = nm,
			value = unname(vec),
			effect = lbl,
			stringsAsFactors = FALSE
		)
		# bootstrap intervals from the per-replicate a / b draws ([n_boot x n_actor],
		# failed replicates are NA rows); none for an effect the bootstrap withholds
		coef_mat <- if (!is.null(bt) &&
				!c(sender = "a", receiver = "b")[[lbl]] %in% .ae_boot_withheld_names(bt, "effects")) {
			if (lbl == "sender") bt$a_coefs else bt$b_coefs
		} else NULL
		if (!is.null(coef_mat) && ncol(coef_mat) == nrow(df)) {
			q <- apply(
				coef_mat,
				2,
				stats::quantile,
				probs = .ae_boot_probs(bt, 0.95),
				na.rm = TRUE
			)
			df$lo <- q[1, ]; df$hi <- q[2, ]
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
	df_all$effect <- c(sender = "Sender", receiver = "Receiver")[df_all$effect]
	ab_subtitle <- if (!all(is.na(df_all$lo))) {
		"With bootstrap 95% intervals"
	} else if (!is.null(bt)) {
		"Point estimates (the bootstrap returns no interval for these effects)"
	} else {
		"Point estimates (refit with bootstrap replicates for intervals)"
	}
	p <- ggplot2::ggplot(
		df_all,
		ggplot2::aes(x = .data$value, y = .data$actor)
	) +
		ggplot2::geom_point()
	if (!all(is.na(df_all$lo))) {
		p <- p + ggplot2::geom_errorbar(
			ggplot2::aes(
				xmin = .data$lo,
				xmax = .data$hi
			),
			orientation = "y",
			width = 0.2,
			alpha = 0.6
		)
	}
	p + ggplot2::geom_vline(xintercept = 0, linetype = 3) +
		ggplot2::facet_wrap(~effect, scales = "free_y") +
		ggplot2::labs(
			title = "Additive effects",
			subtitle = ab_subtitle,
			x = "Effect (linear-predictor scale)",
			y = NULL
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			plot.caption = ggplot2::element_text(hjust = 0),
			axis.ticks = ggplot2::element_blank(),
			legend.position = "top",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y = ggplot2::element_text(color = "white", hjust = 0)
		)
}

#' Print the priors used by an AME / LAME / ame_als fit
#'
#' Bayesian-hygiene helper modelled on \code{rstanarm::prior_summary()}.
#' Prints recorded prior settings and the coefficient-prior specification
#' saved by the fitting build. Older fits without that specification are
#' identified explicitly; their prior is not inferred from the current build. For \code{ame_als} fits it just states that the estimator is a
#' frequentist point estimator and there are no priors.
#'
#' @param object a fitted \code{ame}, \code{lame}, or \code{ame_als} object.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' prior_summary(fit)
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
	.lame_cli_stdout({
		fit_kind <- if (inherits(object, "lame")) "lame" else "ame"
		cli::cli_h2("Recorded prior settings ({fit_kind} fit)")
		cp <- object[["coefficient_prior"]]
		if (!is.null(cp$description)) {
			cli::cli_text("{cp$description}")
		} else {
			cli::cli_text("The coefficient-prior specification was not recorded on this fit. Consult the sampler documentation for the build that produced it; g alone does not specify its covariance.")
		}
		uvp <- object[["static_uv_prior"]]
		if (!is.null(uvp$description)) {
			cli::cli_text("{uvp$description}")
		} else if (!is.null(uvp)) {
			scale_uv <- if (is.matrix(uvp$scale)) uvp$scale[1L, 1L] else 1
			cli::cli_text("Static latent-factor row covariance: inverse-Wishart with {uvp$df} degrees of freedom and scale matrix {round(scale_uv, 4)} times the identity in dimension {uvp$dimension}.")
			# data scale of the additive effects the prior is centred on (R/uv_prior.R)
			if (!is.null(uvp$data_scale)) {
				pm <- uvp$df * uvp$data_scale / (uvp$df - uvp$dimension - 1)
				cli::cli_text("That is Suv0 = {signif(uvp$data_scale, 3)} times the identity, the scale of the additive effects in the data, with the degrees of freedom {uvp$dimension + 2} (2R + 2) times {signif(uvp$df_multiplier, 3)} (the inverse Bernoulli information of the cells; 1 for ordinal, Poisson and normal fits), rounded: a prior mean of {signif(pm, 3)} times the identity.")
			}
		}
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
		# Suv0 and kappa0 are cross-sectional options and do not configure the longitudinal prior
		keys <- c("Sab0", "eta0", "etaab", "s20", "s2u0")
		if (!inherits(object, "lame")) keys <- c(keys, "Suv0", "kappa0")
		for (nm in keys) {
			v <- pr[[nm]]
			if (!is.null(v)) {
				val_str <- if (is.matrix(v)) {
					paste0(
						"matrix(c(",
						paste(round(as.vector(v), 4), collapse = ", "),
						"), ",
						nrow(v),
						", ",
						ncol(v),
						")"
					)
				} else {
					paste(round(unlist(v), 4), collapse = ", ")
				}
				cli::cli_text("  {.code {nm}} = {.val {val_str}}")
			}
		}
		if (!is.null(object$s2_prior_scale)) cli::cli_text("  Resolved residual-variance prior scale = {.val {object$s2_prior_scale}}")
		if (length(pr) == 0L) {
			cli::cli_text(
				cli::col_grey("Other prior hyperparameters: sampler-filled defaults (see ?ame Prior Distributions).")
			)
		}
	})
	invisible(object)
}

#' @rdname prior_summary
#' @method prior_summary lame
#' @export
prior_summary.lame <- function(object, ...) {
	.lame_cli_stdout({
		prior_summary.ame(object, ...)
		pr <- object$prior %||% list()
		if (isTRUE(object$dynamic_uv) || isTRUE(object$dynamic_ab) ||
				isTRUE(object$dynamic_beta)) {
			cli::cli_h3("Dynamic process prior settings")
			dynamic_keys <- c(
				if (isTRUE(object$dynamic_uv)) c(
					"rho_uv_mean", "rho_uv_sd", "sigma_uv_shape", "sigma_uv_scale", "uv_var_cap", "sigma_uv_max", "uv_max_abs"
				),
				if (isTRUE(object$dynamic_ab)) c(
					"rho_ab_mean", "rho_ab_sd", "sigma_ab_shape", "sigma_ab_scale", "ab_var_cap"
				)
			)
			for (nm in dynamic_keys) {
				v <- pr[[nm]]
				if (!is.null(v)) cli::cli_text("  {.code {nm}} = {.val {round(v, 4)}}")
			}
			if (isTRUE(object$dynamic_beta)) {
				for (nm in c(
					"rho_beta_mean", "rho_beta_sd",
					"sigma_beta_shape", "sigma_beta_scale",
					"sigma_beta_init",
					"rho_beta_lower", "rho_beta_upper",
					"beta0_mean", "beta0_var"
				)) {
					v <- pr[[nm]]
					if (!is.null(v)) cli::cli_text("  {.code {nm}} = {.val {round(v, 4)}}")
				}
			}
		}
	})
	invisible(object)
}

#' @rdname prior_summary
#' @method prior_summary ame_als
#' @export
prior_summary.ame_als <- function(object, ...) {
	.lame_cli_stdout({
		cli::cli_h2("Priors (ame_als fit)")
		if (identical(object$non_normal_method, "irls") && object$R > 0L) {
			cli::cli_text("IRLS estimates fixed effects and regularized latent factors. Gaussian factor penalties use empirically updated variance scales; the result is a point estimate, without posterior draws.")
		} else {
			cli::cli_text("The ALS estimator is a frequentist iterative block coordinate descent point estimator. There are no priors.")
		}
		cli::cli_text("For uncertainty, refit with {.code bootstrap = N} (parametric / block bootstrap).")
		cli::cli_text("For Bayesian inference (priors, credible intervals, posterior summaries), use {.fn ame} / {.fn lame} (MCMC).")
	})
	invisible(object)
}

# full Poisson rate, zero on missing cells; an overflowing rate is an error, not a zero count
.als_pois_rate <- function(EZ) {
	lam <- exp(EZ)
	if (any(!is.finite(lam[!is.na(EZ)]))) {
		cli::cli_abort(
			c(
				"Poisson simulation rates exceed the finite numerical range.",
				"i" = "Inspect the fitted predictor and latent variance before simulating."
			),
			class = "lame_als_prediction"
		)
	}
	lam[is.na(EZ)] <- 0
	lam
}
