#' MCMC trace plots and density plots for AME/LAME model parameters
#'
#' Creates diagnostic plots for Markov Chain Monte Carlo (MCMC) samples from
#' AME or LAME models. Displays trace plots to assess convergence and mixing,
#' alongside density plots to visualize posterior distributions.
#'
#' @details
#' This function produces two types of diagnostic plots:
#' \describe{
#'   \item{Trace plots}{Show the evolution of parameter values across MCMC
#'         iterations. Good mixing is indicated by rapid exploration of the
#'         parameter space with no trends or stuck periods.}
#'   \item{Density plots}{Show the posterior distribution of parameters.
#'         Multiple modes may indicate identification issues or convergence problems.}
#' }
#'
#' The plots help diagnose:
#' - Convergence: Has the chain reached the stationary distribution?
#' - Mixing: Is the chain exploring the parameter space efficiently?
#' - Autocorrelation: Are successive samples highly correlated?
#'
#' Parameters displayed include:
#' - Regression coefficients (beta)
#' - Variance components (va, vb, cab, rho, ve)
#' - Dynamic hyperparameters, when the fit has them: the persistence and
#'   innovation of the actor effects (\code{rho_ab}, \code{sigma_ab}) and of
#'   the latent positions (\code{rho_uv}, \code{sigma_uv}), the reset rate of
#'   snap positions (\code{pi_snap}), the per-group persistence and innovation
#'   of dynamic coefficients (\code{rho_beta}, \code{sigma_beta}) and the
#'   per-period dyadic correlation of \code{dynamic_rho} fits
#'
#' A quantity the model fixes draws no panel: the error variance of the
#' binary, ordinal, cbin and frn families (fixed at 1), the dyadic correlation
#' and sender-receiver covariance of bipartite fits, the dyadic correlation of
#' a fit with \code{dcor = FALSE} and the receiver terms of a symmetric fit.
#' Such a quantity is left out only when its draws are constant, so a stuck
#' chain on an estimated quantity always shows. A quantity that
#' \code{summary()} reports as not estimated (for example the effect
#' variances of a sparse panel) is left out as well, with a message naming it.
#'
#' For a fit with several chains, \code{burn.in} and \code{thin} apply within
#' each chain.
#'
#' @param fit An object of class "ame" or "lame" containing MCMC samples
#' @param params Character vector specifying which parameters to plot:
#'        "beta" for regression coefficients, "variance" for variance components,
#'        "dynamic" for the dynamic hyperparameters, or "all" (default) for all three
#' @param include Character vector of specific parameter names to include.
#'        Raw names ("va", "ve", "rho", "rho_uv", "sigma_ab", ...) and their
#'        display labels ("Error variance", "Position persistence (rho_uv)", ...)
#'        are both accepted; display labels match regardless of case.
#'        For dynamic_beta fits, per-period traces are named like
#'        \code{"x_dyad[t1]"}; supplying the base name (e.g. \code{"x_dyad"})
#'        selects all of its periods
#' @param exclude Character vector of specific parameter names to exclude.
#'        Base names match all per-period traces, as with \code{include}
#' @param ncol Number of columns for plot layout (default 3)
#' @param nrow Number of rows for plot layout (default NULL, determined automatically)
#' @param burn.in Number of initial iterations of each chain to exclude as
#'        burn-in (default 0, assumes burn-in already removed)
#' @param thin Thinning interval for display, within each chain (default 1, no thinning)
#' @param title Optional title for the plot
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
#'            nscan = 50, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Basic trace plots for all parameters
#' trace_plot(fit)
#'
#' # Only regression coefficients
#' trace_plot(fit, params = "beta")
#'
#' # Only variance components
#' trace_plot(fit, params = "variance")
#'
#' # The persistence and innovation of dynamic actor effects
#' Y <- lapply(1:3, function(t) YX_nrm$Y[1:20, 1:20])
#' fit_dyn <- lame(Y, R = 0, dynamic_ab = TRUE, nscan = 100, burn = 20, odens = 2,
#'                 verbose = FALSE, chains = 1)
#' trace_plot(fit_dyn, params = "dynamic")
#' }
#' @export
#' @import ggplot2
#' @importFrom stats ave
trace_plot <- function(
	fit,
	params = c("all", "beta", "variance", "dynamic"),
	include = NULL,
	exclude = NULL,
	ncol = 3,
	nrow = NULL,
	burn.in = 0,
	thin = 1,
	title = NULL
) {
	.lame_require_draws(fit, "trace_plot")

	if (!inherits(fit, c("ame", "lame"))) {
		cli::cli_abort("fit must be an object of class 'ame' or 'lame'")
	}

	params <- match.arg(params)
	plot_data <- data.frame()
	# multi-chain fits carry a per-iteration chain_indicator; traces are then
	# coloured by chain
	chain_vec <- fit$chain_indicator
	has_chains <- !is.null(chain_vec) && length(unique(chain_vec)) > 1L

	# display labels for variance components and dynamic hyperparameters
	vc_labels <- c(
		"va" = "Sender variance",
		"vb" = "Receiver variance",
		"cab" = "Sender-receiver covariance",
		"rho" = "Dyadic correlation",
		"ve" = "Error variance"
	)
	dyn_labels <- c(
		"rho_ab" = "Effect persistence (rho_ab)",
		"sigma_ab" = "Effect innovation (sigma_ab)",
		"rho_uv" = "Position persistence (rho_uv)",
		"sigma_uv" = "Position innovation (sigma_uv)",
		"pi_snap" = "Position reset rate (pi_snap)",
		"rho_beta" = "Coefficient persistence (rho_beta)",
		"sigma_beta" = "Coefficient innovation (sigma_beta)"
	)
	all_labels <- c(vc_labels, dyn_labels)

	# draws kept after burn-in and thinning within each chain, with the chain and
	# within-chain iteration of each kept draw
	.keep_draws <- function(n) {
		chain <- if (has_chains && length(chain_vec) == n) chain_vec else rep(1L, n)
		pos <- ave(seq_len(n), chain, FUN = seq_along)
		# a burn-in as long as a chain leaves that chain whole
		len <- ave(seq_len(n), chain, FUN = length)
		b <- ifelse(burn.in > 0 & burn.in < len, burn.in, 0)
		idx <- which(pos > b & ((pos - b - 1) %% thin == 0))
		list(
			idx = idx,
			chain = factor(chain[idx]),
			iter = ave(seq_along(idx), chain[idx], FUN = seq_along)
		)
	}
	.series_frame <- function(kd, value, param_name, type) {
		data.frame(
			iteration = kd$iter,
			value = value,
			parameter = param_name,
			type = type,
			chain = kd$chain,
			stringsAsFactors = FALSE
		)
	}

	# quantities summary() reports as not estimated draw no panel
	hidden <- character(0)
	ab_reason <- .lame_ab_var_unident_reason(fit)
	if (!is.na(ab_reason)) hidden <- .lame_ab_var_hidden_rows(ab_reason)
	if (.lame_uv_hyper_unidentified(fit)) hidden <- c(hidden, "rho_uv", "sigma_uv")
	if (.lame_snap_rate_unidentified(fit)) hidden <- c(hidden, "pi_snap")
	beta_names <- if (!is.null(fit$BETA) && length(dim(fit$BETA)) >= 2L) dimnames(fit$BETA)[[2L]] else NULL
	coef_hidden <- if (length(beta_names)) .lame_ab_coef_hidden_rows(
		.lame_ab_coef_reason(fit),
		beta_names
	) else character(0)
	left_out <- character(0)

	# quantities the model fixes; one is left out only when its draws are
	# constant, so a stuck chain on an estimated quantity still shows
	fixed_vc <- c(
		if (isTRUE(fit$family %in% c("binary", "ordinal", "cbin", "frn"))) "ve",
		if (identical(fit$mode, "bipartite")) c("cab", "rho"),
		if (identical(fit$dcor, FALSE)) "rho",
		if (isTRUE(fit$symmetric)) c("vb", "cab", "rho")
	)
	.is_constant <- function(x) {
		x <- x[is.finite(x)]
		length(x) > 0L && diff(range(x)) <= 1e-10 * max(1, max(abs(x)))
	}


	# regression coefficients
	if (params %in% c("all", "beta") && length(beta_names) && dim(fit$BETA)[1] > 0L) {
		beta_data <- fit$BETA
		beta_is_dyn <- length(dim(beta_data)) == 3L
		kd <- .keep_draws(dim(beta_data)[1])
		nms_p <- beta_names %||% paste0("beta[", seq_len(dim(beta_data)[2]), "]")
		left_out <- c(left_out, intersect(nms_p, coef_hidden))
		if (beta_is_dyn) {
			# dynamic beta: one trace per coefficient and period
			Tt_ <- dim(beta_data)[3]
			nms_t <- dimnames(beta_data)[[3]] %||% paste0("t", seq_len(Tt_))
			for (k in which(!nms_p %in% coef_hidden)) {
				for (t in seq_len(Tt_)) {
					plot_data <- rbind(
						plot_data,
						.series_frame(
							kd,
							beta_data[kd$idx, k, t],
							paste0(nms_p[k], "[", nms_t[t], "]"),
							"Regression"
						)
					)
				}
			}
		} else {
			for (j in which(!nms_p %in% coef_hidden)) {
				plot_data <- rbind(plot_data, .series_frame(kd, beta_data[kd$idx, j], nms_p[j], "Regression"))
			}
		}
	}


	# variance components
	if (params %in% c("all", "variance") && !is.null(fit$VC) && nrow(fit$VC) > 0L) {
		vc_data <- fit$VC
		kd <- .keep_draws(nrow(vc_data))
		vc_names <- colnames(vc_data) %||% paste0("vc[", seq_len(ncol(vc_data)), "]")
		for (j in seq_len(ncol(vc_data))) {
			if (vc_names[j] %in% hidden) {
				left_out <- c(left_out, vc_names[j])
				next
			}
			if (vc_names[j] %in% fixed_vc && .is_constant(vc_data[, j])) next
			plot_data <- rbind(plot_data, .series_frame(kd, vc_data[kd$idx, j], vc_names[j], "Variance"))
		}
	}


	# dynamic hyperparameters: one draw per stored sweep, as BETA and VC
	if (params %in% c("all", "dynamic")) {
		for (v in c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap")) {
			x <- fit[[v]]
			if (!is.numeric(x) || length(x) < 2L) next
			if (v %in% hidden) {
				left_out <- c(left_out, v)
				next
			}
			kd <- .keep_draws(length(x))
			plot_data <- rbind(plot_data, .series_frame(kd, as.numeric(x)[kd$idx], v, "Dynamic"))
		}
		# per-group coefficient hyperparameters and per-period dyadic correlations
		for (v in c("RHO_BETA", "SIGMA_BETA", "RHO")) {
			M <- fit[[v]]
			if (!is.numeric(M) || length(M) < 2L) next
			M <- as.matrix(M)
			kd <- .keep_draws(nrow(M))
			base <- if (v == "RHO") vc_labels[["rho"]] else dyn_labels[[tolower(v)]]
			grp <- colnames(M) %||% paste0("g", seq_len(ncol(M)))
			for (j in seq_len(ncol(M))) {
				param_name <- if (ncol(M) > 1L || v == "RHO") paste0(base, "[", grp[j], "]") else base
				plot_data <- rbind(plot_data, .series_frame(kd, M[kd$idx, j], param_name, "Dynamic"))
			}
		}
	}

	if (length(left_out)) {
		cli::cli_inform(c("i" = "Not drawn: {.val {unique(left_out)}}, which {.fn summary} reports as not estimated."))
	}


	# display labels
	if (nrow(plot_data) > 0L) {
		plot_data$parameter <- ifelse(
			plot_data$parameter %in% names(all_labels),
			all_labels[plot_data$parameter],
			plot_data$parameter
		)
	}

	# include / exclude match the trace name or its base name without a trailing
	# '[..]' suffix (so "x_dyad" selects every period), and accept raw names
	# ("ve", "rho_uv") or display labels in any case
	.map_filter <- function(x) {
		hit <- x %in% names(all_labels)
		x[hit] <- all_labels[x[hit]]
		lab_hit <- match(tolower(x), tolower(all_labels))
		x[!is.na(lab_hit)] <- all_labels[lab_hit[!is.na(lab_hit)]]
		x
	}
	include <- if (!is.null(include)) .map_filter(include) else include
	exclude <- if (!is.null(exclude)) .map_filter(exclude) else exclude
	base_name <- sub("\\[[^]]*\\]$", "", plot_data$parameter)
	if (!is.null(include)) {
		keep <- plot_data$parameter %in% include | base_name %in% include
		plot_data <- plot_data[keep, ]
		base_name <- base_name[keep]
	}
	if (!is.null(exclude)) {
		keep <- !(plot_data$parameter %in% exclude | base_name %in% exclude)
		plot_data <- plot_data[keep, ]
	}

	if (nrow(plot_data) == 0) {
		cli::cli_abort("No parameters to plot after filtering")
	}


	trace_plots <- if (has_chains) {
		ggplot(plot_data, aes(x = iteration, y = value, color = chain, group = chain)) +
			geom_line(alpha = 0.7) +
			facet_wrap(
				~parameter,
				scales = "free_y",
				ncol = ncol,
				nrow = nrow
			) +
			labs(x = "Iteration", y = "Value", color = "Chain") +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing = unit(0.5, "lines")
			)
	} else {
		ggplot(plot_data, aes(x = iteration, y = value)) +
			geom_line(alpha = 0.7) +
			facet_wrap(
				~parameter,
				scales = "free_y",
				ncol = ncol,
				nrow = nrow
			) +
			labs(x = "Iteration", y = "Value") +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing = unit(0.5, "lines")
			)
	}


	density_plots <- if (has_chains) {
		ggplot(plot_data, aes(x = value, color = chain, group = chain)) +
			geom_density() +
			facet_wrap(~parameter, scales = "free", ncol = ncol, nrow = nrow) +
			labs(x = "Value", y = "Density", color = "Chain") +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing = unit(0.5, "lines")
			)
	} else {
		ggplot(plot_data, aes(x = value)) +
			geom_density() +
			facet_wrap(~parameter, scales = "free", ncol = ncol, nrow = nrow) +
			labs(x = "Value", y = "Density") +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing = unit(0.5, "lines")
			)
	}


	if (requireNamespace("patchwork", quietly = TRUE)) {
		# trace 60% height, density 40%
		p <- trace_plots / density_plots + patchwork::plot_layout(heights = c(3, 2))

		if (!is.null(title)) {
			p <- p + patchwork::plot_annotation(title = title)
		} else {
			p <- p + patchwork::plot_annotation(
				title = "MCMC diagnostics: trace plots (top) and posterior densities (bottom)"
			)
		}
	} else {
		cli::cli_warn("Install 'patchwork' for combined trace and density plots. Returning trace plots only.")
		p <- trace_plots
		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle("MCMC trace plots")
		}
	}

	return(p)
}
