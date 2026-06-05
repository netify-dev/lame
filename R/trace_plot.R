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
#' 
#' @param fit An object of class "ame" or "lame" containing MCMC samples
#' @param params Character vector specifying which parameters to plot:
#'        "beta" for regression coefficients, "variance" for variance components,
#'        or "all" (default) for both
#' @param include Character vector of specific parameter names to include
#' @param exclude Character vector of specific parameter names to exclude
#' @param ncol Number of columns for plot layout (default 3)
#' @param nrow Number of rows for plot layout (default NULL, determined automatically)
#' @param burn.in Number of initial iterations to exclude as burn-in when 
#'        calculating statistics (default 0, assumes burn-in already removed)
#' @param thin Thinning interval for display (default 1, no thinning)
#' @param title Optional title for the plot
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
#'
#' # Basic trace plots for all parameters
#' trace_plot(fit)
#'
#' # Only regression coefficients
#' trace_plot(fit, params = "beta")
#'
#' # Only variance components
#' trace_plot(fit, params = "variance")
#' }
#' @export
#' @import ggplot2
#' @import reshape2
#' @importFrom stats ave
trace_plot <- function(
	fit,
	params = c("all", "beta", "variance"),
	include = NULL,
	exclude = NULL,
	ncol = 3,
	nrow = NULL,
	burn.in = 0,
	thin = 1,
	title = NULL
	){
	
	if (!inherits(fit, c("ame", "lame"))) {
		stop("fit must be an object of class 'ame' or 'lame'")
	}

	params <- match.arg(params)
	plot_data <- data.frame()
	# multi-chain fits attach a per-iteration chain_indicator; if present we
	# colour traces by chain so a brms-style user can see chain mixing
	chain_vec <- fit$chain_indicator
	has_chains <- !is.null(chain_vec) && length(unique(chain_vec)) > 1L

	# label mapping for variance components
	vc_labels <- c(
		"va" = "Sender Variance",
		"vb" = "Receiver Variance",
		"cab" = "Sender-Receiver Covariance",
		"rho" = "Dyadic Correlation",
		"ve" = "Error Variance"
	)

	####

	# regression coefficients
	if (params %in% c("all", "beta") && !is.null(fit$BETA) && length(dim(fit$BETA)) >= 2L) {
		beta_data <- fit$BETA
		beta_is_dyn <- length(dim(beta_data)) == 3L
		# burn / thin -- the first dim is iteration for both 2-d and 3-d
		if (burn.in > 0 && burn.in < dim(beta_data)[1]) {
			if (beta_is_dyn) beta_data <- beta_data[-(1:burn.in), , , drop = FALSE]
			else             beta_data <- beta_data[-(1:burn.in), , drop = FALSE]
		}
		if (thin > 1) {
			keep_idx <- seq(1, dim(beta_data)[1], by = thin)
			if (beta_is_dyn) beta_data <- beta_data[keep_idx, , , drop = FALSE]
			else             beta_data <- beta_data[keep_idx, , drop = FALSE]
		}
		# align chain_indicator with any burn-in / thinning we just applied
		beta_chain <- if (has_chains && length(chain_vec) == dim(fit$BETA)[1]) {
			cv <- chain_vec
			if (burn.in > 0 && burn.in < length(cv)) cv <- cv[-(1:burn.in)]
			if (thin > 1) cv <- cv[seq(1, length(cv), by = thin)]
			cv
		} else NULL

		beta_iter <- if (!is.null(beta_chain)) {
			ave(seq_along(beta_chain), beta_chain, FUN = seq_along)
		} else seq_len(dim(beta_data)[1])

		if (beta_is_dyn) {
			# 3-d beta: emit one trace per (coef, period) combination
			p_  <- dim(beta_data)[2]
			Tt_ <- dim(beta_data)[3]
			dn <- dimnames(beta_data)
			nms_p <- if (!is.null(dn)) dn[[2]] else paste0("beta[", seq_len(p_), "]")
			nms_t <- if (!is.null(dn)) dn[[3]] else paste0("t", seq_len(Tt_))
			for (k in seq_len(p_)) {
				for (t in seq_len(Tt_)) {
					param_name <- paste0(nms_p[k], "[", nms_t[t], "]")
					plot_data <- rbind(plot_data, data.frame(
						iteration = beta_iter,
						value = beta_data[, k, t],
						parameter = param_name,
						type = "Regression",
						chain = if (!is.null(beta_chain)) factor(beta_chain) else factor(rep(1L, dim(beta_data)[1])),
						stringsAsFactors = FALSE
					))
				}
			}
		} else if (ncol(beta_data) > 0L) {
			for (j in 1:ncol(beta_data)) {
				param_name <- colnames(beta_data)[j]
				if (is.null(param_name)) param_name <- paste0("beta[", j, "]")
				plot_data <- rbind(plot_data, data.frame(
					iteration = beta_iter,
					value = beta_data[, j],
					parameter = param_name,
					type = "Regression",
					chain = if (!is.null(beta_chain)) factor(beta_chain) else factor(rep(1L, nrow(beta_data))),
					stringsAsFactors = FALSE
				))
			}
		}
	}

	####

	# variance components
	if (params %in% c("all", "variance") && !is.null(fit$VC)) {
		vc_data <- fit$VC
		if (burn.in > 0 && burn.in < nrow(vc_data)) {
			vc_data <- vc_data[-(1:burn.in), , drop = FALSE]
		}

		if (thin > 1) {
			keep_idx <- seq(1, nrow(vc_data), by = thin)
			vc_data <- vc_data[keep_idx, , drop = FALSE]
		}

		vc_chain <- if (has_chains && length(chain_vec) == nrow(fit$VC)) {
			cv <- chain_vec
			if (burn.in > 0 && burn.in < length(cv)) cv <- cv[-(1:burn.in)]
			if (thin > 1) cv <- cv[seq(1, length(cv), by = thin)]
			cv
		} else NULL

		vc_iter <- if (!is.null(vc_chain)) {
			ave(seq_along(vc_chain), vc_chain, FUN = seq_along)
		} else seq_len(nrow(vc_data))

		for (j in 1:ncol(vc_data)) {
			param_name <- colnames(vc_data)[j]
			if (is.null(param_name)) param_name <- paste0("vc[", j, "]")

			plot_data <- rbind(plot_data, data.frame(
				iteration = vc_iter,
				value = vc_data[, j],
				parameter = param_name,
				type = "Variance",
				chain = if (!is.null(vc_chain)) factor(vc_chain) else factor(rep(1L, nrow(vc_data))),
				stringsAsFactors = FALSE
			))
		}
	}

	####

	# apply display labels
	plot_data$parameter <- ifelse(
		plot_data$parameter %in% names(vc_labels),
		vc_labels[plot_data$parameter],
		plot_data$parameter
	)

	if (!is.null(include)) {
		plot_data <- plot_data[plot_data$parameter %in% include, ]
	}
	if (!is.null(exclude)) {
		plot_data <- plot_data[!(plot_data$parameter %in% exclude), ]
	}
	
	if (nrow(plot_data) == 0) {
		stop("No parameters to plot after filtering")
	}

	####

	trace_plots <- if (has_chains) {
		ggplot(plot_data, aes(x = iteration, y = value, color = chain, group = chain)) +
			geom_line(alpha = 0.7) +
			facet_wrap(~ parameter, scales = "free_y",
								ncol = ncol, nrow = nrow) +
			labs(x = "Iteration", y = "Value", color = "Chain") +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				axis.ticks       = element_blank(),
				legend.position  = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing    = unit(0.5, "lines")
			)
	} else {
		ggplot(plot_data, aes(x = iteration, y = value)) +
			geom_line(alpha = 0.7) +
			facet_wrap(~ parameter, scales = "free_y",
								ncol = ncol, nrow = nrow) +
			labs(x = "Iteration", y = "Value") +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				axis.ticks       = element_blank(),
				legend.position  = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing    = unit(0.5, "lines")
			)
	}

	####

	density_plots <- if (has_chains) {
		ggplot(plot_data, aes(x = value, color = chain, group = chain)) +
			geom_density() +
			facet_wrap(~ parameter, scales = "free", ncol = ncol, nrow = nrow) +
			labs(x = "Value", y = "Density", color = "Chain") +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				axis.ticks       = element_blank(),
				legend.position  = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing    = unit(0.5, "lines")
			)
	} else {
		ggplot(plot_data, aes(x = value)) +
			geom_density() +
			facet_wrap(~ parameter, scales = "free", ncol = ncol, nrow = nrow) +
			labs(x = "Value", y = "Density") +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				axis.ticks       = element_blank(),
				legend.position  = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0, size = 9),
				panel.spacing    = unit(0.5, "lines")
			)
	}

	####

	if (requireNamespace("patchwork", quietly = TRUE)) {
		# trace 60% height, density 40%
		p <- trace_plots / density_plots + patchwork::plot_layout(heights = c(3, 2))

		if (!is.null(title)) {
			p <- p + patchwork::plot_annotation(title = title)
		} else {
			p <- p + patchwork::plot_annotation(
				title = "MCMC Diagnostics: Trace Plots (top) and Posterior Densities (bottom)"
			)
		}
	} else {
		warning("Install 'patchwork' for combined trace and density plots. Returning trace plots only.")
		p <- trace_plots
		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle("MCMC Trace Plots")
		}
	}
	
	return(p)
}
