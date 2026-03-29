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
	if (params %in% c("all", "beta") && !is.null(fit$BETA) && ncol(fit$BETA) > 0) {
		beta_data <- fit$BETA
		if (burn.in > 0 && burn.in < nrow(beta_data)) {
			beta_data <- beta_data[-(1:burn.in), , drop = FALSE]
		}
		
		if (thin > 1) {
			keep_idx <- seq(1, nrow(beta_data), by = thin)
			beta_data <- beta_data[keep_idx, , drop = FALSE]
		}
		
		for (j in 1:ncol(beta_data)) {
			param_name <- colnames(beta_data)[j]
			if (is.null(param_name)) param_name <- paste0("beta[", j, "]")
			
			plot_data <- rbind(plot_data, data.frame(
				iteration = 1:nrow(beta_data),
				value = beta_data[, j],
				parameter = param_name,
				type = "Regression",
				stringsAsFactors = FALSE
			))
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
		
		for (j in 1:ncol(vc_data)) {
			param_name <- colnames(vc_data)[j]
			if (is.null(param_name)) param_name <- paste0("vc[", j, "]")
			
			plot_data <- rbind(plot_data, data.frame(
				iteration = 1:nrow(vc_data),
				value = vc_data[, j],
				parameter = param_name,
				type = "Variance",
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

	trace_plots <- ggplot(plot_data, aes(x = iteration, y = value)) +
		geom_line(alpha = 0.7) +
		facet_wrap(~ parameter, scales = "free_y",
							ncol = ncol, nrow = nrow) +
		labs(
			x = "Iteration",
			y = "Value"
		) +
		theme_bw() +
		theme(
			panel.border = element_blank(),
			strip.background = element_rect(fill = "black", color = "black"),
			strip.text = element_text(color = "white", hjust = 0, size = 9),
			panel.spacing = unit(0.5, "lines")
		)

	####

	density_plots <- ggplot(plot_data, aes(x = value)) +
		geom_density() +
		facet_wrap(~ parameter, scales = "free",
							ncol = ncol, nrow = nrow) +
		labs(
			x = "Value",
			y = "Density"
		) +
		theme_bw() +
		theme(
			panel.border = element_blank(),
			strip.background = element_rect(fill = "black", color = "black"),
			strip.text = element_text(color = "white", hjust = 0, size = 9),
			panel.spacing = unit(0.5, "lines")
		)

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