# global variables for r cmd check
if(getRversion() >= "2.15.1") {
	utils::globalVariables("lab")
}

#' Visualize goodness-of-fit statistics for AME and LAME models
#'
#' @description
#' Plots observed network statistics against their posterior predictive
#' distributions to assess model fit.
#' 
#' @details
#' \strong{Overview:}
#' 
#' Goodness-of-fit (GOF) assessment is crucial for network models because standard
#' residual diagnostics are often inadequate for capturing network dependencies.
#' This function implements posterior predictive checking by:
#' 
#' 1. Computing key network statistics from the observed data
#' 2. Generating multiple networks from the model's posterior predictive distribution
#' 3. Computing the same statistics on simulated networks
#' 4. Visualizing the comparison to identify model inadequacies
#' 
#' \strong{Network Statistics Evaluated:}
#' 
#' \emph{For unipartite (square) networks:}
#' \describe{
#'   \item{\code{sd.row} (Out-degree heterogeneity)}{
#'     Standard deviation of row means. High values indicate substantial variation
#'     in how active nodes are as senders/initiators. If the model underestimates
#'     this, it may be missing important sender effects or covariates.
#'   }
#'   \item{\code{sd.col} (In-degree heterogeneity)}{
#'     Standard deviation of column means. High values indicate substantial variation
#'     in node popularity as receivers. Underestimation suggests missing receiver
#'     effects or popularity-related covariates.
#'   }
#'   \item{\code{dyad.dep} (Reciprocity/Mutuality)}{
#'     Correlation between \code{Y[i,j]} and \code{Y[j,i]}. Positive values indicate reciprocity
#'     (mutual ties are more likely). The AME model captures this through the
#'     dyadic correlation parameter rho. Poor fit here points to adjusting
#'     the dcor parameter.
#'   }
#'   \item{\code{triad.dep} (Transitivity/Clustering)}{
#'     Measures tendency for triadic closure (friend of a friend is a friend).
#'     Calculated as correlation between \code{Y[i,j]} and \code{sum(Y[i,k]*Y[k,j])/sqrt(n-2)}.
#'     AME captures this through multiplicative effects (U,V). Poor fit suggests
#'     increasing the latent dimension R.
#'   }
#' }
#' 
#' \emph{For bipartite (rectangular) networks:}
#' \describe{
#'   \item{\code{sd.row} (Type A activity variation)}{
#'     Standard deviation of row means for Type A nodes. Indicates heterogeneity
#'     in how actively Type A nodes connect to Type B nodes.
#'   }
#'   \item{\code{sd.col} (Type B popularity variation)}{
#'     Standard deviation of column means for Type B nodes. Indicates heterogeneity
#'     in how popular Type B nodes are with Type A nodes.
#'   }
#'   \item{\code{four.cycles} (Bipartite clustering)}{
#'     Count of 4-cycles (rectangular paths A1-B1-A2-B2-A1). High values indicate
#'     that pairs of Type A nodes tend to connect to the same Type B nodes.
#'     Captured through bipartite multiplicative effects with appropriate R_row, R_col.
#'   }
#' }
#' 
#' \strong{Interpretation Guide:}
#' 
#' \emph{Histogram plots (static models):}
#' - Dashed orange vertical line (Okabe-Ito \code{#D55E00}): observed statistic
#'   value (dual-encoded on colour and linetype)
#' - Grey histogram: distribution from posterior predictive simulations
#' - Good fit: orange line falls within the bulk of the histogram
#' - Poor fit: orange line in the tail or outside the distribution
#' 
#' \emph{Common model inadequacies and solutions:}
#' \describe{
#'   \item{Observed sd.row/sd.col too high}{
#'     Model underestimates degree heterogeneity.
#'     Solutions: Add row/column covariates (Xrow, Xcol), enable random effects
#'     (rvar=TRUE, cvar=TRUE), or increase their variance.
#'   }
#'   \item{Observed dyad.dep outside distribution}{
#'     Reciprocity not captured well.
#'     Solutions: For positive reciprocity, ensure dcor=TRUE. For negative,
#'     consider transformation or different family.
#'   }
#'   \item{Observed triad.dep too high}{
#'     Clustering/transitivity underestimated.
#'     Solutions: Increase latent dimension R, add network covariates that
#'     capture homophily, or consider including community structure covariates.
#'   }
#'   \item{Multiple statistics showing poor fit}{
#'     Fundamental model misspecification.
#'     Solutions: Change family, add missing covariates, or consider
#'     different model class.
#'   }
#' }
#' 
#' \strong{Mathematical Details:}
#' 
#' The posterior predictive p-value for statistic s is:
#' \deqn{p = P(s(Y_{rep}) \geq s(Y_{obs}) | Y_{obs})}
#' 
#' where Y_rep is drawn from the posterior predictive distribution.
#' Values near 0 or 1 indicate poor fit. The function visualizes the full
#' distribution rather than just p-values for richer diagnostics.
#' 
#' \strong{Customization:}
#' 
#' The function accepts custom GOF statistics through the model's custom_gof
#' argument. These are automatically included in the plot. Custom statistics
#' should capture network features important for your specific application.
#' 
#' \strong{Computational Notes:}
#' 
#' GOF computation involves simulating multiple networks, which can be
#' computationally intensive. The model uses the networks simulated during
#' MCMC if gof=TRUE was specified. For post-hoc GOF, use the gof() function
#' which generates new simulations.
#' 
#' Good model fit is indicated when:
#' - Observed statistics fall within the central 95% of posterior predictive distributions
#' - No systematic patterns of over/underestimation across statistics
#' - Custom statistics (if provided) also show adequate fit
#' - For longitudinal models, observed trajectories track within credible bands
#' 
#' @param fit An object of class "ame" or "lame" containing GOF statistics
#' @param type Character string: "auto" (default), "static", or "longitudinal".
#'        If "auto", determined by model class.
#' @param statistics Character vector of statistics to plot, or \code{NULL}
#'        (default) for the standard panels plus any custom GOF columns.
#'        Unipartite names: "sd.row", "sd.col", "dyad.dep", "triad.dep",
#'        "trans.dep"; bipartite names: "sd.row", "sd.col", "four.cycles".
#'        The internal GOF column names reported by \code{names(fit$GOF)}
#'        ("sd.rowmean", "sd.colmean", "cycle.dep") are accepted as
#'        equivalents; custom GOF columns may also be named. Unknown names
#'        are dropped with a warning.
#' @param credible.level Numeric between 0 and 1; credible interval level for
#'        longitudinal plots (default 0.95)
#' @param ncol Number of columns for faceted plot layout (default 2)
#' @param point.size Size of points in longitudinal plots (default 2)
#' @param line.size Width of lines in plots (default 1)
#' @param title Optional title for the plot
#' @param ... Additional arguments forwarded to the ALS-specific
#'   \code{\link{gof_plot.ame_als}} when \code{fit} inherits from
#'   \code{ame_als} (\code{nsim}, \code{seed}).
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' fit_ame <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, gof = TRUE,
#'                nscan = 100, burn = 10, odens = 1, verbose = FALSE)
#'
#' # Basic GOF plot
#' gof_plot(fit_ame)
#'
#' # Plot only degree-related statistics
#' gof_plot(fit_ame, statistics = c("sd.row", "sd.col"))
#' }
#' @export
#' @import ggplot2
gof_plot <- function(
	fit,
	type = c("auto", "static", "longitudinal"),
	statistics = NULL,
	credible.level = 0.95,
	ncol = 2,
	point.size = 2,
	line.size = 1,
	title = NULL,
	...) {

	# ame_als fits don't carry a posterior-predictive gof chain; delegate to
	# the bootstrap-based gof_plot.ame_als helper so users get something useful
	# instead of a class-check rejection. forward optional als-specific args
	# (nsim, seed) via the dots.
	if (inherits(fit, "ame_als")) {
		dots <- list(...)
		return(do.call(gof_plot.ame_als, c(list(fit = fit), dots)))
	}

	if (!inherits(fit, c("ame", "lame"))) {
		stop("fit must be an object of class 'ame', 'lame', or 'ame_als'")
	}
	
	if (is.null(fit$GOF)) {
		stop("No GOF statistics found. Re-run model with gof = TRUE")
	}
	
	type <- match.arg(type)
	if (type == "auto") {
		type <- ifelse(inherits(fit, "lame"), "longitudinal", "static")
	}
	
	is_bipartite <- FALSE
	if(!is.null(colnames(fit$GOF))) {
		is_bipartite <- "four.cycles" %in% colnames(fit$GOF)
	} else if(!is.null(fit$mode)) {
		is_bipartite <- fit$mode == "bipartite"
	}
	
	# `statistics = null` (the default) means "use the standard panels plus any
	# custom gof columns"; a user-supplied vector is honoured (and may name
	# custom columns). unknown names are dropped with a warning rather than a
	# hard match.arg() error.
	user_stats <- statistics
	if(is_bipartite) {
		stat.names <- c("sd.row" = "sd.rowmean",
										"sd.col" = "sd.colmean",
										"four.cycles" = "four.cycles")
		default_stats <- c("sd.row", "sd.col", "four.cycles")
	} else {
		stat.names <- c("sd.row" = "sd.rowmean",
										"sd.col" = "sd.colmean",
										"dyad.dep" = "dyad.dep",
										"triad.dep" = "cycle.dep",
										"trans.dep" = "trans.dep")
		default_stats <- c("sd.row", "sd.col", "dyad.dep", "triad.dep", "trans.dep")
	}
	# custom gof columns (e.g. from a custom_gof function) become plottable
	custom_cols <- character(0)
	if (!is.null(colnames(fit$GOF))) {
		custom_cols <- setdiff(colnames(fit$GOF), unname(stat.names))
		for (col in custom_cols) stat.names[col] <- col
	}
	if (is.null(user_stats)) {
		statistics <- c(default_stats, custom_cols)
	} else {
		# a statistic may be named either by its documented alias ("sd.row",
		# "triad.dep") or by the internal gof column it maps to ("sd.rowmean",
		# "cycle.dep"), which is what names(fit$GOF) reports. canonicalise the
		# internal spellings onto the alias so both select the same panel once.
		canonical <- names(stat.names)
		names(canonical) <- unname(stat.names)
		resolved <- ifelse(user_stats %in% names(canonical),
												canonical[user_stats], user_stats)
		known <- resolved %in% names(stat.names)
		bad <- unique(user_stats[!known])
		if (length(bad) > 0) {
			avail <- unique(c(names(stat.names), unname(stat.names)))
			cli::cli_warn(c(
				"Unknown gof statistic{?s} dropped: {.val {bad}}.",
				"i" = "Available statistics: {.val {avail}}."
			))
		}
		statistics <- unique(resolved[known])
		if (length(statistics) == 0) statistics <- default_stats
	}
	
	if (type == "static") {
		p <- gof_plot_static(fit, statistics, stat.names, ncol, line.size, title)
	} else {
		p <- gof_plot_longitudinal(fit, statistics, stat.names, credible.level,
															ncol, point.size, line.size, title)
	}
	
	return(p)
}

####

# label mapping for gof statistics
gof_stat_labels <- c(
	"sd.row" = "Sender Degree Heterogeneity",
	"sd.col" = "Receiver Degree Heterogeneity",
	"sd.rowmean" = "Sender Degree Heterogeneity",
	"sd.colmean" = "Receiver Degree Heterogeneity",
	"dyad.dep" = "Dyadic Dependence",
	"triad.dep" = "Triadic Dependence",
	"cycle.dep" = "Triadic Dependence",
	"four.cycles" = "Four Cycles",
	"trans.dep" = "Transitivity"
)

# apply display label to a stat name
gof_label <- function(x) {
	ifelse(x %in% names(gof_stat_labels), gof_stat_labels[x], x)
}

# standard theme for gof plots
gof_theme <- function() {
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

# static gof plots
gof_plot_static <- function(fit, statistics, stat.names, ncol, line.size, title) {

	gof_data <- fit$GOF

	if (is.list(gof_data) && !is.data.frame(gof_data)) {
		stat_names_full <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		n_iter <- nrow(gof_data[[1]])
		gof_matrix <- matrix(NA, n_iter, length(stat_names_full))
		colnames(gof_matrix) <- stat_names_full
		
		for (i in seq_along(stat_names_full)) {
			if (stat_names_full[i] %in% names(gof_data)) {
				if (ncol(gof_data[[stat_names_full[i]]]) > 1) {
					gof_matrix[, i] <- rowMeans(gof_data[[stat_names_full[i]]], na.rm = TRUE)
				} else {
					gof_matrix[, i] <- gof_data[[stat_names_full[i]]][, 1]
				}
			}
		}
		gof_data <- gof_matrix
	} else if (length(dim(gof_data)) == 3) {
		gof_data_avg <- apply(gof_data, c(1, 3), mean, na.rm = TRUE)
		gof_data <- t(gof_data_avg)
	}

	####

	obs_vals <- gof_data[1, ]
	pred_vals <- gof_data[-1, , drop = FALSE]

	plot_list <- list()
	
	for (stat in statistics) {
		if (stat %in% names(stat.names)) {
			stat_col <- stat.names[stat]
		} else {
			stat_col <- stat
		}
		
		if (!stat_col %in% colnames(gof_data)) {
			next
		}
		
		if (nrow(pred_vals) == 0) {
			warning("No posterior predictive samples available for GOF plots")
			next
		}

		pred_stat <- as.numeric(pred_vals[, stat_col])
		pred_stat <- pred_stat[is.finite(pred_stat)]
		if (length(pred_stat) == 0L || !is.finite(obs_vals[stat_col])) {
			next
		}
		
		df <- data.frame(
			value = pred_stat,
			statistic = stat,
			stringsAsFactors = FALSE
		)

		# observed value plotted as a dashed-and-coloured vline so the cue
		# survives both colour-blind viewers (linetype is the secondary
		# encoding) and grayscale printing. shape mapping at zero-extent
		# segments lets the colour-blind-safe okabe-ito red anchor a real
		# legend instead of relying on prose to explain it.
		obs_df <- data.frame(value = obs_vals[stat_col], lab = "Observed")
		p_stat <- ggplot(df, aes(x = value)) +
			geom_histogram(aes(y = after_stat(density), fill = "Posterior predictive"),
			               bins = 30) +
			geom_vline(data = obs_df,
			           aes(xintercept = value, colour = lab, linetype = lab),
			           linewidth = line.size) +
			scale_fill_manual(NULL, values = c("Posterior predictive" = "grey60")) +
			scale_colour_manual(NULL, values = c("Observed" = "#D55E00")) +
			scale_linetype_manual(NULL, values = c("Observed" = "dashed")) +
			labs(
				x = "Statistic Value",
				y = "Density",
				title = gof_label(stat)
			) +
			gof_theme() +
			theme(
				plot.title = element_text(size = 10, hjust = 0.5),
				legend.title = element_blank()
			) +
			guides(fill = guide_legend(order = 1),
			       colour = guide_legend(order = 2),
			       linetype = guide_legend(order = 2))

		plot_list[[stat]] <- p_stat
	}

	####

		if (length(plot_list) == 0L) {
			cli::cli_abort("No finite GOF values are available to plot.")
		} else if (length(plot_list) == 1) {
			p <- plot_list[[1]]
		} else {
		if (requireNamespace("patchwork", quietly = TRUE)) {
			p <- patchwork::wrap_plots(plot_list, ncol = ncol)
		} else if (requireNamespace("gridExtra", quietly = TRUE)) {
			p <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
		} else {
			warning("Install 'patchwork' or 'gridExtra' for better plot layout")
			p <- plot_list[[1]]
		}
	}
	
	if (!is.null(title)) {
		if (requireNamespace("patchwork", quietly = TRUE)) {
			p <- p + patchwork::plot_annotation(title = title)
		}
	} else {
		if (requireNamespace("patchwork", quietly = TRUE)) {
			p <- p + patchwork::plot_annotation(
				title = "Goodness-of-Fit: Observed vs Posterior Predictive"
			)
		}
	}
	
	return(p)
}

####

# longitudinal gof plots
gof_plot_longitudinal <- function(
	fit, statistics, stat.names, credible.level,
	ncol, point.size, line.size, title) {
	
	if (!is.null(fit$GOF_T)) {
		gof_data <- fit$GOF_T
	} else if (!is.null(fit$GOF) && is.list(fit$GOF) && !is.data.frame(fit$GOF)) {
		# use the gof list's actual names -- a hardcoded unipartite set would
		# omit the bipartite "four.cycles" column (and any custom statistic)
		stat_names_full <- names(fit$GOF)

		# gof list: rows = time periods, cols = mcmc iterations
		# first col = observed, rest = simulated
		if (ncol(fit$GOF[[1]]) > 1) {
			gof_data <- list()
			n_time <- nrow(fit$GOF[[1]])  # number of real time periods
			n_iter <- ncol(fit$GOF[[1]])  # number of mcmc iterations (obs + simulated)

			for (t in 1:n_time) {
				# rows = iterations, cols = stats for each time period
				gof_t <- matrix(NA, n_iter, length(stat_names_full))
				colnames(gof_t) <- stat_names_full

				for (i in seq_along(stat_names_full)) {
					if (stat_names_full[i] %in% names(fit$GOF)) {
						gof_t[, i] <- fit$GOF[[stat_names_full[i]]][t, ]
					}
				}
				gof_data[[t]] <- gof_t
			}
		} else {
			warning("No time-indexed GOF found, using static GOF")
			return(gof_plot_static(fit, statistics, stat.names, ncol, line.size, title))
		}
	} else if (!is.null(fit$GOF) && length(dim(fit$GOF)) == 3) {
		gof_data <- list()
		n_time <- dim(fit$GOF)[2]
		stat_names <- dimnames(fit$GOF)[[1]]
		for (t in 1:n_time) {
			gof_t <- t(fit$GOF[, t, ])
			colnames(gof_t) <- stat_names
			gof_data[[t]] <- gof_t
		}
	} else {
		warning("No time-indexed GOF found, using static GOF")
		return(gof_plot_static(fit, statistics, stat.names, ncol, line.size, title))
	}

	####

	alpha <- (1 - credible.level) / 2
	n_time <- length(gof_data)
	plot_data <- data.frame()
	
		for (t in 1:n_time) {
			gof_t <- gof_data[[t]]
			obs_t <- gof_t[1, ]
			pred_t <- gof_t[-1, , drop = FALSE]

			for (stat in statistics) {
				if (stat %in% names(stat.names)) {
					stat_col <- stat.names[stat]
				} else {
					stat_col <- stat
				}
				if (!stat_col %in% colnames(gof_t)) {
					next
				}

				pred_stat <- as.numeric(pred_t[, stat_col])
				pred_stat <- pred_stat[is.finite(pred_stat)]
				obs_val <- as.numeric(obs_t[stat_col])
				if (length(pred_stat) == 0L || !is.finite(obs_val)) {
					next
				}

				plot_data <- rbind(plot_data, data.frame(
					time = t,
					statistic = stat,
					observed = obs_val,
					median = median(pred_stat),
					lower = quantile(pred_stat, alpha, names = FALSE),
					upper = quantile(pred_stat, 1 - alpha, names = FALSE),
					stringsAsFactors = FALSE
				))
			}
		}

		####

		if (nrow(plot_data) == 0L) {
			cli::cli_abort("No finite GOF values are available to plot.")
		}

		# apply display labels
	plot_data$statistic <- gof_label(plot_data$statistic)

	# dual-encode the observed series (okabe-ito red + solid linetype with
	# points; predictive median is grey dashed line, ribbon is grey shading)
	# so the cue survives colour-blind viewers and grayscale printing. the
	# colour / linetype scale provide a real legend instead of the prose-only
	# "red = observed" convention.
	p <- ggplot(plot_data, aes(x = time)) +
		geom_ribbon(aes(ymin = lower, ymax = upper,
		                fill = "Posterior predictive 95%"),
		            alpha = 0.3) +
		geom_line(aes(y = median, colour = "Predictive median",
		              linetype = "Predictive median"),
		          linewidth = line.size) +
		geom_line(aes(y = observed, colour = "Observed",
		              linetype = "Observed"),
		          linewidth = line.size) +
		geom_point(aes(y = observed, colour = "Observed"),
		           size = point.size) +
		scale_fill_manual(NULL, values = c("Posterior predictive 95%" = "grey60")) +
		scale_colour_manual(NULL,
		                    values = c("Observed" = "#D55E00",
		                               "Predictive median" = "grey25")) +
		scale_linetype_manual(NULL,
		                      values = c("Observed" = "solid",
		                                 "Predictive median" = "dashed")) +
		facet_wrap(~ statistic, scales = "free_y", ncol = ncol) +
		labs(
			x = "Time Period",
			y = "Value"
		) +
		gof_theme() +
		theme(
			panel.spacing = unit(1, "lines")
		) +
		guides(fill = guide_legend(order = 1),
		       colour = guide_legend(order = 2),
		       linetype = guide_legend(order = 2))

	if (!is.null(title)) {
		p <- p + ggtitle(title)
	} else {
		p <- p + ggtitle(sprintf(
			"Longitudinal GOF: Observed vs %d%% Credible Intervals",
			round(credible.level * 100)
		))
	}
	
	return(p)
}
