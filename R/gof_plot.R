# global variables for R CMD check
if (getRversion() >= "2.15.1") {
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
#'     Used for families with structural zeros, including binary outcomes.
#'   }
#'   \item{\code{four.cycle.dep} (Weighted bipartite dependence)}{
#'     For Gaussian outcomes, averages products of centered, standardized cell
#'     values around fully observed rectangles. This measures weighted dependence;
#'     a nonzero-edge count would be fixed by the observed-cell pattern.
#'   }
#' }
#'
#' \strong{Interpretation Guide:}
#'
#' \emph{Histogram plots (static models):}
#' - Dashed black vertical line: observed statistic
#'   value (dual-encoded on colour and linetype)
#' - Grey histogram: distribution from posterior predictive simulations
#' - Good fit: dashed line falls within the bulk of the histogram
#' - Poor fit: dashed line in the tail or outside the distribution
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
#'        "trans.dep"; bipartite names: "sd.row", "sd.col", and
#'        "four.cycle.dep" for Gaussian outcomes or "four.cycles" otherwise.
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
#'   \code{ame_als} (\code{nsim}, \code{seed}). Its simulations condition on the
#'   fitted parameter estimates and resample outcome noise; they do not include
#'   parameter uncertainty. See \code{\link{ame_als_bootstrap}} for refit draws.
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit_ame <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, gof = TRUE,
#'                nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
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
	...
) {
	.lame_require_draws(fit, "gof_plot")

	# ALS fits have no posterior-predictive GOF chain: use the point-conditional
	# gof_plot.ame_als, passing extra arguments (nsim, seed) through the dots
	if (inherits(fit, "ame_als")) {
		dots <- list(...)
		return(do.call(gof_plot.ame_als, c(list(fit = fit), dots)))
	}

	if (!inherits(fit, c("ame", "lame"))) {
		cli::cli_abort("fit must be an object of class 'ame', 'lame', or 'ame_als'")
	}

	if (is.null(fit$GOF)) {
		cli::cli_abort("No GOF statistics found. Re-run model with gof = TRUE")
	}

	type <- match.arg(type)
	if (type == "auto") {
		type <- ifelse(inherits(fit, "lame"), "longitudinal", "static")
		# a single-period fit has no trajectory to draw
		if (type == "longitudinal" && is.list(fit$GOF) && is.matrix(fit$GOF[[1]]) &&
				nrow(fit$GOF[[1]]) == 1L) type <- "static"
	}

	is_bipartite <- FALSE
	if (!is.null(colnames(fit$GOF))) {
		is_bipartite <- any(c("four.cycles", "four.cycle.dep") %in% colnames(fit$GOF))
	} else if (!is.null(fit$mode)) {
		is_bipartite <- fit$mode == "bipartite"
	}

	# `statistics = NULL` means the standard panels plus any custom GOF columns;
	# unknown names in a user-supplied vector are dropped with a warning
	user_stats <- statistics
	if (is_bipartite) {
		stat.names <- c(
			"sd.row" = "sd.rowmean",
			"sd.col" = "sd.colmean",
			"four.cycles" = "four.cycles",
			"four.cycle.dep" = "four.cycle.dep"
		)
		default_stats <- c(
			"sd.row",
			"sd.col",
			if (identical(fit$family, "normal"))
				"four.cycle.dep" else "four.cycles"
		)
		if (identical(fit$family, "normal")) {
			if ("four.cycles" %in% user_stats) cli::cli_abort(
				c(
					"A nonzero four-cycle count cannot assess Gaussian predictive fit.",
					"i" = "Use {.val four.cycle.dep}; recompute older stored checks with {.fn gof}."
				),
				class = "lame_gof_uninformative"
			)
			stat.names <- stat.names[names(stat.names) != "four.cycles"]
			stored_names <- if (is.list(fit$GOF)) names(fit$GOF) else colnames(fit$GOF)
			if ("four.cycles" %in% stored_names) cli::cli_inform(c(
				"Stored nonzero four-cycle counts cannot check Gaussian predictive fit.",
				"i" = "Use {.fn gof} to recompute weighted four-cycle dependence."
			))
		}
	} else {
		stat.names <- c(
			"sd.row" = "sd.rowmean",
			"sd.col" = "sd.colmean",
			"dyad.dep" = "dyad.dep",
			"triad.dep" = "cycle.dep",
			"trans.dep" = "trans.dep"
		)
		default_stats <- c("sd.row", "sd.col", "dyad.dep", "triad.dep", "trans.dep")
	}
	# custom GOF columns (e.g. from a custom_gof function) become plottable
	custom_cols <- character(0)
	if (!is.null(colnames(fit$GOF))) {
		custom_cols <- setdiff(colnames(fit$GOF), unname(stat.names))
		if (is_bipartite && identical(fit$family, "normal"))
			custom_cols <- setdiff(custom_cols, "four.cycles")
		for (col in custom_cols) stat.names[col] <- col
	}
	if (is.null(user_stats)) {
		statistics <- c(default_stats, custom_cols)
	} else {
		# a statistic may be named by its alias ("sd.row") or by its internal GOF
		# column ("sd.rowmean"); map both onto the alias so each panel appears once
		canonical <- names(stat.names)
		names(canonical) <- unname(stat.names)
		resolved <- ifelse(
			user_stats %in% names(canonical),
			canonical[user_stats],
			user_stats
		)
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
		p <- gof_plot_longitudinal(
			fit,
			statistics,
			stat.names,
			credible.level,
			ncol,
			point.size,
			line.size,
			title
		)
	}

	return(p)
}


# display labels for GOF statistics
gof_stat_labels <- c(
	"sd.row" = "Sender degree heterogeneity",
	"sd.col" = "Receiver degree heterogeneity",
	"sd.rowmean" = "Sender degree heterogeneity",
	"sd.colmean" = "Receiver degree heterogeneity",
	"dyad.dep" = "Dyadic dependence",
	"triad.dep" = "Triadic dependence",
	"cycle.dep" = "Triadic dependence",
	"four.cycles" = "Four cycles",
	"four.cycle.dep" = "Weighted four-cycle dependence",
	"trans.dep" = "Transitivity"
)

# display label of a statistic name
gof_label <- function(x) {
	ifelse(x %in% names(gof_stat_labels), gof_stat_labels[x], x)
}

# standard theme for GOF plots
gof_theme <- function() {
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

# static GOF plots
gof_plot_static <- function(fit, statistics, stat.names, ncol, line.size, title) {

	gof_data <- fit$GOF

	if (is.list(gof_data) && !is.data.frame(gof_data)) {
		stat_names_full <- names(gof_data)
		n_iter <- ncol(gof_data[[1]])
		gof_matrix <- matrix(NA, n_iter, length(stat_names_full))
		colnames(gof_matrix) <- stat_names_full

		for (i in seq_along(stat_names_full)) {
			if (stat_names_full[i] %in% names(gof_data)) {
				if (ncol(gof_data[[stat_names_full[i]]]) > 1) {
					gof_matrix[, i] <- colMeans(gof_data[[stat_names_full[i]]], na.rm = TRUE)
				} else {
					gof_matrix[, i] <- colMeans(gof_data[[stat_names_full[i]]], na.rm = TRUE)
				}
			}
		}
		gof_data <- gof_matrix
	} else if (length(dim(gof_data)) == 3) {
		gof_data_avg <- apply(gof_data, c(1, 3), mean, na.rm = TRUE)
		gof_data <- t(gof_data_avg)
	}


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
			cli::cli_warn("No posterior predictive samples available for GOF plots")
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

		# observed value as a dashed vline with its own legend entry
		obs_df <- data.frame(value = obs_vals[stat_col], lab = "Observed")
		p_stat <- ggplot(df, aes(x = value)) +
			geom_histogram(
				aes(y = after_stat(density), fill = "Posterior predictive"),
				bins = 30
			) +
			geom_vline(
				data = obs_df,
				aes(xintercept = value, colour = lab, linetype = lab),
				linewidth = line.size
			) +
			scale_fill_manual(NULL, values = c("Posterior predictive" = "grey60")) +
			scale_colour_manual(NULL, values = c("Observed" = "black")) +
			scale_linetype_manual(NULL, values = c("Observed" = "dashed")) +
			labs(
				x = "Statistic value",
				y = "Density",
				title = gof_label(stat)
			) +
			gof_theme() +
			theme(
				plot.title = element_text(size = 10, hjust = 0),
				legend.title = element_blank()
			) +
			guides(
				fill = guide_legend(order = 1),
				colour = guide_legend(order = 2),
				linetype = guide_legend(order = 2)
			)

		plot_list[[stat]] <- p_stat
	}


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
			cli::cli_warn("Install 'patchwork' or 'gridExtra' for better plot layout")
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
				title = "Goodness of fit: observed vs posterior predictive"
			)
		}
	}

	return(p)
}


# longitudinal GOF plots
gof_plot_longitudinal <- function(
	fit,
	statistics,
	stat.names,
	credible.level,
	ncol,
	point.size,
	line.size,
	title
) {

	if (!is.null(fit$GOF_T)) {
		gof_data <- fit$GOF_T
	} else if (!is.null(fit$GOF) && is.list(fit$GOF) && !is.data.frame(fit$GOF)) {
		# take the statistics from the GOF list itself so bipartite and custom
		# columns are included
		stat_names_full <- names(fit$GOF)

		# GOF list: rows are periods, columns are draws (first observed, rest simulated)
		if (ncol(fit$GOF[[1]]) > 1) {
			gof_data <- list()
			n_time <- nrow(fit$GOF[[1]])
			n_iter <- ncol(fit$GOF[[1]]) # observed plus simulated draws

			for (t in 1:n_time) {
				# rows are draws, columns are statistics
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
			cli::cli_warn("No time-indexed GOF found, using static GOF")
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
		cli::cli_warn("No time-indexed GOF found, using static GOF")
		return(gof_plot_static(fit, statistics, stat.names, ncol, line.size, title))
	}


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

			plot_data <- rbind(
				plot_data,
				data.frame(
					time = t,
					statistic = stat,
					observed = obs_val,
					median = median(pred_stat),
					lower = quantile(pred_stat, alpha, names = FALSE),
					upper = quantile(pred_stat, 1 - alpha, names = FALSE),
					stringsAsFactors = FALSE
				)
			)
		}
	}


	if (nrow(plot_data) == 0L) {
		cli::cli_abort("No finite GOF values are available to plot.")
	}

	plot_data$statistic <- gof_label(plot_data$statistic)

	# observed series as a solid line with points, predictive median dashed grey,
	# interval as a grey band
	band_label <- paste0("Posterior predictive ", signif(100 * credible.level, 3), "%")
	p <- ggplot(plot_data, aes(x = time)) +
		geom_ribbon(
			aes(
				ymin = lower,
				ymax = upper,
				fill = band_label
			),
			alpha = 0.3
		) +
		geom_line(
			aes(
				y = median,
				colour = "Predictive median",
				linetype = "Predictive median"
			),
			linewidth = line.size
		) +
		geom_line(
			aes(
				y = observed,
				colour = "Observed",
				linetype = "Observed"
			),
			linewidth = line.size
		) +
		geom_point(
			aes(y = observed, colour = "Observed"),
			size = point.size
		) +
		scale_fill_manual(NULL, values = stats::setNames("grey60", band_label)) +
		scale_colour_manual(
			NULL,
			values = c(
				"Observed" = "black",
				"Predictive median" = "grey25"
			)
		) +
		scale_linetype_manual(
			NULL,
			values = c(
				"Observed" = "solid",
				"Predictive median" = "dashed"
			)
		) +
		facet_wrap(~statistic, scales = "free_y", ncol = ncol) +
		labs(
			x = "Period",
			y = "Value"
		) +
		gof_theme() +
		theme(
			panel.spacing = unit(1, "lines")
		) +
		guides(
			fill = guide_legend(order = 1),
			colour = guide_legend(order = 2),
			linetype = guide_legend(order = 2)
		)

	if (!is.null(title)) {
		p <- p + ggtitle(title)
	} else {
		p <- p + ggtitle(sprintf(
			"Goodness of fit over time: observed vs %d%% credible intervals",
			round(credible.level * 100)
		))
	}

	return(p)
}
