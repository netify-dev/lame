#' Visualize sender and receiver random effects
#'
#' Creates a visualization of the additive sender (row) and receiver (column)
#' random effects from an AME or LAME model. Automatically detects whether effects
#' are static or dynamic and provides appropriate visualization options.
#'
#' @details
#' The additive effects in AME models represent:
#' \describe{
#'   \item{Sender effects (a)}{Actor-specific tendencies to form outgoing ties.
#'         Positive values indicate actors who send more ties than expected;
#'         negative values indicate actors who send fewer ties.}
#'   \item{Receiver effects (b)}{Actor-specific tendencies to receive incoming ties.
#'         Positive values indicate actors who receive more ties than expected;
#'         negative values indicate actors who receive fewer ties.}
#' }
#'
#' For static effects, the plot displays these effects as a dot plot with vertical
#' lines extending from zero to each effect estimate.
#'
#' For dynamic effects (when fit contains a_dynamic/b_dynamic), additional options
#' are available to visualize how effects evolve over time.
#'
#' @param fit An object of class "ame" or "lame" from fitting an AME model
#' @param effect Character string specifying which effect to plot:
#'        "sender" (default), "receiver", or "both". "both" draws the sender
#'        and receiver effects as two panels, each sorted on its own; for
#'        dynamic fits it is available in the snapshot views
#'        (\code{plot_type = "snapshot"}, one period or
#'        \code{time_point = "average"}). Symmetric fits store a single set
#'        of actor effects, so only "sender" applies to them.
#' @param sorted Logical; if TRUE (default), actors are sorted by effect
#'        magnitude. Applies to \code{ame} / \code{lame} fits;
#'        \code{ame_als} fits are always sorted by value.
#' @param labels Logical; if TRUE, actor labels are shown on x-axis (default TRUE
#'        for n <= 50 actors). Applies to \code{ame} / \code{lame} fits.
#' @param title Optional title for the plot (\code{ame} / \code{lame} fits).
#' @param time_point For dynamic effects, which time point to plot (default: last).
#'                   Can be a numeric index, "all" for a faceted plot, or "average" for time-averaged
#' @param plot_type For dynamic effects: "snapshot" (single time), "trajectory" (evolution over time),
#'                  "faceted" (grid of time points), or "ribbon" (effect path with a 95\% posterior
#'                  band per period). For static effects, this parameter is ignored.
#' @param show_actors Character vector of specific actors to highlight (for dynamic trajectory / ribbon plots)
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Visualize sender effects
#' ab_plot(fit, effect = "sender")
#'
#' # Visualize receiver effects without sorting
#' ab_plot(fit, effect = "receiver", sorted = FALSE)
#' }
#' @export
#' @import ggplot2
#' @import cli
#' @importFrom stats reorder
ab_plot <- function(
	fit,
	effect = c("sender", "receiver", "both"),
	sorted = TRUE,
	labels = NULL,
	title = NULL,
	time_point = NULL,
	plot_type = c("snapshot", "trajectory", "faceted", "ribbon"),
	show_actors = NULL
) {
	.lame_require_draws(fit, "ab_plot")

	if (!inherits(fit, c("ame", "lame", "ame_als"))) {
		cli::cli_abort("fit must be an object of class 'ame', 'lame', or 'ame_als'")
	}

	effect <- match.arg(effect)
	plot_type <- match.arg(plot_type)

	is_dynamic <- !is.null(fit$a_dynamic) || !is.null(fit$b_dynamic)

	# static ALS fits have no posterior and use the ame_als lollipop chart;
	# dynamic ALS fits continue below so plot_type and show_actors apply
	if (inherits(fit, "ame_als") && !is_dynamic) {
		return(ab_plot.ame_als(fit, effect = effect))
	}

	# "both" facets the two margins side by side; for dynamic fits only the
	# snapshot views have a single value per actor and margin
	if (effect == "both") {
		if (is_dynamic && (plot_type != "snapshot" || identical(time_point, "all"))) {
			view <- if (identical(time_point, "all")) {
				'time_point = "all"'
			} else {
				paste0('plot_type = "', plot_type, '"')
			}
			cli::cli_abort(c(
				"{.code effect = \"both\"} is only available for static fits and dynamic snapshots (one period, or {.code time_point = \"average\"}).",
				"i" = "Plot {.val sender} and {.val receiver} separately with {.code {view}}."
			))
		}
		return(ab_plot_both_internal(
			fit,
			sorted,
			labels,
			title,
			time_point,
			is_dynamic
		))
	}

	if (is_dynamic) {
		return(ab_plot_dynamic_internal(
			fit,
			effect,
			sorted,
			labels,
			title,
			time_point,
			plot_type,
			show_actors
		))
	}


	# static single-margin plot
	m <- .ab_margin_values(fit, effect)
	mu_df <- data.frame(
		mu = unname(m$mu),
		id = names(m$mu),
		stringsAsFactors = FALSE
	)
	.ab_lollipop(
		mu_df,
		ylabel = m$ylabel,
		sorted = sorted,
		labels = labels,
		title = title,
		default_title = m$default_title
	)
}


# one margin's effect vector with its axis label, panel name and default
# title: the posterior means (APM / BPM) of a static fit, or one period of
# a_dynamic / b_dynamic (default the last) or their average for a dynamic fit
.ab_margin_values <- function(
	fit,
	margin,
	time_point = NULL,
	is_dynamic = FALSE
) {
	is_bip <- identical(fit$mode, "bipartite")
	if (margin == "sender") {
		what <- if (is_bip) "Row actor" else "Sender"
		ylabel <- paste0(what, " effects (a)")
		eff <- if (is_dynamic) fit$a_dynamic else fit$APM
	} else {
		what <- if (is_bip) "Column actor" else "Receiver"
		ylabel <- paste0(what, " effects (b)")
		eff <- if (is_dynamic) fit$b_dynamic else fit$BPM
	}
	if (is.null(eff)) {
		cli::cli_abort(c(
			"This fit has no receiver (b) effects.",
			"i" = "Symmetric models store a single set of actor effects; use {.code effect = \"sender\"}."
		))
	}
	panel <- ylabel
	suffix <- ""
	if (is_dynamic) {
		n_times <- ncol(eff)
		actor_names <- rownames(eff) %||% paste0("Actor", seq_len(nrow(eff)))
		time_labels <- colnames(eff) %||% paste0("T", seq_len(n_times))
		if (identical(time_point, "average")) {
			eff <- rowMeans(eff)
			ylabel <- paste(ylabel, "(time-averaged)")
			suffix <- " - average across time"
		} else {
			if (is.null(time_point)) time_point <- n_times
			if (!is.numeric(time_point)) {
				cli::cli_abort(c(
					"{.arg time_point} must be a number, {.val average}, or {.val all}.",
					"i" = "Got {.val {time_point}}."
				))
			}
			if (time_point > n_times || time_point < 1) {
				cli::cli_abort("{.arg time_point} must be between 1 and {n_times}.")
			}
			suffix <- paste0(" at ", time_labels[time_point])
			eff <- eff[, time_point]
		}
		names(eff) <- actor_names
	} else if (is.null(names(eff))) {
		names(eff) <- paste0("Actor", seq_along(eff))
	}
	default_title <- paste0(
		if (is_dynamic) "Dynamic " else "Additive ",
		tolower(what),
		" effects",
		suffix
	)
	list(
		mu = eff,
		ylabel = ylabel,
		panel = panel,
		suffix = suffix,
		default_title = default_title
	)
}

# lollipop chart of one effect per actor; mu_df has columns mu and id, plus
# margin when both panels are drawn, each panel sorted on its own
.ab_lollipop <- function(mu_df, ylabel, sorted, labels, title, default_title) {
	facet <- "margin" %in% names(mu_df)
	n_actors <- if (facet) max(table(mu_df$margin)) else nrow(mu_df)
	if (is.null(labels)) {
		labels <- n_actors <= 50
	}

	# x is a per-panel key so each panel orders its own actors; the axis maps the
	# keys back to actor names
	mu_df$key <- if (facet) paste(mu_df$margin, mu_df$id, sep = "\r") else mu_df$id
	ord <- if (facet) {
		if (sorted) order(mu_df$margin, mu_df$mu) else order(mu_df$margin)
	} else {
		if (sorted) order(mu_df$mu) else seq_len(nrow(mu_df))
	}
	mu_df$key <- factor(mu_df$key, levels = mu_df$key[ord])
	key_labels <- stats::setNames(mu_df$id, as.character(mu_df$key))

	p <- ggplot(mu_df, aes(x = .data$key, y = .data$mu)) +
		geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
		geom_segment(aes(xend = .data$key, yend = 0)) +
		geom_point(size = 2) +
		scale_x_discrete(labels = key_labels) +
		xlab("") +
		ylab(ylabel) +
		theme_bw() +
		theme(
			panel.border = element_blank(),
			plot.caption = element_text(hjust = 0),
			axis.ticks = element_blank(),
			legend.position = "top"
		)

	if (facet) {
		p <- p + facet_wrap(~margin, scales = "free_x") +
			theme(
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0)
			)
	}

	if (!labels) {
		p <- p + theme(
			axis.text.x = element_blank()
		) + xlab(if (sorted) "Actors (sorted by effect magnitude)" else "Actors")
	} else {
		# smaller, steeper labels as the number of actors grows
		ang <- if (n_actors > 30L) 90 else if (n_actors > 15L) 75 else 45
		size <- if (n_actors > 30L) 7 else if (n_actors > 15L) 8 else 10
		p <- p + theme(
			axis.text.x = element_text(angle = ang, hjust = 1, size = size)
		)
	}

	if (!is.null(title)) {
		p <- p + ggtitle(title)
	} else if (sorted) {
		p <- p + ggtitle(paste(default_title, "(sorted)"))
	} else {
		p <- p + ggtitle(default_title)
	}

	p
}

# sender and receiver effects side by side (effect = "both")
ab_plot_both_internal <- function(
	fit,
	sorted,
	labels,
	title,
	time_point,
	is_dynamic
) {
	a <- .ab_margin_values(fit, "sender", time_point, is_dynamic)
	b <- .ab_margin_values(fit, "receiver", time_point, is_dynamic)
	panels <- c(a$panel, b$panel)
	mu_df <- rbind(
		data.frame(
			mu = unname(a$mu),
			id = names(a$mu),
			margin = panels[1L],
			stringsAsFactors = FALSE
		),
		data.frame(
			mu = unname(b$mu),
			id = names(b$mu),
			margin = panels[2L],
			stringsAsFactors = FALSE
		)
	)
	mu_df$margin <- factor(mu_df$margin, levels = panels)
	is_bip <- identical(fit$mode, "bipartite")
	default_title <- paste0(
		if (is_dynamic) "Dynamic " else "Additive ",
		if (is_bip) "row and column actor" else "sender and receiver",
		" effects",
		a$suffix
	)
	ylabel <- if (identical(time_point, "average")) {
		"Additive effect (time-averaged)"
	} else {
		"Additive effect"
	}
	.ab_lollipop(
		mu_df,
		ylabel = ylabel,
		sorted = sorted,
		labels = labels,
		title = title,
		default_title = default_title
	)
}


#' Internal function for dynamic additive effects plotting
#' @keywords internal
#' @noRd
ab_plot_dynamic_internal <- function(
	fit,
	effect,
	sorted,
	labels,
	title,
	time_point,
	plot_type,
	show_actors
) {

	if (effect == "sender") {
		effects_mat <- fit$a_dynamic
		ylabel <- "Sender effects (a)"
		default_title <- "Dynamic sender effects"
	} else {
		effects_mat <- fit$b_dynamic
		ylabel <- "Receiver effects (b)"
		default_title <- "Dynamic receiver effects"
	}

	n_actors <- nrow(effects_mat)
	n_times <- ncol(effects_mat)
	actor_names <- rownames(effects_mat)
	if (is.null(actor_names)) actor_names <- paste0("Actor", 1:n_actors)

	time_labels <- colnames(effects_mat)
	if (is.null(time_labels)) time_labels <- paste0("T", 1:n_times)

	if (plot_type == "snapshot") {
		if (is.null(time_point)) {
			time_point <- n_times
		} else if (time_point == "average") {
			mu_effect <- rowMeans(effects_mat)
			names(mu_effect) <- actor_names

			mu_df <- data.frame(
				mu = mu_effect,
				id = actor_names,
				stringsAsFactors = FALSE
			)

			if (sorted) {
				mu_df$id <- factor(mu_df$id, levels = mu_df$id[order(mu_df$mu)])
			} else {
				mu_df$id <- factor(mu_df$id, levels = mu_df$id)
			}

			if (is.null(labels)) {
				labels <- n_actors <= 50
			}

			p <- ggplot(mu_df, aes(x = id, y = mu)) +
				geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
				geom_segment(aes(xend = id, yend = 0)) +
				geom_point(size = 2) +
				xlab("") + ylab(paste(ylabel, "(time-averaged)")) +
				theme_bw() +
				theme(
					panel.border = element_blank(),
					plot.caption = element_text(hjust = 0),
					axis.ticks = element_blank(),
					legend.position = "top"
				)

			if (!labels) {
				p <- p + theme(
					axis.text.x = element_blank()
				) + xlab("Actors")
			} else {
				p <- p + theme(
					axis.text.x = element_text(angle = 45, hjust = 1)
				)
			}

			if (!is.null(title)) {
				p <- p + ggtitle(title)
			} else {
				p <- p + ggtitle(paste(default_title, "- average across time"))
			}

			return(p)
		} else if (identical(time_point, "all")) {
			return(ab_plot_dynamic_internal(
				fit,
				effect,
				sorted,
				labels,
				title,
				NULL,
				"faceted",
				show_actors
			))
		}

		if (!is.numeric(time_point)) {
			cli::cli_abort(c(
				"{.arg time_point} must be a number, {.val average}, or {.val all}.",
				"i" = "Got {.val {time_point}}."
			))
		}
		if (time_point > n_times || time_point < 1) {
			cli::cli_abort("{.arg time_point} must be between 1 and {n_times}.")
		}

		mu_effect <- effects_mat[, time_point]
		names(mu_effect) <- actor_names

		mu_df <- data.frame(
			mu = mu_effect,
			id = actor_names,
			stringsAsFactors = FALSE
		)

		if (sorted) {
			mu_df$id <- factor(mu_df$id, levels = mu_df$id[order(mu_df$mu)])
		} else {
			mu_df$id <- factor(mu_df$id, levels = mu_df$id)
		}

		if (is.null(labels)) {
			labels <- n_actors <= 50
		}

		p <- ggplot(mu_df, aes(x = id, y = mu)) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			geom_segment(aes(xend = id, yend = 0)) +
			geom_point(size = 2) +
			xlab("") + ylab(ylabel) +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top"
			)

		if (!labels) {
			p <- p + theme(
				axis.text.x = element_blank()
			) + xlab("Actors")
		} else {
			p <- p + theme(
				axis.text.x = element_text(angle = 45, hjust = 1)
			)
		}

		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "at", time_labels[time_point]))
		}

		return(p)


	} else if (plot_type == "trajectory") {
		traj_data <- data.frame(
			actor = rep(actor_names, times = n_times),
			time = rep(seq_len(n_times), each = n_actors),
			time_label = rep(time_labels, each = n_actors),
			effect = as.vector(effects_mat),
			stringsAsFactors = FALSE
		)

		if (!is.null(show_actors)) {
			traj_data <- traj_data[traj_data$actor %in% show_actors, ]
		} else if (n_actors > 20) {
			# with many actors, show the top and bottom five
			avg_effects <- rowMeans(effects_mat)
			top_actors <- names(sort(avg_effects, decreasing = TRUE)[1:5])
			bottom_actors <- names(sort(avg_effects)[1:5])
			show_actors <- c(top_actors, bottom_actors)
			traj_data <- traj_data[traj_data$actor %in% show_actors, ]

			cli::cli_inform(c(
				"i" = "Showing top 5 and bottom 5 actors by average effect",
				">" = "Use {.arg show_actors} to specify actors to display"
			))
		}

		p <- ggplot(traj_data, aes(x = time, y = effect, group = actor, color = actor)) +
			geom_line(alpha = 0.7, linewidth = 1) +
			geom_point(size = 2) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			scale_x_continuous(breaks = 1:n_times, labels = time_labels) +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top"
			) +
			labs(x = "Time", y = ylabel, color = "Actor")

		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste("Trajectories of", tolower(default_title)))
		}

		return(p)


	} else if (plot_type == "faceted") {
		if (n_times <= 6) {
			show_times <- 1:n_times
		} else {
			# six evenly spaced time points
			show_times <- round(seq(1, n_times, length.out = 6))
		}

		sub <- effects_mat[, show_times, drop = FALSE]
		facet_data <- data.frame(
			actor = rep(actor_names, times = length(show_times)),
			# chronological panel order (facet_wrap would sort the labels alphabetically)
			time = factor(
				rep(time_labels[show_times], each = n_actors),
				levels = time_labels[show_times]
			),
			effect = as.vector(sub),
			stringsAsFactors = FALSE
		)

		p <- ggplot(facet_data, aes(x = reorder(actor, effect), y = effect)) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			geom_segment(aes(xend = actor, yend = 0)) +
			geom_point(size = 1.5) +
			facet_wrap(~time, scales = "free_x") +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text = element_text(color = "white", hjust = 0),
				axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
			) +
			labs(x = "", y = ylabel)

		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "over time"))
		}

		return(p)


	} else if (plot_type == "ribbon") {
		# per-period posterior sd stored by lame() alongside the a_dynamic / b_dynamic means
		sd_mat <- if (effect == "sender") fit$a_dynamic_sd else fit$b_dynamic_sd
		if (is.null(sd_mat)) {
			cli::cli_abort(c(
				"This fit has no per-time uncertainty for the additive effects.",
				"i" = "Ribbon bands need {.code fit$a_dynamic_sd} / {.code fit$b_dynamic_sd}, which this fit does not carry.",
				">" = "Re-fit with the current version, or use {.code plot_type = \"trajectory\"}."
			))
		}

		# 95% normal band, mean +/- 1.96 per-period posterior sd
		ribbon_data <- data.frame(
			actor = rep(actor_names, times = n_times),
			time = rep(seq_len(n_times), each = n_actors),
			time_label = rep(time_labels, each = n_actors),
			mean_effect = as.vector(effects_mat),
			lower = as.vector(effects_mat - 1.96 * sd_mat),
			upper = as.vector(effects_mat + 1.96 * sd_mat),
			stringsAsFactors = FALSE
		)

		if (!is.null(show_actors)) {
			ribbon_data <- ribbon_data[ribbon_data$actor %in% show_actors, ]
		} else if (n_actors > 10) {
			avg_effects <- rowMeans(effects_mat)
			show_actors <- names(sort(abs(avg_effects), decreasing = TRUE)[1:10])
			ribbon_data <- ribbon_data[ribbon_data$actor %in% show_actors, ]

			cli::cli_inform(c(
				"i" = "Showing 10 actors with largest absolute average effects",
				">" = "Use {.arg show_actors} to specify actors to display"
			))
		}

		p <- ggplot(ribbon_data, aes(x = time, y = mean_effect, group = actor)) +
			geom_ribbon(aes(ymin = lower, ymax = upper, fill = actor), alpha = 0.3) +
			geom_line(aes(color = actor), linewidth = 1) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			scale_x_continuous(breaks = seq_len(n_times), labels = time_labels) +
			theme_bw() +
			theme(
				panel.border = element_blank(),
				plot.caption = element_text(hjust = 0),
				axis.ticks = element_blank(),
				legend.position = "top"
			) +
			labs(x = "Time", y = ylabel, color = "Actor", fill = "Actor")

		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "with 95% credible intervals"))
		}

		return(p)
	}
}
