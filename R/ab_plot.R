####

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
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
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
ab_plot <- function(fit,
										effect = c("sender", "receiver", "both"),
										sorted = TRUE,
										labels = NULL,
										title = NULL,
										time_point = NULL,
										plot_type = c("snapshot", "trajectory", "faceted", "ribbon"),
										show_actors = NULL) {

	if (!inherits(fit, c("ame", "lame", "ame_als"))) {
		cli::cli_abort("fit must be an object of class 'ame', 'lame', or 'ame_als'")
	}
	
	effect <- match.arg(effect)
	plot_type <- match.arg(plot_type)

	is_dynamic <- !is.null(fit$a_dynamic) || !is.null(fit$b_dynamic)

	# static ALS fits estimate additive effects but have no posterior; delegate
	# to the ame_als-specific lollipop chart, which facets "both" itself.
	# dynamic ALS fits fall through so plot_type and show_actors are honored.
	if (inherits(fit, "ame_als") && !is_dynamic) {
		return(ab_plot.ame_als(fit, effect = effect))
	}

	# "both" facets the two margins side by side. for dynamic fits only the
	# snapshot views (one period, or the time average) have a single value
	# per actor and margin; the trajectory / faceted / ribbon views already
	# spend their panels and colours on time and actor.
	if (effect == "both") {
		if (is_dynamic && (plot_type != "snapshot" || identical(time_point, "all"))) {
			view <- if (identical(time_point, "all")) {
				'time_point = "all"'
			} else {
				paste0('plot_type = "', plot_type, '"')
			}
			cli::cli_abort(c(
				"{.code effect = \"both\"} is only available for static fits and dynamic snapshots (one period, or {.code time_point = \"average\"}).",
				"i" = "Plot {.val sender} and {.val receiver} separately with {.code {view}}."))
		}
		return(ab_plot_both_internal(fit, sorted, labels, title, time_point,
		                             is_dynamic))
	}

	if (is_dynamic) {
		return(ab_plot_dynamic_internal(fit, effect, sorted, labels, title,
		                                time_point, plot_type, show_actors))
	}

	####

	# static single-margin plot
	m <- .ab_margin_values(fit, effect)
	muDf <- data.frame(mu = unname(m$mu), id = names(m$mu),
	                   stringsAsFactors = FALSE)
	.ab_lollipop(muDf, ylabel = m$ylabel, sorted = sorted, labels = labels,
	             title = title, default_title = m$default_title)
}

####

# one margin's effect vector with its axis label, panel name (the label
# without any time suffix) and default title. static
# fits read the posterior means (apm / bpm); dynamic fits read one column
# of a_dynamic / b_dynamic (default: the last period) or the row means when
# time_point = "average".
.ab_margin_values <- function(fit, margin, time_point = NULL,
                              is_dynamic = FALSE) {
	is_bip <- identical(fit$mode, "bipartite")
	if (margin == "sender") {
		what <- if (is_bip) "Row Actor" else "Sender"
		ylabel <- paste0(what, " Effects (a)")
		eff <- if (is_dynamic) fit$a_dynamic else fit$APM
	} else {
		what <- if (is_bip) "Column Actor" else "Receiver"
		ylabel <- paste0(what, " Effects (b)")
		eff <- if (is_dynamic) fit$b_dynamic else fit$BPM
	}
	if (is.null(eff)) {
		cli::cli_abort(c(
			"This fit has no receiver (b) effects.",
			"i" = "Symmetric models store a single set of actor effects; use {.code effect = \"sender\"}."))
	}
	panel <- ylabel
	suffix <- ""
	if (is_dynamic) {
		n_times <- ncol(eff)
		actor_names <- rownames(eff) %||% paste0("Actor", seq_len(nrow(eff)))
		time_labels <- colnames(eff) %||% paste0("T", seq_len(n_times))
		if (identical(time_point, "average")) {
			eff <- rowMeans(eff)
			ylabel <- paste(ylabel, "(Time-Averaged)")
			suffix <- " - Average Across Time"
		} else {
			if (is.null(time_point)) time_point <- n_times
			if (!is.numeric(time_point)) {
				cli::cli_abort(c(
					"{.arg time_point} must be a number, {.val average}, or {.val all}.",
					"i" = "Got {.val {time_point}}."))
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
	list(mu = eff, ylabel = ylabel, panel = panel, suffix = suffix,
	     default_title = paste0(if (is_dynamic) "Dynamic " else "Additive ",
	                            what, " Effects", suffix))
}

# lollipop chart of one effect per actor. mudf has columns mu and id, plus
# margin (a factor) when the sender and receiver panels are drawn together;
# each panel then sorts on its own.
.ab_lollipop <- function(muDf, ylabel, sorted, labels, title, default_title) {
	facet <- "margin" %in% names(muDf)
	n_actors <- if (facet) max(table(muDf$margin)) else nrow(muDf)
	if (is.null(labels)) {
		labels <- n_actors <= 50
	}

	# x is a per-panel key rather than the actor name so that each panel
	# orders its own actors; the axis maps the keys back to the names
	muDf$key <- if (facet) paste(muDf$margin, muDf$id, sep = "\r") else muDf$id
	ord <- if (facet) {
		if (sorted) order(muDf$margin, muDf$mu) else order(muDf$margin)
	} else {
		if (sorted) order(muDf$mu) else seq_len(nrow(muDf))
	}
	muDf$key <- factor(muDf$key, levels = muDf$key[ord])
	key_labels <- stats::setNames(muDf$id, as.character(muDf$key))

	p <- ggplot(muDf, aes(x = .data$key, y = .data$mu)) +
		geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
		geom_segment(aes(xend = .data$key, yend = 0)) +
		geom_point(size = 2) +
		scale_x_discrete(labels = key_labels) +
		xlab("") +
		ylab(ylabel) +
		theme_bw() +
		theme(
			panel.border     = element_blank(),
			axis.ticks       = element_blank(),
			legend.position  = "top"
		)

	if (facet) {
		p <- p + facet_wrap(~ margin, scales = "free_x") +
			theme(
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0)
			)
	}

	if (!labels) {
		p <- p + theme(
			axis.text.x = element_blank()
		) + xlab(if (sorted) "Actors (Sorted by Effect Magnitude)" else "Actors")
	} else {
		# scale label angle and size with n so n > 15 stays readable.
		ang  <- if (n_actors > 30L) 90 else if (n_actors > 15L) 75 else 45
		size <- if (n_actors > 30L) 7  else if (n_actors > 15L) 8  else 10
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
ab_plot_both_internal <- function(fit, sorted, labels, title, time_point,
                                  is_dynamic) {
	a <- .ab_margin_values(fit, "sender", time_point, is_dynamic)
	b <- .ab_margin_values(fit, "receiver", time_point, is_dynamic)
	panels <- c(a$panel, b$panel)
	muDf <- rbind(
		data.frame(mu = unname(a$mu), id = names(a$mu), margin = panels[1L],
		           stringsAsFactors = FALSE),
		data.frame(mu = unname(b$mu), id = names(b$mu), margin = panels[2L],
		           stringsAsFactors = FALSE))
	muDf$margin <- factor(muDf$margin, levels = panels)
	is_bip <- identical(fit$mode, "bipartite")
	default_title <- paste0(
		if (is_dynamic) "Dynamic " else "Additive ",
		if (is_bip) "Row and Column Actor" else "Sender and Receiver",
		" Effects", a$suffix)
	ylabel <- if (identical(time_point, "average")) {
		"Additive Effect (Time-Averaged)"
	} else {
		"Additive Effect"
	}
	.ab_lollipop(muDf, ylabel = ylabel, sorted = sorted, labels = labels,
	             title = title, default_title = default_title)
}

####

#' Internal function for dynamic additive effects plotting
#' @keywords internal
#' @noRd
ab_plot_dynamic_internal <- function(fit, effect, sorted, labels, title,
																		 time_point, plot_type, show_actors) {
	
	if (effect == "sender") {
		effects_mat <- fit$a_dynamic
		ylabel <- "Sender Effects (a)"
		default_title <- "Dynamic Sender Effects"
	} else {
		effects_mat <- fit$b_dynamic
		ylabel <- "Receiver Effects (b)"
		default_title <- "Dynamic Receiver Effects"
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
			muEff <- rowMeans(effects_mat)
			names(muEff) <- actor_names

			muDf <- data.frame(
				mu = muEff,
				id = actor_names,
				stringsAsFactors = FALSE
			)
			
			if (sorted) {
				muDf$id <- factor(muDf$id, levels = muDf$id[order(muDf$mu)])
			} else {
				muDf$id <- factor(muDf$id, levels = muDf$id)
			}
			
			if (is.null(labels)) {
				labels <- n_actors <= 50
			}
			
			p <- ggplot(muDf, aes(x = id, y = mu)) +
				geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
				geom_segment(aes(xend = id, yend = 0)) +
				geom_point(size = 2) +
				xlab("") + ylab(paste(ylabel, "(Time-Averaged)")) +
				theme_bw() +
				theme(
					panel.border    = element_blank(),
					axis.ticks      = element_blank(),
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
				p <- p + ggtitle(paste(default_title, "- Average Across Time"))
			}
			
			return(p)
		} else if (identical(time_point, "all")) {
			return(ab_plot_dynamic_internal(fit, effect, sorted, labels, title,
			                                NULL, "faceted", show_actors))
		}

		if (!is.numeric(time_point)) {
			cli::cli_abort(c(
				"{.arg time_point} must be a number, {.val average}, or {.val all}.",
				"i" = "Got {.val {time_point}}."))
		}
		if (time_point > n_times || time_point < 1) {
			cli::cli_abort("{.arg time_point} must be between 1 and {n_times}.")
		}
		
		muEff <- effects_mat[, time_point]
		names(muEff) <- actor_names
		
		muDf <- data.frame(
			mu = muEff,
			id = actor_names,
			stringsAsFactors = FALSE
		)
		
		if (sorted) {
			muDf$id <- factor(muDf$id, levels = muDf$id[order(muDf$mu)])
		} else {
			muDf$id <- factor(muDf$id, levels = muDf$id)
		}
		
		# default label visibility
		if (is.null(labels)) {
			labels <- n_actors <= 50
		}
		
		p <- ggplot(muDf, aes(x = id, y = mu)) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			geom_segment(aes(xend = id, yend = 0)) +
			geom_point(size = 2) +
			xlab("") + ylab(ylabel) +
			theme_bw() +
			theme(
				panel.border    = element_blank(),
				axis.ticks      = element_blank(),
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

	####

	} else if (plot_type == "trajectory") {
		traj_data <- data.frame(
			actor      = rep(actor_names, times = n_times),
			time       = rep(seq_len(n_times), each = n_actors),
			time_label = rep(time_labels, each = n_actors),
			effect     = as.vector(effects_mat),
			stringsAsFactors = FALSE
		)
		
		if (!is.null(show_actors)) {
			traj_data <- traj_data[traj_data$actor %in% show_actors, ]
		} else if (n_actors > 20) {
			# too many actors, show top and bottom 5
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
				panel.border    = element_blank(),
				axis.ticks      = element_blank(),
				legend.position = "top"
			) +
			labs(x = "Time", y = ylabel, color = "Actor")
		
		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "Trajectories"))
		}
		
		return(p)

	####

	} else if (plot_type == "faceted") {
		if (n_times <= 6) {
			show_times <- 1:n_times
		} else {
			# select 6 evenly spaced time points
			show_times <- round(seq(1, n_times, length.out = 6))
		}
		
		sub <- effects_mat[, show_times, drop = FALSE]
		facet_data <- data.frame(
			actor  = rep(actor_names, times = length(show_times)),
			# keep panels in chronological order; facet_wrap would otherwise
			# sort the labels alphabetically
			time   = factor(rep(time_labels[show_times], each = n_actors),
			                levels = time_labels[show_times]),
			effect = as.vector(sub),
			stringsAsFactors = FALSE
		)
		
		p <- ggplot(facet_data, aes(x = reorder(actor, effect), y = effect)) +
			geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
			geom_segment(aes(xend = actor, yend = 0)) +
			geom_point(size = 1.5) +
			facet_wrap(~ time, scales = "free_x") +
			theme_bw() +
			theme(
				panel.border     = element_blank(),
				axis.ticks       = element_blank(),
				legend.position  = "top",
				strip.background = element_rect(fill = "black", color = "black"),
				strip.text       = element_text(color = "white", hjust = 0),
				axis.text.x      = element_text(angle = 45, hjust = 1, size = 6)
			) +
			labs(x = "", y = ylabel)
		
		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "Over Time"))
		}
		
		return(p)

	####

	} else if (plot_type == "ribbon") {
		# per-time posterior sd, stored by lame() alongside the a_dynamic /
		# b_dynamic means. older fits (pre-1.3.3) do not carry it.
		sd_mat <- if (effect == "sender") fit$a_dynamic_sd else fit$b_dynamic_sd
		if (is.null(sd_mat)) {
			cli::cli_abort(c(
				"This fit has no per-time uncertainty for the additive effects.",
				"i" = "Ribbon bands need {.code fit$a_dynamic_sd} / {.code fit$b_dynamic_sd}, added in lame 1.3.3.",
				">" = "Re-fit with the current version, or use {.code plot_type = \"trajectory\"}."))
		}

		# 95% normal-approximation band, mean +/- 1.96 * per-time posterior sd
		ribbon_data <- data.frame(
			actor      = rep(actor_names, times = n_times),
			time       = rep(seq_len(n_times), each = n_actors),
			time_label = rep(time_labels, each = n_actors),
			mean_effect = as.vector(effects_mat),
			lower      = as.vector(effects_mat - 1.96 * sd_mat),
			upper      = as.vector(effects_mat + 1.96 * sd_mat),
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
				panel.border    = element_blank(),
				axis.ticks      = element_blank(),
				legend.position = "top"
			) +
			labs(x = "Time", y = ylabel, color = "Actor", fill = "Actor")

		if (!is.null(title)) {
			p <- p + ggtitle(title)
		} else {
			p <- p + ggtitle(paste(default_title, "with 95% Credible Intervals"))
		}

		return(p)
	}
}
