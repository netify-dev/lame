#' Visualize multiplicative effects (latent factors) from AME models
#'
#' Visualizes the multiplicative effects (U and V) of an AME or LAME model:
#' maps of two latent dimensions at a time, or a heatmap that uses every
#' fitted dimension. These latent factors capture network structure beyond
#' what is explained by covariates and additive effects, including
#' clustering, transitivity, and other higher-order dependencies.
#'
#' @details
#' The multiplicative effects in AME models provide a low-rank representation
#' of network structure through latent factors:
#' \describe{
#'   \item{U matrix}{Sender-specific latent positions (row factors)}
#'   \item{V matrix}{Receiver-specific latent positions (column factors)}
#'   \item{UV' product}{Captures dyad-specific effects beyond additive terms}
#' }
#'
#' The visualization can show:
#' \describe{
#'   \item{Circular layout}{Default layout placing nodes on a circle with
#'         latent positions shown as deviations}
#'   \item{Biplot layout}{Shows U and V positions directly in latent space}
#'   \item{Heatmap layout}{Uses every fitted dimension. With
#'         \code{nodes = "row"} or \code{"col"} it shows the similarity of
#'         the nodes of that side (see \code{\link{latent_similarity}}),
#'         ordered by clustering; with \code{nodes = "both"} it shows the
#'         latent product, rows against columns. The latent product is
#'         centred, so positions are relative to the average node of each
#'         side: a node with no structure of its own sits opposite the average
#'         pattern rather than at the origin, and several such nodes look
#'         similar to each other.}
#'   \item{Network overlay}{Optional display of actual network ties}
#' }
#'
#' A map shows two dimensions. When the fit has more than two, the axis labels
#' give each dimension's share of the latent product and a caption names the
#' dimensions shown; three or four values of \code{dims} draw every pair of
#' them as panels of one plot. When a side has more than 50 nodes and
#' \code{show.usernames} is not set, only its 20 nodes with the longest
#' position vectors in the plotted dimensions are labelled and the rest are
#' drawn as small light points. The biplot of a bipartite fit rescales the
#' row positions for display so that their root-mean-square radius equals
#' that of the column positions, and says so in the caption.
#'
#' Interpretation (this is an eigenmodel, not a latent-distance model):
#' - The latent contribution to a tie i -> j is the inner product u_i' v_j:
#'   sender and receiver vectors that align (small angle, large magnitude) raise
#'   the tie propensity, opposed vectors lower it. Euclidean distance between
#'   plotted points is not the operative quantity.
#' - Actors with similar sender (or receiver) vectors have similar connection
#'   patterns; clustering of vectors indicates community structure.
#'
#' @param fit An object of class "ame" or "lame" containing multiplicative effects,
#'        or a network matrix Y if U and V are provided separately
#' @param Y Network matrix (only needed if fit is not provided)
#' @param U Matrix of sender latent factors (extracted from fit if not provided)
#' @param V Matrix of receiver latent factors (extracted from fit if not provided)
#' @param row.names Names for row nodes (defaults to rownames of Y or U)
#' @param col.names Names for column nodes (defaults to colnames of Y or V)
#' @param layout Character string specifying layout: \code{"circle"}
#'        (default), \code{"biplot"} or \code{"heatmap"}
#' @param dims Which latent dimensions to draw in the circle and biplot
#'        layouts: distinct whole numbers within the fitted rank. The default
#'        is the first two (the single one when the fitted rank is 1). Two
#'        values draw one map; three or four values draw every pair of them as
#'        panels. Dynamic fits take exactly two values. The heatmap uses every
#'        dimension and takes no \code{dims}.
#' @param nodes Which side to draw: \code{"both"} (default), \code{"row"}
#'        (row actors; senders in a one-mode network) or \code{"col"} (column
#'        actors; receivers). A symmetric fit has one set of positions: the maps
#'        take only \code{"both"}, and the heatmap takes \code{"row"} or
#'        \code{"col"} (the same similarity of its nodes) or \code{"both"}.
#' @param vscale Scaling factor for V positions relative to U (default 0.8)
#' @param show.edges Logical; if TRUE, show network edges (default FALSE)
#' @param edge.alpha Transparency for edges (default 0.3)
#' @param node.size Node size: a numeric value, \code{"degree"} to scale by degree (needs \code{Y}), or \code{"magnitude"} (the default) to scale by the length of each latent vector
#' @param label.nodes Logical; if TRUE, show node labels (default TRUE)
#' @param label.size Size of node labels (default 3)
#' @param show.usernames Which nodes to label: an integer (that many nodes, ranked by latent magnitude when available, else by degree), \code{TRUE} for all, a character vector of actor names, or \code{NULL} for the default
#' @param sender.color Color for sender/row nodes (default "darkred")
#' @param receiver.color Color for receiver/column nodes (default "darkblue")
#' @param colors Optional vector of colors for nodes (e.g., for communities)
#' @param title Optional title for the plot
#' @param time_point For dynamic UV, which time point to plot (default: last).
#'                   Can be numeric index or "average" for time-averaged
#'                   positions (not for the heatmap, which draws one period)
#' @param plot_type For dynamic UV: "snapshot" (single time), "trajectory" (evolution),
#'                  "faceted" (grid of time points). For static UV, this is ignored.
#' @param show_arrows For trajectory plots, whether to show directional arrows
#' @param highlight Optional character vector of actor names to highlight on a
#'        \code{plot_type = "trajectory"} plot. Highlighted actors are coloured
#'        with the ColorBrewer Dark2 palette; all other actors are
#'        rendered in grey at lower alpha. Ignored for static or snapshot plots.
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{latent_similarity}}, \code{\link{latent_positions}}
#' @examples
#' \donttest{
#' # Fit an AME model with multiplicative effects
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Basic visualization
#' uv_plot(fit)
#'
#' # Use biplot layout
#' uv_plot(fit, layout = "biplot")
#'
#' # three dimensions: every pair as a panel, and the heatmap of all of them
#' fit3 <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 3,
#'             nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#' uv_plot(fit3, layout = "biplot", dims = 1:3)
#' uv_plot(fit3, layout = "heatmap", nodes = "row")
#' }
#' @export
#' @import cli
uv_plot <- function(
	fit = NULL,
	Y = NULL,
	U = NULL,
	V = NULL,
	row.names = NULL,
	col.names = NULL,
	layout = c("circle", "biplot", "heatmap"),
	dims = NULL,
	nodes = c("both", "row", "col"),
	vscale = 0.8,
	show.edges = FALSE,
	edge.alpha = 0.3,
	node.size = "magnitude",
	label.nodes = TRUE,
	label.size = 3,
	show.usernames = NULL,
	sender.color = "darkred",
	receiver.color = "darkblue",
	colors = NULL,
	title = NULL,
	time_point = NULL,
	plot_type = c("snapshot", "trajectory", "faceted"),
	show_arrows = TRUE,
	highlight = NULL
) {
	.lame_require_draws(fit, "uv_plot")

	layout <- match.arg(layout)
	nodes <- match.arg(nodes)
	plot_type <- match.arg(plot_type)
	if (layout == "heatmap" && !is.null(dims)) {
		cli::cli_abort(c(
			"{.code layout = \"heatmap\"} uses every fitted dimension, so it takes no {.arg dims}.",
			"i" = "Drop {.arg dims}, or use {.code layout = \"circle\"} or {.code \"biplot\"} to draw chosen dimensions."
		))
	}
	user_uv <- !is.null(U) || !is.null(V)
	is_sym <- FALSE
	is_dynamic <- FALSE
	st <- NULL
	rank_txt <- NULL

	if (!is.null(fit)) {
		if (!inherits(fit, c("ame", "lame", "ame_als"))) {
			cli::cli_abort("fit must be an object of class 'ame', 'lame' or 'ame_als'")
		}
		# an ALS fit with a withheld latent product has no positions to plot
		if (is.null(U) && .ae_product_is_withheld(fit)) {
			pw <- fit$eb$product_withheld
			cli::cli_abort(
				c(
					"This fit has no latent positions to plot: the latent product is not identified by the data at rank {pw$rank} ({signif(pw$per_loading, 2)} noise variances of signal per loading, below {pw$threshold}).",
					"i" = "Compare ranks with {.code fit$eb$product_withheld$log_marginal}, or fit {.fn ame} / {.fn lame} for the posterior of the positions."
				),
				class = "lame_product_withheld"
			)
		}
		# positions factored from the posterior mean products
		fit <- .lame_positions(fit)
		is_sym <- isTRUE(fit$symmetric)
		# the similarity heatmap is defined for symmetric fits; the maps are not
		if (is_sym && nodes != "both" && layout != "heatmap") .uv_abort_symmetric_nodes(nodes)
		if (is.null(U)) U <- fit$U
		if (is.null(V)) V <- fit[["V"]]
		# symmetric fits store U and L only; plot one latent space with V = U
		if (is_sym && is.null(V) && !is.null(U)) V <- U
		if (is.null(U)) {
			cli::cli_abort(c(
				"This fit has no latent positions to plot.",
				"i" = "It was fitted with {.code R = 0}, or with a one-sided bipartite rank, which drops the multiplicative term."
			))
		}
		if (layout == "heatmap" && NCOL(U) > 0L) {
			h_fit <- if (user_uv) .uv_matrix_fit(U, V, row.names, col.names) else fit
			return(.uv_heatmap(h_fit, nodes, time_point, plot_type, title))
		}
		rank_txt <- .lame_rank_text(fit)
		if (length(dim(U)) == 3) {
			is_dynamic <- TRUE
			k <- .uv_rank(U, V, is_sym)
			dims <- .uv_check_dims(dims, k, rank_txt, dynamic = TRUE)
			# snapshots use the static layout below; trajectory and faceted plots
			# use the dynamic renderer
			if (plot_type != "snapshot") {
				return(uv_plot_dynamic_internal(
					fit,
					U,
					V,
					time_point,
					plot_type,
					show_arrows,
					title,
					label.nodes,
					label.size,
					colors,
					highlight = highlight,
					dims = dims,
					nodes = nodes,
					k = k
				))
			}
			n_times <- dim(U)[3]
			if (is.null(time_point)) {
				time_point <- n_times
				cli::cli_inform(c(
					"i" = "Detected dynamic UV effects",
					"*" = "Using last time point for visualization"
				))
			}
			if (time_point == "average") {
				U <- apply(U, c(1, 2), mean)
				if (!is.null(V) && length(dim(V)) == 3) V <- apply(V, c(1, 2), mean)
				cli::cli_inform("Using time-averaged positions")
			} else {
				# keep the n x R shape and actor names when slicing a period
				dU <- dim(U)
				U <- matrix(
					U[, , time_point],
					nrow = dU[1],
					ncol = dU[2],
					dimnames = list(dimnames(U)[[1]], NULL)
				)
				if (!is.null(V) && length(dim(V)) == 3) {
					dV <- dim(V)
					V <- matrix(
						V[, , time_point],
						nrow = dV[1],
						ncol = dV[2],
						dimnames = list(dimnames(V)[[1]], NULL)
					)
				}
			}
		} else if (!user_uv) {
			st <- .lame_dim_strength(fit)
		}
		if (is.null(Y)) {
			if (!is.null(fit$Y)) {
				Y <- fit$Y
			} else if (!is.null(fit$YPM)) {
				Y <- fit$YPM
			}
			if (is.list(Y)) {
				if (length(Y) > 0) {
					Y_avg <- Y[[1]]
					if (length(Y) > 1) {
						for (t in 2:length(Y)) {
							Y_avg <- Y_avg + Y[[t]]
						}
						Y_avg <- Y_avg / length(Y)
					}
					Y <- Y_avg
				}
			} else if (length(dim(Y)) == 3) {
				Y <- apply(Y, c(1, 2), mean, na.rm = TRUE)
			}
		}
	}


	# null or zero-size factors mean no latent space
	.uv_empty <- function(M) is.null(M) || NCOL(M) == 0L || NROW(M) == 0L
	if (.uv_empty(U) && .uv_empty(V)) {
		if (!is.null(fit)) {
			# an R = 0 fit has no latent space, so no SVD of the data is used
			cli::cli_abort(paste0(
				"This fit has no multiplicative effects (R = 0); there is no ",
				"estimated latent space to plot. Refit with R > 0, or use ",
				"ab_plot() for the additive effects."
			))
		}
		if (is.null(Y)) {
			cli::cli_abort("Must provide either a fit object with U and V, or a Y matrix to compute factors")
		}
		Y_centered <- Y - mean(Y, na.rm = TRUE)
		Y_centered[is.na(Y_centered)] <- 0
		svd_Y <- svd(Y_centered)
		U <- svd_Y$u
		V <- svd_Y$v * vscale
		if (is.null(rownames(U))) rownames(U) <- rownames(Y)
		if (is.null(rownames(V))) rownames(V) <- colnames(Y)
	}
	if (is.null(fit) && layout == "heatmap") {
		return(.uv_heatmap(.uv_matrix_fit(U, V, row.names, col.names), nodes, NULL, plot_type, title))
	}

	k <- .uv_rank(U, V, is_sym)
	if (!is_dynamic) dims <- .uv_check_dims(dims, k, rank_txt, dynamic = FALSE)

	n_row <- nrow(U)
	n_col <- nrow(V)
	is_bip <- !is.null(fit) && identical(fit$mode, "bipartite")
	row_type_label <- if (is_bip) "Row actor" else "Sender"
	col_type_label <- if (is_bip) "Column actor" else "Receiver"
	if (is.null(row.names)) {
		row.names <- rownames(U)
		if (is.null(row.names)) row.names <- paste0(if (is_bip) "R" else "S", 1:n_row)
	}
	if (is.null(col.names)) {
		col.names <- rownames(V)
		if (is.null(col.names)) col.names <- paste0(if (is_bip) "C" else "R", 1:n_col)
	}
	show_row <- nodes != "col"
	show_col <- nodes != "row"
	if (show.edges && !is.null(Y) && nodes != "both") {
		cli::cli_abort(c(
			"Edges join row and column nodes, so {.arg show.edges} needs both sides.",
			"i" = "Use {.code nodes = \"both\"}, or drop {.arg show.edges}."
		))
	}

	# one panel per pair of the shown dimensions
	panels <- if (length(dims) <= 2L) list(dims) else utils::combn(dims, 2L, simplify = FALSE)
	multi <- length(panels) > 1L
	two_cols <- function(A, pr) {
		if (length(pr) == 1L) cbind(A[, pr, drop = FALSE], 0) else A[, pr, drop = FALSE]
	}

	# biplot coordinates: each dimension scaled by the larger sd of its two sides
	rescaled <- FALSE
	if (layout == "biplot") {
		U_s <- U[, dims, drop = FALSE]
		V_s <- V[, dims, drop = FALSE] * vscale
		V_r <- V[, dims, drop = FALSE]
		for (j in seq_along(dims)) {
			d <- dims[j]
			u_sd <- sd(U[, d])
			v_sd <- sd(V[, d])
			if (u_sd > 1e-10 && v_sd > 1e-10) {
				scale_factor <- max(u_sd, v_sd)
				U_s[, j] <- U[, d] / scale_factor
				V_s[, j] <- V[, d] * vscale / scale_factor
				V_r[, j] <- V[, d] / scale_factor
			}
		}
		# bipartite sides differ in size, so the row cloud is rescaled for display
		if (is_bip && show_row && show_col) {
			rms_u <- sqrt(mean(rowSums(U_s^2)))
			rms_v <- sqrt(mean(rowSums(V_r^2)))
			if (rms_u > 1e-10 && rms_v > 1e-10) {
				U_s <- U_s * (rms_v / rms_u)
				rescaled <- TRUE
			}
		}
	}

	# per-node attributes shared by all panels
	keep <- c(rep(show_row, n_row), rep(show_col, n_col))
	side <- c(rep("row", n_row), rep("col", n_col))
	node_base <- data.frame(
		name = c(row.names, col.names),
		type = factor(c(rep(row_type_label, n_row), rep(col_type_label, n_col))),
		label = "",
		stringsAsFactors = FALSE
	)
	strength <- c(
		sqrt(rowSums(U[, dims, drop = FALSE]^2)),
		sqrt(rowSums(V[, dims, drop = FALSE]^2))
	)
	degree <- if (!is.null(Y)) c(rowSums(Y > 0, na.rm = TRUE), colSums(Y > 0, na.rm = TRUE))
	node_base <- node_base[keep, , drop = FALSE]
	strength <- strength[keep]
	side <- side[keep]
	if (!is.null(degree)) degree <- degree[keep]

	if (!is.null(colors)) {
		if (length(colors) != nrow(node_base)) {
			cli::cli_warn("Length of colors doesn't match number of nodes, ignoring colors")
			colors <- NULL
		} else {
			node_base$color <- colors
		}
	}

	# a side over 50 nodes shows its 20 strongest; the rest are light points
	side_n <- table(factor(side, levels = c("row", "col")))
	large <- names(side_n)[side_n > 50]
	faint <- rep(FALSE, nrow(node_base))
	auto_large <- is.null(show.usernames) && length(large) > 0L
	if (auto_large) {
		for (s in large) {
			idx <- which(side == s)
			top <- idx[order(strength[idx], decreasing = TRUE)[seq_len(20L)]]
			faint[setdiff(idx, top)] <- TRUE
		}
	}

	if (label.nodes || !is.null(show.usernames)) {
		if (!is.null(show.usernames)) {
			if (is.logical(show.usernames) && show.usernames) {
				node_base$label <- node_base$name
			} else if (is.character(show.usernames)) {
				node_base$label <- ifelse(
					node_base$name %in% show.usernames,
					node_base$name,
					""
				)
			} else if (is.numeric(show.usernames)) {
				n_top <- min(show.usernames, nrow(node_base))
				if (layout == "circle") {
					top_indices <- order(strength, decreasing = TRUE)[1:n_top]
				} else if (!is.null(degree)) {
					top_indices <- order(degree, decreasing = TRUE)[1:n_top]
				} else {
					top_indices <- 1:n_top
				}
				node_base$label <- ifelse(
					seq_len(nrow(node_base)) %in% top_indices,
					node_base$name,
					""
				)
			}
		} else if (auto_large) {
			node_base$label <- ifelse(faint, "", node_base$name)
		} else if (label.nodes) {
			if (nrow(node_base) <= 50) {
				node_base$label <- node_base$name
			} else {
				if (!is.null(degree)) {
					top_nodes <- degree >= quantile(degree, 0.9)
					node_base$label <- ifelse(top_nodes, node_base$name, "")
				} else {
					node_base$label <- ""
				}
			}
		}
	}

	# coordinates of every panel
	panel_names <- vapply(panels, function(pr) .uv_panel_name(pr, st), character(1))
	node_data <- do.call(rbind, lapply(seq_along(panels), function(i) {
		pr <- panels[[i]]
		if (layout == "circle") {
			xy <- .uv_circle_xy(two_cols(U, pr), two_cols(V, pr), vscale)
			pos <- rbind(xy$u, xy$v)
			mag <- c(xy$mu, xy$mv)
		} else {
			cols <- match(pr, dims)
			pos <- rbind(two_cols(U_s, cols), two_cols(V_s, cols))
			mag <- NULL
		}
		out <- data.frame(
			x = pos[keep, 1],
			y = pos[keep, 2],
			node_base,
			stringsAsFactors = FALSE
		)
		if (!is.null(mag)) out$magnitude <- mag[keep]
		if (multi) out$panel <- factor(panel_names[i], levels = panel_names)
		out
	}))
	rownames(node_data) <- NULL
	faint_all <- rep(faint, length(panels))

	if (is.character(node.size)) {
		if (node.size == "degree" && !is.null(Y)) {
			sz <- degree
			sz <- 1 + 4 * (sz - min(sz)) / (max(sz) - min(sz) + 1)
			node_data$size <- rep(sz, length(panels))
		} else if (node.size == "magnitude" && layout == "circle") {
			# sizes rank the magnitudes within each panel
			grp <- if (multi) node_data$panel else rep(1L, nrow(node_data))
			node_data$size <- stats::ave(
				node_data$magnitude,
				grp,
				FUN = function(m) 1 + 4 * (m - min(m)) / (max(m) - min(m) + 0.01)
			)
		} else {
			node_data$size <- 3
		}
	} else {
		node_data$size <- node.size
	}


	p <- ggplot2::ggplot(node_data, ggplot2::aes(x = x, y = y))

	if (layout == "circle") {
		circle_df <- data.frame(
			r = c(vscale, 1), # inner circle for receivers, outer for senders
			type = c(paste(col_type_label, "circle"), paste(row_type_label, "circle"))
		)
		circle_df <- circle_df[c(show_col, show_row), , drop = FALSE]
		p <- p +
			ggforce::geom_circle(
				data = circle_df,
				ggplot2::aes(x0 = 0, y0 = 0, r = r),
				color = "gray80",
				fill = NA,
				linetype = "dotted",
				inherit.aes = FALSE
			)
	}

	if (show.edges && !is.null(Y)) {
		ij <- which(!is.na(Y) & Y > 0, arr.ind = TRUE)
		ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
		if (nrow(ij) > 0) {
			n_nodes <- n_row + n_col
			edge_data <- do.call(rbind, lapply(seq_along(panels), function(i) {
				off <- (i - 1L) * n_nodes
				out <- data.frame(
					x = node_data$x[off + ij[, 1]],
					y = node_data$y[off + ij[, 1]],
					xend = node_data$x[off + n_row + ij[, 2]],
					yend = node_data$y[off + n_row + ij[, 2]],
					weight = Y[ij]
				)
				if (multi) out$panel <- factor(panel_names[i], levels = panel_names)
				out
			}))
			p <- p + ggplot2::geom_segment(
				data = edge_data,
				ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
				alpha = edge.alpha,
				color = "gray50"
			)
		}
	}

	point_data <- if (any(faint_all)) node_data[!faint_all, , drop = FALSE]
	if (any(faint_all)) {
		faint_data <- node_data[faint_all, , drop = FALSE]
		p <- p + ggplot2::geom_point(
			data = faint_data,
			mapping = if (!is.null(colors)) {
				ggplot2::aes(color = color, shape = type)
			} else {
				ggplot2::aes(color = type, shape = type)
			},
			size = 0.6,
			alpha = 0.3,
			show.legend = FALSE
		)
	}
	if (!is.null(colors)) {
		p <- p + ggplot2::geom_point(
			data = point_data,
			ggplot2::aes(color = color, shape = type, size = size)
		)
	} else {
		# one legend name so colour and shape merge into a single legend
		p <- p + ggplot2::geom_point(
			data = point_data,
			ggplot2::aes(color = type, shape = type, size = size)
		) +
			ggplot2::scale_color_manual(
				values = setNames(
					c(sender.color, receiver.color),
					c(row_type_label, col_type_label)
				),
				name = "Node type"
			) +
			ggplot2::scale_shape_discrete(name = "Node type")
	}

	if (any(node_data$label != "")) {
		if (is.null(colors)) {
			p <- p + ggrepel::geom_text_repel(
				data = node_data[node_data$label != "", ],
				ggplot2::aes(x = x, y = y, label = label, color = type),
				size = label.size,
				max.overlaps = Inf,
				show.legend = FALSE
			)
		} else {
			p <- p + ggrepel::geom_text_repel(
				data = node_data[node_data$label != "", ],
				ggplot2::aes(x = x, y = y, label = label),
				size = label.size,
				max.overlaps = Inf
			)
		}
	}

	p <- p +
		ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			plot.caption = ggplot2::element_text(hjust = 0),
			axis.ticks = ggplot2::element_blank(),
			axis.text = ggplot2::element_blank(),
			legend.position = "top"
		)
	if (multi) {
		# the strips name the two dimensions of each panel
		p <- p +
			ggplot2::facet_wrap(~panel, ncol = 3L) +
			ggplot2::theme(
				axis.title = ggplot2::element_blank(),
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text = ggplot2::element_text(color = "white", hjust = 0)
			)
	} else if (k > 2L) {
		p <- p + ggplot2::labs(
			x = .uv_dim_label(dims[1], st),
			y = if (length(dims) > 1L) .uv_dim_label(dims[2], st)
		)
		if (length(dims) == 1L) p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
	} else if (identical(layout, "biplot")) {
		# axis titles only for the biplot; circle coordinates are arbitrary
		d2 <- if (length(dims) > 1L) dims[2] else 2L
		p <- p + ggplot2::labs(x = paste("Latent dimension", dims[1]), y = paste("Latent dimension", d2))
	} else {
		p <- p + ggplot2::theme(axis.title = ggplot2::element_blank())
	}
	p <- p + ggplot2::coord_fixed(clip = "off")

	caption <- c(
		.uv_dims_caption(dims, k, st),
		if (rescaled) "Row positions rescaled for display."
	)
	if (length(caption) > 0L) p <- p + ggplot2::labs(caption = .uv_wrap(caption))

	if (!is.null(title)) {
		p <- p + ggplot2::ggtitle(title)
	} else if (layout == "circle") {
		p <- p + ggplot2::ggtitle("Multiplicative effects (circular layout)")
	} else {
		p <- p + ggplot2::ggtitle("Multiplicative effects (biplot)")
	}

	return(p)
}

# circle coordinates of two latent columns: directions on the unit circle for
# rows and the vscale circle for columns, radii jittered by magnitude rank
.uv_circle_xy <- function(U2, V2, vscale) {
	n_row <- nrow(U2)
	n_col <- nrow(V2)
	mu <- sqrt(rowSums(U2^2))
	mv <- sqrt(rowSums(V2^2))
	u_norm <- U2
	v_norm <- V2
	ok_u <- mu > 1e-10
	ok_v <- mv > 1e-10
	u_norm[ok_u, ] <- U2[ok_u, , drop = FALSE] / mu[ok_u]
	v_norm[ok_v, ] <- V2[ok_v, , drop = FALSE] / mv[ok_v]
	v_norm <- v_norm * vscale
	jitter_factor <- 0.15
	u_radii <- 1 + jitter_factor * (rank(mu) / (n_row + 1) - 0.5)
	u_plot <- u_norm * u_radii
	v_radii <- vscale + jitter_factor * vscale * (rank(mv) / (n_col + 1) - 0.5)
	angle <- atan2(v_norm[, 2], v_norm[, 1])
	v_plot <- cbind(v_radii * cos(angle), v_radii * sin(angle))
	list(u = unname(u_plot), v = unname(v_plot), mu = mu, mv = mv)
}

# fitted rank of the positions: the number of dimensions of the product
.uv_rank <- function(U, V, symmetric) {
	r_u <- dim(U)[2L]
	if (symmetric || is.null(V)) return(r_u)
	min(r_u, dim(V)[2L])
}

# the dimensions to draw: the first two by default (one at rank 1); a dynamic
# fit is drawn one plane at a time
.uv_check_dims <- function(dims, k, rank_txt = NULL, dynamic = FALSE) {
	rank_say <- if (is.null(rank_txt)) paste("rank", k) else rank_txt
	if (is.null(dims)) return(seq_len(min(2L, k)))
	ok <- is.numeric(dims) && length(dims) >= 1L && all(is.finite(dims)) &&
		all(dims == round(dims)) && !anyDuplicated(dims)
	if (!ok || any(dims < 1) || any(dims > k)) {
		cli::cli_abort(c(
			"{.arg dims} must be distinct whole numbers between 1 and {k}: this fit has {rank_say}.",
			"x" = "You supplied {.code dims = {deparse(dims)}}."
		))
	}
	dims <- as.integer(dims)
	if (length(dims) > 4L) {
		cli::cli_abort(c(
			"{.fn uv_plot} draws at most four dimensions at once (every pair as a panel); you asked for {length(dims)}.",
			"i" = "Pick up to four, e.g. {.code dims = 1:4}, or use {.code layout = \"heatmap\"}, which uses all {k} dimensions."
		))
	}
	if (dynamic && length(dims) != min(2L, k)) {
		cli::cli_abort(c(
			"A fit with dynamic positions is drawn one plane at a time: give {.arg dims} exactly two values, e.g. {.code dims = c(1, 3)}.",
			"i" = "This fit has {rank_say}."
		))
	}
	dims
}

# axis label of a dimension, with its share of the latent product when known
.uv_dim_label <- function(d, st) {
	if (is.null(st) || !is.finite(st$share[d])) return(paste("Dimension", d))
	sprintf("Dimension %d (%s%% of latent product)", d, formatC(100 * st$share[d], format = "f", digits = 0))
}

# strip label of a panel: the across and up dimensions
.uv_panel_name <- function(pr, st) {
	short <- function(d) {
		if (is.null(st) || !is.finite(st$share[d])) return(paste("dimension", d))
		sprintf("dimension %d (%s%%)", d, formatC(100 * st$share[d], format = "f", digits = 0))
	}
	if (length(pr) == 1L) return(paste0("D", substring(short(pr[1]), 2L)))
	paste0("Across: ", short(pr[1]), "; up: ", short(pr[2]))
}

# "1 and 2", "1, 2 and 3"
.uv_and <- function(x) {
	if (length(x) == 1L) return(as.character(x))
	paste(paste(utils::head(x, -1L), collapse = ", "), "and", utils::tail(x, 1L))
}

# caption naming the dimensions shown when the fit has more than two, with their
# share of the latent product when known, and the tie sentence for a shown
# dimension tied with one that is not shown
.uv_dims_caption <- function(dims, k, st) {
	if (k <= 2L) return(NULL)
	head_txt <- sprintf(
		"%s %s of %d shown",
		if (length(dims) == 1L) "Dimension" else "Dimensions",
		.uv_and(dims),
		k
	)
	line <- if (!is.null(st) && all(is.finite(st$share[dims]))) {
		sprintf("%s: %s%% of the latent product.", head_txt, formatC(100 * sum(st$share[dims]), format = "f", digits = 0))
	} else {
		paste0(head_txt, ".")
	}
	runs <- Filter(function(r) any(r %in% dims) && any(!r %in% dims), .lame_tie_runs(st))
	c(line, .lame_tie_sentence(st, runs))
}

# caption lines wrapped for the plot width
.uv_wrap <- function(lines, width = 85L) {
	paste(vapply(lines, function(l) paste(strwrap(l, width), collapse = "\n"), character(1)), collapse = "\n")
}

# a symmetric fit has one set of positions
.uv_abort_symmetric_nodes <- function(nodes) {
	cli::cli_abort(c(
		"A symmetric fit has one set of positions, so {.arg nodes} takes only {.val both}; you supplied {.val {nodes}}.",
		"i" = "Drop {.arg nodes}, or use {.code layout = \"heatmap\"} with {.code nodes = \"row\"} for the similarity of its nodes."
	))
}

# a fit-like object for positions passed as matrices: identity middle matrix,
# no posterior sds
.uv_matrix_fit <- function(U, V, row.names = NULL, col.names = NULL) {
	U <- as.matrix(U)
	V <- as.matrix(V)
	if (!is.null(row.names)) rownames(U) <- row.names
	if (!is.null(col.names)) rownames(V) <- col.names
	if (is.null(rownames(U))) rownames(U) <- paste0("S", seq_len(nrow(U)))
	if (is.null(rownames(V))) rownames(V) <- paste0("R", seq_len(nrow(V)))
	k <- min(ncol(U), ncol(V))
	structure(
		list(
			U = U[, seq_len(k), drop = FALSE],
			V = V[, seq_len(k), drop = FALSE],
			mode = "unipartite",
			symmetric = FALSE,
			uv_positions = "posterior_mean_product"
		),
		class = "ame"
	)
}

# heatmap of a fit's positions over every dimension: similarity of one side's
# nodes (latent_similarity), or the latent product of rows against columns
.uv_heatmap <- function(fit, nodes, time_point, plot_type, title) {
	if (plot_type != "snapshot") {
		cli::cli_abort(c(
			"{.code layout = \"heatmap\"} draws one period, so it takes only {.code plot_type = \"snapshot\"}.",
			"i" = "Choose the period with {.arg time_point}."
		))
	}
	is_dyn <- length(dim(fit$U)) == 3L
	tp <- if (is_dyn) time_point
	parts <- .uv_similarity_parts(fit, tp)
	is_bip <- identical(fit$mode, "bipartite")
	is_sym <- isTRUE(fit$symmetric)
	side_lab <- c(
		row = if (is_bip) "row actors" else if (is_sym) "nodes" else "senders",
		col = if (is_bip) "column actors" else if (is_sym) "nodes" else "receivers"
	)
	k <- parts$k
	dims_txt <- if (k == 1L) "1 dimension" else sprintf("all %d dimensions", k)
	at_txt <- if (!is.null(parts$period)) paste(" at", parts$period) else ""
	div_low <- "#2166AC"
	div_high <- "#B2182B"
	cap <- 2000L

	if (nodes != "both") {
		n_side <- if (nodes == "row") nrow(parts$U) else nrow(parts$V)
		if (n_side > cap) {
			other <- if (nodes == "row") "col" else "row"
			n_other <- if (nodes == "row") nrow(parts$V) else nrow(parts$U)
			cli::cli_abort(c(
				"The {nodes} side ({side_lab[[nodes]]}) has {format(n_side, big.mark = ',')} nodes; the similarity heatmap draws at most {format(cap, big.mark = ',')}.",
				"i" = "Use {.code nodes = \"{other}\"} ({format(n_other, big.mark = ',')} {side_lab[[other]]}), or {.code nodes = \"both\"} for the latent product of rows against columns."
			))
		}
		sim <- latent_similarity(fit, nodes, tp)
		ord <- attr(sim, "order")
		weak <- attr(sim, "weak")
		nm <- rownames(sim)
		n <- nrow(sim)
		S_plot <- unclass(sim)
		attr(S_plot, "order") <- NULL
		attr(S_plot, "weak") <- NULL
		is_weak <- nm %in% weak
		S_plot[is_weak, ] <- NA
		S_plot[, is_weak] <- NA
		pos <- integer(n)
		pos[ord] <- seq_len(n)
		heat_data <- data.frame(
			x = pos[rep(seq_len(n), each = n)],
			y = pos[rep(seq_len(n), times = n)],
			fill_value = as.vector(S_plot)
		)
		p <- ggplot2::ggplot(heat_data, ggplot2::aes(x = x, y = y, fill = fill_value))
		p <- p + ggplot2::geom_raster()
		p <- p + ggplot2::scale_fill_gradient2(
			low = div_low,
			mid = "white",
			high = div_high,
			midpoint = 0,
			limits = c(-1, 1),
			na.value = "grey75",
			name = paste0("Similarity (", dims_txt, ")")
		)
		axis_lab <- paste0(toupper(substring(side_lab[[nodes]], 1, 1)), substring(side_lab[[nodes]], 2))
		long_lab <- paste(axis_lab, "(ordered by clustering)")
		p <- .uv_heat_axes(
			p,
			n,
			nm[ord],
			n,
			nm[ord],
			long_lab,
			long_lab,
			axis_lab,
			axis_lab
		)
		p <- p + ggplot2::coord_fixed()
		caption <- NULL
		if (length(weak) > 0L) {
			shown <- if (length(weak) > 10L) {
				paste0(paste(utils::head(weak, 10L), collapse = ", "), " and ", length(weak) - 10L, " more")
			} else .uv_and(weak)
			caption <- sprintf("Not distinguishable from zero (grey, placed last): %s.", shown)
		}
		p <- p + ggplot2::labs(
			title = title %||% paste0("Latent similarity of ", side_lab[[nodes]], at_txt),
			caption = if (!is.null(caption)) .uv_wrap(caption)
		)
		return(.uv_heat_theme(p, n))
	}

	# latent product of rows against columns
	U <- parts$U
	V <- parts$V
	P <- (U %*% parts$M) %*% t(V)
	side_order <- function(A, which_side) {
		if (nrow(A) <= cap) {
			sim <- latent_similarity(fit, which_side, tp)
			return(list(
				ord = match(rownames(sim)[attr(sim, "order")], rownames(A)),
				how = "ordered by clustering",
				n_weak = length(attr(sim, "weak"))
			))
		}
		ang <- if (ncol(A) >= 2L) atan2(A[, 2], A[, 1]) else A[, 1]
		list(
			ord = order(ang),
			how = if (ncol(A) >= 2L) "ordered by angle in dimensions 1 and 2" else "ordered by dimension 1",
			n_weak = 0L
		)
	}
	ro <- side_order(U, "row")
	co <- side_order(V, "col")
	n_r <- nrow(P)
	n_c <- ncol(P)
	r_pos <- integer(n_r)
	r_pos[ro$ord] <- seq_len(n_r)
	c_pos <- integer(n_c)
	c_pos[co$ord] <- seq_len(n_c)
	heat_data <- data.frame(
		x = c_pos[rep(seq_len(n_c), each = n_r)],
		y = r_pos[rep(seq_len(n_r), times = n_c)],
		fill_value = as.vector(P)
	)
	lim <- max(abs(P), na.rm = TRUE)
	if (!is.finite(lim) || lim == 0) lim <- 1
	is_sym <- isTRUE(fit$symmetric)
	row_lab <- if (is_bip) "Row actors" else if (is_sym) "Actors" else "Senders"
	col_lab <- if (is_bip) "Column actors" else if (is_sym) "Actors" else "Receivers"
	p <- ggplot2::ggplot(heat_data, ggplot2::aes(x = x, y = y, fill = fill_value))
	p <- p + ggplot2::geom_raster()
	p <- p + ggplot2::scale_fill_gradient2(
		low = div_low,
		mid = "white",
		high = div_high,
		midpoint = 0,
		limits = c(-lim, lim),
		name = paste0("Latent product (", dims_txt, ")")
	)
	p <- .uv_heat_axes(
		p,
		n_c,
		colnames(P)[co$ord] %||% rownames(V)[co$ord],
		n_r,
		rownames(P)[ro$ord] %||% rownames(U)[ro$ord],
		paste0(col_lab, " (", co$how, ")"),
		paste0(row_lab, " (", ro$how, ")"),
		col_lab,
		row_lab
	)
	# nodes without a defined direction sit at the end of their side
	weak_txt <- c(
		if (ro$n_weak > 0L) sprintf("%d %s", ro$n_weak, tolower(row_lab)),
		if (co$n_weak > 0L && !is_sym) sprintf("%d %s", co$n_weak, tolower(col_lab))
	)
	p <- p + ggplot2::labs(
		title = title %||% paste0("Latent product of ", if (is_sym) "actors" else paste(tolower(row_lab), "and", tolower(col_lab)), at_txt),
		caption = if (length(weak_txt) > 0L) {
			.uv_wrap(sprintf("Placed last: %s not distinguishable from zero.", .uv_and(weak_txt)))
		}
	)
	.uv_heat_theme(p, n_c)
}

# heatmap axes: names on a side with at most 60 nodes, else only its order
.uv_heat_axes <- function(p, n_x, x_names, n_y, y_names, x_title_long, y_title_long, x_title, y_title) {
	if (n_x <= 60L) {
		p <- p + ggplot2::scale_x_continuous(breaks = seq_len(n_x), labels = x_names, expand = c(0, 0))
	} else {
		p <- p + ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0))
	}
	if (n_y <= 60L) {
		p <- p + ggplot2::scale_y_reverse(breaks = seq_len(n_y), labels = y_names, expand = c(0, 0))
	} else {
		p <- p + ggplot2::scale_y_reverse(breaks = NULL, expand = c(0, 0))
	}
	p + ggplot2::labs(
		x = if (n_x <= 60L) x_title else x_title_long,
		y = if (n_y <= 60L) y_title else y_title_long
	)
}

.uv_heat_theme <- function(p, n_x) {
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::theme(
		panel.border = ggplot2::element_blank(),
		axis.line = ggplot2::element_blank(),
		axis.ticks = ggplot2::element_blank(),
		legend.position = "top",
		plot.title = ggplot2::element_text(hjust = 0),
		plot.caption = ggplot2::element_text(hjust = 0)
	)
	# many column names need rotating to fit
	if (n_x > 10L && n_x <= 60L) {
		p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
	}
	p
}


#' Internal function for dynamic UV plotting
#' @keywords internal
#' @noRd
uv_plot_dynamic_internal <- function(
	fit,
	U,
	V,
	time_point,
	plot_type,
	show_arrows,
	title,
	label.nodes,
	label.size,
	colors,
	highlight = NULL,
	dims = 1:2,
	nodes = "both",
	k = dim(U)[2]
) {

	# colorbrewer Dark2 for small actor sets
	dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

	n_U <- dim(U)[1]
	n_V <- if (!is.null(V)) dim(V)[1] else n_U
	n <- n_U
	n_times <- dim(U)[3]
	is_bip <- !is.null(fit) && identical(fit$mode, "bipartite")

	u_names <- rownames(U)
	if (is.null(u_names)) u_names <- paste0(if (is_bip) "R" else "S", 1:n_U)
	v_names <- if (!is.null(V)) rownames(V) else u_names
	if (is.null(v_names)) v_names <- paste0(if (is_bip) "C" else "R", 1:n_V)

	time_labels <- dimnames(U)[[3]]
	if (is.null(time_labels)) time_labels <- paste0("T", 1:n_times)

	# the two shown dimensions of each side; rank 1 gets a zero second coordinate
	.uv_dim_cols <- function(A) {
		if (is.null(A)) return(NULL)
		if (length(dims) == 2L) return(A[, dims, , drop = FALSE])
		A_pad <- array(
			0,
			dim = c(dim(A)[1], 2, dim(A)[3]),
			dimnames = list(dimnames(A)[[1]], NULL, dimnames(A)[[3]])
		)
		A_pad[, 1, ] <- A[, dims[1], ]
		A_pad
	}
	U <- .uv_dim_cols(U)
	V <- .uv_dim_cols(V)
	show_u <- nodes != "col"
	show_v <- nodes != "row" && !is.null(V)
	x_lab <- paste("Dimension", dims[1])
	y_lab <- paste("Dimension", if (length(dims) == 2L) dims[2] else 2L)
	caption <- .uv_dims_caption(dims, k, NULL)
	caption <- if (length(caption) > 0L) .uv_wrap(caption)

	# one row per node and period, nodes outer, periods inner
	.uv_paths <- function(A, nm, type) {
		n_a <- dim(A)[1]
		data.frame(
			x = as.vector(t(matrix(A[, 1, ], n_a, n_times))),
			y = as.vector(t(matrix(A[, 2, ], n_a, n_times))),
			time = rep(seq_len(n_times), n_a),
			time_label = rep(time_labels, n_a),
			node = rep(nm, each = n_times),
			type = type,
			stringsAsFactors = FALSE
		)
	}

	if (plot_type == "trajectory") {
		traj_data <- rbind(
			if (show_u) .uv_paths(U, u_names, if (is_bip) "Row" else "Sender"),
			if (show_v) .uv_paths(V, v_names, if (is_bip) "Col" else "Receiver")
		)

		# colour by actor: the highlighted subset, else Dark2 for up to 8 actors
		n_actors <- length(unique(traj_data$node))
		if (!is.null(highlight)) {
			traj_data$highlight <- ifelse(
				traj_data$node %in% highlight,
				as.character(traj_data$node),
				"Other"
			)
			# highlighted actors first in the legend, "Other" last
			lvl <- c(intersect(unique(traj_data$node), highlight), "Other")
			traj_data$highlight <- factor(traj_data$highlight, levels = lvl)
			n_h <- length(highlight)
			h_cols <- dark2[seq_len(min(n_h, length(dark2)))]
			names(h_cols) <- highlight[seq_len(min(n_h, length(dark2)))]
			h_cols <- c(h_cols, Other = "grey75")
		}
		p <- ggplot2::ggplot(traj_data, ggplot2::aes(x = x, y = y, group = interaction(node, type))) +
			ggplot2::geom_path(
				ggplot2::aes(
					color = if (!is.null(highlight)) highlight else node,
					linetype = type
				),
				alpha = if (!is.null(highlight)) 0.4 else 0.5
			) +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				plot.caption = ggplot2::element_text(hjust = 0),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "top"
			) +
			ggplot2::coord_fixed(clip = "off") +
			ggplot2::labs(
				x = x_lab,
				y = y_lab,
				colour = if (!is.null(highlight)) "Actor" else NULL,
				linetype = "Node type",
				shape = "Node type"
			)
		if (!is.null(caption)) p <- p + ggplot2::labs(caption = caption)
		if (!is.null(highlight)) {
			p <- p + ggplot2::scale_colour_manual(values = h_cols)
		} else if (n_actors <= 8L) {
			p <- p + ggplot2::scale_colour_manual(values = dark2[seq_len(n_actors)])
		}

		if (show_arrows) {
			arrow_data <- traj_data[traj_data$time == n_times, ]
			prev_data <- traj_data[traj_data$time == n_times - 1, ]
			arrow_data$x_prev <- prev_data$x
			arrow_data$y_prev <- prev_data$y

			p <- p + ggplot2::geom_segment(
				data = arrow_data,
				ggplot2::aes(
					x = x_prev,
					y = y_prev,
					xend = x,
					yend = y,
					color = node
				),
				arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
				linewidth = 1
			)
		}

		start_data <- traj_data[traj_data$time == 1, ]
		end_data <- traj_data[traj_data$time == n_times, ]

		p <- p +
			ggplot2::geom_point(
				data = start_data,
				ggplot2::aes(color = node, shape = type),
				size = 3,
				alpha = 0.5
			) +
			ggplot2::geom_point(
				data = end_data,
				ggplot2::aes(color = node, shape = type),
				size = 4
			)

		if (label.nodes && n <= 20) {
			p <- p + ggrepel::geom_text_repel(
				data = end_data,
				ggplot2::aes(label = node, color = node),
				size = label.size,
				max.overlaps = Inf
			)
		}

		if (!is.null(title)) {
			p <- p + ggplot2::ggtitle(title)
		} else {
			p <- p + ggplot2::ggtitle("UV trajectories over time")
		}

		return(p)


	} else if (plot_type == "faceted") {
		if (n_times <= 6) {
			show_times <- 1:n_times
		} else {
			show_times <- round(seq(1, n_times, length.out = 6))
		}

		facet_data <- do.call(rbind, lapply(show_times, function(t) {
			rbind(
				if (show_u) {
					data.frame(
						x = U[, 1, t],
						y = U[, 2, t],
						node = u_names,
						time = time_labels[t],
						type = if (is_bip) "Row" else "Sender",
						stringsAsFactors = FALSE
					)
				},
				if (show_v) {
					data.frame(
						x = V[, 1, t],
						y = V[, 2, t],
						node = v_names,
						time = time_labels[t],
						type = if (is_bip) "Col" else "Receiver",
						stringsAsFactors = FALSE
					)
				}
			)
		}))
		rownames(facet_data) <- NULL

		p <- ggplot2::ggplot(facet_data, ggplot2::aes(x = x, y = y)) +
			ggplot2::geom_point(
				ggplot2::aes(
					color = if (!is.null(colors)) colors else type,
					shape = type
				),
				size = 2
			) +
			# coord_fixed() needs fixed facet scales
			ggplot2::facet_wrap(~time, scales = "fixed") +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border = ggplot2::element_blank(),
				plot.caption = ggplot2::element_text(hjust = 0),
				axis.ticks = ggplot2::element_blank(),
				legend.position = "top",
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text = ggplot2::element_text(color = "white", hjust = 0)
			) +
			ggplot2::coord_fixed(clip = "off") +
			ggplot2::labs(
				x = x_lab,
				y = y_lab,
				colour = if (is.null(colors)) "Node type" else NULL,
				shape = "Node type"
			)
		if (!is.null(caption)) p <- p + ggplot2::labs(caption = caption)

		if (label.nodes && n <= 10) {
			p <- p + ggrepel::geom_text_repel(
				ggplot2::aes(label = node),
				size = label.size * 0.8,
				max.overlaps = Inf
			)
		}

		if (!is.null(title)) {
			p <- p + ggplot2::ggtitle(title)
		} else {
			p <- p + ggplot2::ggtitle("UV positions over time")
		}

		return(p)
	}
}
