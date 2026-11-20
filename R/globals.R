# ggplot2 non-standard evaluation variables, declared for R CMD check
utils::globalVariables(c(
	"mu", "ymax", "ymin", "id", "mean_effect", # ab_plot
	"odmax", # get_start_vals
	"X1", "X2", "eff", "actor", "tPch", # uv_plot
	"Var2", "value", "x", "..density..", "y", "q90", "q95", "actual", # gof_plot
	"time", # gof_plot
	"Var1", "lo90", "hi90", "lo95", "hi95", "chain", "line", # trace_plot
	"iteration", "parameter", "type", "effect", "name", "index", # plot methods
	"observed", "median", "lower", "upper", "statistic", # plot methods
	"U1", "U2", "V1", "V2", "from", "to", "weight", # plot methods
	"x_from", "y_from", "x_to", "y_to", "node", # plot methods
	"xend", "yend", "color", "size", "label", # uv_plot aesthetics
	"x_prev", "y_prev", # uv_plot trajectory arrows
	"fill_value", # uv_plot heatmap
	"r", # ggforce::geom_circle aes
	".data" # aes pronoun, supplied by the ggplot2 data mask
))

# make a covariate-slice name end in `_dyad`, `_row` or `_col` exactly once
# (a `.<sfx>` ending becomes `_<sfx>`), so repeated calls are harmless
.lame_apply_suffix <- function(nms, sfx) {
	if (is.null(nms) || length(nms) == 0L) return(nms)
	has_us <- grepl(paste0("_", sfx, "$"), nms)
	has_dot <- grepl(paste0("\\.", sfx, "$"), nms)
	out <- nms
	out[has_dot] <- sub(
		paste0("\\.", sfx, "$"),
		paste0("_", sfx),
		out[has_dot]
	)
	bare <- !has_us & !has_dot
	if (any(bare)) out[bare] <- paste0(out[bare], "_", sfx)
	out
}

# the pooled node-coefficient name of a symmetric fit for a row (or column) slice name:
# only the generated markers change, the `_row` / `_col` suffix and the `X<row|col>` of a
# default name, so `growth_row` becomes `growth_node` and `Xrow1_row` becomes `Xnode1_node`
.lame_node_name <- function(nms) {
	out <- sub("_(row|col)$", "_node", nms)
	sub("^X(row|col)([0-9]+)_node$", "Xnode\\2_node", out)
}

# position in `table` (all coefficient names of one fit) of the coefficient each design-slice
# name in `x` feeds. the rule is chosen for the whole fit: a slice takes its own name, else its
# pooled node name (.lame_node_name); when that leaves a `_node` coefficient of the fit without a
# slice, the fit was saved by 1.4.2, whose samplers replaced every "row" in the row-slice name
# (`growth_row` became `gnodeth_node`, its twin `growth_col` feeds the same coefficient)
.lame_node_match <- function(x, table) {
	own <- match(x, table)
	pick <- function(key) {
		hit <- own
		hit[is.na(own)] <- match(key[is.na(own)], table)
		hit
	}
	hit <- pick(.lame_node_name(x))
	node <- grep("_node$", table)
	if (all(node %in% hit)) return(hit)
	row_form <- sub("^Xcol([0-9]+)_row$", "Xrow\\1_row", sub("_col$", "_row", x))
	pick(ifelse(grepl("_row$", row_form), gsub("row", "node", row_form, fixed = TRUE), NA_character_))
}
