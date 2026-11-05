# ggplot2 non-standard evaluation variables, declared for R CMD check
utils::globalVariables(c(
	"mu", "ymax", "ymin", "id", "mean_effect", # ab_plot
	"odmax", # get_start_vals
	"X1", "X2", "eff", "actor", "tPch", # uv_plot
	"Var2", "value", "x", "..density..", "y", "q90", "q95", "actual", # gof_plot
	"time", # gof_plot
	"Var1", "lo90", "hi90", "lo95", "hi95", "chain", # trace_plot
	"iteration", "parameter", "type", "effect", "name", "index", # plot methods
	"observed", "median", "lower", "upper", "statistic", # plot methods
	"U1", "U2", "V1", "V2", "from", "to", "weight", # plot methods
	"x_from", "y_from", "x_to", "y_to", "node", # plot methods
	"xend", "yend", "color", "size", "label", # uv_plot aesthetics
	"x_prev", "y_prev", # uv_plot trajectory arrows
	"r" # ggforce::geom_circle aes
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
