# Global variables used in ggplot2 functions
# These are column names used in non-standard evaluation
utils::globalVariables(c(
  "mu", "ymax", "ymin", "id",  # ab_plot
  "odmax",  # get_start_vals
  "X1", "X2", "eff", "actor", "tPch",  # uv_plot
  "Var2", "value", "x", "..density..", "y", "q90", "q95", "actual",  # gof_plot
  "time",  # gof_plot
  "Var1", "lo90", "hi90", "lo95", "hi95",  # trace_plot
  "iteration", "parameter", "type", "effect", "name", "index",  # plot methods
  "observed", "median", "lower", "upper", "statistic",  # plot methods
  "U1", "U2", "V1", "V2", "from", "to", "weight",  # plot methods
  "x_from", "y_from", "x_to", "y_to", "node",  # plot methods
  "xend", "yend", "color", "size", "label"  # uv_plot aesthetics
))