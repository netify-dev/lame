# Global variables used in ggplot2 functions
# These are column names used in non-standard evaluation
utils::globalVariables(c(
  "mu", "ymax", "ymin", "id",  # abPlot
  "odmax",  # getStartVals
  "X1", "X2", "eff", "actor", "tPch",  # ggCirc
  "Var2", "value", "x", "..density..", "y", "q90", "q95", "actual",  # gofPlot
  "time",  # gofPlot_long
  "Var1", "lo90", "hi90", "lo95", "hi95"  # paramPlot
))