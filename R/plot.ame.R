#' Simple diagnostic plot for AME model fit
#' 
#' Provides a quick visual summary of an AME model fit, focusing on key 
#' convergence diagnostics. For comprehensive diagnostics, use the specific
#' plotting functions: trace_plot(), gof_plot(), ab_plot(), and uv_plot().
#' 
#' @details
#' By default, this function simply calls trace_plot() to show MCMC trace
#' plots and posterior distributions for key parameters. This provides a
#' quick check of model convergence and mixing.
#' 
#' For more detailed visualizations, use the specialized functions:
#' \describe{
#'   \item{trace_plot()}{MCMC diagnostics and posterior distributions}
#'   \item{gof_plot()}{Goodness-of-fit assessment}
#'   \item{ab_plot()}{Additive sender/receiver effects}
#'   \item{uv_plot()}{Multiplicative latent factors}
#' }
#' 
#' @param x an object of class "ame" from fitting an AME model
#' @param ... additional arguments passed to trace_plot()
#' @return A ggplot2 object from trace_plot() (invisibly)
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{ame}}, \code{\link{trace_plot}}, \code{\link{gof_plot}}, 
#'          \code{\link{ab_plot}}, \code{\link{uv_plot}}
#' @examples
#' \donttest{
#' # Fit an AME model
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = TRUE,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
#'
#' # Quick diagnostic plot (shows trace plots)
#' plot(fit)
#' }
#' @method plot ame
#' @export
plot.ame <- function(x, ...) {
	p <- trace_plot(x, ...)
	
	invisible(p)
}