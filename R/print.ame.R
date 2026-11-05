#' Print method for AME model objects
#'
#' @description
#' Displays a formatted summary of a fitted AME (Additive and Multiplicative Effects) model.
#' This method provides a concise overview of model structure, parameter estimates,
#' and goodness-of-fit statistics without generating new data.
#'
#' @param x an object of class "ame" from fitting an AME model
#' @param ... additional arguments (not currently used)
#'
#' @return Invisibly returns the input object (for method chaining)
#'
#' @details
#' The print method displays:
#' \itemize{
#'   \item Model type (unipartite/bipartite, symmetric/asymmetric)
#'   \item Network dimensions and observation count
#'   \item Family and link function used
#'   \item Number of MCMC iterations
#'   \item Parameter counts for regression coefficients and latent factors
#'   \item Basic convergence diagnostics if available
#' }
#'
#' Unlike \code{simulate}, this method only formats existing results
#' for display and does not perform any new computations or data generation.
#'
#' @examples
#' \donttest{
#' # Fit model
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE, chains = 1)
#'
#' # Display summary
#' print(fit)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @seealso
#' \code{\link{summary.ame}} for detailed summaries,
#' \code{\link{simulate.ame}} for generating new networks,
#' \code{\link{predict.ame}} for predictions
#' @method print ame
#' @export
print.ame <- function(x, ...) {
	.lame_cli_stdout({
		cli::cli_h1("Additive and multiplicative effects (AME) model")

		# number of actors
		if (!is.null(x$n)) {
			n <- x$n
		} else if (!is.null(x$APM)) {
			n <- length(x$APM)
		} else if (!is.null(x[["Y", exact = TRUE]])) {
			n <- nrow(x[["Y", exact = TRUE]])
		} else {
			n <- NA
		}

		# bipartite dimensions; identical() handles a fit without $mode, and the
		# exact lookup keeps $Y from partial-matching $YPM
		if (identical(x$mode, "bipartite") && !is.null(x[["Y", exact = TRUE]])) {
			n_row <- nrow(x[["Y", exact = TRUE]])
			n_col <- ncol(x[["Y", exact = TRUE]])
		} else {
			n_row <- n
			n_col <- n
		}
		n_stored <- nrow(x$BETA)
		cli::cli_h3("Network structure")
		cli::cli_bullets(c(
			"*" = "Network dimensions: {.val {n_row}} x {.val {n_col}}",
			"*" = "Stored posterior samples: {.val {n_stored}}"
		))

		config_info <- character()
		if (!is.null(x$model.name)) {
			config_info <- c(config_info, "*" = "Model name: {.field {x$model.name}}")
		}
		if (!is.null(x$family)) {
			config_info <- c(config_info, "*" = "Family: {.field {x$family}}")
		}
		if (!is.null(x$mode)) {
			config_info <- c(config_info, "*" = "Mode: {.field {x$mode}}")
		}
		if (!is.null(x$symmetric)) {
			config_info <- c(config_info, "*" = "Symmetric: {.val {x$symmetric}}")
		}
		if (isTRUE(x$dynamic_uv)) {
			config_info <- c(config_info, "*" = "Dynamic latent positions: {.code dynamic_uv = TRUE}")
		}
		if (isTRUE(x$dynamic_ab)) {
			config_info <- c(config_info, "*" = "Dynamic additive effects: {.code dynamic_ab = TRUE}")
		}
		if (length(config_info)) {
			cli::cli_h3("Model configuration")
			cli::cli_bullets(config_info)
		}

		cli::cli_h3("Model parameters")
		param_info <- c("*" = "Regression coefficients: {.val {ncol(x$BETA)}}")

		if (!is.null(x$rvar) && x$rvar) {
			param_info <- c(param_info, "*" = "Row/sender effects: enabled")
		}
		if (!is.null(x$cvar) && x$cvar) {
			param_info <- c(param_info, "*" = "Column/receiver effects: enabled")
		}
		if (!is.null(x$dcor) && x$dcor) {
			param_info <- c(param_info, "*" = "Dyadic correlation: enabled")
		}

		if (!is.null(x$R) && x$R > 0) {
			param_info <- c(param_info, "*" = "Multiplicative effects dimension: {.val {x$R}}")
		}
		cli::cli_bullets(param_info)

		# statements on what the fit does not report (.lame_ab_statements)
		for (st_ in .lame_ab_statements(x)) cli::cli_alert_info("{st_}")

		if (!is.null(.chain_line <- .lame_chain_print_line(x$auto_extend, x))) cli::cli_alert_info(.chain_line)

		cli::cli_alert_info("Use {.code summary(object)} for detailed results")
	})
	invisible(x)
}
