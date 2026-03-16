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
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
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
	cat("\nAdditive and Multiplicative Effects (AME) Model\n")
	cat("================================================\n")
	
	# node count
	if(!is.null(x$n)) {
		n <- x$n
	} else if(!is.null(x$APM)) {
		n <- length(x$APM)
	} else if(!is.null(x$Y)) {
		n <- nrow(x$Y)
	} else {
		n <- NA
	}
	
	# bipartite dimensions
	if(x$mode == "bipartite" && !is.null(x$Y)) {
		cat("\nNetwork dimensions: ", nrow(x$Y), "x", ncol(x$Y), "\n")
	} else {
		cat("\nNetwork dimensions: ", n, "x", n, "\n")  
	}
	
	nscan <- nrow(x$BETA)
	cat("MCMC iterations: ", nscan, "\n")
	
	# model configuration
	if (!is.null(x$model.name)) {
		cat("Model name: ", x$model.name, "\n")
	}
	if (!is.null(x$family)) {
		cat("Family: ", x$family, "\n")
	}
	if (!is.null(x$mode)) {
		cat("Mode: ", x$mode, "\n")
	}
	if (!is.null(x$symmetric)) {
		cat("Symmetric: ", x$symmetric, "\n")
	}
	if (isTRUE(x$dynamic_uv)) {
		cat("Dynamic UV: TRUE\n")
	}
	if (isTRUE(x$dynamic_ab)) {
		cat("Dynamic AB: TRUE\n")
	}

	# parameter summary
	cat("\nNumber of parameters:\n")
	cat("  Regression coefficients: ", ncol(x$BETA), "\n")
	
	# additive effects
	if(!is.null(x$rvar) && x$rvar) {
		cat("  Row/sender effects: enabled\n")
	}
	if(!is.null(x$cvar) && x$cvar) {
		cat("  Column/receiver effects: enabled\n") 
	}
	if(!is.null(x$dcor) && x$dcor) {
		cat("  Dyadic correlation: enabled\n")
	}
	
	# multiplicative effects dimension
	if (!is.null(x$R) && x$R > 0) {
		cat("  Multiplicative effects dimension: ", x$R, "\n")
	}
	
	
	cat("\nUse summary(object) for detailed results\n")
	
	invisible(x)
}