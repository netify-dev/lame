# average standardized edge product over fully observed rectangles
# scaling before centering avoids overflow or underflow in the standard deviation
.lame_bipartite_four_cycle_dependence <- function(Y) {
	observed <- is.finite(Y)
	common <- tcrossprod(1 * observed)
	pairs <- upper.tri(common)
	n_rectangles <- sum(common[pairs] * (common[pairs] - 1) / 2)
	if (n_rectangles == 0) return(NA_real_)
	y <- Y[observed]
	scale <- max(abs(y))
	if (!is.finite(scale) || scale == 0) return(NA_real_)
	y <- y / scale
	spread <- stats::sd(y)
	if (!is.finite(spread) || spread == 0) return(NA_real_)
	E <- matrix(0, nrow(Y), ncol(Y))
	E[observed] <- (y - mean(y)) / spread
	cross <- tcrossprod(E)
	squares <- tcrossprod(E^2)
	sum((cross[pairs]^2 - squares[pairs]) / 2) / n_rectangles
}

.lame_bipartite_gof_names <- function(family = NULL) {
	c(
		"sd.rowmean",
		"sd.colmean",
		if (identical(family, "normal"))
			"four.cycle.dep" else "four.cycles"
	)
}

#' Goodness of fit statistics for bipartite networks
#'
#' Calculates goodness of fit statistics specifically designed for bipartite
#' networks, evaluating degree heterogeneity and higher-order dependencies.
#'
#' @param Y a bipartite relational data matrix (nA x nB rectangular matrix)
#'   where Y\\[i,j\\] represents the relationship from node i in set A to
#'   node j in set B. Missing values (NA) are allowed and will be handled
#'   appropriately.
#' @param family optional outcome family. For \code{"normal"}, use weighted
#'   four-cycle dependence instead of a nonzero-edge count. Model fitting
#'   and post-hoc diagnostics supply this automatically.
#' @param warn_square logical; if \code{TRUE} (default) a warning is issued
#'   when \code{Y} is square (a possible unipartite matrix passed by mistake).
#'   Set \code{FALSE} for a genuinely square bipartite network.
#' @return A named numeric vector containing bipartite-specific goodness-of-fit
#'   statistics:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means. Measures the
#'       heterogeneity in out-degree from set A nodes (sender effects).
#'       Higher values indicate more variation in how active A nodes are.}
#'     \item{sd.colmean}{Standard deviation of column means. Measures the
#'       heterogeneity in in-degree to set B nodes (receiver effects).
#'       Higher values indicate more variation in how popular B nodes are.}
#'     \item{four.cycle.dep}{For Gaussian outcomes, weighted dependence over
#'       complete observed rectangles, as defined below. Replaces the count.}
#'     \item{four.cycles}{Count of four-cycles (also called 4-paths or squares) in
#'       the bipartite network. A four-cycle occurs when two nodes from set A
#'       (e.g., i and k) both connect to the same two nodes in set B (e.g., j and l),
#'       forming a closed path: i->j->k->l->i. This measures the tendency for pairs
#'       of A-nodes to share multiple common B-node connections, capturing a form of
#'       clustering specific to bipartite networks. High four-cycle counts indicate
#'       that connections are not random but show patterns of shared preferences or
#'       co-occurrence. For example, in a user-item network, many four-cycles suggest
#'       that users who like one item tend to also like other items that co-occur
#'       with it.}
#'   }
#' @details
#' For Gaussian outcomes, \code{four.cycle.dep} is the mean product of
#' four globally centered, sample-standardized values around each fully
#' observed rectangle. It measures dependence among edge magnitudes.
#' It returns \code{NA} when no complete rectangle is observed or the
#' observed values have zero variance. Missing cells are not zeros.
#' A nonzero-edge count is unsuitable for Gaussian predictive checks:
#' continuous draws are nonzero almost surely, so fixing the observation
#' mask fixes their four-cycle count regardless of model fit.
#'
#' For bipartite networks, reciprocity and triadic closure are not meaningful
#' concepts since edges only exist between the two node sets. Instead, this
#' function focuses on:
#' \itemize{
#'   \item Degree heterogeneity in both node sets
#'   \item Four-cycles as the simplest higher-order dependence pattern
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
#' @examples
#' \donttest{
#' # Create a random bipartite network
#' Y <- matrix(rnorm(10*12), 10, 12)
#'
#' # Calculate GOF statistics
#' gof_stats_bipartite(Y, family = "normal")
#' }
#' @export
gof_stats_bipartite <- function(Y, warn_square = TRUE, family = NULL) {

	if (warn_square && nrow(Y) == ncol(Y)) {
		cli::cli_warn("Y appears to be square. Use gof_stats() for unipartite networks.")
	}

	gof <- numeric(3)
	gof[1] <- sd(rowMeans(Y, na.rm = TRUE), na.rm = TRUE)
	gof[2] <- sd(colMeans(Y, na.rm = TRUE), na.rm = TRUE)
	if (identical(family, "normal")) {
		gof[3] <- .lame_bipartite_four_cycle_dependence(Y)
	} else {
		Y_arr <- array(Y, dim = c(nrow(Y), ncol(Y), 1))
		gof[3] <- as.numeric(count_four_cycles_bip_cpp(Y_arr))
	}
	names(gof) <- .lame_bipartite_gof_names(family)

	return(gof)
}
