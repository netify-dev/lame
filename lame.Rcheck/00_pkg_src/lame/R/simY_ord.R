#' Simulate an ordinal relational matrix
#'
#' Simulates an ordinal relational matrix whose marginal category
#' frequencies match those of an observed ordinal matrix.
#'
#' A latent Gaussian relational array is drawn from the social relations
#' model with mean \code{EZ} and within-dyad correlation \code{rho} via
#' \code{simZ}. The latent values are then discretized into ordinal
#' categories by cutting them at empirical quantiles chosen so that the
#' proportion of entries falling in each category reproduces the marginal
#' distribution of the observed categories of \code{Y}.
#'
#' @usage simY_ord(EZ, rho, Y)
#' @param EZ square matrix giving the expected value of the latent Z matrix
#' @param rho scalar giving the within-dyad correlation
#' @param Y ordinal relational data matrix
#' @return a square matrix
#' @author lame authors
#' @keywords internal
#' @export simY_ord
simY_ord <-
function(EZ, rho, Y) {
	# observed ordinal categories, in increasing order (NA dropped)
	observed <- as.vector(Y)
	categories <- sort(unique(observed))
	n_cat <- length(categories)

	# target marginal: cumulative share of mass at or below each category
	counts <- tabulate(match(observed, categories), nbins = n_cat)
	cum_share <- cumsum(counts / sum(counts))

	# latent Gaussian draw under the social relations model
	z_latent <- simZ(EZ, rho)
	if (nrow(z_latent) == ncol(z_latent)) diag(z_latent) <- NA

	# interior cutpoints: quantiles of the latent field at the target
	# cumulative shares (the final share of 1 needs no upper cut)
	interior_probs <- cum_share[-n_cat]
	thresholds <- quantile(z_latent, probs = interior_probs, na.rm = TRUE)

	# category index = 1 + (number of thresholds a latent value exceeds);
	# left.open matches the "z <= cut" boundary used to bucket the values
	z_vec <- as.vector(z_latent)
	bucket <- findInterval(z_vec, thresholds, left.open = TRUE) + 1L

	simulated <- matrix(categories[bucket],
		nrow = nrow(z_latent), ncol = ncol(z_latent),
		dimnames = dimnames(z_latent))
	simulated
}