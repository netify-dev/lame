# scale of residuals(type = "pearson"): the family's standard deviation at the
# fitted value mu (YPM, which carries a period's exposure)
#   normal        ve
#   binary, cbin  mu (1 - mu)
#   poisson       mu + mu^2 (exp(ve) - 1), for y ~ Poisson(exp(z)), z ~ N(eta, ve)
#   ordinal       category variance at the posterior-mean predictor, using the
#                 posterior-mean cutpoints, or for data-induced cutpoints those
#                 that reproduce the observed category shares (as simY_ord does)
#   frn           none: the response is a rank within its row
# a cell with zero variance has no residual (NA)

# variance of the ordinal levels at predictor eta (probit, unit latent
# variance) when P(y <= levels[j]) = pnorm(alpha[j] - eta)
.ordinal_response_variance <- function(eta, alpha, levels) {
	m1 <- m2 <- eta * 0
	lower <- eta * 0
	for (j in seq_along(levels)) {
		upper <- if (j < length(levels)) stats::pnorm(alpha[j] - eta) else eta * 0 + 1
		p <- upper - lower
		m1 <- m1 + levels[j] * p
		m2 <- m2 + levels[j]^2 * p
		lower <- upper
	}
	pmax(m2 - m1^2, 0)
}

# cutpoints at which a latent field with means eta and unit variance meets the
# cumulative category shares: the large-sample form of simY_ord's quantiles
.ordinal_induced_cutpoints <- function(eta, shares) {
	eta <- eta[is.finite(eta)]
	if (!length(eta)) {
		cli::cli_abort("The ordinal predictor has no finite cell to cut.", class = "lame_pearson_residuals")
	}
	lo <- min(eta) - 10
	hi <- max(eta) + 10
	vapply(shares[-length(shares)], function(s) {
		stats::uniroot(function(cut) mean(stats::pnorm(cut - eta)) - s, c(lo, hi), tol = 1e-12)$root
	}, numeric(1))
}

# the posterior-mean residual variance ve of a fit
.pearson_ve <- function(fit) {
	VC <- fit[["VC", exact = TRUE]]
	if (!is.matrix(VC) || nrow(VC) < 1L) {
		cli::cli_abort(
			c(
				"This fit stores no draws of the residual variance, so Pearson residuals cannot be formed.",
				"i" = "Use {.code type = \"response\"}, or refit with {.fn lame::ame} or {.fn lame::lame}."
			),
			class = "lame_pearson_residuals"
		)
	}
	mean(VC[, if ("ve" %in% colnames(VC)) "ve" else ncol(VC)], na.rm = TRUE)
}

# variance of one period's responses: mu is the fitted value; eta (the
# posterior-mean predictor) and y (observed) are used by the ordinal family only
.pearson_variance <- function(fit, mu, eta = NULL, y = NULL) {
	family <- fit[["family", exact = TRUE]] %||% "normal"
	switch(
		family,
		normal = mu * 0 + .pearson_ve(fit),
		binary = ,
		cbin = mu * (1 - mu),
		poisson = mu + mu^2 * expm1(.pearson_ve(fit)),
		ordinal = {
			if (is.null(eta)) {
				cli::cli_abort(
					"Ordinal Pearson residuals need the fitted linear predictor.", class = "lame_pearson_residuals"
				)
			}
			if (isTRUE(fit[["symmetric", exact = TRUE]])) eta <- eta / 2 + t(eta) / 2
			if (!identical(fit[["mode", exact = TRUE]], "bipartite") && nrow(eta) == ncol(eta)) diag(eta) <- NA_real_
			if (identical(fit[["ordinal_cutpoints", exact = TRUE]], "explicit")) {
				levels <- fit[["ordinal_levels", exact = TRUE]]
				ALPHA <- fit[["ALPHA", exact = TRUE]]
				if (length(levels) > 2L && (!is.matrix(ALPHA) || ncol(ALPHA) != length(levels) - 2L)) {
					cli::cli_abort(
						"Explicit ordinal Pearson residuals need the saved cutpoint draws {.code fit$ALPHA}.",
						class = "lame_pearson_residuals"
					)
				}
				alpha <- c(0, if (length(levels) > 2L) colMeans(ALPHA))
				.ordinal_validate_cutpoints(alpha, levels)
			} else {
				seen <- y[is.finite(y)]
				levels <- sort(unique(seen))
				shares <- cumsum(tabulate(match(seen, levels), nbins = length(levels))) / length(seen)
				alpha <- .ordinal_induced_cutpoints(eta, shares)
			}
			.ordinal_response_variance(eta, alpha, levels)
		},
		cli::cli_abort(
			c(
				"Pearson residuals are not defined for {.code family = \"{family}\"}.",
				"i" = "A fixed-rank nomination is a rank within its row, so the variance of a cell depends on the other cells of the row.",
				"i" = "Use {.code type = \"response\"}."
			),
			class = "lame_pearson_residuals"
		)
	)
}

# (y - mu) / sd with sd the family's standard deviation at the fitted value;
# NA where the variance is zero or not finite
.pearson_residual <- function(fit, y, mu, eta = NULL) {
	v <- .pearson_variance(fit, mu, eta, y)
	v[!is.finite(v) | v <= 0] <- NA_real_
	(y - mu) / sqrt(v)
}
