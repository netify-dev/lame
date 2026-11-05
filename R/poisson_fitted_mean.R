# fitted mean count of a Poisson fit, e exp(eta + s2 / 2) for y ~ Poisson(e exp(z)),
# z ~ N(eta, s2), at the posterior means of eta and s2; replicated counts are
# not averaged because they have no finite mean under the inverse-gamma prior on s2
.poisson_fitted_mean <- function(EZ, s2, exposure = 1, drop_diagonal = FALSE) {
	if (length(s2) != 1L || !is.finite(s2) || s2 < 0 ||
			length(exposure) != 1L || !is.finite(exposure) || exposure < 0) {
		cli::cli_abort(
			"The Poisson fitted mean needs a finite non-negative residual variance and exposure.",
			class = "lame_poisson_fitted"
		)
	}
	out <- if (exposure == 0) EZ * 0 else exp(EZ + s2 / 2 + log(exposure))
	if (any(!is.finite(out[!is.na(EZ)]))) {
		cli::cli_abort(
			"The Poisson fitted mean overflowed; no clipped mean was substituted.",
			class = "lame_poisson_fitted"
		)
	}
	if (isTRUE(drop_diagonal) && nrow(out) == ncol(out)) diag(out) <- NA_real_
	out
}

# exposure of period t, or one for a fit without exposure
.poisson_exposure_at <- function(exposure, t) {
	if (is.null(exposure) || !length(exposure)) 1 else exposure[min(t, length(exposure))]
}

# posterior mean residual variance of a Poisson fit
.poisson_residual_variance <- function(fit) {
	VC <- fit[["VC", exact = TRUE]]
	if (!is.matrix(VC) || !"ve" %in% colnames(VC) || !any(is.finite(VC[, "ve"]))) {
		cli::cli_abort(
			c(
				"This Poisson fit stores no draws of the residual variance, so a mean count cannot be formed.",
				"i" = "Refit with {.fn lame::ame} or {.fn lame::lame}, which keep {.code fit$VC}."
			),
			class = "lame_poisson_fitted"
		)
	}
	mean(VC[, "ve"], na.rm = TRUE)
}

# fitted mean counts of a Poisson fit: a matrix, or a list of period matrices
# for a panel (from its stored EZ)
.poisson_fitted_from_fit <- function(fit, EZ = NULL, exposure = fit[["period_exposure", exact = TRUE]]) {
	s2 <- .poisson_residual_variance(fit)
	unipartite <- !identical(fit[["mode", exact = TRUE]], "bipartite")
	EZ <- EZ %||% reconstruct_EZ(fit)
	if (is.list(EZ)) {
		out <- lapply(
			seq_along(EZ),
			function(t)
				.poisson_fitted_mean(EZ[[t]], s2, .poisson_exposure_at(exposure, t), unipartite)
		)
		names(out) <- names(EZ)
		return(out)
	}
	out <- .poisson_fitted_mean(EZ, s2, 1, unipartite)
	template <- fit[["YPM", exact = TRUE]]
	if (!is.null(dimnames(template)) && identical(dim(template), dim(out)))
		dimnames(out) <- dimnames(template)
	out
}

# response-scale value of a linear predictor eta in period t; only Poisson
# fits read the fit, for s2 (unless supplied) and the period exposure
.fit_to_response <- function(fit, eta, t = 1L, s2 = NULL) {
	family <- fit[["family", exact = TRUE]] %||% "normal"
	if (!identical(family, "poisson")) return(transform_to_response(eta, family))
	transform_to_response(
		eta,
		family,
		s2 = s2 %||% .poisson_residual_variance(fit),
		exposure = .poisson_exposure_at(fit[["period_exposure", exact = TRUE]], t)
	)
}
