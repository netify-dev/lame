# saved scores distinguish observed marginal predictions from likelihood factors.
.lame_validate_score_request <- function(family, cutpoints, method, enabled) {
	if (!enabled || !identical(method, "observed_exact")) return(invisible(NULL))
	unsupported <- family %in% c("cbin", "frn") ||
		(identical(family, "ordinal") && !identical(cutpoints, "explicit"))
	if (unsupported) {
		reason <- if (family == "ordinal")
			"The data-induced ordinal model is a rank likelihood, not an empirical-cutpoint cumulative-probit likelihood." else
			"Nomination caps and rank constraints are not independent Bernoulli observations."
		guidance <- if (family == "ordinal")
			"Refit with ordinal_cutpoints = 'explicit' to obtain supported observed marginal scores." else
			"Use log_lik_method = 'observed_ghk' in lame() or 'observed_approx' in ame() for explicitly approximate diagnostics, or disable save_log_lik."
		cli::cli_abort(
			c(
				"Observed exact scoring is unavailable for family {.val {family}}.",
				"i" = reason,
				"i" = guidance
			),
			class = "lame_observed_score_unsupported"
		)
	}
	invisible(NULL)
}

.lame_score_index <- function(Y, bipartite = FALSE, symmetric = FALSE) {
	d <- dim(Y)
	if (length(d) == 2L) Y <- array(Y, c(d, 1L))
	mask <- is.finite(Y)
	if (!isTRUE(bipartite)) {
		for (t in seq_len(dim(Y)[3L])) {
			if (isTRUE(symmetric)) mask[, , t] <- mask[, , t] & upper.tri(Y[, , t])
			else diag(mask[, , t]) <- FALSE
		}
	}
	which(mask, arr.ind = TRUE)
}

.lame_explicit_ordinal_loglik <- function(y, eta, alpha, levels) {
	if (length(levels) < 2L || any(!is.finite(levels)) ||
			any(diff(levels) <= 0) || length(alpha) != length(levels) - 1L ||
			any(!is.finite(alpha)) || alpha[1L] != 0 || any(diff(alpha) <= 0) ||
			length(y) != length(eta) || any(!is.finite(eta)))
		cli::cli_abort("Explicit ordinal scoring requires finite matched predictors, sorted levels and anchored ordered cutpoints.")
	code <- match(y, levels)
	if (anyNA(code)) cli::cli_abort("Observed ordinal categories must match the fitted ordinal_levels.")
	bounds <- c(-Inf, alpha, Inf)
	log_phi_diff(
		bounds[code + 1L] - eta,
		bounds[code] - eta,
		width = bounds[code + 1L] - bounds[code]
	)
}

# rebuild the saved predictor after every parameter block has moved (no RNG)
.lame_saved_predictor <- function(
	Xlist,
	beta_path,
	a_path,
	b_path,
	U,
	V,
	G = NULL,
	rank_active = TRUE,
	offset = NULL,
	symmetric = FALSE
) {
	nt <- length(Xlist); na <- nrow(a_path); nb <- nrow(b_path)
	out <- array(0, c(na, nb, nt))
	at <- function(x, t) {
		if (length(dim(x)) != 3L) return(x)
		k <- if (dim(x)[3L] == 1L) 1L else t
		matrix(x[, , k], dim(x)[1L], dim(x)[2L])
	}
	for (t in seq_len(nt)) {
		X <- Xlist[[t]]; p <- dim(X)[3L]
		eta <- outer(a_path[, t], b_path[, t], "+")
		if (p > 0L) eta <- eta + matrix(matrix(X, na * nb, p) %*% beta_path[t, ], na, nb)
		if (rank_active) {
			u <- at(U, t); v <- at(V, t); g <- at(G, t)
			eta <- eta + if (is.null(g)) tcrossprod(u, v) else u %*% g %*% t(v)
		}
		if (!is.null(offset)) eta <- eta + offset[, , t]
		if (symmetric) eta <- (eta + t(eta)) / 2
		out[, , t] <- eta
	}
	if (any(!is.finite(out))) cli::cli_abort("The complete saved predictor is not finite.")
	out
}

.lame_score_contract <- function(
	family,
	method,
	cutpoints,
	bipartite,
	symmetric,
	dcor,
	index,
	rho = NULL
) {
	latent <- identical(method, "augmented")
	exact_family <- family %in% c("normal", "binary", "poisson") ||
		(family == "ordinal" && identical(cutpoints, "explicit"))
	approximate <- !latent && !exact_family
	independent <- isTRUE(bipartite) || isTRUE(symmetric) ||
		(!isTRUE(dcor) && (is.null(rho) || all(is.finite(rho) & rho == 0)))
	valid <- !latent && !approximate && independent
	reason <- if (latent) "Augmented Z scores are latent diagnostics, not observed-data likelihood factors." else
	if (approximate) "This rank/cap/empirical-threshold score is an approximate diagnostic, not the fitted observed-data likelihood." else
	if (!independent) "Directed cell marginals do not factorize when residual dyad correlation is estimated or nonzero; joint dyad likelihood factors are required." else
		"Observed conditional likelihood factors use unique cells (unordered dyads for symmetric fits)."
	list(
		version = 1L,
		family = family,
		method = method,
		response_scale = if (latent) "latent" else if (approximate) "approximate_proxy" else "observed",
		estimand = if (latent) "latent_gaussian_diagnostic" else if (approximate) "approximate_diagnostic" else "observed_cell_marginal",
		numerical_method = if (family == "poisson" && !latent) "adaptive_log_integral_with_concavity_tail_bounds" else
		if (approximate && method == "observed_ghk") "legacy_approximate_GHK_rank_estimator" else
		if (approximate) "legacy_approximate_response_proxy" else
		if (family == "ordinal") "normal_CDF_with_scaled_short_interval_integral" else "closed_form",
		numerical_accuracy = if (family == "poisson" && !latent)
			list(
				integration_relative_error_criterion = 1e-9, omitted_tail_relative_bound_per_side = 1e-12,
				reason = "The numerical integral must meet its error criterion and concavity tail bounds, or fail explicitly."
			) else NULL,
		observation_unit = if (bipartite) "bipartite_cell" else if (symmetric) "unordered_dyad" else "directed_cell",
		index = index,
		factorization = list(valid_for_loo = valid, reason = reason)
	)
}

.lame_require_loo_contract <- function(x) {
	contract <- x[["log_lik_contract"]]
	if (is.null(contract) || !identical(contract$version, 1L))
		cli::cli_abort(
			c(
				"The saved scores have no verified observed likelihood factorization contract.",
				"i" = "Refit with the current package and save_log_lik = TRUE. Existing scores remain readable with read_log_lik()."
			),
			class = "lame_loo_contract_missing"
		)
	if (!isTRUE(contract$factorization$valid_for_loo))
		cli::cli_abort(
			c(
				"LOO and WAIC require observed-data likelihood factors.",
				"i" = contract$factorization$reason,
				"i" = "The saved marginal or diagnostic scores remain available with read_log_lik()."
			),
			class = "lame_loo_factorization_unsupported"
		)
	invisible(contract)
}
