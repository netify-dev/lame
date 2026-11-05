# diagnostics use identifiable products, never raw factor coordinates
.lame_robust_identified_draws <- function(object) {
	z <- object$draws
	D <- object$data
	ns <- nrow(z$hyper)
	label <- function(x, prefix, names) {
		x <- as.matrix(x)
		if (!ncol(x)) return(x)
		colnames(x) <- paste0(prefix, "[", names, "]")
		x
	}
	H <- z$hyper
	if (!D$rank) H <- H[, !colnames(H) %in% c("rho_uv", "var_uv"), drop = FALSE]
	M <- cbind(
		H,
		label(z$beta, "beta", colnames(D$X)),
		label(z$mu, "intercept", D$years)
	)
	if (!is.null(z$tau_ab)) M <- cbind(
		M,
		label(
			t(matrix(z$tau_ab, ncol = ns)),
			"actor_precision",
			as.vector(outer(c("sender", "receiver"), D$years, paste, sep = "@"))
		)
	)
	if (!is.null(z$omega)) {
		ijt <- which(D$observed, arr.ind = TRUE)
		nm <- paste0(D$actor_names[ijt[, 1L]], "->", D$actor_names[ijt[, 2L]], "@", D$years[ijt[, 3L]])
		M <- cbind(
			M,
			label(
				t(matrix(z$omega, ncol = ns)[which(D$observed), , drop = FALSE]),
				"observation_precision",
				nm
			)
		)
	}
	actor_period <- as.vector(outer(D$actor_names, D$years, paste, sep = "@"))
	for (side in c("A", "B"))
		M <- cbind(M, label(t(matrix(z[[side]], ncol = ns)), side, actor_period))
	for (side in c("a0", "b0"))
		M <- cbind(M, label(t(z[[side]]), side, D$actor_names))
	ij <- which(D$pairs, arr.ind = TRUE)
	pair_names <- paste(D$actor_names[ij[, 1L]], D$actor_names[ij[, 2L]], sep = "->")
	if (D$rank) {
		uv <- matrix(NA_real_, ns, nrow(ij) * D$tt)
		for (s in seq_len(ns)) for (t in seq_len(D$tt)) {
			U <- matrix(z$U[, , t, s], D$n, D$rank)
			V <- matrix(z$V[, , t, s], D$n, D$rank)
			uv[s, seq_len(nrow(ij)) + (t - 1L) * nrow(ij)] <- tcrossprod(U, V)[D$pairs]
		}
		M <- cbind(M, label(uv, "UV", as.vector(outer(pair_names, D$years, paste, sep = "@"))))
	}
	if (D$persistent) M <- cbind(
		M,
		label(t(matrix(z$lambda, ncol = ns)[which(D$pairs), , drop = FALSE]), "precision", pair_names)
	)
	M
}

.lame_robust_diagnostics <- function(chain_list) {
	if (!is.list(chain_list) || length(chain_list) < 2L ||
			!all(vapply(chain_list, inherits, logical(1), "lame_robust")))
		cli::cli_abort("Supply at least two independently seeded lame_robust fits.", call = NULL)
	first <- chain_list[[1L]]
	if (!all(vapply(
		chain_list,
		function(f) identical(f$data, first$data) &&
			identical(f$prior, first$prior),
		logical(1)
	)))
		cli::cli_abort("Robust diagnostic chains must use the same data, model, and priors.", call = NULL)
	seeds <- vapply(chain_list, function(f) as.integer(f$seed), integer(1))
	if (anyDuplicated(seeds))
		cli::cli_abort("Robust diagnostic chains must have different sampler seeds.", call = NULL)
	ns <- vapply(chain_list, function(f) nrow(f$draws$hyper), integer(1))
	if (length(unique(ns)) != 1L || ns[1L] < 4L)
		cli::cli_abort(
			"Robust diagnostic chains need equal lengths and at least four retained draws each.", call = NULL
		)
	if (!requireNamespace("posterior", quietly = TRUE))
		cli::cli_abort(
			"Install the 'posterior' package for robust-model rank-split R-hat and ESS diagnostics.", call = NULL
		)
	chains <- lapply(chain_list, .lame_robust_identified_draws)
	if (!all(vapply(chains, function(x) all(is.finite(x)), logical(1))))
		cli::cli_abort("Robust diagnostic draws must all be finite.", call = NULL)
	if (any(vapply(chains[-1L], identical, logical(1), chains[[1L]])))
		cli::cli_abort("Identical posterior draws cannot establish between-chain convergence.", call = NULL)
	names <- colnames(chains[[1L]])
	result <- matrix(
		NA_real_,
		length(names),
		5L,
		dimnames = list(names, c("rhat", "rhat_classic", "ess", "ess_bulk", "ess_tail"))
	)
	for (j in seq_along(names)) {
		x <- vapply(chains, function(z) z[, j], numeric(ns[1L]))
		W <- mean(apply(x, 2L, stats::var))
		B <- ns[1L] * stats::var(colMeans(x))
		if (is.finite(W) && W > 0) {
			result[j, "rhat_classic"] <- sqrt(((ns[1L] - 1) * W + B) / (ns[1L] * W))
			result[j, "ess"] <- min(as.numeric(coda::effectiveSize(coda::mcmc(as.vector(x)))), length(x))
		}
		result[j, "rhat"] <- posterior::rhat(x)
		result[j, "ess_bulk"] <- posterior::ess_bulk(x)
		result[j, "ess_tail"] <- posterior::ess_tail(x)
	}
	pass <- apply(
		result,
		1L,
		function(x) all(is.finite(x)) &&
			x["rhat"] < 1.05 && x["rhat_classic"] < 1.05 &&
			all(x[c("ess", "ess_bulk", "ess_tail")] > 0)
	)
	list(
		rhat = result[, "rhat"],
		rhat_classic = result[, "rhat_classic"],
		ess = result[, "ess"],
		ess_bulk = result[, "ess_bulk"],
		ess_tail = result[, "ess_tail"],
		param_names = names,
		n_chains = length(chains),
		n_samples = ns,
		seeds = seeds,
		diagnostics = as.data.frame(result),
		parameter_pass = pass,
		pass = all(pass),
		rhat_threshold = 1.05,
		scope = "Coefficients, time intercepts, actor paths and levels, estimated hyperparameters, all off-diagonal latent products, persistent directed-pair precisions, any role-period actor precisions, and temporary precisions at observed cells. The screen reports split R-hat and effective sample sizes for these parameters."
	)
}
