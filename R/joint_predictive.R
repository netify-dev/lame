# identify the MCMC iterations shared by every retained predictor component
.joint_predictive_state <- function(fit, burn_in = 0L, thin = 1L) {
	S <- dim(fit$BETA)[1L]
	fail <- function() cli::cli_abort(
		c(
			"Joint posterior prediction requires aligned saved additive and latent draws.",
			"i" = "Refit with posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)."
		),
		class = "lame_predictive_draws"
	)
	if (!is.null(fit$THETA_ACTOR) || !is.null(fit$dynamic_beta_per_actor)) cli::cli_abort(
		"Joint predictive simulation for actor-specific coefficient paths is not implemented.", class = "lame_predictive_model"
	)
	if (length(S) != 1L || S < 1L || length(burn_in) != 1L || !is.finite(burn_in) ||
			burn_in < 0L || burn_in >= S || burn_in != as.integer(burn_in) ||
			length(thin) != 1L || !is.finite(thin) || thin < 1L || thin != as.integer(thin))
		cli::cli_abort("burn_in must leave at least one draw and thin must be a positive integer.")
	bip <- identical(fit$mode, "bipartite")
	sym <- isTRUE(fit$symmetric) && !bip
	state <- list(
		a = fit$a_draws %||% fit$a_samples,
		b = fit$b_draws %||% fit$b_samples,
		U = fit$U_draws %||% fit$U_samples,
		V = fit$V_draws %||% fit$V_samples,
		G = fit$G_samples
	)
	if (sym) state$b <- state$a
	index <- function(x, supplied = NULL) {
		if (is.null(x)) return(NULL)
		if (!is.numeric(x) || is.null(dim(x)) || any(!is.finite(x))) fail()
		n <- utils::tail(dim(x), 1L)
		idx <- supplied %||% if (n == S) seq_len(S) else NULL
		if (is.null(idx) || length(idx) != n || anyNA(idx) || any(!is.finite(idx)) ||
				any(idx != as.integer(idx) | idx < 1L | idx > S) || anyDuplicated(idx)) fail()
		as.integer(idx)
	}
	state$ab_index <- index(state$a %||% state$b, fit$ab_draw_index)
	state$uv_index <- index(state$U, fit$uv_draw_index)
	need <- c(
		if (!identical(fit$rvar, FALSE)) "a",
		if (!identical(fit$cvar, FALSE) || (sym && !identical(fit$rvar, FALSE))) "b"
	)
	# explicit rank metadata takes precedence over placeholder coordinate arrays
	ranks <- c(fit[["R"]], fit[["R_row"]], fit[["R_col"]])
	rank <- if (length(ranks)) any(ranks > 0) else
		any(if (!is.null(fit[["U"]])) dim(fit[["U"]])[2L] > 0 else FALSE)
	if (rank) need <- c(need, "U", "V", if (bip) "G")
	if (any(vapply(state[need], is.null, logical(1)))) fail()
	for (nm in intersect(need, c("a", "b"))) if (!identical(
		index(state[[nm]], fit$ab_draw_index),
		state$ab_index
	)) fail()
	for (nm in intersect(need, c("U", "V", "G"))) if (!identical(
		index(state[[nm]], fit$uv_draw_index),
		state$uv_index
	)) fail()
	idx <- seq.int(as.integer(burn_in) + 1L, S, by = as.integer(thin))
	if (length(intersect(need, c("a", "b")))) idx <- intersect(idx, state$ab_index)
	if (rank) idx <- intersect(idx, state$uv_index)
	if (!length(idx) || !is.matrix(fit$VC) || nrow(fit$VC) != S) fail()
	state$index <- idx
	state
}

# aligned joint draws for simulate(type = "posterior"); a fit without them is refused with
# an error naming the approximate route, and other failures are raised as they are
.simulate_joint_state <- function(fit, burn_in = 0L, thin = 1L) {
	tryCatch(
		.joint_predictive_state(fit, burn_in, thin),
		lame_predictive_draws = function(e) {
			cli::cli_abort(
				c(
					"{.code simulate(type = \"posterior\")} replicates the network from aligned saved draws of the coefficients, additive effects and positions, and this fit does not carry them.",
					"i" = "Refit with {.code posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)}.",
					"i" = "{.code type = \"approximate\"} simulates from the posterior means and standard deviations of the effects instead; it does not propagate their joint posterior uncertainty."
				),
				class = c("lame_simulate_draws", "lame_predictive_draws")
			)
		}
	)
}

# contract one complete, matched posterior state on the stored actor grid
.joint_predictor <- function(fit, state, s, t = 1L, design = NULL) {
	B <- fit$BETA
	beta <- if (length(dim(B)) == 3L) B[s, , t] else B[s, ]
	names(beta) <- dimnames(B)[[2L]]
	X <- design %||% if (inherits(fit, "lame")) fit$Xlist[[t]] else fit[["X"]]
	grid <- .fit_actor_order(fit)
	base <- if (is.list(fit$YPM)) fit$YPM[[t]] else fit$YPM
	n <- if (!is.null(state$a)) dim(state$a)[1L] else if (!is.null(state$U)) dim(state$U)[1L] else length(grid$rows) %||% nrow(base)
	m <- if (!is.null(state$b)) dim(state$b)[1L] else if (!is.null(state$V)) dim(state$V)[1L] else length(grid$cols) %||% ncol(base)
	if (!n) n <- nrow(base)
	if (!m) m <- ncol(base)
	X <- .lame_forecast_design(X, names(beta), grid$rows, grid$cols, n, m)
	if (is.null(X) && length(beta)) cli::cli_abort("Joint prediction requires the original covariate design.")
	eta <- matrix(0, n, m, dimnames = list(grid$rows, grid$cols))
	for (j in seq_along(beta)) eta <- eta + beta[j] * X[, , j]
	ab <- match(s, state$ab_index); uv <- match(s, state$uv_index)
	at <- function(x, k, dynamic_dim) {
		if (is.null(x)) return(NULL)
		d <- dim(x); ix <- lapply(d, seq_len); ix[[length(d)]] <- k
		if (length(d) == dynamic_dim) ix[[length(d) - 1L]] <- t
		z <- do.call(`[`, c(list(x), ix, list(drop = FALSE)))
		if (dynamic_dim == 3L) as.numeric(z) else matrix(z, d[1L], d[2L])
	}
	a <- at(state$a, ab, 3L); b <- at(state$b, ab, 3L)
	if (!is.null(a)) eta <- eta + a
	if (!is.null(b)) eta <- sweep(eta, 2L, b, "+")
	if (!is.null(state$U)) {
		u <- at(state$U, uv, 4L); v <- at(state$V, uv, 4L)
		if (identical(fit$mode, "bipartite")) {
			g <- matrix(state$G[, , uv], dim(state$G)[1L], dim(state$G)[2L])
			eta <- eta + u %*% g %*% t(v)
		} else eta <- eta + u %*% t(v)
	}
	eta
}

# outcome and returned latent field share the same observation-noise draw
.predictive_observation <- function(eta, fit, s, t = 1L) {
	family <- fit$family; bip <- identical(fit$mode, "bipartite")
	sym <- isTRUE(fit$symmetric) && !bip
	if (sym) eta <- (eta + t(eta)) / 2
	vc <- fit$VC[s, ]
	ve <- if (family %in% c("normal", "poisson")) unname(vc["ve"]) else 1
	rho <- if (bip || sym || !"rho" %in% names(vc)) 0 else unname(vc["rho"])
	if (!bip && !sym && isTRUE(fit$dynamic_rho)) {
		if (!is.matrix(fit$RHO) || nrow(fit$RHO) != dim(fit$BETA)[1L] ||
				ncol(fit$RHO) < t)
			cli::cli_abort(
				"Period-specific predictive reciprocity requires aligned saved RHO draws.",
				class = "lame_predictive_draws"
			)
		rho <- unname(fit$RHO[s, t])
	}
	if (length(ve) != 1L || !is.finite(ve) || ve < 0 || !is.finite(rho) || abs(rho) > 1)
		cli::cli_abort("Predictive simulation requires finite valid observation variance and correlation draws.")
	alpha <- if (identical(family, "ordinal") && identical(fit$ordinal_cutpoints, "explicit"))
		.ordinal_saved_cutpoints(fit, s) else NULL
	Z <- if (bip) eta + matrix(stats::rnorm(length(eta), sd = sqrt(ve)), nrow(eta), ncol(eta)) else simZ(eta, rho, ve)
	if (!bip) diag(Z) <- NA_real_
	if (sym) Z[lower.tri(Z)] <- t(Z)[lower.tri(Z)]
	Y <- switch(
		family,
		normal = Z,
		binary = 1 * (Z > 0),
		poisson = {
			exposure <- (fit$period_exposure %||% rep(1, t))[t]
			if (!is.finite(exposure) || exposure < 0) cli::cli_abort("Prediction requires finite nonnegative period exposure.")
			lambda <- if (exposure == 0) matrix(0, nrow(Z), ncol(Z)) else exp(Z + log(exposure))
			if (any(!is.finite(lambda[!is.na(Z)]))) cli::cli_abort("The simulated Poisson rate overflowed; no clipped count draw was substituted.")
			out <- matrix(NA_real_, nrow(Z), ncol(Z)); ok <- is.finite(lambda)
			out[ok] <- stats::rpois(sum(ok), lambda[ok]); out
		},
		ordinal = {
			if (identical(fit$ordinal_cutpoints, "explicit")) {
				.ordinal_categories(Z, alpha, fit$ordinal_levels)
			} else {
				obs <- fit$Y_obs %||% fit$Y
				if (is.list(obs)) obs <- obs[[t]] else if (length(dim(obs)) == 3L) obs <- obs[, , t]
				levels <- fit$ordinal_levels %||% sort(unique(as.vector(obs[!is.na(obs)])))
				if (!length(levels)) cli::cli_abort("Ordinal simulation requires the observed category levels.")
				count <- tabulate(match(as.vector(obs), levels), nbins = length(levels))
				cuts <- stats::quantile(Z, cumsum(count / sum(count))[-length(levels)], na.rm = TRUE)
				matrix(levels[findInterval(as.vector(Z), cuts, left.open = TRUE) + 1L], nrow(Z), ncol(Z))
			}
		},
		cbin = NULL,
		frn = NULL,
		cli::cli_abort("Unsupported predictive observation family.")
	)
	if (family %in% c("cbin", "frn")) {
		zz <- Z; zz[is.na(zz)] <- -Inf
		yo <- .sim_absent_yo(fit, t, nrow(Z), bip)
		if (!is.null(yo)) zz[is.na(yo)] <- -Inf
		cap <- rep_len(fit$odmax %||% Inf, nrow(Z)); Y <- matrix(0, nrow(Z), ncol(Z))
		for (i in seq_len(nrow(Z))) {
			eligible <- which(is.finite(zz[i, ]) & zz[i, ] > 0)
			eligible <- eligible[order(zz[i, eligible], decreasing = TRUE)]
			if (is.finite(cap[i])) eligible <- head(eligible, cap[i])
			if (length(eligible)) Y[i, eligible] <- if (family == "cbin") 1 else rank(zz[i, eligible])
		}
		Y[is.na(Z)] <- NA_real_
		if (!is.null(yo)) Y[is.na(yo)] <- NA_real_
	}
	if (!bip) diag(Y) <- NA_real_
	if (sym) Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]
	dimnames(Y) <- dimnames(Z) <- dimnames(eta)
	list(Y = Y, Z = Z)
}

.simulate_joint <- function(fit, state, nsim, n_time = 1L, return_latent = FALSE, designs = NULL) {
	if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim))
		cli::cli_abort("nsim must be a positive integer.")
	idx <- state$index[sample.int(length(state$index), nsim, replace = TRUE)]
	out <- list(
		Y = vector("list", nsim),
		family = fit$family,
		mode = fit$mode,
		method = "joint posterior",
		draw_index = idx
	)
	if (return_latent) out$Z <- vector("list", nsim)
	long <- inherits(fit, "lame")
	for (i in seq_len(nsim)) {
		ys <- zs <- vector("list", n_time)
		for (t in seq_len(n_time)) {
			eta <- .joint_predictor(fit, state, idx[i], t, designs[[t]])
			obs <- .predictive_observation(eta, fit, idx[i], t)
			ys[[t]] <- .sim_period_grid(obs$Y, fit, t, identical(fit$mode, "bipartite"))
			if (return_latent) zs[[t]] <- .sim_period_grid(obs$Z, fit, t, identical(fit$mode, "bipartite"))
		}
		out$Y[[i]] <- if (long) ys else ys[[1L]]
		if (return_latent) out$Z[[i]] <- if (long) zs else zs[[1L]]
	}
	if (long) out$n_time <- n_time
	class(out) <- if (long) "lame.sim" else "ame.sim"
	out
}

# names for unnamed new covariates: the fitted names of that role when the counts match
.joint_role_names <- function(beta_names, role, k, prefix) {
	suffix <- paste0("_", role, "$")
	fitted <- beta_names[grepl(suffix, beta_names)]
	if (length(fitted) == k) return(sub(suffix, "", fitted))
	paste0(prefix, seq_len(k))
}

.joint_new_design <- function(fit, newdata, t = NULL) {
	if (is.null(newdata)) return(NULL)
	ord <- .fit_actor_order(fit)
	base <- if (is.list(fit$YPM)) fit$YPM[[1L]] else fit$YPM
	n <- length(ord$rows); m <- length(ord$cols)
	if (!n) n <- nrow(base)
	if (!m) m <- ncol(base)
	bn <- dimnames(fit$BETA)[[2L]]
	layers <- list(); labels <- character()
	if ("intercept" %in% bn) { layers <- list(matrix(1, n, m)); labels <- "intercept" }
	for (role in c("row", "col", "dyad")) {
		x <- newdata[[paste0("X", role)]]
		if (!is.null(t) && is.list(x) && !is.data.frame(x)) x <- x[[t]]
		if (is.null(x)) next
		if (role == "dyad") {
			if (is.matrix(x)) x <- array(x, c(dim(x), 1L))
			x <- .reorder_newdata_actors(x, ord$rows, ord$cols)
			if (length(dim(x)) != 3L || !identical(as.integer(dim(x)[1:2]), as.integer(c(n, m))))
				cli::cli_abort("New dyadic covariates do not match the fitted actor grid.")
			nm <- dimnames(x)[[3L]] %||% .joint_role_names(bn, "dyad", dim(x)[3L], "Xdyad")
			layers <- c(layers, lapply(seq_len(dim(x)[3L]), function(j) x[, , j]))
		} else {
			x <- as.matrix(x); actors <- if (role == "row") ord$rows else ord$cols
			if (!is.null(rownames(x)) && !is.null(actors)) {
				if (!all(actors %in% rownames(x))) cli::cli_abort("New nodal covariates omit fitted actors.")
				x <- x[actors, , drop = FALSE]
			}
			if (nrow(x) != if (role == "row") n else m) cli::cli_abort("New nodal covariates do not match the fitted actor grid.")
			nm <- colnames(x) %||% .joint_role_names(bn, role, ncol(x), paste0("X", role))
			layers <- c(layers, lapply(seq_len(ncol(x)), function(j) matrix(x[, j], n, m, byrow = role == "col")))
		}
		labels <- c(labels, .lame_apply_suffix(nm, role))
	}
	X <- array(0, c(n, m, length(layers)), dimnames = list(ord$rows, ord$cols, labels))
	for (j in seq_along(layers)) X[, , j] <- layers[[j]]
	X <- .lame_forecast_design(X, bn, ord$rows, ord$cols, n, m)
	.lame_tempered_prediction_design(fit, X)
}
