# compact terminal states retain the dependence needed for posterior forecasts
.lame_forecast_trim <- function(state, n_draws) {
	out <- lapply(state, function(x) {
		d <- dim(x)
		idx <- lapply(d, seq_len)
		idx[[length(d)]] <- seq_len(n_draws)
		do.call(`[`, c(list(x), idx, list(drop = FALSE)))
	})
	out$draw_index <- seq_len(n_draws)
	out
}

.lame_forecast_state <- function(object) {
	S <- dim(object$BETA)[1L]
	state <- object[["forecast_state"]]
	if (is.null(state)) {
		# fits without forecast_state take the terminal states from full posterior storage
		last <- function(x, dynamic) {
			if (is.null(x)) return(NULL)
			d <- dim(x)
			if (!dynamic) return(x)
			idx <- lapply(d, seq_len)
			idx[[length(d) - 1L]] <- d[length(d) - 1L]
			y <- do.call(`[`, c(list(x), idx, list(drop = FALSE)))
			dim(y) <- d[-(length(d) - 1L)]
			y
		}
		state <- list(
			a = last(object$a_draws %||% object$a_samples, !is.null(object$a_draws)),
			b = last(object$b_draws %||% object$b_samples, !is.null(object$b_draws)),
			a_level = object$a_level_samples,
			b_level = object$b_level_samples,
			U = last(object$U_draws %||% object$U_samples, !is.null(object$U_draws)),
			V = last(object$V_draws %||% object$V_samples, !is.null(object$V_draws)),
			G = object$G_samples,
			draw_index = object$ab_draw_index %||% seq_len(S)
		)
		# longitudinal latent storage is aligned with beta, while effects may be thinned
		for (nm in c("U", "V", "G")) if (!is.null(state[[nm]])) {
			if (dim(state[[nm]])[3L] != S) {
				cli::cli_abort(
					"Saved latent draws are not aligned with the coefficient draws.",
					class = "lame_forecast_draws"
				)
			}
			state[[nm]] <- state[[nm]][, , state$draw_index, drop = FALSE]
		}
	}
	needed <- c("a", "b")
	if (isTRUE(object$dynamic_ab)) needed <- c(needed, "a_level", "b_level")
	if (any(c(object[["R"]], object$R_row, object$R_col) > 0)) {
		needed <- c(needed, "U", "V")
		if (identical(object$mode, "bipartite")) needed <- c(needed, "G")
	}
	idx <- state$draw_index
	valid_index <- is.numeric(idx) && length(idx) > 0L && all(is.finite(idx)) &&
		all(idx == as.integer(idx) & idx >= 1L & idx <= S) && !anyDuplicated(idx)
	valid <- valid_index && all(vapply(needed, function(nm) {
		x <- state[[nm]]
		!is.null(x) && is.numeric(x) && length(dim(x)) >= 2L &&
			utils::tail(dim(x), 1L) == length(idx) && all(is.finite(x))
	}, logical(1)))
	if (!valid) {
		cli::cli_abort(
			c(
				"This fit lacks aligned joint posterior states for forecasting.",
				"i" = "Refit with the current version of {.fn lame}; terminal states are retained automatically."
			),
			class = "lame_forecast_draws"
		)
	}
	state
}

# one prior transition per actor; heavy-tailed scales are shared across dimensions
.lame_forecast_uv_step <- function(
	state,
	rho,
	sigma,
	kind = "ar1",
	pi = NULL,
	kappa = NULL,
	nu = 4
) {
	n <- nrow(state)
	mu <- rho * state
	sd <- rep(sigma, n)
	if (identical(kind, "snap")) {
		if (is.null(pi) || !is.finite(pi) || pi < 0 || pi > 1) {
			cli::cli_abort(
				"Snap forecasts require a saved reset probability.",
				class = "lame_forecast_draws"
			)
		}
		reset <- stats::rbinom(n, 1L, pi) == 1L
		mu[reset, ] <- 0
		sd[reset] <- .lame_snap_kappa(kappa, rho, sigma)
	} else if (identical(kind, "t")) {
		sd <- sigma / sqrt(stats::rgamma(n, shape = nu / 2, rate = nu / 2))
	} else if (!identical(kind, "ar1")) {
		cli::cli_abort("Unsupported position transition in this forecast.")
	}
	mu + matrix(stats::rnorm(length(state)), n, ncol(state)) * sd
}

# pool terminal states across chains; NULL when any chain lacks them
.lame_pool_forecast_states <- function(chains, offsets) {
	if (!any(vapply(chains, function(f) !is.null(f$forecast_state), logical(1)))) return(NULL)
	states <- lapply(chains, function(f) tryCatch(.lame_forecast_state(f), error = function(e) NULL))
	if (any(vapply(states, is.null, logical(1)))) return(NULL)
	fields <- setdiff(unique(unlist(lapply(states, names))), "draw_index")
	out <- list()
	for (nm in fields) {
		xs <- lapply(states, `[[`, nm)
		if (all(vapply(xs, is.null, logical(1)))) next
		if (any(vapply(xs, is.null, logical(1)))) return(NULL)
		out[[nm]] <- abind::abind(xs, along = length(dim(xs[[1L]])))
	}
	out$draw_index <- unlist(Map(function(s, o) s$draw_index + o, states, offsets), use.names = FALSE)
	out
}

# put a complete stored design on the forecast grid and coefficient layout
.lame_forecast_design <- function(X, beta_names, rows, cols, n, m) {
	if (is.null(X)) return(NULL)
	if (length(dim(X)) == 2L) X <- array(X, c(dim(X), 1L))
	dn <- dimnames(X)
	if (!is.null(rows) && !is.null(dn[[1L]])) {
		if (!all(rows %in% dn[[1L]])) cli::cli_abort("Future covariates omit fitted row actors.")
		X <- X[match(rows, dn[[1L]]), , , drop = FALSE]
	}
	if (!is.null(cols) && !is.null(dn[[2L]])) {
		if (!all(cols %in% dn[[2L]])) cli::cli_abort("Future covariates omit fitted column actors.")
		X <- X[, match(cols, dn[[2L]]), , drop = FALSE]
	}
	if (!identical(as.integer(dim(X)[1:2]), as.integer(c(n, m))))
		cli::cli_abort("Future covariates must cover the fitted actor grid.")
	if (is.null(beta_names)) return(X)
	xn <- dimnames(X)[[3L]]
	if (is.null(xn)) {
		if (dim(X)[3L] != length(beta_names)) cli::cli_abort("The forecast design does not match the coefficients.")
		return(X)
	}
	out <- array(0, c(n, m, length(beta_names)), dimnames = list(rows, cols, beta_names))
	for (j in seq_along(beta_names)) {
		hit <- which(xn == beta_names[j])
		if (!length(hit)) hit <- which(sub("_(row|col)$", "_node", xn) == beta_names[j])
		if (!length(hit)) cli::cli_abort("The forecast design lacks coefficient {.val {beta_names[j]}}.")
		out[, , j] <- Reduce(`+`, lapply(hit, function(k) X[, , k]))
	}
	out
}

# forecast law of rw2 and matern32 coefficient paths: the Gaussian conditional of the h future periods
# given the fitted ones under the prior of the path (its level plus its change since the first period),
# on the fitted time positions extended by the common gap. the change has precision K / (sigma^2 lam),
# so the conditional mean A beta does not depend on the scale and the conditional covariance is
# sigma^2 lam C. returns, per kind present, A (h x T), the lower Cholesky factor L of C and the
# coefficients that follow the kind
.lame_beta_smooth_forecast_map <- function(object, h) {
	B <- object$BETA
	nt <- dim(B)[3L]
	names_beta <- dimnames(B)[[2L]]
	mask <- object$beta_dynamic_mask
	groups <- object$beta_dynamic_groups
	kinds <- .lame_beta_block_kinds(object)
	coef_kind <- rep("", length(names_beta))
	coef_kind[mask] <- unname(kinds[groups[mask]])
	present <- intersect(c("rw2", "matern32"), coef_kind)
	if (!length(present)) return(NULL)
	ti <- object$time_index
	tpos <- if (is.null(ti)) seq_len(nt) else as.numeric(ti)
	gap <- if (nt > 1L) tpos[nt] - tpos[nt - 1L] else 1
	all_pos <- c(tpos, tpos[nt] + gap * seq_len(h))
	past <- seq_len(nt)
	fut <- nt + seq_len(h)
	# the length scale of the fit: the prior's, or the default of the fitted periods
	ls <- object$prior$matern32_length_scale %||% max(2, nt / 4)
	out <- list()
	for (kind in present) {
		K <- .lame_beta_change_precision(kind, nt + h, ls, all_pos)
		C <- solve(K[fut, fut, drop = FALSE])
		C <- (C + t(C)) / 2
		A <- -C %*% K[fut, past, drop = FALSE]
		out[[kind]] <- list(A = A, L = t(chol(C)), coefs = which(coef_kind == kind))
	}
	out
}

# one forecast draw of the rw2 and matern32 coefficients: `path` the fitted p x T path of the draw,
# `sd` the per-coefficient scale sigma sqrt(lam), `beta_path` the h x p forecast to fill in
.lame_beta_smooth_forecast_draw <- function(beta_path, path, sd, map) {
	h <- nrow(beta_path)
	for (m in map) {
		for (k in m$coefs) {
			beta_path[, k] <- drop(m$A %*% path[k, ]) + sd[k] * drop(m$L %*% stats::rnorm(h))
		}
	}
	beta_path
}
