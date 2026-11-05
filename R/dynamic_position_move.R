# whole-path moves of an actor's positions under dynamic_uv for the Poisson
# family (a pointwise likelihood factor) and the binary family (a probit
# factor), and two exact draws for the symmetric unit-variance gauge. levels,
# coefficients and static positions use R/latent_translation_move.R; the
# ordinal family has no path move.
#
# the Gibbs sweep moves a position by about sqrt(s2 / cells) per sweep however
# wide its posterior, so on sparse panels whole paths mix slowly.
#
# Poisson: shifting u_it by d_t and the cells z_ijt by d_t' w_jt leaves every
# residual z - E[z], hence the Gaussian density of the latent array, unchanged,
# so the path's conditional is its AR(1) prior times the Poisson factor
#     sum_{j, t} [ y_ijt delta_ijt - exp(z_ijt + offset_t) (exp(delta_ijt) - 1) ],  delta_ijt = d_t' w_jt,
# drawn by an elliptical slice step (Murray, Adams and MacKay 2010).
# binary: the latent array is integrated out as in R/binary_actor_collapsed_move.R
# (shared dyad factor drawn once per visit, a probit factor per cell) and is
# redrawn given the factor at the end of the visit.
# partner weights are v_jt (directed), G v_jt or u_it' G (bipartite) and l * u_jt
# (symmetric, where the diagonal cell moves by sum_r l_r (u'^2 - u^2)). sender
# and receiver paths are visited in turn; rho_uv and sigma_uv keep their own
# updates, and the rail prior$uv_max_abs enters as the indicator of the box.

# the covariance of an AR(1) path of length nt, per coordinate: stationary start N(0, sigma^2 / (1 - rho^2)),
# transitions N(rho u_{t-1}, sigma^2 / lambda_t) (lambda NULL: all 1)
.lame_dpm_cov <- function(nt, rho, sigma, lambda = NULL) {
	v <- numeric(nt)
	v[1L] <- sigma^2 / (1 - rho^2)
	cv <- matrix(0, nt, nt)
	cv[1L, 1L] <- v[1L]
	if (nt > 1L) {
		for (t in 2:nt) {
			v[t] <- rho^2 * v[t - 1L] + sigma^2 / (if (is.null(lambda)) 1 else lambda[t])
			cv[t, t] <- v[t]
			for (s in seq_len(t - 1L)) {
				cv[s, t] <- cv[t, s] <- rho^(t - s) * v[s]
			}
		}
	}
	cv
}

# the root of the prior of a path (rk coordinates by nt periods, stored column by column): kron(L, I_rk)
.lame_dpm_root <- function(nt, rk, rho, sigma, lambda = NULL) {
	ch <- tryCatch(chol(.lame_dpm_cov(nt, rho, sigma, lambda)), error = function(e) NULL)
	if (is.null(ch)) return(NULL)
	kronecker(t(ch), diag(rk))
}

# shift of the predictor of a line's cells when its actor's position moves by
# `d` (rk by nt); `w` the partner weights (cells by rk by nt)
.lame_dpm_delta <- function(w, d) {
	nl <- dim(w)[1L]
	nt <- dim(w)[3L]
	out <- matrix(0, nl, nt)
	for (r in seq_len(dim(w)[2L])) out <- out + matrix(w[, r, ], nl, nt) * rep(d[r, ], each = nl)
	out
}

# Poisson factor of a line's cells when its actor's position moves by `d`: `yl`
# the counts (zero without a datum), `el` exp(z + offset) on observed cells
# (zero elsewhere), `w` the partner weights
.lame_dpm_poisson_ll <- function(yl, el, w, d) {
	dl <- .lame_dpm_delta(w, d)
	sum(yl * dl - el * expm1(dl))
}

# probit factor of a line's observed cells when its actor's position moves by
# `d`: `b` their predictor (with any shared dyad factor), `sg` their signs (+1
# event, -1 non-event), `idx` their positions in the line (cells by periods),
# `w` the partner weights, `sd_cell` the cell sd given the shared factor
.lame_dpm_probit_ll <- function(b, sg, idx, w, d, sd_cell = 1) {
	sum(pnorm(sg * ((b + .lame_dpm_delta(w, d)[idx]) / sd_cell), log.p = TRUE))
}

# visit every path; shared by the families. `mk_ll(i, row, w, ezl, zl, ezc)`
# returns the log likelihood of actor i's line as a function of the path change
# `d` (rk by nt): `w` the partner weights, `ezl` and `zl` the predictor and
# latent lines (NULL when Z is not carried), `ezc` the column a symmetric
# network's pair mean reads. `carry` moves Z with the predictor (Poisson).
# `U` (nr x RA x N), `V` (nc x RB x N) the paths (symmetric: `V` NULL, `Ld` the
# weights; bipartite: `G` RA by RB); `max_abs` the rail; `lambda_u`, `lambda_v`
# the local innovation precisions of the t kind
.lame_dpm_drive <- function(cx, Z, EZ, U, V, G, Ld, rho_uv, sigma_uv, max_abs, lambda_u, lambda_v, mk_ll, carry) {
	nr <- cx$nr
	nc <- cx$nc
	nt <- cx$N
	sym <- identical(cx$mode, "sym")
	bip <- identical(cx$mode, "bip")
	stats <- c(dpm_calls = 0, dpm_evals = 0, dpm_stuck = 0)
	state <- function() list(Z = Z, EZ = EZ, U = U, V = V, stats = stats)
	if (!(is.finite(rho_uv) && abs(rho_uv) < 1 && is.finite(sigma_uv) && sigma_uv > 0)) return(state())
	rail <- is.finite(max_abs) && max_abs > 0
	local_scales <- !is.null(lambda_u) || !is.null(lambda_v)
	root_ar <- NULL
	rk_u <- dim(U)[2L]
	rk_v <- if (sym) rk_u else dim(V)[2L]
	if (!local_scales) {
		root_ar <- lapply(unique(c(rk_u, rk_v)), function(k) .lame_dpm_root(nt, k, rho_uv, sigma_uv))
		names(root_ar) <- as.character(unique(c(rk_u, rk_v)))
		if (any(vapply(root_ar, is.null, logical(1)))) return(state())
	}
	for (row in if (sym) TRUE else c(TRUE, FALSE)) {
		n_act <- if (row) nr else nc
		rk <- if (row) rk_u else rk_v
		# partner weights (line cells by coordinates by periods), fixed over the phase
		# unless the network is symmetric
		wfull <- NULL
		if (bip) {
			wfull <- if (row) array(
				vapply(seq_len(nt), function(t) matrix(V[, , t], nc, rk_v) %*% t(G), numeric(nc * rk_u)),
				c(nc, rk_u, nt)
			)
			else array(
				vapply(seq_len(nt), function(t) matrix(U[, , t], nr, rk_u) %*% G, numeric(nr * rk_v)),
				c(nr, rk_v, nt)
			)
		} else if (!sym) wfull <- if (row) V else U
		for (i in sample.int(n_act)) {
			cur <- matrix(if (row) U[i, , ] else V[i, , ], rk, nt)
			w <- wfull
			if (sym) {
				# pairs {i, j}: weights l * u_j of the other actors (the self cell has none)
				w <- U
				for (r in seq_len(rk)) w[, r, ] <- Ld[r] * U[, r, ]
				w[i, , ] <- 0
			}
			lam <- if (row) lambda_u else lambda_v
			root <- if (local_scales) .lame_dpm_root(
				nt,
				rk,
				rho_uv,
				sigma_uv,
				if (is.null(lam)) NULL else lam[i, ]
			) else root_ar[[as.character(rk)]]
			if (is.null(root)) next
			n_line <- if (row) nc else nr
			ll <- mk_ll(
				i,
				row,
				w,
				matrix(if (row) EZ[i, , ] else EZ[, i, ], n_line, nt),
				if (carry) matrix(if (row) Z[i, , ] else Z[, i, ], n_line, nt),
				if (sym) matrix(EZ[, i, ], n_line, nt)
			)
			loglik <- function(x) {
				up <- matrix(x, rk, nt)
				if (rail && any(abs(up) > max_abs)) return(-Inf)
				ll(up - cur)
			}
			ans <- .lame_ess_step(as.vector(cur), numeric(rk * nt), root, loglik, ll(matrix(0, rk, nt)))
			stats["dpm_calls"] <- stats["dpm_calls"] + 1
			stats["dpm_evals"] <- stats["dpm_evals"] + ans$evals
			stats["dpm_stuck"] <- stats["dpm_stuck"] + ans$stuck
			new <- matrix(ans$x, rk, nt)
			d <- new - cur
			if (!any(d != 0)) next
			dl <- .lame_dpm_delta(w, d)
			if (sym) dl[i, ] <- colSums(Ld * (new^2 - cur^2))
			# move the line by dl; a symmetric network moves the pairs {i, j} in both
			# triangles and the diagonal cell once (dl[i, ], not repeated by the column)
			if (sym) {
				dc <- dl
				dc[i, ] <- 0
				EZ[i, , ] <- EZ[i, , ] + dl
				if (carry) Z[i, , ] <- Z[i, , ] + dl
				for (t in seq_len(nt)) {
					EZ[, i, t] <- EZ[, i, t] + dc[, t]
					if (carry) Z[, i, t] <- Z[, i, t] + dc[, t]
				}
			} else if (row) {
				EZ[i, , ] <- EZ[i, , ] + dl
				if (carry) Z[i, , ] <- Z[i, , ] + dl
			} else {
				EZ[, i, ] <- EZ[, i, ] + dl
				if (carry) Z[, i, ] <- Z[, i, ] + dl
			}
			if (row) U[i, , ] <- new else V[i, , ] <- new
		}
	}
	state()
}

# Poisson visit. `cx` from .lame_ltm_context(); `Z`, `EZ` the latent array and
# its mean (nr x nc x N); other arguments as for .lame_dpm_drive()
.lame_dpm_visit <- function(
	cx,
	Z,
	EZ,
	U,
	V = NULL,
	G = NULL,
	Ld = NULL,
	rho_uv,
	sigma_uv,
	max_abs = Inf,
	lambda_u = NULL,
	lambda_v = NULL
) {
	nr <- cx$nr
	nc <- cx$nc
	nt <- cx$N
	off <- cx$offset
	obs <- cx$obs3
	yv <- array(ifelse(obs, cx$Y, 0), dim(obs))
	mk_ll <- function(i, row, w, ezl, zl, ezc) {
		n_line <- if (row) nc else nr
		yl <- matrix(if (row) yv[i, , ] else yv[, i, ], n_line, nt)
		ol <- matrix(if (row) obs[i, , ] else obs[, i, ], n_line, nt)
		el <- exp(zl + rep(off, each = n_line))
		el[!ol] <- 0
		function(d) .lame_dpm_poisson_ll(yl, el, w, d)
	}
	.lame_dpm_drive(cx, Z, EZ, U, V, G, Ld, rho_uv, sigma_uv, max_abs, lambda_u, lambda_v, mk_ll, carry = TRUE)
}

# index sets of the binary collapsed update for the network's mode (observed
# cells of every line with their signs) plus the fields the driver reads;
# `base` reuses an existing context
.lame_dpm_binary_context <- function(Y, mode, base = NULL) {
	cx <- base
	if (is.null(cx)) cx <- if (identical(mode, "bip")) .lame_actor_move_context_bip(Y)
	else if (identical(mode, "sym")) .lame_actor_move_context_sym(Y)
	else .lame_actor_move_context(Y)
	cx$mode <- mode
	cx$nr <- if (identical(mode, "bip")) cx$nA else cx$n
	cx$nc <- if (identical(mode, "bip")) cx$nB else cx$n
	cx
}

# binary visit. `cx` from .lame_dpm_binary_context(); `rho` the dyadic
# correlation (read by the symmetric latent kernel; unused when bipartite). the
# latent array is integrated out of the visit and redrawn at the end
.lame_dpm_visit_binary <- function(
	cx,
	Z,
	EZ,
	Y,
	rho,
	U,
	V = NULL,
	G = NULL,
	Ld = NULL,
	rho_uv,
	sigma_uv,
	max_abs = Inf,
	lambda_u = NULL,
	lambda_v = NULL
) {
	nr <- cx$nr
	nc <- cx$nc
	nt <- cx$N
	sym <- identical(cx$mode, "sym")
	bip <- identical(cx$mode, "bip")
	off <- if (!sym && !bip) .lame_dyad_factor_draw(Z, EZ, rho, cx$upper) else NULL
	sd_cell <- if (!sym && !bip) sqrt(1 - abs(rho)) else 1
	mk_ll <- function(i, row, w, ezl, zl, ezc) {
		n_line <- if (row) nc else nr
		sd_ <- if (sym) cx$rows[[i]] else if (row) cx$rows[[i]] else cx$cols[[i]]
		base <- if (sym) (ezl + ezc) / 2
		else if (bip) ezl
		else ezl + matrix(if (row) off[i, , ] else off[, i, ], n_line, nt)
		b <- base[sd_$idx]
		function(d) .lame_dpm_probit_ll(b, sd_$sg, sd_$idx, w, d, sd_cell)
	}
	o <- .lame_dpm_drive(cx, Z, EZ, U, V, G, Ld, rho_uv, sigma_uv, max_abs, lambda_u, lambda_v, mk_ll, carry = FALSE)
	EZ <- o$EZ
	if (sym) {
		for (t in seq_len(nt)) Z[, , t] <- rZ_bin_sym_fc(Z[, , t], EZ[, , t], Y[, , t], rho)
	} else if (bip) {
		Z <- rZ_bin_bip_batch_cpp(Z, EZ, Y)
	} else {
		Z <- .lame_binary_z_given_factor(Z, EZ, off, rho, cx$i1, cx$i0, cx$ina, cx$idg)
	}
	o$Z <- Z
	o
}

# symmetric gauge: coordinate paths have unit stationary variance and the
# weights l (l_r ~ N(0, 1)) carry the scale and sign of each direction of
# u_t diag(l) u_t'. two slow directions get exact moves. the scale ridge:
# scaling coordinate r by c and l_r by 1 / c^2 leaves the predictor unchanged,
# so a Metropolis step on log c uses the prior ratio and the Jacobian
# c^(n N - 2), rejecting paths outside the rail. the persistence: random-walk
# Metropolis on rho_uv given the transitions, with its proposal scale from the
# information at the least-squares persistence (a function of the paths alone,
# so the walk stays symmetric)

# Metropolis draw of the gauge persistence; U_cube n x R x N paths. target as in .lame_rho_uv_unit_mh()
.lame_dpm_rho_unit_mh <- function(U_cube, rho_current, prior_mean = 0.9, prior_sd = 0.1, scales = c(2.4, 0.6, 6)) {
	N <- dim(U_cube)[3L]
	if (is.null(N) || N < 2L) return(rho_current)
	n <- dim(U_cube)[1L]
	R <- dim(U_cube)[2L]
	s11 <- 0; s12 <- 0; s22 <- 0
	for (t in 2:N) {
		U1 <- matrix(U_cube[, , t - 1L], n, R)
		U2 <- matrix(U_cube[, , t], n, R)
		s11 <- s11 + sum(U1 * U1)
		s12 <- s12 + sum(U1 * U2)
		s22 <- s22 + sum(U2 * U2)
	}
	m <- n * R * (N - 1L)
	if (!is.finite(s11) || s11 <= 0) return(rho_current)
	log_post <- function(r) {
		v <- 1 - r * r
		-0.5 * m * log(v) - (s22 - 2 * r * s12 + r * r * s11) / (2 * v) + stats::dnorm(
			r,
			prior_mean,
			prior_sd,
			log = TRUE
		)
	}
	lo <- -0.99
	hi <- 0.99
	rc <- min(max(rho_current, lo), hi)
	rhat <- min(max(s12 / s11, -0.98), 0.98)
	info <- m * (1 + rhat^2) / (1 - rhat^2)^2
	for (sc in scales) {
		prop <- rc + stats::rnorm(1, 0, min(max(sc / sqrt(info), 1e-3), 0.2))
		if (prop <= lo || prop >= hi) next
		la <- log_post(prop) - log_post(rc)
		if (is.finite(la) && log(stats::runif(1)) < la) rc <- prop
	}
	rc
}

# scale ridge move. `U` the n x R x N paths, `Ld` the R weights (prior variance
# `l_var`), `max_abs` the rail; the product is unchanged
.lame_dpm_sym_scale <- function(U, Ld, rho_uv, max_abs = Inf, taus = c(0.02, 0.1, 0.4), l_var = 1) {
	nr <- dim(U)[1L]
	rk <- dim(U)[2L]
	nt <- dim(U)[3L]
	stats <- c(dpm_scale_attempts = 0, dpm_scale_accepts = 0)
	if (!(is.finite(rho_uv) && abs(rho_uv) < 1)) return(list(U = U, Ld = Ld, stats = stats))
	for (r in seq_len(rk)) {
		u <- matrix(U[, r, ], nr, nt)
		q <- sum(u[, 1L]^2)
		if (nt > 1L) q <- q + sum((u[, -1L, drop = FALSE] - rho_uv * u[, -nt, drop = FALSE])^2) / (1 - rho_uv^2)
		for (tau in taus) {
			cc <- exp(rnorm(1, 0, tau))
			stats["dpm_scale_attempts"] <- stats["dpm_scale_attempts"] + 1
			if (is.finite(max_abs) && cc * max(abs(U[, r, ])) > max_abs) next
			la <- -(cc^2 - 1) / 2 * q - 0.5 * Ld[r]^2 / l_var * (cc^-4 - 1) + (nr * nt - 2) * log(cc)
			if (log(runif(1)) < la) {
				U[, r, ] <- cc * U[, r, ]
				Ld[r] <- Ld[r] / cc^2
				q <- q * cc^2
				stats["dpm_scale_accepts"] <- stats["dpm_scale_accepts"] + 1
			}
		}
	}
	list(U = U, Ld = Ld, stats = stats)
}
