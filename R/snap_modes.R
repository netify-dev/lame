# the snap kind of dynamic_uv on bipartite and symmetric panels. the reset-or-drift prior of each path,
#     m(y | x) = (1 - pi) N(y; rho x, sigma^2 I) + pi N(y; 0, kappa^2 I),
# is the one of directed panels: a bipartite panel has its row paths (nA x RA x N) and its column paths
# (nB x RB x N), each with its own resets, under one (rho_uv, sigma_uv, pi_snap) and the product
# U_t G V_t'; a symmetric panel has one set of unit-variance paths (sigma_uv^2 = 1 - rho_uv^2, so the
# linked reset scale is 1) and the product U_t diag(L) U_t'. given the resets the path prior is Gaussian,
# exp(-Q / (2 sigma^2)) with Q the transitions weighted by w = sigma^2 / kappa^2 and k = 0 at a reset
# (.lame_uv_dyn_kind_weights), so every gauge move of the ar1 kind that acts on the paths alone keeps its
# derivation with those weights:
# - the bipartite scale split and common gauge (R/uv_split_dynamic.R, R/uv_gauge_dynamic.R) read Q through
#   the weighted moments;
# - the period signs (.lame_uv_sign_ar1) and the signs along the singular directions of G
#   (.lame_uv_svd_sign_ar1) see only the cross products of drifting transitions, a reset position being
#   symmetric about zero;
# - the symmetric split of a column against its weight reads the weighted moment, and the period signs of
#   one coordinate (.lame_uv_sym_sign_draw) the cross products of drifting transitions;
# - the symmetric persistence sees the drifting transitions only, as the first period and a reset position
#   (scale 1 when linked, prior$snap_kappa otherwise) do not involve rho_uv; its non-centred draw
#   (R/uv_nc_step.R) holds the reset positions, whose standardised values are N(0, 1) whatever rho_uv;
# - the joint scale of bipartite paths, sigma_uv and G is a scale family under the linked reset scale, as
#   under ar1; a fixed scale adds the prior ratio of the reset positions (.lame_snap_scale_extra).
# the mirror of the persistence (R/uv_rho_flip.R) holds under the snap prior as it stands

# signs (+1 event, -1 non-event, 0 no datum) of the cells of a bipartite or symmetric binary context
# (R/binary_actor_collapsed_move_bip.R, R/binary_actor_collapsed_move.R), as an array of dimension `d`
.lame_snap_sign_mask <- function(cx, d) {
	mask <- array(0, d)
	n_line <- d[2L]
	for (i in seq_along(cx$rows)) {
		r <- cx$rows[[i]]
		if (!length(r$idx)) next
		mask[cbind(i, (r$idx - 1L) %% n_line + 1L, (r$idx - 1L) %/% n_line + 1L)] <- r$sg
	}
	mask
}

# path move of a bipartite or symmetric panel; arguments as for .lame_dpm_snap_visit(). a symmetric visit
# returns V_t = U_t diag(Ld) and the resets on both sides
.lame_dpm_snap_visit_mode <- function(cx, Z, EZ, rho, U, V, rho_uv, sigma_uv, snap, max_abs, family, s2, mode, G, Ld, Y) {
	sym <- identical(mode, "sym")
	bip <- identical(mode, "bip")
	still <- list(
		Z = Z,
		EZ = EZ,
		U = U,
		V = V,
		delta_u = snap$delta_u,
		delta_v = if (sym) snap$delta_u else snap$delta_v,
		stats = c(dpm_calls = 0, dpm_evals = 0, dpm_stuck = 0, snap_scales = 0, snap_scales_accepted = 0)
	)
	R_u <- dim(U)[2L]
	ok <- (bip || sym) && is.finite(rho_uv) && abs(rho_uv) < 1 && is.finite(sigma_uv) && sigma_uv > 0 &&
		is.finite(snap$pi) && snap$pi > 0 && snap$pi < 1 && is.finite(snap$kappa) && snap$kappa > 0 &&
		is.finite(s2) && s2 > 0 &&
		(!sym || (length(Ld) == R_u && all(is.finite(Ld)))) &&
		(!bip || (is.matrix(G) && identical(dim(G), c(R_u, dim(V)[2L])) && all(is.finite(G))))
	if (!ok) return(still)
	if (!is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	N <- dim(Z)[3L]
	counts <- numeric(0)
	period_off <- numeric(0)
	if (identical(family, "poisson")) {
		kind <- 1L
		mask <- array(as.numeric(cx$obs3), dim(Z))
		counts <- array(ifelse(cx$obs3, cx$Y, 0), dim(Z))
		period_off <- as.numeric(cx$offset)
	} else if (identical(family, "gaussian")) {
		kind <- 2L
		mask <- array(as.numeric(is.finite(Z)), dim(Z))
		if (sym) for (t in seq_len(N)) diag(mask[, , t]) <- 0
	} else {
		kind <- 0L
		mask <- .lame_snap_sign_mask(cx, dim(Z))
	}
	o <- snap_path_visit_mode_cpp(
		EZ,
		Z,
		numeric(0),
		mask,
		counts,
		period_off,
		U,
		if (sym) U else V,
		if (bip) G else matrix(0, 0L, 0L),
		if (sym) as.numeric(Ld) else numeric(0),
		rho_uv,
		sigma_uv,
		snap$kappa,
		snap$pi,
		snap$delta_u,
		if (sym) snap$delta_u else snap$delta_v,
		1,
		max_abs,
		s2,
		numeric(N),
		kind,
		if (bip) 1L else 2L
	)
	if (sym) {
		o$delta_v <- o$delta_u
		o$V <- o$U
		for (t in seq_len(N)) o$V[, , t] <- matrix(o$U[, , t], dim(U)[1L], R_u) %*% diag(Ld, R_u)
	}
	if (kind == 0L) {
		if (bip) {
			o$Z <- rZ_bin_bip_batch_cpp(Z, o$EZ, Y)
		} else {
			for (t in seq_len(N)) o$Z[, , t] <- rZ_bin_sym_fc(Z[, , t], o$EZ[, , t], Y[, , t], rho)
		}
	}
	o
}

# gibbs sweep of the positions and resets of a bipartite snap fit (src/snap_modes.cpp), then the period
# signs and the signs along the singular directions of G, each seeing the drifting transitions only
.lame_snap_bip_kernel <- function(U_cube, V_cube, E, G, rho_uv, sigma_uv, s2, kappa, pi_snap, delta_u, delta_v, max_abs) {
	o <- rUV_dynamic_snap_bip_fc_cpp(
		U_cube,
		V_cube,
		E,
		G,
		rho_uv,
		sigma_uv,
		s2,
		kappa,
		pi_snap,
		delta_u,
		delta_v,
		max_abs = max_abs
	)
	keep_u <- 1 - o$delta_u
	keep_v <- 1 - o$delta_v
	o <- .lame_uv_svd_sign_ar1(o, G, rho_uv, sigma_uv, max_abs = max_abs, weight_u = keep_u, weight_v = keep_v)
	.lame_uv_sign_ar1(o, rho_uv, sigma_uv, keep_u, keep_v)
}

# draw of the reset rate given the resets of every eligible transition (both sides of a bipartite panel)
.lame_snap_pi_draw <- function(delta_u, delta_v, pi_a, pi_b) {
	N <- ncol(delta_u)
	n_snap <- sum(delta_u[, -1L]) + sum(delta_v[, -1L])
	n_elig <- (nrow(delta_u) + nrow(delta_v)) * (N - 1L)
	stats::rbeta(1, pi_a + n_snap, pi_b + (n_elig - n_snap))
}

# persistence and innovation scale of bipartite snap paths given the resets, as for directed snap paths, under
# the variance prior of `prior` (R/uv_var_prior.R); returns list(rho, sigma, u) with the recycled uniform of the
# persistence draw for the mirror move
.lame_snap_hyper_step <- function(U_cube, V_cube, rho_uv, sigma_uv, delta_u, delta_v, prior, stat_var_cap) {
	linked <- is.null(prior$snap_kappa)
	rho_draw <- .lame_sample_rho_uv_drift_only(
		U_cube,
		V_cube,
		sigma_uv,
		rho_uv,
		FALSE,
		delta_u,
		delta_v,
		prior_mean = prior$rho_uv_mean,
		prior_sd = prior$rho_uv_sd,
		reset_stationary = linked,
		stat_var_cap = stat_var_cap,
		with_u = TRUE,
		stationary_prior = if (.lame_uv_var_stationary(prior)) c(prior$uv_var_shape, prior$uv_var_scale)
	)
	sigma <- .lame_sample_sigma_uv_drift_only(
		U_cube,
		V_cube,
		rho_draw$rho,
		sigma_uv,
		FALSE,
		delta_u,
		delta_v,
		.lame_uv_ig_shape(prior),
		.lame_uv_ig_scale(prior, rho_draw$rho),
		sigma_max = prior$sigma_uv_max %||% Inf,
		stat_var_cap = stat_var_cap,
		reset_stationary = linked
	)
	list(rho = rho_draw$rho, sigma = sigma, u = rho_draw$u)
}

# the joint scale step of bipartite paths, sigma_uv and G (R/uv_nc_step.R) holds under the snap prior: with the
# linked reset scale every factor of a path is a scale family in sigma_uv, as under ar1. a fixed reset scale
# kappa does not scale, so each reset position x adds -(c^2 - 1) ||x||^2 / (2 kappa^2) + R log c. NULL when
# nothing is added (ar1, or snap with the linked scale)
.lame_snap_scale_extra <- function(kappa_fixed, U_cube, V_cube, delta_u, delta_v) {
	if (is.null(kappa_fixed) || is.null(delta_u)) return(NULL)
	st <- .lame_snap_reset_ss(U_cube, V_cube, FALSE, delta_u, delta_v)
	k2 <- kappa_fixed * kappa_fixed
	function(x) -(exp(2 * x) - 1) * st$ss / (2 * k2) + st$cnt * x
}

# orientation of symmetric snap paths. with v = u |L|^(1/2) the product is sum_r sign(l_r) v_ir v_jr, so for two
# coordinates r, q the rotations (signs alike) or the boosts (signs unlike) B of the pair leave it unchanged; in the
# coordinates of the paths the map is M = D B D^-1, D = diag(|l_r|, |l_q|)^(1/2), with det M = 1. applied to every
# period from t0 on it changes the path prior only, through the weighted moments of the transitions (the first
# period's too when t0 = 1), so the angle or rapidity is drawn by a slice step along its one-parameter group (Haar
# measure, unit Jacobian; Liu and Sabatti 2000), cut at max_abs. a reset position enters the moments at weight
# sigma^2 / kappa^2 with no transition from the period before; the conditional draws turn a period's orientation
# only one actor at a time. `U_cube` the paths (n x R x N), `Ld` the weights, `delta` the resets, `snap_w` the
# reset weight; returns the paths and the number of steps that moved
.lame_snap_sym_orient_draw <- function(U_cube, Ld, rho, sigma, delta, snap_w, max_abs = Inf, width = 0.25) {
	d <- dim(U_cube)
	n <- d[1L]; R <- d[2L]; N <- d[3L]
	still <- list(U = U_cube, moved = 0L)
	if (R < 2L || !all(is.finite(U_cube)) || !all(is.finite(Ld)) || any(Ld == 0)) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	s2 <- sigma * sigma
	k <- 1 * (delta < 0.5)
	k[, 1L] <- 1
	w <- ifelse(k > 0.5, 1, snap_w)
	moved <- 0L
	for (r in seq_len(R - 1L)) for (q in (r + 1L):R) {
		cols <- c(r, q)
		boost <- sign(Ld[r]) != sign(Ld[q])
		dg <- sqrt(abs(Ld[cols]))
		map <- function(a) {
			B <- if (boost) matrix(c(cosh(a), sinh(a), sinh(a), cosh(a)), 2L) else matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2L)
			diag(dg) %*% B %*% diag(1 / dg)
		}
		for (t0 in seq_len(N)) {
			X <- U_cube[, cols, , drop = FALSE]
			# weighted moments of the transitions inside the tail and of the one into it
			S <- matrix(0, 2L, 2L)
			if (t0 == 1L) S <- S + (1 - rho * rho) * crossprod(X[, , 1L], w[, 1L] * X[, , 1L])
			if (t0 < N) for (t in (t0 + 1L):N) {
				e <- X[, , t] - rho * k[, t] * X[, , t - 1L]
				S <- S + crossprod(e, w[, t] * e)
			}
			A <- C <- matrix(0, 2L, 2L)
			if (t0 > 1L) {
				A <- crossprod(X[, , t0], w[, t0] * X[, , t0])
				C <- rho * crossprod(X[, , t0], (w[, t0] * k[, t0]) * X[, , t0 - 1L])
			}
			xmax <- max(abs(X[, , t0:N]))
			lg <- function(a) {
				M <- map(a)
				if (is.finite(max_abs) && xmax * max(abs(M)) * sqrt(2) > max_abs &&
						max(abs(apply(X[, , t0:N, drop = FALSE], 3L, function(x) x %*% M))) > max_abs) return(-Inf)
				-(sum(M * ((S + A) %*% M)) - 2 * sum(M * C)) / (2 * s2)
			}
			lim <- if (boost) c(-Inf, Inf) else c(-pi, pi)
			a <- .lame_nc_slice(lg, 0, lim[1L], lim[2L], w = width)$x
			if (!is.finite(a) || a == 0) next
			M <- map(a)
			for (t in t0:N) U_cube[, cols, t] <- matrix(U_cube[, cols, t], n, 2L) %*% M
			moved <- moved + 1L
		}
	}
	list(U = U_cube, moved = moved)
}
