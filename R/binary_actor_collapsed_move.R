# collapsed per-actor update of the binary sampler. the Z-conditional draw of an actor whose cells
# carry no event is as narrow as its data precision, so the augmented chain random-walks across a
# wide posterior; here each actor's blocks are drawn with Z integrated out (dyad errors split into a
# shared factor f drawn given Z, the rest integrated out, leaving one probit factor per observed cell;
# Z is redrawn given f afterwards) by elliptical slice steps (Murray, Adams and MacKay 2010).
# per actor, in turn: (a) sender and (b) receiver level translations, the centring compensated through
# the intercept (a whole intercept path under dynamic_beta), then (c) sender and (d) receiver positions
# under a per-visit draw of Suv (tau_r when symmetric). positions are left to the sampler in the diffuse
# phase of lame() burn-in; missing cells carry no likelihood factor. interweaving: Yu and Meng (2011)

# one elliptical slice step for x with prior N(mu, root root'), current log likelihood ll0 and log
# likelihood function `loglik`
.lame_ess_step <- function(x, mu, root, loglik, ll0) {
	nu <- drop(root %*% rnorm(length(x)))
	xc <- x - mu
	thr <- ll0 + log(runif(1))
	theta <- runif(1, 0, 2 * pi)
	lo <- theta - 2 * pi
	hi <- theta
	evals <- 0L
	repeat {
		xp <- mu + xc * cos(theta) + nu * sin(theta)
		llp <- loglik(xp)
		evals <- evals + 1L
		if (is.finite(llp) && llp > thr) return(list(x = xp, ll = llp, evals = evals, stuck = FALSE))
		if (theta < 0) lo <- theta else hi <- theta
		if (evals >= 100L) return(list(x = x, ll = ll0, evals = evals, stuck = TRUE))
		theta <- runif(1, lo, hi)
	}
}

# index sets for the per-actor likelihood and latent refresh: for row i's (column j's) n x N slice,
# the positions of its observed off-diagonal cells, the partner of each and its sign (+1 event, -1 non-event)
.lame_actor_move_context <- function(Y) {
	n <- dim(Y)[1L]
	N <- dim(Y)[3L]
	obs <- is.finite(Y)
	dg <- array(FALSE, dim = dim(Y))
	dg[cbind(rep(seq_len(n), N), rep(seq_len(n), N), rep(seq_len(N), each = n))] <- TRUE
	obs <- obs & !dg
	sgn <- array(0, dim = dim(Y))
	sgn[obs] <- 2 * Y[obs] - 1
	side <- function(i, row) {
		o <- if (row) matrix(obs[i, , ], n, N) else matrix(obs[, i, ], n, N)
		s <- if (row) matrix(sgn[i, , ], n, N) else matrix(sgn[, i, ], n, N)
		idx <- which(o)
		list(idx = idx, partner = ((idx - 1L) %% n) + 1L, sg = s[idx])
	}
	rows <- lapply(seq_len(n), side, row = TRUE)
	cols <- lapply(seq_len(n), side, row = FALSE)
	list(
		n = n,
		N = N,
		rows = rows,
		cols = cols,
		i1 = which(Y == 1 & !dg),
		i0 = which(Y == 0 & !dg),
		ina = which(is.na(Y) & !dg),
		idg = which(dg),
		upper = upper.tri(matrix(0, n, n))
	)
}

# factor covariance for the positions: an inverse-Wishart draw given the current rows, as rUV_rep_fc
# (lame, `prior` NULL; diffuse n I in the first half of burn-in) or rSuv_fc (ame, `prior` =
# list(df, scale)) draw it
.lame_actor_move_suv <- function(U, V, n, R, shrink, prior = NULL) {
	if (is.null(prior) && !isTRUE(shrink)) return(diag(n, 2L * R))
	W <- cbind(U, V)
	if (is.null(prior)) prior <- .lame_uv_prior_unit(R)
	df0 <- prior$df
	scale0 <- prior$scale
	tryCatch(
		solve(stats::rWishart(1L, n + df0, solve(scale0 + crossprod(W)))[, , 1L]),
		error = function(e) NULL
	)
}
# column precisions of the symmetric positions: the host's Gamma conditional, or 1 / n in the diffuse phase
.lame_actor_move_tau <- function(U, n, R, shrink) {
	if (!isTRUE(shrink)) return(rep(1 / n, R))
	vapply(seq_len(R), function(r) stats::rgamma(1L, (2 + n) / 2, rate = (1 + sum(U[, r]^2)) / 2), numeric(1))
}
# the intercept of a centred fit moves with each level translation, so a shift c adds
# -(q c^2 + 2 l c) / 2 to its log prior; `ic` is NULL for static coefficients (l = (iV0 beta)[1],
# q = iV0[1, 1]), kind "static" for a static intercept with precision prec, or kind "path" for an
# AR(1) intercept path, where q = 1 / v0 + (N - 1) (1 - rho)^2 / sd2 and
# l = (beta_1 - m0) / v0 + (1 - rho) sum_t e_t / sd2
.lame_actor_move_icpt <- function(ic, iV0, beta) {
	if (is.null(ic)) return(c(l = unname(drop(iV0 %*% beta)[1L]), q = unname(iV0[1L, 1L])))
	if (identical(ic$kind, "static")) return(c(l = unname(ic$prec * ic$value), q = unname(ic$prec)))
	path <- ic$path
	k <- length(path)
	e <- if (k > 1L) path[-1L] - ic$rho * path[-k] else numeric(0)
	c(
		l = unname((path[1L] - ic$m0) / ic$v0 + (1 - ic$rho) * sum(e) / ic$sd2),
		q = unname(1 / ic$v0 + (k - 1L) * (1 - ic$rho)^2 / ic$sd2)
	)
}
# intercept after a shift by `by` (beta[1] for static coefficients, else the value or path in
# `ic`); returns list(beta, ic)
.lame_actor_move_icpt_shift <- function(ic, beta, by) {
	if (is.null(ic)) beta[1L] <- beta[1L] + by
	else if (identical(ic$kind, "static")) ic$value <- ic$value + by
	else ic$path <- ic$path + by
	list(beta = beta, ic = ic)
}
# log prior of a level translation of actor i in its shift delta, -(a2 delta^2 + 2 a1 delta) / 2 + const,
# from (a_k, b_k) ~ N(0, Sab) and, for centred lame levels, the intercept prior; the centred shift is
# a + delta (e_i - 1/n) with beta_0 + delta / n, the unconstrained (ame) shift a + delta
.lame_actor_move_level_prior <- function(i, row, la, lb, iS, iV0, beta, n, centred = TRUE, ic = NULL) {
	# unnamed elements, so labelled start values do not rename the result
	li <- unname(la[i])
	bi <- unname(lb[i])
	a1 <- if (row) iS[1L, 1L] * li + iS[1L, 2L] * bi else iS[2L, 2L] * bi + iS[1L, 2L] * li
	a2 <- if (row) iS[1L, 1L] else iS[2L, 2L]
	if (centred) {
		lq <- .lame_actor_move_icpt(ic, iV0, beta)
		a1 <- a1 + lq[["l"]] / n
		a2 <- a2 * (n - 1) / n + lq[["q"]] / n^2
	}
	c(a1 = a1, a2 = a2)
}
# the same for a symmetric network: one effect per actor with variance va; the centred shift pairs
# with beta_0 + 2 delta / n
.lame_actor_move_level_prior_sym <- function(i, la, va, iV0, beta, n, centred = TRUE, ic = NULL) {
	a1 <- unname(la[i]) / va
	a2 <- 1 / va
	if (centred) {
		lq <- .lame_actor_move_icpt(ic, iV0, beta)
		a1 <- a1 + 2 * lq[["l"]] / n
		a2 <- a2 * (n - 1) / n + 4 * lq[["q"]] / n^2
	}
	c(a1 = a1, a2 = a2)
}
# one visit: `la`, `lb` the additive levels (static a, b or dynamic a_lev, b_lev; centred in lame,
# free in ame), `a_mat` / `b_mat` the dynamic paths (NULL when static), U and V the static positions,
# `centred` the lame constraint, `uv_prior` the ame factor prior, `ic` the dynamic_beta intercept
# prior and `lev` which margins (sender, receiver) are in the model; returns the state and counters
.lame_actor_move_visit <- function(
	Z,
	EZ,
	rho,
	beta,
	la,
	lb,
	a_mat,
	b_mat,
	U,
	V,
	Sab,
	iV0,
	shrink,
	cx,
	centred = TRUE,
	uv_prior = NULL,
	ic = NULL,
	lev = c(TRUE, TRUE)
) {
	n <- cx$n
	N <- cx$N
	R <- if (is.null(U)) 0L else ncol(U)
	dyn <- !is.null(a_mat)
	sd_cell <- sqrt(1 - abs(rho))
	stats <- c(level_calls = 0, level_evals = 0, uv_calls = 0, uv_evals = 0, stuck = 0)
	iS <- if (any(lev)) tryCatch(solve(Sab), error = function(e) NULL) else diag(2L)
	# an unusable covariance leaves the state to the conditional updates this visit
	if (is.null(iS) || any(!is.finite(iS)))
		return(list(
			Z = Z,
			EZ = EZ,
			beta = beta,
			ic = ic,
			la = la,
			lb = lb,
			a_mat = a_mat,
			b_mat = b_mat,
			U = U,
			V = V,
			stats = stats
		))
	off <- .lame_dyad_factor_draw(Z, EZ, rho, cx$upper)
	# prior of the positions given the factor covariance; positions are left to the sampler during the
	# diffuse first half of lame() burn-in, where uninformed positions can drift to any scale
	uv <- R > 0L && (isTRUE(shrink) || !is.null(uv_prior))
	if (uv) {
		pri <- tryCatch({
			Suv <- .lame_actor_move_suv(U, V, n, R, shrink, uv_prior)
			P <- solve(Suv)
			iu <- seq_len(R); iv <- R + seq_len(R)
			cov_u <- solve(P[iu, iu, drop = FALSE]); cov_v <- solve(P[iv, iv, drop = FALSE])
			cov_u <- (cov_u + t(cov_u)) / 2; cov_v <- (cov_v + t(cov_v)) / 2
			list(
				root_u = t(chol(cov_u)),
				root_v = t(chol(cov_v)),
				reg_u = -cov_u %*% P[iu, iv, drop = FALSE], # mean of u given v: reg_u %*% v
				reg_v = -cov_v %*% P[iv, iu, drop = FALSE]
			)
		}, error = function(e) NULL)
		if (is.null(pri) || any(!is.finite(unlist(pri)))) uv <- FALSE
	}
	ez_off <- function(i, row) {
		if (row) matrix(EZ[i, , ], n, N) + matrix(off[i, , ], n, N)
		else matrix(EZ[, i, ], n, N) + matrix(off[, i, ], n, N)
	}
	w_of <- function(i) {
		if (centred) { w <- rep(-1 / n, n); w[i] <- 1 - 1 / n } else { w <- numeric(n); w[i] <- 1 }
		w
	}
	for (i in sample.int(n)) {
		# (a) sender level, (b) receiver level
		for (row in c(TRUE, FALSE)) {
			if (!lev[if (row) 1L else 2L]) next
			sd_ <- if (row) cx$rows[[i]] else cx$cols[[i]]
			base <- ez_off(i, row)[sd_$idx]
			ll_level <- function(d) sum(pnorm(sd_$sg * ((base + d) / sd_cell), log.p = TRUE))
			pr <- .lame_actor_move_level_prior(i, row, la, lb, iS, iV0, beta, n, centred, ic)
			ans <- .lame_ess_step(0, -pr[["a1"]] / pr[["a2"]], matrix(1 / sqrt(pr[["a2"]])), ll_level, ll_level(0))
			stats["level_calls"] <- stats["level_calls"] + 1
			stats["level_evals"] <- stats["level_evals"] + ans$evals
			stats["stuck"] <- stats["stuck"] + ans$stuck
			d <- ans$x
			if (d != 0) {
				w <- w_of(i) * d
				if (row) {
					la <- la + w
					if (dyn) a_mat <- a_mat + w
					EZ[i, , ] <- EZ[i, , ] + d
				} else {
					lb <- lb + w
					if (dyn) b_mat <- b_mat + w
					EZ[, i, ] <- EZ[, i, ] + d
				}
				if (centred) {
					sh <- .lame_actor_move_icpt_shift(ic, beta, d / n)
					beta <- sh$beta; ic <- sh$ic
				}
			}
		}
		if (uv) {
			# (c) sender position u_i, (d) receiver position v_i
			for (row in c(TRUE, FALSE)) {
				sd_ <- if (row) cx$rows[[i]] else cx$cols[[i]]
				cur <- if (row) U[i, ] else V[i, ]
				opp <- if (row) V else U
				base <- ez_off(i, row)[sd_$idx] - drop(opp %*% cur)[sd_$partner]
				ll_uv <- function(x) sum(pnorm(
					sd_$sg * ((base + drop(opp %*% x)[sd_$partner]) / sd_cell),
					log.p = TRUE
				))
				other <- if (row) V[i, ] else U[i, ]
				mu <- drop((if (row) pri$reg_u else pri$reg_v) %*% other)
				ans <- .lame_ess_step(cur, mu, if (row) pri$root_u else pri$root_v, ll_uv, ll_uv(cur))
				stats["uv_calls"] <- stats["uv_calls"] + 1
				stats["uv_evals"] <- stats["uv_evals"] + ans$evals
				stats["stuck"] <- stats["stuck"] + ans$stuck
				dx <- ans$x - cur
				if (any(dx != 0)) {
					shift <- drop(opp %*% dx)
					if (row) {
						U[i, ] <- ans$x
						EZ[i, , ] <- EZ[i, , ] + shift
					} else {
						V[i, ] <- ans$x
						EZ[, i, ] <- EZ[, i, ] + shift
					}
				}
			}
		}
	}
	Z <- .lame_binary_z_given_factor(Z, EZ, off, rho, cx$i1, cx$i0, cx$ina, cx$idg)
	list(
		Z = Z,
		EZ = EZ,
		beta = beta,
		ic = ic,
		la = la,
		lb = lb,
		a_mat = a_mat,
		b_mat = b_mat,
		U = U,
		V = V,
		stats = stats
	)
}

# index sets of a symmetric panel: the observed value of pair {i, j} is Y_ij, or Y_ji when that is
# missing; rows[[i]] holds the positions of the observed pairs in row i's slice, partners and signs
.lame_actor_move_context_sym <- function(Y) {
	n <- dim(Y)[1L]
	N <- dim(Y)[3L]
	up <- upper.tri(matrix(0, n, n))
	ysym <- array(NA_real_, dim = dim(Y))
	for (t in seq_len(N)) {
		y <- Y[, , t]
		pair <- y[up]
		alt <- t(y)[up]
		pair[is.na(pair)] <- alt[is.na(pair)]
		m <- matrix(NA_real_, n, n)
		m[up] <- pair
		m[lower.tri(m)] <- t(m)[lower.tri(m)]
		ysym[, , t] <- m
	}
	side <- function(i) {
		y <- matrix(ysym[i, , ], n, N)
		y[i, ] <- NA_real_
		idx <- which(is.finite(y))
		list(idx = idx, partner = ((idx - 1L) %% n) + 1L, sg = 2 * y[idx] - 1)
	}
	list(n = n, N = N, rows = lapply(seq_len(n), side))
}

# one visit for a symmetric binary network: `la` the level (static a or dynamic a_lev), `a_mat` the
# dynamic path (NULL when static), U the positions with weights `Ld` (V = U diag(Ld)), `va` the effect
# variance, `ic` the dynamic_beta intercept prior and `lev` whether the actor effect is in the model;
# Z is refreshed with the host kernel at diagonal variance 1 + rho
.lame_actor_move_visit_sym <- function(
	Z,
	EZ,
	Y,
	rho,
	beta,
	la,
	a_mat,
	U,
	Ld,
	va,
	iV0,
	shrink,
	cx,
	centred = TRUE,
	ic = NULL,
	lev = TRUE
) {
	n <- cx$n
	N <- cx$N
	R <- if (is.null(U)) 0L else ncol(U)
	dyn <- !is.null(a_mat)
	stats <- c(level_calls = 0, level_evals = 0, uv_calls = 0, uv_evals = 0, stuck = 0)
	if (lev && (!is.finite(va) || va <= 0))
		return(list(Z = Z, EZ = EZ, beta = beta, ic = ic, la = la, a_mat = a_mat, U = U, stats = stats))
	# positions are left to the sampler during the diffuse first half of burn-in (see .lame_actor_move_visit)
	uv <- R > 0L && all(is.finite(Ld)) && isTRUE(shrink)
	if (uv) tau <- .lame_actor_move_tau(U, n, R, shrink)
	pm_row <- function(i) (matrix(EZ[i, , ], n, N) + matrix(EZ[, i, ], n, N)) / 2
	for (i in sample.int(n)) {
		sd_ <- cx$rows[[i]]
		# (a) level
		d <- 0
		if (lev) {
			base <- pm_row(i)[sd_$idx]
			ll_level <- function(d) sum(pnorm(sd_$sg * (base + d), log.p = TRUE))
			pr <- .lame_actor_move_level_prior_sym(i, la, va, iV0, beta, n, centred, ic)
			ans <- .lame_ess_step(0, -pr[["a1"]] / pr[["a2"]], matrix(1 / sqrt(pr[["a2"]])), ll_level, ll_level(0))
			stats["level_calls"] <- stats["level_calls"] + 1
			stats["level_evals"] <- stats["level_evals"] + ans$evals
			stats["stuck"] <- stats["stuck"] + ans$stuck
			d <- ans$x
		}
		if (d != 0) {
			w <- numeric(n)
			if (centred) { w[] <- -d / n; w[i] <- d * (1 - 1 / n) } else w[i] <- d
			la <- la + w
			if (dyn) a_mat <- a_mat + w
			dv <- rep(d, n)
			dv[i] <- 0
			EZ[i, , ] <- EZ[i, , ] + dv
			EZ[, i, ] <- EZ[, i, ] + dv
			EZ[i, i, ] <- EZ[i, i, ] + 2 * d
			if (centred) {
				sh <- .lame_actor_move_icpt_shift(ic, beta, 2 * d / n)
				beta <- sh$beta; ic <- sh$ic
			}
		}
		# (b) position u_i
		if (uv) {
			cur <- U[i, ]
			base <- pm_row(i)[sd_$idx] - drop(U %*% (Ld * cur))[sd_$partner]
			ll_uv <- function(x) sum(pnorm(sd_$sg * (base + drop(U %*% (Ld * x))[sd_$partner]), log.p = TRUE))
			ans <- .lame_ess_step(cur, numeric(R), diag(1 / sqrt(tau), R), ll_uv, ll_uv(cur))
			stats["uv_calls"] <- stats["uv_calls"] + 1
			stats["uv_evals"] <- stats["uv_evals"] + ans$evals
			stats["stuck"] <- stats["stuck"] + ans$stuck
			dx <- ans$x - cur
			if (any(dx != 0)) {
				shift <- drop(U %*% (Ld * dx))
				shift[i] <- 0
				EZ[i, , ] <- EZ[i, , ] + shift
				EZ[, i, ] <- EZ[, i, ] + shift
				EZ[i, i, ] <- EZ[i, i, ] + sum(Ld * (ans$x^2 - cur^2))
				U[i, ] <- ans$x
			}
		}
	}
	for (t in seq_len(N)) Z[, , t] <- rZ_bin_sym_fc(Z[, , t], EZ[, , t], Y[, , t], rho)
	list(Z = Z, EZ = EZ, beta = beta, ic = ic, la = la, a_mat = a_mat, U = U, stats = stats)
}
