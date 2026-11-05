# translation moves for the Poisson and ordinal families: an actor's level (or a
# coefficient, scale or position) shifts together with the latent cells it enters,
# so every residual is unchanged and the shift is drawn from its prior times the
# data factor. these break the slow (level, latent row) ridge of sparse panels.
# - actor translation: elliptical slice step for Poisson
#   (Murray, Adams and MacKay 2010), exact truncated normal for ordinal; with an
#   intercept the levels stay centred and the intercept absorbs delta / n
# - ordinal location: shift of the levels' means with the whole latent array
# - Poisson coefficients, dyadic (s2, rho) scale, and static positions
# missing cells carry no data factor and simply move with their row or column

# integer category codes of an ordinal array (1 for the lowest value, NA if missing);
# `Y_int` holds the codes of an explicit-cutpoint fit
.lame_ltm_codes <- function(Y, Y_int = NULL) {
	if (!is.null(Y_int)) {
		W <- array(as.integer(Y_int), dim = dim(Y))
	} else {
		u <- sort(unique(Y[is.finite(Y)]))
		W <- array(match(Y, u), dim = dim(Y))
	}
	W[!is.finite(Y)] <- NA_integer_
	W
}

# data of a fit as the update reads them; a symmetric ordinal panel whose
# triangles differ is not served (ok = FALSE, the host sampler runs alone)
.lame_ltm_context <- function(Y, family, mode, offset = NULL, Y_int = NULL) {
	d <- dim(Y)
	nr <- d[1L]
	nc <- d[2L]
	N <- if (length(d) >= 3L) d[3L] else 1L
	Y <- array(Y, dim = c(nr, nc, N))
	sym <- identical(mode, "sym")
	ok <- TRUE
	if (!identical(mode, "bip")) {
		dg <- cbind(rep(seq_len(nr), N), rep(seq_len(nr), N), rep(seq_len(N), each = nr))
		Y[dg] <- NA_real_
	}
	if (sym) {
		up <- upper.tri(matrix(0, nr, nc))
		for (tp in seq_len(N)) {
			y <- Y[, , tp]
			pair <- y[up]
			alt <- t(y)[up]
			if (identical(family, "ordinal") && (!identical(is.na(y), is.na(t(y))) ||
					any(y != t(y), na.rm = TRUE))) ok <- FALSE
			pair[is.na(pair)] <- alt[is.na(pair)]
			m <- matrix(NA_real_, nr, nc)
			m[up] <- pair
			m[lower.tri(m)] <- t(m)[lower.tri(m)]
			Y[, , tp] <- m
		}
	}
	cx <- list(
		family = family,
		mode = mode,
		nr = nr,
		nc = nc,
		N = N,
		ok = ok,
		offset = if (is.null(offset)) numeric(N) else as.numeric(offset)
	)
	if (identical(family, "ordinal")) {
		cx$W <- .lame_ltm_codes(Y, Y_int)
		cx$K <- if (any(is.finite(Y))) max(cx$W, na.rm = TRUE) else 0L
		# cells of every category in every period (linear positions within the period's matrix)
		cx$cells <- lapply(seq_len(N), function(tp) {
			w <- cx$W[, , tp]
			lapply(seq_len(cx$K), function(k) which(w == k))
		})
	} else {
		cx$Y <- Y
		# cells with a Poisson factor: the observed ones, a symmetric pair once (its upper cell)
		obs <- is.finite(Y)
		if (sym) for (tp in seq_len(N)) obs[, , tp][lower.tri(obs[, , tp])] <- FALSE
		cx$ok_idx <- which(obs)
		cx$y_ok <- Y[cx$ok_idx]
		cx$off_ok <- cx$offset[arrayInd(cx$ok_idx, c(nr, nc, N))[, 3L]]
		# observed cells of each actor's lines, the log exposure of every cell and the
		# counts summed over periods
		cx$obs3 <- is.finite(Y)
		cx$offarr <- array(rep(cx$offset, each = nr * nc), c(nr, nc, N))
		cx$ysj <- rowSums(ifelse(cx$obs3, Y, 0), dims = 2L)
		cx$ys_row <- rowSums(cx$ysj)
		cx$ys_col <- colSums(cx$ysj)
	}
	if (!identical(mode, "bip")) {
		# upper cell (i < j) of every dyad and its mirror, in every period
		up <- which(upper.tri(matrix(0, nr, nc)))
		ij <- arrayInd(up, c(nr, nc))
		lo <- ij[, 2L] + (ij[, 1L] - 1L) * nr
		step <- (seq_len(N) - 1L) * nr * nc
		cx$pu <- rep(up, N) + rep(step, each = length(up))
		cx$pl <- rep(lo, N) + rep(step, each = length(lo))
		cx$pd <- rep(seq.int(1L, by = nr + 1L, length.out = nr), N) + rep(step, each = nr)
		if (!identical(family, "ordinal")) {
			# dyad (and triangle) of every cell with a Poisson factor
			k_u <- match(cx$ok_idx, cx$pu)
			cx$ok_dyad <- ifelse(is.na(k_u), match(cx$ok_idx, cx$pl), k_u)
			cx$ok_up <- !is.na(k_u)
		}
	}
	cx
}

# Poisson sums of a whole phase: E = exp(z + offset) on observed cells, S its sum over
# each actor's line, ys the line's counts
.lame_ltm_pois_phase <- function(cx, Z, row) {
	E <- exp(Z + cx$offarr)
	E[!cx$obs3] <- 0
	if (row || identical(cx$mode, "sym")) {
		list(E = E, S = rowSums(E, dims = 1L), ys = cx$ys_row)
	} else {
		list(E = E, S = rowSums(matrix(colSums(matrix(E, nrow = cx$nr)), cx$nc, cx$N)), ys = cx$ys_col)
	}
}

# interval of shifts delta that keep an ordinal actor's latent cells where the data
# allow them (rank bounds, or cutpoint intervals when alpha is given)
.lame_ltm_ord_bounds <- function(cx, Z, i, row, alpha = NULL) {
	W <- cx$W
	sym <- identical(cx$mode, "sym")
	lo <- -Inf
	hi <- Inf
	af <- if (is.null(alpha)) NULL else c(-Inf, alpha, Inf)
	for (t in seq_len(cx$N)) {
		if (sym || row) {
			zm <- Z[i, , t]
			wm <- W[i, , t]
		} else {
			zm <- Z[, i, t]
			wm <- W[, i, t]
		}
		ok <- !is.na(wm)
		if (!any(ok)) next
		zm <- zm[ok]
		wm <- wm[ok]
		if (!is.null(af)) {
			lo <- max(lo, max(af[wm] - zm))
			hi <- min(hi, min(af[wm + 1L] - zm))
			next
		}
		zt <- Z[, , t]
		wt <- W[, , t]
		if (sym) {
			zo <- zt[-i, -i, drop = FALSE]
			wo <- wt[-i, -i, drop = FALSE]
		} else if (row) {
			zo <- zt[-i, , drop = FALSE]
			wo <- wt[-i, , drop = FALSE]
		} else {
			zo <- zt[, -i, drop = FALSE]
			wo <- wt[, -i, drop = FALSE]
		}
		ko <- !is.na(wo)
		zo <- zo[ko]
		wo <- wo[ko]
		if (!length(zo)) next
		for (k in unique(wm)) {
			mk <- zm[wm == k]
			below <- wo < k
			above <- wo > k
			if (any(below)) lo <- max(lo, max(zo[below]) - min(mk))
			if (any(above)) hi <- min(hi, min(zo[above]) - max(mk))
		}
	}
	c(lo, hi)
}

# cache of each category's extremes per period, so rank bounds need no full scan;
# rebuilt at the start of each phase (rows, then columns)
.lame_ltm_ord_cache <- function(cx, Z) {
	K <- cx$K
	N <- cx$N
	gmax <- gmin <- matrix(NA_real_, K, N)
	lmax <- lmin <- matrix(NA_integer_, K, N)
	for (t in seq_len(N)) {
		zt <- Z[, , t]
		for (k in seq_len(K)) {
			idx <- cx$cells[[t]][[k]]
			if (!length(idx)) next
			zc <- zt[idx]
			a <- which.max(zc)
			b <- which.min(zc)
			gmax[k, t] <- zc[a]
			lmax[k, t] <- idx[a]
			gmin[k, t] <- zc[b]
			lmin[k, t] <- idx[b]
		}
	}
	list(gmax = gmax, gmin = gmin, lmax = lmax, lmin = lmin)
}

# bounds for actor i from the cached extremes; `upd` holds what
# .lame_ltm_ord_update() needs
.lame_ltm_ord_bounds_cached <- function(cx, Z, cache, i, row) {
	nr <- cx$nr
	N <- cx$N
	K <- cx$K
	sym <- identical(cx$mode, "sym")
	mover <- sym || row
	ZM <- matrix(if (mover) Z[i, , ] else Z[, i, ], ncol = N)
	WM <- matrix(if (mover) cx$W[i, , ] else cx$W[, i, ], ncol = N)
	cols_n <- seq_len(N)
	mmax <- mmin <- matrix(NA_real_, K, N)
	lm <- lmn <- matrix(NA_integer_, K, N)
	obs <- !is.na(WM)
	for (k in seq_len(K)) {
		mk <- obs & WM == k
		has <- colSums(mk) > 0
		if (!any(has)) next
		A <- ZM
		A[!mk] <- -Inf
		a <- max.col(t(A), ties.method = "first")
		B <- ZM
		B[!mk] <- Inf
		b <- max.col(t(-B), ties.method = "first")
		mmax[k, has] <- A[cbind(a, cols_n)][has]
		mmin[k, has] <- B[cbind(b, cols_n)][has]
		lm[k, has] <- as.integer(if (mover) i + (a - 1L) * nr else a + (i - 1L) * nr)[has]
		lmn[k, has] <- as.integer(if (mover) i + (b - 1L) * nr else b + (i - 1L) * nr)[has]
	}
	inside <- function(loc) {
		r <- ((loc - 1L) %% nr) + 1L
		cc <- ((loc - 1L) %/% nr) + 1L
		if (sym) r == i | cc == i else if (row) r == i else cc == i
	}
	omax <- cache$gmax
	omin <- cache$gmin
	olmax <- cache$lmax
	olmin <- cache$lmin
	for (q in which(!is.na(omax) & inside(olmax))) {
		k <- ((q - 1L) %% K) + 1L
		tp <- ((q - 1L) %/% K) + 1L
		idx <- cx$cells[[tp]][[k]]
		idx <- idx[!inside(idx)]
		if (length(idx)) {
			zc <- Z[, , tp][idx]
			a <- which.max(zc)
			omax[q] <- zc[a]
			olmax[q] <- idx[a]
		} else {
			omax[q] <- NA_real_
			olmax[q] <- NA_integer_
		}
	}
	for (q in which(!is.na(omin) & inside(olmin))) {
		k <- ((q - 1L) %% K) + 1L
		tp <- ((q - 1L) %/% K) + 1L
		idx <- cx$cells[[tp]][[k]]
		idx <- idx[!inside(idx)]
		if (length(idx)) {
			zc <- Z[, , tp][idx]
			a <- which.min(zc)
			omin[q] <- zc[a]
			olmin[q] <- idx[a]
		} else {
			omin[q] <- NA_real_
			olmin[q] <- NA_integer_
		}
	}
	om <- omax
	om[is.na(om)] <- -Inf
	on <- omin
	on[is.na(on)] <- Inf
	lowmax <- matrix(-Inf, K, N)
	highmin <- matrix(Inf, K, N)
	for (k in seq_len(K)[-1L]) lowmax[k, ] <- pmax(lowmax[k - 1L, ], om[k - 1L, ])
	for (k in rev(seq_len(K))[-1L]) highmin[k, ] <- pmin(highmin[k + 1L, ], on[k + 1L, ])
	lo <- -Inf
	hi <- Inf
	if (any(!is.na(mmin))) {
		lo <- max(lowmax - mmin, na.rm = TRUE)
		hi <- min(highmin - mmax, na.rm = TRUE)
	}
	list(
		bounds = c(lo, hi),
		upd = list(
			omax = omax,
			omin = omin,
			olmax = olmax,
			olmin = olmin,
			mmax = mmax,
			mmin = mmin,
			lm = lm,
			lmn = lmn
		)
	)
}

# cached extremes after the actor's cells moved by delta
.lame_ltm_ord_update <- function(cache, upd, delta) {
	has <- !is.na(upd$mmax)
	vmax <- upd$mmax + delta
	win <- has & (is.na(upd$omax) | vmax > upd$omax)
	cache$gmax[win] <- vmax[win]
	cache$lmax[win] <- upd$lm[win]
	lose <- has & !win
	cache$gmax[lose] <- upd$omax[lose]
	cache$lmax[lose] <- upd$olmax[lose]
	vmin <- upd$mmin + delta
	win <- has & (is.na(upd$omin) | vmin < upd$omin)
	cache$gmin[win] <- vmin[win]
	cache$lmin[win] <- upd$lmn[win]
	lose <- has & !win
	cache$gmin[lose] <- upd$omin[lose]
	cache$lmin[lose] <- upd$olmin[lose]
	cache
}

# one visit of every level (rows, then columns, each in random order; a symmetric
# network's actors once); returns the updated state and the counters
.lame_ltm_visit <- function(
	cx,
	Z,
	EZ,
	la,
	lb,
	a_mat,
	b_mat,
	beta,
	Sab,
	iV0,
	centred = TRUE,
	ic = NULL,
	lev = c(TRUE, TRUE),
	alpha = NULL
) {
	nr <- cx$nr
	nc <- cx$nc
	sym <- identical(cx$mode, "sym")
	bip <- identical(cx$mode, "bip")
	pois <- identical(cx$family, "poisson")
	dyn <- !is.null(a_mat)
	stats <- c(level_calls = 0, level_evals = 0, stuck = 0, skipped = 0)
	state <- function() list(
		Z = Z,
		EZ = EZ,
		la = la,
		lb = lb,
		a_mat = a_mat,
		b_mat = b_mat,
		beta = beta,
		ic = ic,
		stats = stats
	)
	iS <- NULL
	if (sym) {
		va <- Sab[1L, 1L]
		if (!(is.finite(va) && va > 0)) return(state())
	} else if (bip) {
		va <- Sab[1L, 1L]
		vb <- Sab[2L, 2L]
	} else {
		iS <- if (any(lev)) tryCatch(solve(Sab), error = function(e) NULL) else diag(2L)
		if (is.null(iS) || any(!is.finite(iS))) return(state())
	}
	for (row in if (sym) TRUE else c(TRUE, FALSE)) {
		if (!lev[if (row) 1L else 2L]) next
		n_side <- if (row) nr else nc
		cache <- if (!pois && is.null(alpha)) .lame_ltm_ord_cache(cx, Z)
		ph <- if (pois) .lame_ltm_pois_phase(cx, Z, row)
		# latent shifts of a Poisson phase are applied at the end of the phase
		dsh <- numeric(n_side)
		# Poisson margin with static, uncentred levels runs as one compiled phase, drawing
		# random numbers in the same order as the loop below
		cpp_phase <- pois && !centred && !dyn && !identical(getOption("lame.ltm_engine"), "r")
		if (cpp_phase) {
			mode_ <- if (sym) 1L else if (bip) 2L else 0L
			x_ <- if (row) la else lb
			o_ <- if (row) lb else la
			res_ <- ltm_level_phase_cpp(
				as.numeric(ph$S),
				as.numeric(ph$ys),
				if (sym) as.numeric(ph$E) else numeric(0),
				if (sym) as.integer(dim(ph$E)) else integer(0),
				as.numeric(x_),
				as.numeric(o_),
				if (is.null(iS)) matrix(0, 2L, 2L) else iS,
				if (sym || (bip && row)) va else if (bip) vb else 0,
				mode_,
				row,
				as.integer(sample.int(n_side))
			)
			if (row) la[] <- res_$x else lb[] <- res_$x
			dsh <- res_$dsh
			stats[c(
				"level_calls", "level_evals", "stuck", "skipped"
			)] <- stats[c("level_calls", "level_evals", "stuck", "skipped")] + res_$stats
		}
		for (i in if (cpp_phase) integer(0) else sample.int(n_side)) {
			pr <- if (sym) .lame_actor_move_level_prior_sym(i, la, va, iV0, beta, nr, centred, ic)
			else if (bip) .lame_actor_move_level_prior_bip(
				i,
				if (row) la else lb,
				if (row) va else vb,
				iV0,
				beta,
				centred,
				ic
			)
			else .lame_actor_move_level_prior(i, row, la, lb, iS, iV0, beta, nr, centred, ic)
			if (!all(is.finite(pr)) || pr[["a2"]] <= 0) {
				stats["skipped"] <- stats["skipped"] + 1
				next
			}
			if (pois) {
				ys <- ph$ys[i]
				s <- ph$S[i]
				ans <- if (identical(getOption("lame.ltm_engine"), "r")) {
					ll <- function(dl) ys * dl - s * expm1(dl)
					.lame_ess_step(0, -pr[["a1"]] / pr[["a2"]], matrix(1 / sqrt(pr[["a2"]])), ll, 0)
				} else ess_pois_scalar_cpp(ys, s, -pr[["a1"]] / pr[["a2"]], 1 / sqrt(pr[["a2"]]))
				stats["level_calls"] <- stats["level_calls"] + 1
				stats["level_evals"] <- stats["level_evals"] + ans$evals
				stats["stuck"] <- stats["stuck"] + ans$stuck
				d <- ans$x
			} else {
				cb <- if (is.null(alpha)) .lame_ltm_ord_bounds_cached(cx, Z, cache, i, row)
				bd <- if (is.null(alpha)) cb$bounds else .lame_ltm_ord_bounds(cx, Z, i, row, alpha)
				d <- if (bd[1L] < 0 && bd[2L] > 0)
					tryCatch(
						.lame_checked_truncated_gaussian(
							-pr[["a1"]] / pr[["a2"]],
							1 / sqrt(pr[["a2"]]),
							bd[1L],
							bd[2L]
						),
						error = function(e) NULL
					)
				else NULL
				if (is.null(d)) {
					stats["skipped"] <- stats["skipped"] + 1
					next
				}
				stats["level_calls"] <- stats["level_calls"] + 1
			}
			if (!pois && is.null(alpha)) cache <- .lame_ltm_ord_update(cache, cb$upd, d)
			if (d == 0) next
			if (pois) {
				ed <- exp(d)
				if (sym) {
					# pairs {i, j} are shared with actor j, whose sum gains the change of the pair
					Ei <- matrix(ph$E[i, , ], ncol = cx$N)
					ph$S <- ph$S + (ed - 1) * rowSums(Ei)
					ph$E[i, , ] <- Ei * ed
					ph$E[, i, ] <- Ei * ed
				}
				ph$S[i] <- ph$S[i] * ed
			}
			n_lev <- if (row) nr else nc
			w <- numeric(n_lev)
			if (centred) {
				w[] <- -d / n_lev
				w[i] <- d * (1 - 1 / n_lev)
			} else w[i] <- d
			if (sym) {
				# pairs {i, j} gain d; the diagonal cell (i, i) gains 2 d
				if (pois) dsh[i] <- d else {
					Z[i, , ] <- Z[i, , ] + d
					Z[, i, ] <- Z[, i, ] + d
					EZ[i, , ] <- EZ[i, , ] + d
					EZ[, i, ] <- EZ[, i, ] + d
				}
				la <- la + w
				if (dyn) a_mat <- a_mat + w
			} else if (row) {
				if (pois) dsh[i] <- d else {
					Z[i, , ] <- Z[i, , ] + d
					EZ[i, , ] <- EZ[i, , ] + d
				}
				la <- la + w
				if (dyn) a_mat <- a_mat + w
			} else {
				if (pois) dsh[i] <- d else {
					Z[, i, ] <- Z[, i, ] + d
					EZ[, i, ] <- EZ[, i, ] + d
				}
				lb <- lb + w
				if (dyn) b_mat <- b_mat + w
			}
			if (centred) {
				sh <- .lame_actor_move_icpt_shift(ic, beta, (if (sym) 2 else 1) * d / n_lev)
				beta <- sh$beta
				ic <- sh$ic
			}
		}
		if (pois && any(dsh != 0)) {
			# a row's cells gain its shift, a column's cells theirs, a symmetric pair both actors' shifts
			add <- if (sym) dsh + rep(dsh, each = nr) else if (row) dsh else rep(dsh, each = nr)
			Z <- Z + add
			EZ <- EZ + add
		}
	}
	state()
}

# location move of an ordinal latent array: the intercept (or the levels' means)
# and every latent value shift together, drawn exactly from the prior, truncated
# to the strip allowed by explicit cutpoints
.lame_ltm_location <- function(
	cx,
	Z,
	EZ,
	la,
	lb,
	a_mat,
	b_mat,
	beta,
	Sab,
	iV0,
	centred = FALSE,
	ic = NULL,
	lev = c(TRUE, TRUE),
	alpha = NULL
) {
	sym <- identical(cx$mode, "sym")
	bip <- identical(cx$mode, "bip")
	same <- list(Z = Z, EZ = EZ, la = la, lb = lb, a_mat = a_mat, b_mat = b_mat, beta = beta, ic = ic, shift = 0)
	# interval of total shifts that keep every categorised value in its cutpoint interval
	bd <- c(-Inf, Inf)
	if (!is.null(alpha)) {
		ok <- !is.na(cx$W)
		af <- c(-Inf, alpha, Inf)
		bd <- c(max(af[cx$W[ok]] - Z[ok]), min(af[cx$W[ok] + 1L] - Z[ok]))
		if (!(bd[1L] < 0 && bd[2L] > 0)) return(same)
	}
	tn <- function(m, v) {
		if (!is.finite(bd[1L]) && !is.finite(bd[2L])) return(stats::rnorm(1L, m, sqrt(v)))
		tryCatch(.lame_checked_truncated_gaussian(m, sqrt(v), bd[1L], bd[2L]), error = function(e) NULL)
	}
	if (centred) {
		lq <- .lame_actor_move_icpt(ic, iV0, beta)
		if (!all(is.finite(lq)) || lq[["q"]] <= 0) return(same)
		shift <- tn(-lq[["l"]] / lq[["q"]], 1 / lq[["q"]])
		if (is.null(shift)) return(same)
		sh <- .lame_actor_move_icpt_shift(ic, beta, shift)
		return(list(
			Z = Z + shift,
			EZ = EZ + shift,
			la = la,
			lb = lb,
			a_mat = a_mat,
			b_mat = b_mat,
			beta = sh$beta,
			ic = sh$ic,
			shift = shift
		))
	}
	if (!any(lev)) return(same)
	if (sym) {
		va <- Sab[1L, 1L]
		if (!(is.finite(va) && va > 0)) return(same)
		m <- -mean(la)
		S <- matrix(va / length(la))
		w <- 2
		idx <- 1L
	} else if (bip) {
		idx <- which(lev)
		dv <- c(Sab[1L, 1L] / length(la), Sab[2L, 2L] / length(lb))[idx]
		if (any(!is.finite(dv)) || any(dv <= 0)) return(same)
		m <- -c(mean(la), mean(lb))[idx]
		S <- diag(dv, length(idx))
		w <- rep(1, length(idx))
	} else {
		iS <- tryCatch(solve(Sab), error = function(e) NULL)
		if (is.null(iS) || any(!is.finite(iS))) return(same)
		idx <- which(lev)
		n <- length(la)
		P <- n * iS[idx, idx, drop = FALSE]
		S <- tryCatch(solve(P), error = function(e) NULL)
		if (is.null(S) || any(!is.finite(S))) return(same)
		S <- (S + t(S)) / 2
		m <- -drop(S %*% (n * iS[idx, , drop = FALSE] %*% c(mean(la), mean(lb))))
		w <- rep(1, length(idx))
	}
	ch <- tryCatch(chol(S), error = function(e) NULL)
	if (is.null(ch)) return(same)
	cc <- m + drop(crossprod(ch, stats::rnorm(length(m))))
	if (is.finite(bd[1L]) || is.finite(bd[2L])) {
		Sw <- drop(S %*% w)
		vs <- sum(w * Sw)
		tot <- tn(sum(w * m), vs)
		if (is.null(tot)) return(same)
		cc <- cc + Sw * (tot - sum(w * cc)) / vs
	}
	c1 <- c2 <- 0
	if (sym) {
		c1 <- cc[1L]
	} else {
		if (lev[1L]) c1 <- cc[match(1L, idx)]
		if (lev[2L]) c2 <- cc[match(2L, idx)]
	}
	shift <- if (sym) 2 * c1 else c1 + c2
	la <- la + c1
	if (!is.null(a_mat)) a_mat <- a_mat + c1
	if (!sym) {
		lb <- lb + c2
		if (!is.null(b_mat)) b_mat <- b_mat + c2
	}
	list(
		Z = Z + shift,
		EZ = EZ + shift,
		la = la,
		lb = lb,
		a_mat = a_mat,
		b_mat = b_mat,
		beta = beta,
		ic = ic,
		shift = shift
	)
}

# Poisson coefficient translations: coefficient k and the latent array move by
# delta x_k together, drawn by an elliptical slice step
.lame_ltm_coef <- function(cx, Z, EZ, beta, Xk, iV0) {
	stats <- c(coef_calls = 0, coef_evals = 0, coef_stuck = 0)
	for (k in sample.int(length(beta))) {
		x <- Xk[[k]]
		xo <- x[cx$ok_idx]
		if (!all(is.finite(xo)) || !any(xo != 0)) next
		ez <- exp(Z[cx$ok_idx] + cx$off_ok)
		ysx <- sum(cx$y_ok * xo)
		ll <- function(dl) dl * ysx - sum(ez * expm1(dl * xo))
		a1 <- drop(iV0 %*% beta)[k]
		a2 <- iV0[k, k]
		if (!is.finite(a1) || !is.finite(a2) || a2 <= 0) next
		ans <- .lame_ess_step(0, -a1 / a2, matrix(1 / sqrt(a2)), ll, 0)
		stats["coef_calls"] <- stats["coef_calls"] + 1
		stats["coef_evals"] <- stats["coef_evals"] + ans$evals
		stats["coef_stuck"] <- stats["coef_stuck"] + ans$stuck
		d <- ans$x
		if (d == 0) next
		Z <- Z + d * x
		EZ <- EZ + d * x
		beta[k] <- beta[k] + d
	}
	list(Z = Z, EZ = EZ, beta = beta, stats = stats)
}

# joint Metropolis scale move of s2 and rho with the latent residuals (Poisson),
# scaling the dyad contrasts p and q; `s2_prior` is c(shape0, rate0) of the
# inverse-gamma prior, and `diag_in_target` scales the diagonal when the host's s2
# conditional counts it (ame() does)
.lame_ltm_scale <- function(
	cx,
	Z,
	EZ,
	s2,
	rho,
	s2_prior,
	two = FALSE,
	rho_logprior = NULL,
	steps = c(0.03, 0.12, 0.45),
	tries = 6L,
	diag_in_target = FALSE
) {
	mode <- cx$mode
	bip <- identical(mode, "bip")
	sym <- identical(mode, "sym")
	stats <- c(scale_attempts = 0, scale_accepts = 0, rho_attempts = 0, rho_accepts = 0)
	# a symmetric pair carries one latent value around the (symmetric) predictor
	if (sym && max(abs(EZ[cx$pu] - EZ[cx$pl])) > 1e-8) return(list(Z = Z, s2 = s2, rho = rho, stats = stats))
	lpi_s <- function(s) -(s2_prior[["shape0"]] + 1) * log(s) - s2_prior[["rate0"]] / s
	# proposals are judged on the observed cells only; carried residuals are scaled by
	# the accepted factors and the array is rebuilt once at the end
	ok <- cx$ok_idx
	base <- EZ[ok] + cx$off_ok
	e0 <- Z - EZ
	eo <- e0[ok]
	if (two) {
		p0 <- (e0[cx$pu] + e0[cx$pl]) / 2
		q0 <- (e0[cx$pu] - e0[cx$pl]) / 2
		po <- p0[cx$ok_dyad]
		qo <- ifelse(cx$ok_up, 1, -1) * q0[cx$ok_dyad]
	}
	sp <- sq <- 1
	ll_of <- function(zz) sum(cx$y_ok * zz - exp(zz))
	ll0 <- ll_of(Z[ok] + cx$off_ok)
	moved <- FALSE
	for (k in seq_len(tries)) {
		tau <- steps[sample.int(length(steps), 1L)]
		type <- if (two) sample.int(3L, 1L) else 1L
		rho_only <- type > 1L
		if (two) {
			if (!rho_only) {
				# both contrasts scale: log c_plus and log c_minus normal
				lc <- stats::rnorm(2L, 0, tau)
				cp <- exp(lc[1L])
				cm <- exp(lc[2L])
				vp <- s2 * (1 + rho)
				vm <- s2 * (1 - rho)
				s2n <- (cp * vp + cm * vm) / 2
				rhon <- (cp * vp - cm * vm) / (cp * vp + cm * vm)
			} else {
				# correlation alone moves (s2 held); a chain at |rho| >= 0.99 is left alone
				if (type == 3L && abs(rho) >= 0.99) next
				rhon <- if (type == 2L) rho + tau * stats::rnorm(1L) else stats::runif(1L, -0.99, 0.99)
				stats["rho_attempts"] <- stats["rho_attempts"] + 1
				if (abs(rhon) >= 1) next
				s2n <- s2
				cp <- (1 + rhon) / (1 + rho)
				cm <- (1 - rhon) / (1 - rho)
			}
			zn <- base + sp * sqrt(cp) * po + sq * sqrt(cm) * qo
			# Jacobian and prior terms of the scale move; the rho-only move is a rotation
			# that leaves only the prior ratio of rho
			extra <- if (!rho_only) log(cp) + log(cm) + log(s2 / s2n) + lpi_s(s2n) - lpi_s(s2) + rho_logprior(rhon) - rho_logprior(rho)
			else rho_logprior(rhon) - rho_logprior(rho)
		} else {
			cc <- exp(stats::rnorm(1L, 0, tau))
			s2n <- cc * s2
			rhon <- rho
			zn <- base + sp * sqrt(cc) * eo
			extra <- log(cc) + lpi_s(s2n) - lpi_s(s2)
		}
		stats["scale_attempts"] <- stats["scale_attempts"] + 1
		ll1 <- ll_of(zn)
		la <- ll1 - ll0 + extra
		if (is.finite(la) && log(stats::runif(1L)) < la) {
			if (two) {
				sp <- sp * sqrt(cp)
				sq <- sq * sqrt(cm)
			} else sp <- sp * sqrt(cc)
			ll0 <- ll1
			s2 <- s2n
			rho <- rhon
			moved <- TRUE
			stats["scale_accepts"] <- stats["scale_accepts"] + 1
			if (rho_only) stats["rho_accepts"] <- stats["rho_accepts"] + 1
		}
	}
	if (moved) {
		if (bip) {
			Z <- EZ + sp * e0
		} else if (two) {
			Z[cx$pu] <- EZ[cx$pu] + sp * p0 + sq * q0
			Z[cx$pl] <- EZ[cx$pl] + sp * p0 - sq * q0
			# diagonal cell in the target has variance s2 (1 + rho) = v_plus
			if (diag_in_target) Z[cx$pd] <- EZ[cx$pd] + sp * e0[cx$pd]
		} else if (sym) {
			# one latent value per unordered pair, mirrored
			Z[cx$pu] <- EZ[cx$pu] + sp * e0[cx$pu]
			Z[cx$pl] <- Z[cx$pu]
		} else {
			Z[cx$pu] <- EZ[cx$pu] + sp * e0[cx$pu]
			Z[cx$pl] <- EZ[cx$pl] + sp * e0[cx$pl]
			if (diag_in_target) Z[cx$pd] <- EZ[cx$pd] + sp * e0[cx$pd]
		}
	}
	list(Z = Z, s2 = s2, rho = rho, stats = stats)
}

# static positions of the Poisson family: u_i (or v_j) moves with the latent values
# of its line, drawn by an elliptical slice step from the host's factor prior times
# the line's Poisson factor; `shrink` says whether the hierarchical factor prior is
# in force, `Suv` and `tau` fix the prior (tests)
.lame_ltm_uv_visit <- function(
	cx,
	Z,
	EZ,
	U,
	V = NULL,
	G = NULL,
	Ld = NULL,
	shrink = TRUE,
	uv_prior = NULL,
	Suv = NULL,
	tau = NULL
) {
	mode <- cx$mode
	sym <- identical(mode, "sym")
	bip <- identical(mode, "bip")
	nr <- cx$nr
	nc <- cx$nc
	R <- ncol(U)
	stats <- c(uv_calls = 0, uv_evals = 0, uv_stuck = 0)
	state <- function() list(Z = Z, EZ = EZ, U = U, V = V, stats = stats)
	if (R == 0L) return(state())
	pri <- NULL
	if (sym) {
		if (!(all(is.finite(Ld)) && isTRUE(shrink))) return(state())
		if (is.null(tau)) tau <- .lame_actor_move_tau(U, nr, R, shrink)
	} else if (!bip) {
		# early burn-in uses a diffuse prior, so the positions are left to the host
		if (!isTRUE(shrink) && is.null(uv_prior) && is.null(Suv)) return(state())
		pri <- tryCatch({
			S <- if (is.null(Suv)) .lame_actor_move_suv(U, V, nr, R, shrink, uv_prior) else Suv
			P <- solve(S)
			iu <- seq_len(R)
			iv <- R + seq_len(R)
			cov_u <- solve(P[iu, iu, drop = FALSE])
			cov_v <- solve(P[iv, iv, drop = FALSE])
			cov_u <- (cov_u + t(cov_u)) / 2
			cov_v <- (cov_v + t(cov_v)) / 2
			list(
				root_u = t(chol(cov_u)),
				root_v = t(chol(cov_v)),
				reg_u = -cov_u %*% P[iu, iv, drop = FALSE],
				reg_v = -cov_v %*% P[iv, iu, drop = FALSE]
			)
		}, error = function(e) NULL)
		if (is.null(pri) || any(!is.finite(unlist(pri)))) return(state())
	}
	# Poisson factor of each line: counts and exp(z + offset) summed over periods,
	# formed once per phase
	sides <- if (sym) TRUE else c(TRUE, FALSE)
	for (row in sides) {
		sj <- rowSums(.lame_ltm_pois_phase(cx, Z, TRUE)$E, dims = 2L)
		# latent shifts are applied at the end of the phase
		Dsh <- matrix(0, nr, nc)
		# compiled phase of a margin, drawing random numbers in the same order as the loop below
		if (!identical(getOption("lame.ltm_engine"), "r")) {
			Rr <- ncol(U)
			pm <- if (sym) list(
				mov = U,
				oth = U,
				opp = U,
				reg = matrix(0, 0L, 0L),
				root = diag(1 / sqrt(tau), Rr),
				w = Ld,
				mode = 1L
			)
			else if (bip && row) list(
				mov = U,
				oth = U,
				opp = tcrossprod(V, G),
				reg = matrix(0, 0L, 0L),
				root = diag(Rr),
				w = rep(1, Rr),
				mode = 2L
			)
			else if (bip) list(
				mov = V,
				oth = V,
				opp = U %*% G,
				reg = matrix(0, 0L, 0L),
				root = diag(ncol(V)),
				w = rep(1, ncol(V)),
				mode = 2L
			)
			else if (row) list(mov = U, oth = V, opp = V, reg = pri$reg_u, root = pri$root_u, w = rep(1, Rr), mode = 0L)
			else list(mov = V, oth = U, opp = U, reg = pri$reg_v, root = pri$root_v, w = rep(1, ncol(V)), mode = 0L)
			res_ <- ltm_uv_phase_cpp(
				pm$mov,
				pm$oth,
				pm$opp,
				pm$reg,
				pm$root,
				as.numeric(pm$w),
				sj,
				cx$ysj,
				pm$mode,
				row,
				as.integer(sample.int(if (row) nr else nc))
			)
			if (row) U[] <- res_$M else V[] <- res_$M
			sj <- res_$sj
			Dsh <- res_$Dsh
			stats[c("uv_calls", "uv_evals", "uv_stuck")] <- stats[c("uv_calls", "uv_evals", "uv_stuck")] + res_$stats
			Z <- Z + as.vector(Dsh)
			EZ <- EZ + as.vector(Dsh)
			next
		}
		for (i in sample.int(if (row) nr else nc)) {
			if (row && !sym && !bip) {
				cur <- U[i, ]
				opp <- V
				mu <- drop(pri$reg_u %*% V[i, ])
				root <- pri$root_u
			} else if (!row && !bip) {
				cur <- V[i, ]
				opp <- U
				mu <- drop(pri$reg_v %*% U[i, ])
				root <- pri$root_v
			} else if (sym) {
				cur <- U[i, ]
				opp <- U
				mu <- numeric(R)
				root <- diag(1 / sqrt(tau), R)
			} else if (row) {
				cur <- U[i, ]
				opp <- tcrossprod(V, G)
				mu <- numeric(R)
				root <- diag(R)
			} else {
				cur <- V[i, ]
				opp <- U %*% G
				mu <- numeric(ncol(V))
				root <- diag(ncol(V))
			}
			w <- if (sym) Ld else rep(1, length(cur))
			lsum <- if (row || sym) list(ys = cx$ysj[i, ], s = sj[i, ]) else list(ys = cx$ysj[, i], s = sj[, i])
			ans <- if (identical(getOption("lame.ltm_engine"), "r")) {
				if (!any(lsum$ys > 0) && !any(lsum$s > 0)) {
					ll <- function(x) 0
				} else {
					ll <- function(x) {
						d <- drop(opp %*% (w * (x - cur)))
						if (sym) d[i] <- 0
						sum(lsum$ys * d - lsum$s * expm1(d))
					}
				}
				.lame_ess_step(cur, mu, root, ll, 0)
			} else ess_pois_uv_cpp(cur, mu, root, opp, w, lsum$ys, lsum$s, if (sym) i else 0L)
			stats["uv_calls"] <- stats["uv_calls"] + 1
			stats["uv_evals"] <- stats["uv_evals"] + ans$evals
			stats["uv_stuck"] <- stats["uv_stuck"] + ans$stuck
			dx <- ans$x - cur
			if (!any(dx != 0)) next
			shift <- drop(opp %*% (w * dx))
			if (sym) {
				dg <- sum(Ld * (ans$x^2 - cur^2))
				shift[i] <- 0
				sj[i, ] <- sj[i, ] * exp(shift)
				sj[, i] <- sj[, i] * exp(shift)
				Dsh[i, ] <- Dsh[i, ] + shift
				Dsh[, i] <- Dsh[, i] + shift
				Dsh[i, i] <- Dsh[i, i] + dg
				U[i, ] <- ans$x
			} else if (row) {
				sj[i, ] <- sj[i, ] * exp(shift)
				Dsh[i, ] <- shift
				U[i, ] <- ans$x
			} else {
				sj[, i] <- sj[, i] * exp(shift)
				Dsh[, i] <- shift
				V[i, ] <- ans$x
			}
		}
		Z <- Z + as.vector(Dsh)
		EZ <- EZ + as.vector(Dsh)
	}
	state()
}
