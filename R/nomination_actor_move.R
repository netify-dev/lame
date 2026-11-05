# receiver-level move of the cbin and frn nomination samplers: b_j and its latent
# column Z[, j] shift by a common delta drawn exactly from its truncated normal
# conditional. nom_col_shift_cpp is the compiled sweep, .lame_nom_col_shift_r the
# R reference used in the tests

# row maxima of a matrix (-Inf for a row without a finite entry), ties by position
.lame_nom_rowmax <- function(M) {
	if (!nrow(M)) return(numeric(0))
	if (!ncol(M)) return(rep(-Inf, nrow(M)))
	M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
}

# cell classes of one slice and what the bounds read besides the latent values
#   cls: -1 free, 0 non-nominee, cbin 1 nominee; frn_c the rank y >= 1;
#        frn_h the order position q >= 1 of an observed nominee, -2 a hidden nominee
#   kind "cbin", "frn_c" (frn, no hidden cells) or "frn_h" (frn with hidden cells)
#   room / full: whether the row is below its cap (cbin, frn_c) or at it (frn_h)
#   YL / P: the column holding each rank (frn_c) or each order position (frn_h)
.lame_nom_col_prep <- function(Zt, Yt, family, YLt, cap, bip) {
	n_r <- nrow(Zt)
	n_c <- ncol(Zt)
	hid <- .lame_nom_hidden_cells(Yt, bip)
	any_hid <- any(hid)
	pr <- list()
	if (identical(family, "cbin")) {
		Ya <- if (any_hid) .lame_cbin_implied_labels(Zt, Yt, cap, hid) else Yt
		cls <- matrix(-1L, n_r, n_c)
		cls[!is.na(Ya) & Ya > 0] <- 1L
		cls[!is.na(Ya) & Ya == 0] <- 0L
		if (!bip) diag(cls) <- -1L
		pr$kind <- "cbin"
		pr$room <- rowSums(cls == 1L) < cap
		isnom <- cls == 1L
	} else if (!any_hid) {
		Yr <- Yt
		Yr[is.na(Yr)] <- -1
		cls <- matrix(as.integer(Yr), n_r, n_c)
		cls[cls < 0L] <- -1L
		if (!bip) diag(cls) <- -1L
		pr$kind <- "frn_c"
		pr$YL <- YLt
		pr$max_rank <- ncol(YLt)
		pr$room <- rowSums(!is.na(Yt) & Yt > 0) < cap
		isnom <- cls >= 1L
	} else {
		Ya <- .lame_cbin_implied_labels(Zt, Yt, cap, hid)
		obs <- !is.na(Yt)
		if (!bip) diag(obs) <- FALSE
		o1 <- obs & Yt > 0
		h1 <- hid & !is.na(Ya) & Ya > 0
		non <- (obs & Yt == 0) | (hid & !is.na(Ya) & Ya == 0)
		k_row <- rowSums(o1)
		k_max <- max(c(0L, k_row))
		P <- matrix(NA_integer_, n_r, max(k_max, 1L))
		cls <- matrix(-1L, n_r, n_c)
		cls[non] <- 0L
		cls[h1] <- -2L
		if (k_max > 0L) {
			idx <- which(o1)
			rr <- row(o1)[idx]
			cc <- col(o1)[idx]
			o <- order(rr, Yt[idx])
			rr <- rr[o]
			cc <- cc[o]
			qpos <- sequence(tabulate(rr, n_r))
			P[cbind(rr, qpos)] <- cc
			cls[cbind(rr, cc)] <- qpos
		}
		pr$kind <- "frn_h"
		pr$P <- P
		pr$k_max <- k_max
		pr$full <- k_row + rowSums(h1) >= cap
		isnom <- cls >= 1L | cls == -2L
	}
	pr$cls <- cls
	pr$isnom <- isnom
	pr
}

# state of one slice for the R reference sweep: cell classes and row aggregates
.lame_nom_col_slice <- function(Zt, Yt, family, YLt, cap, bip) {
	pr <- .lame_nom_col_prep(Zt, Yt, family, YLt, cap, bip)
	st <- new.env(parent = emptyenv())
	for (nm in names(pr)) assign(nm, pr[[nm]], envir = st)
	st$Z <- Zt
	cls <- pr$cls
	# row maximum over non-nominees and negated row minimum over nominees
	st$Zm0 <- Zt
	st$Zm0[cls != 0L] <- -Inf
	st$Zmn <- -Zt
	st$Zmn[!pr$isnom] <- -Inf
	st$M0 <- .lame_nom_rowmax(st$Zm0)
	st$Nq <- .lame_nom_rowmax(st$Zmn)
	st
}

# interval (lo, hi) each cell of column j may take given the other cells of its row
.lame_nom_col_limits <- function(st, j) {
	n_r <- nrow(st$Z)
	cl <- st$cls[, j]
	lo <- rep(-Inf, n_r)
	hi <- rep(Inf, n_r)
	Zt <- st$Z
	i0 <- which(cl == 0L)
	if (st$kind == "cbin") {
		i1 <- which(cl == 1L)
		if (length(i1)) lo[i1] <- pmax(0, st$M0[i1])
		if (length(i0)) hi[i0] <- ifelse(st$room[i0], 0, -st$Nq[i0])
	} else if (st$kind == "frn_c") {
		YL <- st$YL
		mr <- st$max_rank
		if (length(i0)) {
			h <- Zt[cbind(i0, YL[i0, 1L])]
			h[is.na(h)] <- Inf
			h[st$room[i0]] <- 0
			hi[i0] <- h
		}
		iy <- which(cl >= 1L)
		if (length(iy)) {
			y <- cl[iy]
			l <- rep(-Inf, length(iy))
			u <- rep(Inf, length(iy))
			f1 <- y == 1L
			if (any(f1)) l[f1] <- pmax(0, st$M0[iy[f1]])
			f2 <- !f1
			if (any(f2)) {
				v <- Zt[cbind(iy[f2], YL[cbind(iy[f2], y[f2] - 1L)])]
				v[is.na(v)] <- -Inf
				l[f2] <- v
			}
			fu <- y < mr
			if (any(fu)) {
				v <- Zt[cbind(iy[fu], YL[cbind(iy[fu], y[fu] + 1L)])]
				v[is.na(v)] <- Inf
				u[fu] <- v
			}
			lo[iy] <- l
			hi[iy] <- u
		}
	} else {
		P <- st$P
		if (length(i0)) {
			hi[i0] <- ifelse(st$full[i0], -st$Nq[i0], 0)
		}
		ih <- which(cl == -2L)
		if (length(ih)) lo[ih] <- pmax(0, st$M0[ih])
		iy <- which(cl >= 1L)
		if (length(iy)) {
			q <- cl[iy]
			l <- rep(-Inf, length(iy))
			u <- rep(Inf, length(iy))
			f1 <- q == 1L
			if (any(f1)) l[f1] <- pmax(0, st$M0[iy[f1]])
			f2 <- !f1
			if (any(f2)) {
				v <- Zt[cbind(iy[f2], P[cbind(iy[f2], q[f2] - 1L)])]
				v[is.na(v)] <- -Inf
				l[f2] <- v
			}
			fu <- q < st$k_max
			if (any(fu)) {
				v <- Zt[cbind(iy[fu], P[cbind(iy[fu], q[fu] + 1L)])]
				v[is.na(v)] <- Inf
				u[fu] <- v
			}
			lo[iy] <- l
			hi[iy] <- u
		}
	}
	list(lo = lo, hi = hi)
}

# translate column j of a slice by delta and keep the row aggregates current
.lame_nom_col_apply <- function(st, j, delta) {
	old <- st$Z[, j]
	new <- old + delta
	st$Z[, j] <- new
	cl <- st$cls[, j]
	r0 <- which(cl == 0L)
	if (length(r0)) {
		st$Zm0[r0, j] <- new[r0]
		nz <- new[r0]
		m <- st$M0[r0]
		inc <- nz >= m
		st$M0[r0[inc]] <- nz[inc]
		rc <- r0[!inc & old[r0] >= m]
		if (length(rc)) st$M0[rc] <- .lame_nom_rowmax(st$Zm0[rc, , drop = FALSE])
	}
	rn <- which(st$isnom[, j])
	if (length(rn)) {
		st$Zmn[rn, j] <- -new[rn]
		nq <- -new[rn]
		m <- st$Nq[rn]
		inc <- nq >= m
		st$Nq[rn[inc]] <- nq[inc]
		rc <- rn[!inc & -old[rn] >= m]
		if (length(rc)) st$Nq[rc] <- .lame_nom_rowmax(st$Zmn[rc, , drop = FALSE])
	}
	invisible(NULL)
}

# one draw of N(ez, 1) on (lb, ub) from a supplied uniform, so the R and
# compiled sweeps can consume identical random numbers
.lame_nom_tn_draw <- function(ez, lb, ub, u) {
	u <- min(max(u, .Machine$double.xmin), 1 - 1e-16)
	if (!is.na(ub - ez) && (ub - ez) <= 0) {
		la <- stats::pnorm(lb - ez, log.p = TRUE)
		lc <- stats::pnorm(ub - ez, log.p = TRUE)
		ez + stats::qnorm(lc + log(u + (1 - u) * exp(la - lc)), log.p = TRUE)
	} else {
		lqa <- stats::pnorm(lb - ez, lower.tail = FALSE, log.p = TRUE)
		lqb <- stats::pnorm(ub - ez, lower.tail = FALSE, log.p = TRUE)
		ez + stats::qnorm(lqa + log((1 - u) + u * exp(lqb - lqa)), lower.tail = FALSE, log.p = TRUE)
	}
}

# conditional prior of the receiver levels given the sender levels: mean m and
# variance s2 (no covariance for a bipartite panel)
.lame_nom_col_prior <- function(a, Sab, bip, n_c) {
	if (bip) {
		list(m = rep(0, n_c), s2 = max(Sab[2L, 2L], 1e-12))
	} else {
		r_ba <- if (is.finite(Sab[1L, 1L]) && Sab[1L, 1L] > 0) Sab[1L, 2L] / Sab[1L, 1L] else 0
		list(m = rep_len(as.numeric(a) * r_ba, n_c), s2 = max(Sab[2L, 2L] - r_ba * Sab[1L, 2L], 1e-12))
	}
}

# R reference sweep over the columns; ord the visiting order, u one uniform per
# visited column
.lame_nom_col_shift_r <- function(Z, Y, family, YL, cap, b, bip, center, beta, iV0, pri, ord, u) {
	n_c <- dim(Z)[2L]
	N <- dim(Z)[3L]
	m <- pri$m
	s2 <- pri$s2
	sts <- lapply(
		seq_len(N),
		function(t) .lame_nom_col_slice(
			Z[, , t],
			Y[, , t],
			family,
			if (identical(family, "frn")) YL[[t]] else NULL,
			cap,
			bip
		)
	)
	calls <- 0
	skips <- 0
	for (k in seq_along(ord)) {
		j <- ord[k]
		dlo <- -Inf
		dhi <- Inf
		for (st in sts) {
			lim <- .lame_nom_col_limits(st, j)
			zj <- st$Z[, j]
			fl <- is.finite(lim$lo)
			if (any(fl)) dlo <- max(dlo, lim$lo[fl] - zj[fl])
			fh <- is.finite(lim$hi)
			if (any(fh)) dhi <- min(dhi, lim$hi[fh] - zj[fh])
		}
		if (dlo > 0 || dhi < 0) {
			skips <- skips + 1
			next
		}
		if (center) {
			a2 <- (n_c - 1) / n_c / s2 + iV0[1L, 1L] / n_c^2
			a1 <- ((b[j] - m[j]) - mean(b - m)) / s2 + drop(iV0 %*% beta)[1L] / n_c
			sd_d <- 1 / sqrt(a2)
			delta <- sd_d * .lame_nom_tn_draw((-a1 / a2) / sd_d, dlo / sd_d, dhi / sd_d, u[k])
		} else {
			sd_d <- sqrt(s2)
			val <- sd_d * .lame_nom_tn_draw(m[j] / sd_d, (b[j] + dlo) / sd_d, (b[j] + dhi) / sd_d, u[k])
			delta <- val - b[j]
		}
		if (!is.finite(delta)) {
			skips <- skips + 1
			next
		}
		calls <- calls + 1
		for (st in sts) .lame_nom_col_apply(st, j, delta)
		if (center) {
			b <- b - delta / n_c
			b[j] <- b[j] + delta
			beta[1L] <- beta[1L] + delta / n_c
		} else {
			b[j] <- b[j] + delta
		}
	}
	for (t in seq_len(N)) Z[, , t] <- sts[[t]]$Z
	list(Z = Z, b = b, beta = beta, calls = calls, skips = skips)
}

# one sweep of receiver moves over the columns of a cbin / frn panel
#   YL      frn: per-period matrices of the column holding each rank
#   odmax   the nomination cap (a scalar or one per row)
#   cols_ok the columns that may move (a function of the data alone)
#   center  levels centred over the columns with an intercept (lame()): the
#           other levels and beta[1] take the compensation
#   ic      intercept prior of a dynamic_beta fit (.lame_actor_move_icpt)
#   engine  "cpp" (nom_col_shift_cpp) or "r" (the reference)
# random numbers are drawn here, the visiting order first and then one uniform
# per visited column, whichever engine runs the sweep
.lame_nom_col_shift <- function(
	Z,
	Y,
	family,
	YL,
	odmax,
	a,
	b,
	Sab,
	bip,
	cols_ok,
	center = FALSE,
	beta = NULL,
	iV0 = NULL,
	engine = "cpp",
	ic = NULL
) {
	n_r <- dim(Z)[1L]
	n_c <- dim(Z)[2L]
	N <- dim(Z)[3L]
	cap <- rep_len(as.numeric(odmax), n_r)
	pri <- .lame_nom_col_prior(a, Sab, bip, n_c)
	cols <- which(cols_ok)
	ord <- if (length(cols) > 1L) cols[sample.int(length(cols))] else cols
	u <- stats::runif(length(ord))
	# the intercept prior enters only through its curvature q and linear term l, so a
	# one-coefficient stand-in with value l / q and precision q carries it through
	# either engine; the total intercept shift is read back from it
	if (isTRUE(center) && !is.null(ic)) {
		lq <- .lame_actor_move_icpt(ic, iV0, beta)
		stand <- unname(lq[["l"]] / lq[["q"]])
		if (!is.finite(stand) || !is.finite(lq[["q"]]) || lq[["q"]] <= 0)
			return(list(Z = Z, b = b, beta = beta, ic = ic, calls = 0, skips = length(ord)))
		out <- .lame_nom_col_shift_engine(
			Z,
			Y,
			family,
			YL,
			cap,
			b,
			bip,
			TRUE,
			stand,
			matrix(unname(lq[["q"]]), 1L, 1L),
			pri,
			ord,
			u,
			engine,
			n_r,
			n_c,
			N
		)
		sh <- .lame_actor_move_icpt_shift(ic, beta, out$beta[1L] - stand)
		return(list(Z = out$Z, b = out$b, beta = sh$beta, ic = sh$ic, calls = out$calls, skips = out$skips))
	}
	out <- .lame_nom_col_shift_engine(
		Z,
		Y,
		family,
		YL,
		cap,
		b,
		bip,
		center,
		beta,
		iV0,
		pri,
		ord,
		u,
		engine,
		n_r,
		n_c,
		N
	)
	out$ic <- ic
	out
}
# the sweep itself, in the chosen engine
.lame_nom_col_shift_engine <- function(
	Z,
	Y,
	family,
	YL,
	cap,
	b,
	bip,
	center,
	beta,
	iV0,
	pri,
	ord,
	u,
	engine,
	n_r,
	n_c,
	N
) {
	if (identical(engine, "r"))
		return(.lame_nom_col_shift_r(Z, Y, family, YL, cap, b, bip, center, beta, iV0, pri, ord, u))
	prs <- lapply(
		seq_len(N),
		function(t) .lame_nom_col_prep(
			Z[, , t],
			Y[, , t],
			family,
			if (identical(family, "frn")) YL[[t]] else NULL,
			cap,
			bip
		)
	)
	kind <- match(vapply(prs, function(p) p$kind, ""), c("cbin", "frn_c", "frn_h"))
	flag <- matrix(unlist(lapply(prs, function(p) if (p$kind == "frn_h") p$full else p$room)), n_r, N)
	mr <- vapply(prs, function(p) as.integer(switch(p$kind, frn_c = p$max_rank, frn_h = p$k_max, 0L)), 1L)
	mw <- max(1L, vapply(prs, function(p) as.integer(switch(p$kind, frn_c = ncol(p$YL), frn_h = ncol(p$P), 1L)), 1L))
	map <- array(NA_integer_, c(n_r, mw, N))
	for (t in seq_len(N)) {
		mp <- switch(prs[[t]]$kind, frn_c = prs[[t]]$YL, frn_h = prs[[t]]$P, NULL)
		if (!is.null(mp) && length(mp)) map[, seq_len(ncol(mp)), t] <- as.integer(mp)
	}
	cls <- array(unlist(lapply(prs, function(p) p$cls)), c(n_r, n_c, N))
	storage.mode(Z) <- "double"
	out <- nom_col_shift_cpp(
		Z,
		cls,
		flag,
		map,
		as.integer(kind),
		as.integer(mr),
		as.integer(ord),
		u,
		as.numeric(b),
		as.numeric(pri$m),
		pri$s2,
		isTRUE(center),
		if (center) as.numeric(beta) else numeric(0),
		if (center) as.numeric(iV0[1L, ]) else numeric(0)
	)
	list(Z = out$Z, b = out$b, beta = if (center) out$beta else beta, calls = out$calls, skips = out$skips)
}

# level move of a centred cbin / frn panel (lame()): the free rows C shift by
# delta with their latent rows, the pinned rows R (rows with room) stay, and the
# centred levels and intercept are compensated
#     rows of C   a_i + delta n_R / n,   latent row + delta
#     rows of R   a_i - delta n_C / n,   latent row unchanged
#     intercept   + delta n_C / n
# residuals are unchanged, so delta is an exact truncated normal draw; the row
# effects a_i are required to tell C from R
.lame_nom_level_shift <- function(
	Z,
	a,
	b,
	Sab,
	odmax,
	masks,
	bip,
	rows_ok,
	beta,
	iV0,
	ic = NULL,
	center = TRUE
) {
	n_r <- dim(Z)[1L]
	none <- list(Z = Z, a = a, beta = beta, ic = ic, calls = 0, free = 0L)
	bd <- .lame_nom_shift_bounds(Z, odmax, masks)
	free <- rows_ok & bd$lo < 0 & is.infinite(bd$hi)
	n_c <- sum(free)
	if (n_c == 0L) return(none)
	n_p <- n_r - n_c
	if (center) {
		w <- ifelse(free, n_p / n_r, -n_c / n_r)
		w_b <- n_c / n_r
	} else {
		# levels free of centring (ame()): only the rows with room are compensated and
		# the intercept takes the whole shift
		w <- ifelse(free, 0, -1)
		w_b <- 1
	}
	if (bip) {
		m <- rep(0, n_r)
		s2 <- max(Sab[1L, 1L], 1e-12)
	} else {
		r_ab <- if (is.finite(Sab[2L, 2L]) && Sab[2L, 2L] > 0) Sab[1L, 2L] / Sab[2L, 2L] else 0
		m <- b * r_ab
		s2 <- max(Sab[1L, 1L] - r_ab * Sab[1L, 2L], 1e-12)
	}
	lq <- .lame_actor_move_icpt(ic, iV0, beta)
	Q <- sum(w^2) / s2 + lq[["q"]] * w_b^2
	lin <- sum(w * (a - m)) / s2 + lq[["l"]] * w_b
	if (!is.finite(Q) || Q <= 0 || !is.finite(lin)) return(none)
	sd_d <- 1 / sqrt(Q)
	lower <- max(bd$lo[free])
	delta <- sd_d * rtnorm_interval_logp((-lin / Q) / sd_d, lower / sd_d, Inf)
	if (!is.finite(delta)) return(none)
	Z[free, , ] <- Z[free, , ] + delta
	a <- a + w * delta
	sh <- .lame_actor_move_icpt_shift(ic, beta, w_b * delta)
	list(Z = Z, a = a, beta = sh$beta, ic = sh$ic, calls = 1, free = n_c)
}

# prior state of the intercept of a dynamic_beta fit, in the form
# .lame_actor_move_icpt reads: its path (first state N(beta0_mean, beta0_var),
# AR(1) transitions) when it is dynamic, otherwise its static value with prior
# N(0, g s2). `k` is the column of the path or the position among the static
# coefficients that the move writes back to.
.lame_nom_icpt_state <- function(
	beta_dyn,
	beta_dyn_path,
	beta_static,
	rho_by_group,
	sigma_by_group,
	Lambda_full,
	group_id,
	beta0_mean,
	beta0_cov,
	g,
	s2
) {
	kd <- match(1L, beta_dyn$dynamic_idx)
	if (!is.na(kd)) {
		return(list(
			kind = "path",
			k = kd,
			path = beta_dyn_path[, kd],
			rho = rho_by_group[group_id[kd]],
			sd2 = sigma_by_group[group_id[kd]]^2 * Lambda_full[kd, kd],
			m0 = beta0_mean[kd],
			v0 = beta0_cov[kd, kd]
		))
	}
	ks <- match(1L, beta_dyn$static_idx)
	list(
		kind = "static",
		k = ks,
		value = beta_static[ks],
		prec = 1 / (.lame_static_g_scale(g, beta_dyn$static_idx)[ks] * s2)
	)
}

# intercept-means ridge of a bipartite fit with uncentred levels (ame()):
# beta_0 - c1 - c2, a + c1, b + c2 leaves the likelihood unchanged, so (c1, c2)
# is an exact Gaussian draw from the conditional prior of the shift. a level not
# in the model is left out.
.lame_nom_ridge_bip <- function(beta, a, b, iV0, va, vb, rvar = TRUE, cvar = TRUE) {
	none <- list(beta = beta, a = a, b = b, calls = 0)
	act <- c(isTRUE(rvar) && is.finite(va) && va > 0, isTRUE(cvar) && is.finite(vb) && vb > 0)
	if (!any(act) || !length(beta)) return(none)
	m1 <- drop(iV0 %*% beta)[1L]
	v0 <- iV0[1L, 1L]
	P <- matrix(v0, 2L, 2L) + diag(c(length(a) / va, length(b) / vb), 2L)
	h <- c(m1 - sum(a) / va, m1 - sum(b) / vb)
	P <- P[act, act, drop = FALSE]
	h <- h[act]
	ch <- tryCatch(chol(P), error = function(e) NULL)
	if (is.null(ch) || !all(is.finite(h))) return(none)
	cc <- solve(P, h) + backsolve(ch, stats::rnorm(length(h)))
	if (act[1L]) a <- a + cc[1L]
	if (act[2L]) b <- b + cc[sum(act)]
	beta[1L] <- beta[1L] - sum(cc)
	list(beta = beta, a = a, b = b, calls = 1)
}

# collapsed reciprocity update of the cbin / frn samplers (directed, static rho):
# non-nominations and free cells are integrated out, r takes Metropolis steps
# against the dyad terms (bivariate normal density, Phi, Phi2), and when r moves
# the integrated cells are redrawn given the nominations. `pw` is the exponent of
# (1 - r^2) in the rho prior (.lame_rho_prior_power)

# ceiling of the non-nominations of each row of a slice (see .lame_nom_col_limits)
.lame_nom_row_ceiling <- function(st) {
	weakest <- -st$Nq
	if (identical(st$kind, "frn_h")) ifelse(st$full, weakest, 0) else ifelse(st$room, 0, weakest)
}
# dyads of every slice, by kind, with the quantities the terms read
.lame_nom_rho_dyads <- function(Z, EZ, Y, family, YL, odmax) {
	n <- dim(Z)[1L]
	N <- dim(Z)[3L]
	cap <- rep_len(as.numeric(odmax), n)
	up <- which(upper.tri(matrix(0, n, n)))
	ij <- arrayInd(up, c(n, n))
	lo <- ij[, 2L] + (ij[, 1L] - 1L) * n
	out <- list(
		slices = vector("list", N),
		bb = matrix(0, 0, 2),
		bf_b = numeric(0),
		bf_c = numeric(0),
		ff = matrix(0, 0, 2)
	)
	for (t in seq_len(N)) {
		st <- .lame_nom_col_slice(
			Z[, , t],
			Y[, , t],
			family,
			if (identical(family, "frn")) YL[[t]] else NULL,
			cap,
			FALSE
		)
		eta <- EZ[, , t]
		e <- Z[, , t] - eta
		slack <- matrix(.lame_nom_row_ceiling(st), n, n) - eta
		cl <- st$cls
		nm <- st$isnom
		n1 <- nm[up]; n2 <- nm[lo]
		f1 <- cl[up] == 0L; f2 <- cl[lo] == 0L
		g1 <- cl[up] == -1L; g2 <- cl[lo] == -1L
		bb <- n1 & n2
		b1f2 <- n1 & f2
		f1b2 <- f1 & n2
		ff <- f1 & f2
		out$bb <- rbind(out$bb, cbind(e[up][bb], e[lo][bb]))
		out$bf_b <- c(out$bf_b, e[up][b1f2], e[lo][f1b2])
		out$bf_c <- c(out$bf_c, slack[lo][b1f2], slack[up][f1b2])
		out$ff <- rbind(out$ff, cbind(slack[up][ff], slack[lo][ff]))
		out$slices[[t]] <- list(
			up = up,
			lo = lo,
			e = e,
			eta = eta,
			slack = slack,
			b1f2 = b1f2,
			f1b2 = f1b2,
			ff = ff,
			gg = g1 & g2,
			g1x = g1 & !g2,
			g2x = g2 & !g1,
			dg = (1:n) + (0:(n - 1L)) * n
		)
	}
	out
}
# the log density of r (up to a constant)
.lame_nom_rho_logf <- function(r, d, pw) {
	ld <- log1p(r) + log1p(-r)
	v <- pw * ld
	if (nrow(d$bb)) v <- v - 0.5 * nrow(d$bb) * ld - sum(d$bb^2) / (2 * (1 - r^2)) + r * sum(d$bb[, 1L] * d$bb[, 2L]) / (1 - r^2)
	if (length(d$bf_b)) v <- v + sum(stats::pnorm((d$bf_c - r * d$bf_b) / sqrt(1 - r^2), log.p = TRUE))
	if (nrow(d$ff)) {
		p <- pbivnorm::pbivnorm(d$ff[, 1L], d$ff[, 2L], r)
		v <- v + sum(log(pmax(p, 1e-300)))
	}
	v
}
# redraw every integrated-out cell from its conditional given the nominations, at r
.lame_nom_rho_redraw <- function(Z, d, r) {
	sr <- sqrt(1 - r^2)
	for (t in seq_along(d$slices)) {
		s <- d$slices[[t]]
		z <- Z[, , t]
		e <- s$e
		eta <- s$eta
		# a non-nomination of a nomination's dyad
		for (k in 1:2) {
			sel <- if (k == 1L) s$b1f2 else s$f1b2
			if (!any(sel)) next
			fc <- if (k == 1L) s$lo[sel] else s$up[sel]
			bc <- if (k == 1L) s$up[sel] else s$lo[sel]
			w <- rtnorm_interval_logp(rep(0, length(fc)), -Inf, (s$slack[fc] - r * e[bc]) / sr)
			ef <- r * e[bc] + sr * w
			e[fc] <- ef
			z[fc] <- eta[fc] + ef
		}
		# two non-nominations: by rejection from the sequential proposal
		if (any(s$ff)) {
			up <- s$up[s$ff]; lo <- s$lo[s$ff]
			pend <- seq_along(up)
			e1 <- numeric(length(up)); e2 <- numeric(length(up))
			for (pass in 1:5000) {
				if (!length(pend)) break
				x1 <- rtnorm_interval_logp(rep(0, length(pend)), -Inf, s$slack[up[pend]])
				b2 <- (s$slack[lo[pend]] - r * x1) / sr
				ok <- stats::runif(length(pend)) < stats::pnorm(b2)
				if (any(ok)) {
					w <- rtnorm_interval_logp(rep(0, sum(ok)), -Inf, b2[ok])
					e1[pend[ok]] <- x1[ok]
					e2[pend[ok]] <- r * x1[ok] + sr * w
					pend <- pend[!ok]
				}
			}
			if (length(pend)) {
				# a pair the rejection did not reach keeps its value until the next sweep
				e1[pend] <- e[up[pend]]; e2[pend] <- e[lo[pend]]
			}
			e[up] <- e1; e[lo] <- e2
			z[up] <- eta[up] + e1; z[lo] <- eta[lo] + e2
		}
		# free cells: jointly when both are free, otherwise given the partner
		if (any(s$gg)) {
			up <- s$up[s$gg]; lo <- s$lo[s$gg]
			x1 <- stats::rnorm(length(up))
			e[up] <- x1; e[lo] <- r * x1 + sr * stats::rnorm(length(up))
			z[up] <- eta[up] + e[up]; z[lo] <- eta[lo] + e[lo]
		}
		if (any(s$g1x)) {
			up <- s$up[s$g1x]; lo <- s$lo[s$g1x]
			e[up] <- r * e[lo] + sr * stats::rnorm(length(up)); z[up] <- eta[up] + e[up]
		}
		if (any(s$g2x)) {
			up <- s$up[s$g2x]; lo <- s$lo[s$g2x]
			e[lo] <- r * e[up] + sr * stats::rnorm(length(lo)); z[lo] <- eta[lo] + e[lo]
		}
		z[s$dg] <- eta[s$dg] + sqrt(1 + r) * stats::rnorm(length(s$dg))
		Z[, , t] <- z
	}
	Z
}
.lame_nom_rho_collapsed <- function(
	Z,
	EZ,
	rho,
	Y,
	family,
	YL,
	odmax,
	pw = .lame_rho_prior_power(),
	tries = 4L,
	step = 0.12
) {
	none <- list(Z = Z, rho = rho, accepts = 0L)
	if (!is.finite(rho) || abs(rho) >= 1) return(none)
	d <- .lame_nom_rho_dyads(Z, EZ, Y, family, YL, odmax)
	cur <- .lame_nom_rho_logf(rho, d, pw)
	if (!is.finite(cur)) return(none)
	r <- rho
	acc <- 0L
	for (k in seq_len(tries)) {
		if (stats::runif(1L) < 0.1) {
			# independent uniform draw on (-0.99, 0.99)
			if (abs(r) >= 0.99) next
			prop <- stats::runif(1L, -0.99, 0.99)
			lq <- 0
		} else {
			prop <- r + step * stats::rnorm(1L)
			# reflected walk, symmetric on (-1, 1)
			prop <- ((prop + 1) %% 4) - 1
			if (prop > 1) prop <- 2 - prop
			if (abs(prop) >= 1) next
			lq <- 0
		}
		new <- .lame_nom_rho_logf(prop, d, pw)
		if (is.finite(new) && log(stats::runif(1L)) < new - cur + lq) {
			r <- prop
			cur <- new
			acc <- acc + 1L
		}
	}
	if (r == rho) return(none)
	list(Z = .lame_nom_rho_redraw(Z, d, r), rho = r, accepts = acc)
}
