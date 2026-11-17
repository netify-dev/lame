# R references of the compiled kernels of the Poisson translation moves, the symmetric actor-period
# move, the nomination feasibility check and the tempered transition. each takes the arguments of its
# kernel, and local_move_r_engine() runs the moves on the references in place of the kernels

# R reference of the transition, with the same steps, sums (long double, as in C++) and random numbers;
# `z` the mirrored latent array, `root_n` the lower Cholesky factor of a path's AR(1) correlation
.lame_uv_sym_temper_r = function(U_cube, L, z, dt, r, root_n, l_var, bs, cycles, max_abs) {
	n = dt$n; N = dt$N
	K = length(bs) - 1L
	U = U_cube
	Lr = L[r]
	okm = dt$ok
	Y0 = dt$Y
	Y0[!is.finite(Y0)] = 0
	E = array(0, c(n, n, N))
	for (t in seq_len(N)) E[, , t] = ifelse(okm[, , t], exp(z[, , t] + dt$off[t]), 0)
	up = upper.tri(matrix(0, n, n))
	ll = sum(
		unlist(
			lapply(
				seq_len(N),
				function(t) { o = okm[, , t] & up; (Y0[, , t] * z[, , t] - exp(z[, , t] + dt$off[t]))[o] }
			)
		)
	)
	rail = is.finite(max_abs)
	actor = function(a, b) {
		ui = U[a, r, ]
		uj = matrix(U[, r, ], n, N)
		diffv = function(x) {
			du = x - ui
			dz = Lr * uj * matrix(rep(du, each = n), n, N)
			em = expm1(dz)
			dz[a, ] = 0
			em[a, ] = 0
			cells = matrix(okm[a, , ], n, N)
			cells[a, ] = FALSE
			list(dz = dz, em = em, val = sum((matrix(Y0[a, , ], n, N) * dz - matrix(E[a, , ], n, N) * em)[cells]))
		}
		e = stats::rnorm(N)
		nu = vapply(seq_len(N), function(t) { m = 0; for (q in seq_len(t)) m = m + root_n[t, q] * e[q]; m }, 0)
		thr = log(stats::runif(1L))
		th = stats::runif(1L, 0, 2 * pi)
		lo = th - 2 * pi
		hi = th
		ev = 0L
		repeat {
			xp = ui * cos(th) + nu * sin(th)
			ev = ev + 1L
			feasible = !(rail && any(abs(xp) > max_abs))
			if (feasible) {
				d = diffv(xp)
				if (is.finite(d$val) && b * d$val > thr) break
			}
			if (th < 0) lo = th else hi = th
			if (ev >= 60L) return(invisible(NULL))
			th = stats::runif(1L, lo, hi)
		}
		for (t in seq_len(N)) {
			for (j in seq_len(n)) {
				if (j == a) next
				z[a, j, t] <<- z[a, j, t] + d$dz[j, t]
				z[j, a, t] <<- z[j, a, t] + d$dz[j, t]
				if (okm[a, j, t]) { f = 1 + d$em[j, t]; E[a, j, t] <<- E[a, j, t] * f; E[j, a, t] <<- E[j, a, t] * f }
			}
			z[a, a, t] <<- z[a, a, t] + Lr * (xp[t] * xp[t] - ui[t] * ui[t])
			U[a, r, t] <<- xp[t]
		}
		ll <<- ll + d$val
		invisible(NULL)
	}
	weight = function(b) {
		L0 = Lr
		base = array(0, c(n, n, N))
		for (t in seq_len(N)) base[, , t] = tcrossprod(U[, r, t])
		diffv = function(v) sum(
			unlist(
				lapply(
					seq_len(N),
					function(t) { o = okm[, , t] & up; dz = (v - L0) * base[, , t]; (Y0[, , t] * dz - E[, , t] * expm1(dz))[o] }
				)
			)
		)
		nuL = sqrt(l_var) * stats::rnorm(1L)
		thr = log(stats::runif(1L))
		th = stats::runif(1L, 0, 2 * pi)
		lo = th - 2 * pi
		hi = th
		ev = 0L
		repeat {
			xp = L0 * cos(th) + nuL * sin(th)
			ev = ev + 1L
			val = diffv(xp)
			if (is.finite(val) && b * val > thr) {
				dv = xp - L0
				for (t in seq_len(N)) {
					dz = dv * base[, , t]
					z[, , t] <<- z[, , t] + dz
					E[, , t] <<- ifelse(okm[, , t], E[, , t] * (1 + expm1(dz)), E[, , t])
				}
				Lr <<- xp
				ll <<- ll + val
				return(invisible(NULL))
			}
			if (th < 0) lo = th else hi = th
			if (ev >= 60L) return(invisible(NULL))
			th = stats::runif(1L, lo, hi)
		}
	}
	pass = function(b, rev, with_weight) {
		steps = c(sample.int(n), if (with_weight) 0L)
		if (rev) steps = rev(steps)
		for (st in steps) if (st == 0L) weight(b) else actor(st, b)
	}
	logacc = 0
	for (j in seq_len(K)) {
		logacc = logacc + (bs[j + 1L] - bs[j]) * ll
		for (cc in seq_len(cycles)) pass(bs[j + 1L], FALSE, cc == 1L)
	}
	for (j in K:1) {
		logacc = logacc + (bs[j] - bs[j + 1L]) * ll
		for (cc in seq_len(cycles)) pass(bs[j + 1L], TRUE, cc == cycles)
	}
	Lo = L
	Lo[r] = Lr
	list(z = z, U = U, L = Lo, logacc = logacc, ll = ll)
}

# the rows that break the constraints, in R (the reference of nom_bad_rows_cpp)
.lame_nom_bad_rows_r = function(Z, Y, YL, room, bip) {
	nom = !is.na(Y) & Y > 0
	nr = nrow(Z)
	if (!any(nom)) return(rep(FALSE, nr))
	non = !is.na(Y) & Y == 0
	Zn = Z
	Zn[!nom] = Inf
	mn = apply(Zn, 1L, min)
	Zo = Z
	Zo[!non] = -Inf
	mx = apply(Zo, 1L, max)
	hrow = rowSums(.lame_nom_hidden_cells(Y, bip)) > 0L
	room_chk = room & !hrow
	bad = is.finite(mn) & (mn <= 0 | (is.finite(mx) & mx >= mn) | (room_chk & is.finite(mx) & mx >= 0))
	if (!is.null(YL)) {
		for (i in which(is.finite(mn) & !bad)) {
			cols = YL[i, ]
			cols = cols[!is.na(cols)]
			if (length(cols) > 1L && any(diff(Z[i, cols]) <= 0)) bad[i] = TRUE
		}
	}
	bad
}

# the tempered transition with the arguments of uv_temper_transition_cpp (`r` zero-based, the data as vectors)
uv_temper_transition_r = function(z, U, L, r, Y, ok, off, root, lvar, bs, cycles, max_abs) {
	N = length(off); n = as.integer(round(sqrt(length(z) / N))); R = length(L)
	dt = list(Y = array(Y, c(n, n, N)), ok = array(ok == 1L, c(n, n, N)), off = off, n = n, N = N)
	.lame_uv_sym_temper_r(array(U, c(n, R, N)), L, array(z, c(n, n, N)), dt, r + 1L, root, lvar, bs, cycles, max_abs)
}

# the rows that break the constraints with the arguments of nom_bad_rows_cpp (an empty `YL` for cbin)
nom_bad_rows_r = function(Z, Y, YL, room, bip) {
	.lame_nom_bad_rows_r(Z, Y, if (length(YL)) YL else NULL, room, bip)
}

# elliptical slice step of a Poisson translation d of total count `ys` and rate sum `s` (ess_pois_scalar_cpp)
ess_pois_scalar_r = function(ys, s, mu, sd) {
	.lame_ess_step(0, mu, matrix(sd), function(d) ys * d - s * expm1(d), 0)
}

# the same for cells with rates `ez` moving by d x (ess_pois_cells_cpp)
ess_pois_cells_r = function(ez, x, ysx, mu, sd) {
	.lame_ess_step(0, mu, matrix(sd), function(d) d * ysx - ppm_shift_sum_cpp(ez, x, d), 0)
}

# elliptical slice step of a position row `cur` whose line of cells moves by opp (w (x - cur)); `skip` the
# one-based diagonal cell of a symmetric row, 0 otherwise
ess_pois_uv_r = function(cur, mu, root, opp, w, ys, s, skip) {
	ll = if (!any(ys > 0) && !any(s > 0)) function(x) 0 else function(x) {
		d = drop(opp %*% (w * (x - cur)))
		if (skip > 0L) d[skip] = 0
		sum(ys * d - s * expm1(d))
	}
	.lame_ess_step(cur, mu, root, ll, 0)
}

# one margin of the uncentred, static Poisson level translations in the order `perm` (ltm_level_phase_cpp):
# `S`, `ys` the line sums, `E` the rates (a symmetric network), `x` the moving levels, `o` the others, `iS`
# the inverse level covariance (directed), `v` the level variance (symmetric, bipartite), `mode` 0 directed,
# 1 symmetric, 2 bipartite
ltm_level_phase_r = function(S, ys, E, dimE, x, o, iS, v, mode, row, perm) {
	sym = mode == 1L
	if (sym) E = array(E, dimE)
	dsh = numeric(length(x))
	st = c(0, 0, 0, 0)
	for (i in perm) {
		if (mode == 0L) {
			a1 = if (row) iS[1L, 1L] * x[i] + iS[1L, 2L] * o[i] else iS[2L, 2L] * x[i] + iS[1L, 2L] * o[i]
			a2 = if (row) iS[1L, 1L] else iS[2L, 2L]
		} else {
			a1 = x[i] / v
			a2 = 1 / v
		}
		if (!all(is.finite(c(a1, a2))) || a2 <= 0) {
			st[4L] = st[4L] + 1
			next
		}
		ans = ess_pois_scalar_r(ys[i], S[i], -a1 / a2, 1 / sqrt(a2))
		st[1:3] = st[1:3] + c(1, ans$evals, ans$stuck)
		d = ans$x
		if (d == 0) next
		ed = exp(d)
		if (sym) {
			Ei = matrix(E[i, , ], ncol = dimE[3L])
			S = S + (ed - 1) * rowSums(Ei)
			E[i, , ] = Ei * ed
			E[, i, ] = Ei * ed
		}
		S[i] = S[i] * ed
		dsh[i] = d
		x[i] = x[i] + d
	}
	list(x = x, dsh = dsh, stats = st)
}

# one margin of the Poisson position translations in the order `perm` (ltm_uv_phase_cpp): `mov` the moving
# positions, `oth` the other side's (the regression of a directed row on its partner), `opp` the positions the
# cells of a line move by, `sj` the rate sums and `ysj` the counts of the dyads summed over periods; `mode` 0
# directed, 1 symmetric (`w` the weights), 2 bipartite
ltm_uv_phase_r = function(mov, oth, opp, reg, root, w, sj, ysj, mode, row, perm) {
	sym = mode == 1L
	Dsh = matrix(0, nrow(sj), ncol(sj))
	st = c(0, 0, 0)
	for (i in perm) {
		cur = mov[i, ]
		op = if (sym) mov else opp
		mu = if (mode == 0L) drop(reg %*% oth[i, ]) else numeric(length(cur))
		ys = if (row || sym) ysj[i, ] else ysj[, i]
		s = if (row || sym) sj[i, ] else sj[, i]
		ans = ess_pois_uv_r(cur, mu, root, op, w, ys, s, if (sym) i else 0L)
		st = st + c(1, ans$evals, ans$stuck)
		dx = ans$x - cur
		if (!any(dx != 0)) next
		shift = drop(op %*% (w * dx))
		if (sym) {
			dg = sum(w * (ans$x^2 - cur^2))
			shift[i] = 0
			sj[i, ] = sj[i, ] * exp(shift)
			sj[, i] = sj[, i] * exp(shift)
			Dsh[i, ] = Dsh[i, ] + shift
			Dsh[, i] = Dsh[, i] + shift
			Dsh[i, i] = Dsh[i, i] + dg
		} else if (row) {
			sj[i, ] = sj[i, ] * exp(shift)
			Dsh[i, ] = shift
		} else {
			sj[, i] = sj[, i] * exp(shift)
			Dsh[, i] = shift
		}
		mov[i, ] = ans$x
	}
	list(M = mov, sj = sj, Dsh = Dsh, stats = st)
}

# the symmetric actor-period translations of one period in the order `perm` (pab_sym_period_cpp)
pab_sym_period_r = function(Zt, Et, y, off, eff, m, v, centred, perm) {
	n = length(eff)
	pairs = 0L
	evals = 0L
	obs = is.finite(y)
	y0 = ifelse(obs, y, 0)
	apply_shift = function(i, d) {
		Zt[i, ] <<- Zt[i, ] + d
		Zt[, i] <<- Zt[, i] + d
		Et[i, ] <<- Et[i, ] + d
		Et[, i] <<- Et[, i] + d
	}
	stats_of = function(i, k) {
		j = setdiff(seq_len(n), c(i, k))
		o = obs[i, j] & is.finite(Zt[i, j])
		c(sum(y0[i, j][o]), sum(exp(Zt[i, j][o] + off)))
	}
	if (centred) {
		p = perm
		np = length(p) %/% 2L
		for (q in seq_len(np)) {
			i = p[q]
			k = p[np + q]
			si = stats_of(i, k)
			sk = stats_of(k, i)
			prec = 1 / v[i] + 1 / v[k]
			mu = ((m[i] - eff[i]) / v[i] + (eff[k] - m[k]) / v[k]) / prec
			sl = .lame_pab_slice(mu, prec, si[1L] - sk[1L], si[2L], sk[2L])
			d = sl$delta
			if (d != 0) {
				eff[i] = eff[i] + d
				eff[k] = eff[k] - d
				apply_shift(i, d)
				apply_shift(k, -d)
			}
			pairs = pairs + 1L
			evals = evals + sl$evals
		}
	} else {
		for (i in perm) {
			si = stats_of(i, 0L)
			sl = .lame_pab_slice(m[i] - eff[i], 1 / v[i], si[1L], si[2L], 0)
			d = sl$delta
			eff[i] = eff[i] + d
			apply_shift(i, d)
			pairs = pairs + 1L
			evals = evals + sl$evals
		}
	}
	list(Z = Zt, E = Et, eff = eff, pairs = pairs, evals = evals)
}

# runs the Poisson translation moves, the symmetric actor-period move, the nomination feasibility check and the
# tempered transition of the calling frame on the R references in place of their compiled kernels
local_move_r_engine = function(env = parent.frame()) {
	testthat::local_mocked_bindings(
		ess_pois_scalar_cpp = ess_pois_scalar_r,
		ess_pois_cells_cpp = ess_pois_cells_r,
		ltm_level_phase_cpp = ltm_level_phase_r,
		ltm_uv_phase_cpp = ltm_uv_phase_r,
		pab_sym_period_cpp = pab_sym_period_r,
		nom_bad_rows_cpp = nom_bad_rows_r,
		uv_temper_transition_cpp = uv_temper_transition_r,
		.package = "lame",
		.env = env
	)
}
