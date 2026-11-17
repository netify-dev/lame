# exactness of the cbin / frn latent samplers with hidden cells (R/rZ_nomination_hidden.R): one sweep from the joint
# law of latent array and data must keep that law, so paired differences of any statistic before and after have mean zero

# the rule: row i nominates its min(cap, number of positive eligible cells) strongest cells; frn records ranks
hx_rule = function(Z, cap, elig, frn) {
	Y = matrix(0, nrow(Z), ncol(Z))
	for (i in seq_len(nrow(Z))) {
		z = Z[i, ]; z[!elig[i, ]] = -Inf
		pos = which(z > 0); k = min(cap, length(pos))
		if (k > 0) {
			top = pos[order(z[pos], decreasing = TRUE)[seq_len(k)]]
			Y[i, top] = if (frn) rank(z[top]) else 1
		}
	}
	Y
}

# one sweep from the joint law: z values returned for the paired differences of the statistics
hx_joint_law = function(family, bip, rho, absent, reps = 2500, seed = 1, n = 5L, cap = 2, share = 0.25) {
	set.seed(seed)
	frn = family == 'frn'
	nc = if (bip) n - 1L else n
	eta = -0.3 + outer(rnorm(n, 0, .5), rnorm(nc, 0, .5), '+')
	elig = matrix(1L, n, nc)
	if (!bip) diag(elig) = 0L
	if (absent) { elig[n, ] = 0L; elig[, nc] = 0L }
	cells = which(elig == 1L, arr.ind = TRUE)
	# statistics of the array and of its hidden cells relative to the observed cells of their row; the hidden cells are
	# read from the data, not the package, so the check does not depend on what it tests
	stat_of = function(Z, hid, n1) {
		v = unlist(lapply(seq_len(nrow(cells)), function(q) { z = Z[cells[q, 1], cells[q, 2]]; c(z, z^2, z > 0) }))
		lab = vapply(seq_len(n), function(i) {
			el = which(elig[i, ] == 1L)
			sum(vapply(el[hid[i, el]], function(j) Z[i, j] > 0 && sum(Z[i, el] > Z[i, j]) < cap, NA))
		}, 1)
		K = vapply(seq_len(n), function(i) sum(Z[i, elig[i, ] == 1L] > 0), 1)
		g = c(full = 1, room = 2, none = 3)
		cls = ifelse(n1 >= cap, 1, ifelse(n1 == 0, 3, 2))
		hz = vapply(1:3, function(c) sum(Z[hid & rep(cls == c, nc)]), 1)
		hp = vapply(1:3, function(c) sum(Z[hid & rep(cls == c, nc)] > 0), 1)
		# the products of the two cells of every dyad, those of an absent actor's dyads included (cells no nomination can sit in)
		pr = if (!bip && rho != 0) vapply(
			which(upper.tri(Z)),
			function(p) { ij = arrayInd(p, dim(Z)); (Z[ij[1], ij[2]] - eta[ij[1], ij[2]]) * (Z[ij[2], ij[1]] - eta[ij[2], ij[1]]) },
			1
		) else NULL
		c(v, sum(lab), lab, K, hz, hp, pr)
	}
	d = NULL
	for (r in seq_len(reps)) {
		e1 = matrix(rnorm(n * nc), n, nc); e2 = matrix(rnorm(n * nc), n, nc)
		Z = eta + e1
		if (!bip) {
			low = row(Z) > col(Z)
			Z[low] = (eta + rho * t(e1) + sqrt(1 - rho^2) * e2)[low]
			diag(Z) = diag(eta) + sqrt(1 + rho) * rnorm(n)
		}
		Y = hx_rule(Z, rep(cap, n), elig, frn)
		Y[elig == 0L] = NA
		Y[sample(which(elig == 1L), round(share * sum(elig == 1L)))] = NA
		hid = is.na(Y) & elig == 1L
		if (frn) for (i in seq_len(n)) { k = which(!is.na(Y[i, ]) & Y[i, ] > 0); if (length(k)) Y[i, k] = rank(Y[i, k]) }
		odobs = rowSums(!is.na(Y) & Y > 0)
		n1 = odobs
		YL = NULL
		if (frn) {
			YL = matrix(NA_integer_, n, max(c(1, Y), na.rm = TRUE))
			for (i in seq_len(n)) for (q in seq_len(ncol(YL))) { j = which(Y[i, ] == q); if (length(j)) YL[i, q] = j }
		}
		Z1 = if (bip) { if (frn) lame:::rZ_frn_bip_fc(
			Z,
			eta,
			Y,
			YL,
			rep(cap, n),
			odobs
		) else lame:::rZ_cbin_bip_fc(Z, eta, Y, rep(cap, n), odobs) }
		else { if (frn) lame:::rZ_frn_fc(
			Z,
			eta,
			rho,
			Y,
			YL,
			rep(cap, n),
			odobs
		) else lame:::rZ_cbin_fc(Z, eta, rho, Y, rep(cap, n), odobs) }
		d = rbind(d, stat_of(Z1, hid, n1) - stat_of(Z, hid, n1))
	}
	keep = apply(d, 2, sd) > 1e-12
	colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrow(d)))
}

test_that('one sweep of the hidden-cell samplers leaves the joint law of latent values and data unchanged', {
	skip_on_cran()
	# family, bipartite, dyadic correlation, an actor absent from the slice
	cfg = list(
		list('cbin', FALSE, 0.5, TRUE),
		list('frn', FALSE, 0.4, TRUE),
		list('cbin', TRUE, 0, FALSE),
		list('frn', TRUE, 0, TRUE),
		list('cbin', FALSE, 0, FALSE),
		list('cbin', FALSE, 0.5, FALSE)
	)
	for (k in seq_along(cfg)) {
		z = hx_joint_law(cfg[[k]][[1]], cfg[[k]][[2]], rho = cfg[[k]][[3]], absent = cfg[[k]][[4]], seed = 1100 + k)
		expect_lt(max(abs(z)), 4.6)
	}
	# the check has power: a sampler that draws the hidden cells from their prior, as cells no nomination can sit in, fails it
	testthat::local_mocked_bindings(
		.package = 'lame',
		.lame_nom_hidden_cells = function(Y, bip) matrix(FALSE, nrow(Y), ncol(Y))
	)
	expect_gt(max(abs(hx_joint_law('cbin', FALSE, 0, FALSE, reps = 1500, seed = 1112))), 8)
})

# cbin / frn fits with hidden cells: the latent draws, the predictive nominations, the row-shift move and the reports

# a panel of cbin nominations with a fifth of the cells hidden
.hidden_cbin_panel = function(n = 12, cap = 3, mu = -0.2, seed = 41, bip = FALSE) {
	set.seed(seed)
	nb = if (bip) 9 else n
	z = mu + outer(rnorm(n, 0, 0.6), rnorm(nb, 0, 0.4), "+") + matrix(rnorm(n * nb), n, nb)
	if (!bip) diag(z) = -Inf
	y = matrix(0, n, nb)
	for (i in seq_len(n)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		y[i, top[z[i, top] > 0]] = 1
	}
	if (!bip) diag(y) = NA
	off = which(!is.na(y))
	y[sample(off, round(0.2 * length(off)))] = NA
	y
}

# the nomination rule applied to a latent row set: the top min(cap, #positive)
# positive cells of each row
.nominate = function(Z, cap, sq) {
	if (sq) diag(Z) = -Inf
	t(apply(Z, 1, function(z) {
		k = min(cap, sum(z > 0)); y = numeric(length(z))
		if (k > 0) y[order(z, decreasing = TRUE)[seq_len(k)]] = 1
		y
	}))
}

test_that("cbin latent draws stay consistent with the observed cells when cells are hidden", {
	skip_on_cran()
	for (bip in c(FALSE, TRUE)) {
		Y = .hidden_cbin_panel(bip = bip)
		n = nrow(Y); nb = ncol(Y); cap = 3
		EZ = matrix(-0.2, n, nb)
		Z = matrix(rnorm(n * nb), n, nb)
		na_off = is.na(Y); if (!bip) diag(na_off) = FALSE
		odobs = rowSums(Y > 0, na.rm = TRUE)
		obs = !is.na(Y)
		ok = logical(200)
		for (s in 1:200) {
			room = odobs + rowSums(na_off & Z > 0)
			Z = if (bip) lame:::rZ_cbin_bip_fc(Z, EZ, Y, rep(cap, n), room)
			else lame:::rZ_cbin_fc(Z, EZ, 0.2, Y, rep(cap, n), room)
			ok[s] = all(.nominate(Z, cap, !bip)[obs] == Y[obs])
		}
		# after every sweep the nominations the latent values imply match
		# every observed cell
		expect_true(all(ok))
	}
})

test_that("a hidden cell of a row whose observed nominations fill the cap stays below them", {
	skip_on_cran()
	n = 10; cap = 3
	Y = matrix(0, n, n); diag(Y) = NA
	Y[1, 2:4] = 1; Y[1, 8:10] = NA
	EZ = matrix(0.5, n, n)
	Z = matrix(rnorm(n * n), n, n)
	set.seed(3)
	hid_max = nom_min = numeric(300)
	for (s in 1:300) {
		Z = lame:::rZ_cbin_fc(
			Z,
			EZ,
			0,
			Y,
			rep(cap, n),
			rowSums(Y > 0, na.rm = TRUE) + rowSums(is.na(Y) & row(Y) != col(Y) & Z > 0)
		)
		hid_max[s] = max(Z[1, 8:10]); nom_min[s] = min(Z[1, 2:4])
	}
	# the cap is spent on the observed nominations, so no hidden cell of the
	# row can outrank them
	expect_true(all(hid_max < nom_min))
})

test_that("complete cbin panels take the unchanged draw", {
	skip_on_cran()
	set.seed(5); n = 12
	Y = matrix(0, n, n); for (i in 1:n) Y[i, sample(setdiff(1:n, i), sample(1:4, 1))] = 1; diag(Y) = NA
	Z = matrix(rnorm(n * n), n, n); EZ = matrix(rnorm(n * n, -0.5), n, n)
	# the sweep first repairs a start outside the nomination constraints, so the reference starts from the repaired state
	Z = lame:::.lame_nom_feasible(Z, Y, NULL, rep(3, n), rowSums(Y > 0, na.rm = TRUE), FALSE)
	set.seed(9); z_new = lame:::rZ_cbin_fc(Z, EZ, 0.3, Y, rep(3, n), rowSums(Y > 0, na.rm = TRUE))
	set.seed(9); z_old = lame:::.rZ_cbin_sweep(Z, EZ, 0.3, Y, rep(3, n), rowSums(Y > 0, na.rm = TRUE))
	expect_identical(z_new, z_old)
})

test_that("the likelihood-ratio saturation flag needs an estimate near the threshold", {
	skip_on_cran()
	thr = lame:::.lame_odmax_sat_threshold
	# half the rows at a cap of 5, a fifth of the cells hidden: the estimate
	# sits well below the threshold, and the flag must not be set whatever
	# the likelihood-ratio test says
	n_checked = 0
	for (seed in 1:6) {
		set.seed(seed)
		Ys = lapply(1:4, function(t) {
			m = matrix(0, 30, 30)
			for (i in 1:30) m[i, sample(setdiff(1:30, i), sample(c(2, 3, 4, 5, 5, 5), 1))] = 1
			diag(m) = NA
			off = which(!is.na(m)); m[sample(off, round(0.2 * length(off)))] = NA
			m
		})
		sh = lame:::.lame_odmax_saturated_share(Ys, 5, "cbin")
		if (sh < thr - 0.08) {
			n_checked = n_checked + 1
			expect_false(isTRUE(attr(sh, "not_ruled_out")))
			expect_true(is.na(lame:::.lame_ab_var_unident_reason(list(
				family = "cbin",
				mode = "unipartite",
				odmax_saturated_share = sh
			))))
		}
	}
	expect_gt(n_checked, 3)
})

test_that("the compiled hidden-cell draw gives the R reference's draws", {
	skip_on_cran()
	for (bip in c(FALSE, TRUE)) {
		Y = .hidden_cbin_panel(bip = bip, seed = 43)
		n = nrow(Y); nb = ncol(Y)
		hid = lame:::.lame_nom_hidden_cells(Y, bip)
		set.seed(8); Z = matrix(rnorm(n * nb), n, nb); EZ = matrix(rnorm(n * nb, -0.3), n, nb)
		expect_identical(
			lame:::.lame_cbin_implied_labels(Z, Y, 3, hid),
			.lame_cbin_implied_labels_r(Z, Y, 3, hid)
		)
		rho = if (bip) 0 else 0.3
		set.seed(9); z_c = lame:::.lame_cbin_hidden_moves(Z, EZ, rho, Y, 3, hid, bip)
		set.seed(9); z_r = .lame_cbin_hidden_moves_r(Z, EZ, rho, Y, 3, hid, bip)
		expect_equal(z_c, z_r, tolerance = 1e-10)
	}
})

test_that("frn latent draws keep the observed order and membership when cells are hidden", {
	skip_on_cran()
	set.seed(44)
	n = 12; cap = 4
	z = -0.2 + outer(rnorm(n, 0, 0.6), rnorm(n, 0, 0.4), "+") + matrix(rnorm(n * n), n, n)
	diag(z) = -Inf
	Y = matrix(0, n, n)
	for (i in seq_len(n)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		top = top[z[i, top] > 0]
		if (length(top)) Y[i, top] = rank(z[i, top])
	}
	diag(Y) = NA
	off = which(!is.na(Y)); Y[sample(off, round(0.2 * length(off)))] = NA
	YL = t(apply(Y, 1, function(y) match(seq_len(cap), y)))
	EZ = matrix(-0.2, n, n); Z = matrix(rnorm(n * n), n, n)
	na_off = is.na(Y); diag(na_off) = FALSE
	obs = !is.na(Y)
	# from a random start the order constraints take a few sweeps to be met
	# (as with complete panels); from then on every sweep must meet them
	ok = logical(150)
	for (s in 1:160) {
		Z = lame:::rZ_frn_fc(Z, EZ, 0.2, Y, YL, rep(cap, n), rowSums(Y > 0, na.rm = TRUE) + rowSums(na_off & Z > 0))
		if (s <= 10) next
		nom = .nominate(Z, cap, TRUE)
		ord_ok = all(vapply(seq_len(n), function(i) {
			j = which(obs[i, ] & Y[i, ] > 0)
			length(j) < 2 || all(diff(Z[i, j[order(Y[i, j])]]) > 0)
		}, logical(1)))
		ok[s - 10] = all((nom[obs] > 0) == (Y[obs] > 0)) && ord_ok
	}
	# after every sweep the nominations the latent values imply match every
	# observed cell and the observed nominations keep their observed order
	expect_true(all(ok))
})

test_that("the posterior-predictive nomination draw lets hidden cells hold nominations", {
	skip_on_cran()
	Y = .hidden_cbin_panel()
	n = nrow(Y)
	Y[5, ] = NA; Y[, 5] = NA
	hid = lame:::.lame_nom_hidden_cells(Y, FALSE)
	YO = lame:::.lame_nom_predictive_yo(Y, FALSE)
	# hidden cells stay in contention; the diagonal and the absent actor do not
	expect_true(all(YO[hid] == 0))
	expect_true(all(is.na(diag(YO))) && all(is.na(YO[5, ])) && all(is.na(YO[, 5])))
	expect_identical(is.na(YO), !(!is.na(Y) | hid))
	EZ = matrix(0.3, n, n)
	set.seed(2); d1 = simY_frn(EZ, 0, 3, YO = YO)
	YO2 = matrix(0, n, n); YO2[5, ] = NA; YO2[, 5] = NA
	set.seed(2); d2 = simY_frn(EZ, 0, 3, YO = YO2)
	# the same draw as a panel with no hidden cell: every row fills its cap
	# over the cells a nomination can sit in
	expect_identical(d1, d2)
	expect_true(all(rowSums(d1 > 0, na.rm = TRUE)[-5] <= 3))
	expect_true(sum(d1[hid] > 0) > 0)
})

test_that("YPM of a cbin fit with hidden cells matches the observed nomination rate", {
	skip_on_cran()
	set.seed(61)
	n = 20; TT = 2; cap = 4
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Y = lapply(seq_len(TT), function(t) {
		z = 0.3 + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
		diag(z) = -Inf
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			y[i, top[z[i, top] > 0]] = 1
		}
		diag(y) = NA
		off = which(!is.na(y)); y[sample(off, round(0.2 * length(off)))] = NA
		dimnames(y) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n))
		y
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "cbin",
		odmax = cap,
		R = 0,
		nscan = 400,
		burn = 100,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 3
	)))
	obs = lapply(Y, function(y) !is.na(y))
	yp = unlist(lapply(seq_len(TT), function(t) fit$YPM[[t]][obs[[t]]]))
	yo = unlist(lapply(seq_len(TT), function(t) Y[[t]][obs[[t]]]))
	# on the observed cells the predictive mean matches the observed mean
	expect_lt(abs(mean(yp) / mean(yo) - 1), 0.08)
	hid = lame:::.lame_nom_hidden_cells(Y[[1]], FALSE)
	expect_true(all(is.finite(fit$YPM[[1]][hid])))
	expect_true(all(is.na(diag(fit$YPM[[1]]))))
	# the settling check of hidden-cell fits is stored
	expect_named(fit$hidden_cell_ess, c("intercept", "va", "draws"))
})

# brute-force check of the shifts that keep a latent row consistent with its
# observed cells, over a grid of shifts
.shift_feasible = function(Z, Y, cap, i, deltas) {
	TT = dim(Z)[3]
	rows = lapply(seq_len(TT), function(t) {
		y = Y[, , t]
		hid = lame:::.lame_nom_hidden_cells(y, nrow(y) != ncol(y))
		elig = !is.na(y) | hid
		if (nrow(y) == ncol(y)) diag(elig) = FALSE
		o = !is.na(y[i, ]); if (nrow(y) == ncol(y)) o[i] = FALSE
		zi = Z[i, , t]; zi[!elig[i, ]] = -Inf
		list(z = zi, o = o, y = y[i, ])
	})
	vapply(deltas, function(d) {
		all(vapply(rows, function(r) {
			zi = r$z + d
			k = min(cap, sum(zi > 0)); nom = numeric(length(zi))
			if (k > 0) nom[order(zi, decreasing = TRUE)[seq_len(k)]] = 1
			all((nom[r$o] > 0) == (r$y[r$o] > 0))
		}, logical(1)))
	}, logical(1))
}

test_that("the row-shift bounds are exactly the shifts that keep a row consistent", {
	skip_on_cran()
	for (bip in c(FALSE, TRUE)) {
		Y1 = .hidden_cbin_panel(bip = bip, seed = 41); Y2 = .hidden_cbin_panel(bip = bip, seed = 42)
		if (!bip) { Y2[7, ] = NA; Y2[, 7] = NA }
		Y = array(c(Y1, Y2), c(dim(Y1), 2))
		n = nrow(Y1); nb = ncol(Y1); cap = 3
		Z = array(rnorm(n * nb * 2), c(n, nb, 2))
		# a consistent state from the latent draw
		for (s in 1:40) for (t in 1:2) {
			na_off = is.na(Y[, , t]); if (!bip) diag(na_off) = FALSE
			room = rowSums(Y[, , t] > 0, na.rm = TRUE) + rowSums(na_off & Z[, , t] > 0)
			Z[, , t] = if (bip) lame:::rZ_cbin_bip_fc(Z[, , t], matrix(-0.2, n, nb), Y[, , t], rep(cap, n), room)
			else lame:::rZ_cbin_fc(Z[, , t], matrix(-0.2, n, nb), 0, Y[, , t], rep(cap, n), room)
		}
		bd = lame:::.lame_nom_shift_bounds(Z, cap, lame:::.lame_nom_shift_masks(Y, bip))
		deltas = seq(-3, 3, by = 0.005)
		for (i in seq_len(n)) {
			f = .shift_feasible(Z, Y, cap, i, deltas)
			expect_true(f[deltas == 0])
			# feasible exactly on (lo, hi), up to the grid
			inside = deltas > bd$lo[i] + 0.005 & deltas < bd$hi[i] - 0.005
			outside = deltas < bd$lo[i] - 0.005 | deltas > bd$hi[i] + 0.005
			expect_true(all(f[inside]) && !any(f[outside]))
		}
	}
})

test_that("the row-shift move keeps the panel consistent and draws the row effect from its truncated prior", {
	skip_on_cran()
	Y1 = .hidden_cbin_panel(seed = 41); Y2 = .hidden_cbin_panel(seed = 43)
	Y = array(c(Y1, Y2), c(dim(Y1), 2))
	n = nrow(Y1); cap = 3
	masks = lame:::.lame_nom_shift_masks(Y, FALSE)
	Z = array(rnorm(n * n * 2), c(n, n, 2))
	a = rnorm(n, 0, 0.3); b = rnorm(n, 0, 0.3)
	Sab = matrix(c(0.3, 0.05, 0.05, 0.2), 2, 2)
	obs = !is.na(Y)
	for (t in 1:2) { o = obs[, , t]; diag(o) = FALSE; obs[, , t] = o }
	ok = logical(100)
	for (s in 1:140) {
		for (t in 1:2) {
			na_off = is.na(Y[, , t]); diag(na_off) = FALSE
			room = rowSums(Y[, , t] > 0, na.rm = TRUE) + rowSums(na_off & Z[, , t] > 0)
			Z[, , t] = lame:::rZ_cbin_fc(Z[, , t], matrix(-0.2, n, n) + a, 0, Y[, , t], rep(cap, n), room)
		}
		sh = lame:::.lame_nom_row_shift(Z, a, b, Sab, cap, masks, FALSE, rep(TRUE, n))
		Z = sh$Z; a = sh$a
		if (s > 40) ok[s - 40] = all(vapply(1:2, function(t) {
			all((.nominate(Z[, , t], cap, TRUE)[obs[, , t]] > 0) == (Y[, , t][obs[, , t]] > 0))
		}, logical(1)))
	}
	expect_true(all(ok))
	# from one state, many draws of the move: the row effect follows its
	# conditional prior truncated to (a + lo, a + hi)
	bd = lame:::.lame_nom_shift_bounds(Z, cap, masks)
	r = Sab[1, 2] / Sab[2, 2]; s = sqrt(Sab[1, 1] - r * Sab[1, 2])
	m = b * r
	draws = replicate(4000, lame:::.lame_nom_row_shift(Z, a, b, Sab, cap, masks, FALSE, rep(TRUE, n))$a)
	for (i in seq_len(n)) {
		lo = (a[i] + bd$lo[i] - m[i]) / s; hi = (a[i] + bd$hi[i] - m[i]) / s
		mu = m[i] + s * (dnorm(lo) - dnorm(hi)) / (pnorm(hi) - pnorm(lo))
		expect_lt(abs(mean(draws[i, ]) - mu), 4 * sd(draws[i, ]) / sqrt(4000) + 1e-8)
	}
})

test_that("simulate() leaves the actors absent from a period out of the nomination draw", {
	skip_on_cran()
	set.seed(71)
	n = 14; cap = 3
	nm = sprintf("a%02d", 1:n)
	Y = lapply(1:2, function(t) {
		keep = if (t == 1) nm else nm[-(1:3)]
		k = length(keep)
		z = 1.5 + matrix(rnorm(k * k), k, k); diag(z) = -Inf
		y = matrix(0, k, k)
		for (i in seq_len(k)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			y[i, top[z[i, top] > 0]] = 1
		}
		diag(y) = NA; dimnames(y) = list(keep, keep)
		y
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "cbin",
		odmax = cap,
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5
	)))
	s = simulate(fit, nsim = 20, seed = 3)
	y2 = s$Y[[1]][[2]]
	expect_identical(dim(y2), c(11L, 11L))
	rs = unlist(lapply(s$Y, function(traj) rowSums(traj[[2]] > 0, na.rm = TRUE)))
	# the rows of the second period fill their cap over the actors present
	expect_gt(mean(rs == cap), 0.85)
})

test_that("a bipartite ame() cbin fit with hidden cells predicts the observed nomination rate", {
	skip_on_cran()
	set.seed(63)
	nA = 24; nB = 18; cap = 6
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	z = 0.2 + outer(a, b, "+") + matrix(rnorm(nA * nB), nA, nB)
	Y = matrix(0, nA, nB)
	for (i in seq_len(nA)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		Y[i, top[z[i, top] > 0]] = 1
	}
	dimnames(Y) = list(sprintf("r%02d", 1:nA), sprintf("c%02d", 1:nB))
	Y[sample(length(Y), round(0.2 * length(Y)))] = NA
	fit = suppressWarnings(suppressMessages(ame(
		Y,
		family = "cbin",
		mode = "bipartite",
		odmax = cap,
		nscan = 400,
		burn = 100,
		odens = 4,
		verbose = FALSE,
		gof = TRUE,
		seed = 4
	)))
	obs = !is.na(Y)
	expect_lt(abs(mean(fit$YPM[obs]) / mean(Y[obs]) - 1), 0.08)
	expect_true(all(is.finite(fit$YPM[!obs])))
	expect_true(all(is.finite(fit$VC[, "va"])))
	expect_true(all(is.finite(fit$GOF)))
})

test_that("a hidden-cell cbin / frn fit records the effective sample size of its intercept and va", {
	skip_on_cran()
	set.seed(9)
	slow = cumsum(rnorm(200, 0, 0.05)); fast = rnorm(200)
	fake = list(
		family = "cbin",
		BETA = cbind(intercept = slow, x1_dyad = fast),
		VC = cbind(va = fast + 3, vb = fast + 3)
	)
	he = lame:::.lame_nom_hidden_ess(fake, TRUE)
	expect_lt(he[["intercept"]], 20)
	expect_gt(he[["va"]], 100)
	expect_equal(he[["draws"]], 200)
	expect_null(lame:::.lame_nom_hidden_ess(fake, FALSE))
	expect_null(lame:::.lame_nom_hidden_ess(modifyList(fake, list(family = "binary")), TRUE))
	# a fit with hidden cells carries the record and reports the intercept and va whatever the chain's effective size
	Y = .hidden_cbin_panel(n = 14, seed = 47)
	fit = suppressWarnings(suppressMessages(ame(
		Y,
		family = "cbin",
		odmax = 3,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)))
	expect_named(fit$hidden_cell_ess, c("intercept", "va", "draws"))
	fit$hidden_cell_ess[["intercept"]] = 5
	expect_false(grepl("not usable", paste(capture.output(print(summary(fit))), collapse = " "), ignore.case = TRUE))
	expect_false(any(is.na(summary(fit)$variance["va", c("Estimate", "StdError")])))
	# a complete panel carries none
	Yc = Y; Yc[is.na(Yc)] = 0; diag(Yc) = NA
	fit_c = suppressWarnings(suppressMessages(ame(
		Yc,
		family = "cbin",
		odmax = 3,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)))
	expect_null(fit_c$hidden_cell_ess)
})
