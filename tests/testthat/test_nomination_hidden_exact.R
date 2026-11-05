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
