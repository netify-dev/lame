# a latent slice outside the constraints of the nomination samplers is moved into them before the sweep
# (R/nomination_feasible.R); a slice inside is returned as it is
skip_on_cran()

nf_check = function(Z, Y, YL, odmax, odobs) {
	ok = TRUE
	for (i in seq_len(nrow(Z))) {
		nom = which(!is.na(Y[i, ]) & Y[i, ] > 0); non = which(!is.na(Y[i, ]) & Y[i, ] == 0)
		if (length(nom)) {
			if (min(Z[i, nom]) <= 0) ok = FALSE
			if (length(non) && max(Z[i, non]) >= min(Z[i, nom])) ok = FALSE
			if (length(non) && isTRUE(odobs[i] < odmax[i]) && max(Z[i, non]) >= 0) ok = FALSE
			if (!is.null(YL)) { cc = YL[i, ]; cc = cc[!is.na(cc)]; if (length(cc) > 1L && any(diff(Z[i, cc]) <= 0)) ok = FALSE }
		}
	}
	ok
}

nf_panel = function(fam, nr = 12L, nc = 9L, k = 3L, seed = 4L) {
	set.seed(seed)
	z = matrix(rnorm(nr * nc, -0.5), nr, nc); Y = matrix(0L, nr, nc); YL = matrix(NA_integer_, nr, k)
	for (i in seq_len(nr)) { top = order(z[i, ], decreasing = TRUE)[1:k]; if (fam == 'frn') { Y[i, top] = k:1; YL[i, ] = rev(top) } else Y[i, top] = 1L }
	cand = which(Y == 0); Y[sample(cand, 6)] = NA
	list(Y = Y, YL = if (fam == 'frn') YL else NULL, odmax = rep(k, nr), odobs = rowSums(!is.na(Y) & Y > 0))
}

test_that('a latent slice inside the constraints is returned as it is, and one outside is moved inside', {
	ns = asNamespace('lame')
	for (fam in c('cbin', 'frn')) {
		d = nf_panel(fam)
		set.seed(2)
		# a feasible slice: nominees ordered above zero and above the non-nominees, non-nominees below zero when the row has room
		Z = matrix(rnorm(length(d$Y), -1), nrow(d$Y), ncol(d$Y))
		for (i in seq_len(nrow(Z))) {
			nom = which(!is.na(d$Y[i, ]) & d$Y[i, ] > 0); if (fam == 'frn') nom = nom[order(d$Y[i, nom])]
			Z[i, nom] = 1 + 0.5 * seq_along(nom); non = which(
				!is.na(d$Y[i, ]) & d$Y[i, ] == 0
			); Z[i, non] = pmin(Z[i, non], -0.1)
		}
		expect_true(nf_check(Z, d$Y, d$YL, d$odmax, d$odobs))
		expect_identical(ns$.lame_nom_feasible(Z, d$Y, d$YL, d$odmax, d$odobs, TRUE), Z)
		# breaking rows: ranks reversed, a non-nominee above a nominee, a nominee below zero, a row with room and a positive non-nominee
		Zb = Z; Zb[1, which(!is.na(d$Y[1, ]) & d$Y[1, ] > 0)] = rev(Zb[1, which(!is.na(d$Y[1, ]) & d$Y[1, ] > 0)]) - 5
		Zb[2, which(d$Y[2, ] == 0)[1]] = 9; Zb[3, which(d$Y[3, ] > 0)[1]] = -2
		expect_false(nf_check(Zb, d$Y, d$YL, d$odmax, d$odobs))
		Zr = ns$.lame_nom_feasible(Zb, d$Y, d$YL, d$odmax, d$odobs, TRUE)
		expect_true(nf_check(Zr, d$Y, d$YL, d$odmax, d$odobs))
		# rows that satisfied the constraints are untouched
		keep = setdiff(seq_len(nrow(Z)), 1:3)
		expect_identical(Zr[keep, ], Zb[keep, ])
	}
})

test_that('the sweep leaves a slice that starts outside the constraints inside them from the first draw, as it did not before', {
	ns = asNamespace('lame')
	for (fam in c('cbin', 'frn')) {
		d = nf_panel(fam, seed = 6L)
		set.seed(3)
		Z = matrix(rnorm(length(d$Y), -1), nrow(d$Y), ncol(d$Y))
		# the stuck state: the rank order of a row reversed and a non-nominee above the weakest nominee
		i = 4L; nom = which(d$Y[i, ] > 0); nom = nom[order(d$Y[i, nom])]
		Z[i, nom] = rev(c(0.03, 0.7, 3.1)[seq_along(nom)]); Z[i, which(d$Y[i, ] == 0)[1]] = 2.4
		EZ = matrix(-0.5, nrow(Z), ncol(Z))
		f = if (fam == 'frn') function(Z) ns$rZ_frn_bip_fc(
			Z,
			EZ,
			d$Y,
			d$YL,
			d$odmax,
			d$odobs
		) else function(Z) ns$rZ_cbin_bip_fc(Z, EZ, d$Y, d$odmax, d$odobs)
		for (s in 1:5) Z = f(Z)
		expect_true(nf_check(Z, d$Y, d$YL, d$odmax, d$odobs))
	}
})

test_that('the repair returns every state of the joint law of the latent array and the data unchanged, hidden cells and rooms included', {
	ns = asNamespace('lame')
	set.seed(21)
	for (bip in c(FALSE, TRUE)) for (frn in c(FALSE, TRUE)) {
		n = 8L; nc = if (bip) 6L else n; cap = 2
		moved = 0L
		for (r in 1:300) {
			Z = -0.3 + matrix(rnorm(n * nc), n, nc) + outer(rnorm(n, 0, .5), rnorm(nc, 0, .5), '+')
			el = matrix(TRUE, n, nc); if (!bip) diag(el) = FALSE
			if (r %% 3 == 0L) { el[n, ] = FALSE; el[, nc] = FALSE }
			Y = matrix(0, n, nc)
			for (i in 1:n) { z = Z[i, ]; z[!el[i, ]] = -Inf; pos = which(z > 0); k = min(cap, length(pos)); if (k) { top = pos[order(z[pos], decreasing = TRUE)[seq_len(k)]]; Y[i, top] = if (frn) rank(z[top]) else 1 } }
			Y[!el] = NA; Y[sample(which(el), round(0.3 * sum(el)))] = NA
			if (frn) for (i in 1:n) { k = which(!is.na(Y[i, ]) & Y[i, ] > 0); if (length(k)) Y[i, k] = rank(Y[i, k]) }
			odobs = rowSums(!is.na(Y) & Y > 0)
			YL = NULL
			if (frn) { YL = matrix(NA_integer_, n, cap); for (i in 1:n) for (q in 1:cap) { j = which(Y[i, ] == q); if (length(j)) YL[i, q] = j[1] } }
			if (!identical(ns$.lame_nom_feasible(Z, Y, YL, rep(cap, n), odobs, bip), Z)) moved = moved + 1L
		}
		expect_identical(moved, 0L)
	}
})
