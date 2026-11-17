skip_on_cran()

library(lame)

# the fast route returns a panel's actors as the sampler does, sorted by name, and fits the panel exactly as it
# fits the same panel sorted by hand

ao_panel = function(kind = "directed", n = 12, tt = 3, seed = 1, x_named = TRUE) {
	set.seed(seed)
	bip = kind == "bipartite"
	nc = if (bip) 9 else n
	rn = sample(sprintf("r%02d", seq_len(n)))
	cn = if (bip) sample(sprintf("c%02d", seq_len(nc))) else rn
	xr = matrix(rnorm(n), n, 1, dimnames = list(if (x_named) rn, "xr"))
	Y = list(); Xd = list()
	for (t in seq_len(tt)) {
		x = matrix(rnorm(n * nc), n, nc)
		y = 0.2 + 0.6 * x + 0.4 * xr[, 1] + matrix(rnorm(n * nc), n, nc)
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		Xd[[t]] = array(x, c(n, nc, 1), dimnames = if (x_named) list(rn, cn, "x") else list(NULL, NULL, "x"))
	}
	names(Y) = names(Xd) = paste0("t", seq_len(tt))
	list(Y = Y, Xd = Xd, Xr = rep(list(xr), tt), rn = rn, cn = cn, bip = bip)
}
ao_by_hand = function(d) {
	ro = order(d$rn); co = order(d$cn)
	d$Y = lapply(d$Y, function(y) y[ro, co])
	d$Xd = lapply(d$Xd, function(x) x[ro, co, , drop = FALSE])
	d$Xr = lapply(d$Xr, function(x) x[ro, , drop = FALSE])
	d
}
ao_fit = function(d, ...) suppressWarnings(suppressMessages(lame(
	d$Y, Xdyad = d$Xd, Xrow = d$Xr, R = 1, mode = if (d$bip) "bipartite" else "unipartite", verbose = FALSE, ...
)))
ao_same = function(f, g, lab) {
	expect_identical(coef(f), coef(g), label = lab)
	expect_identical(fitted(f), fitted(g), label = lab)
	expect_identical(f$a, g$a, label = lab)
	expect_identical(f$b, g$b, label = lab)
	expect_identical(f$U, g$U, label = lab)
}

test_that("a panel with unsorted actor names is fitted as the same panel sorted by hand", {
	for (kind in c("directed", "bipartite")) {
		for (x_named in c(TRUE, FALSE)) {
			d = ao_panel(kind, seed = 2, x_named = x_named)
			for (dyn in c(FALSE, TRUE)) {
				lab = paste(kind, x_named, dyn)
				f = ao_fit(d, method = "als", dynamic_uv = dyn)
				ao_same(f, ao_fit(ao_by_hand(d), method = "als", dynamic_uv = dyn), lab)
				expect_identical(rownames(fitted(f)[[1]]), sort(d$rn), label = lab)
				expect_identical(colnames(fitted(f)[[1]]), sort(d$cn), label = lab)
			}
		}
	}
})

test_that("the fast route and the sampler return the same actor order", {
	d = ao_panel("directed", seed = 3)
	f = ao_fit(d, method = "als")
	m = ao_fit(d, burn = 20, nscan = 40, odens = 4, chains = 1, gof = FALSE)
	expect_identical(dimnames(fitted(f)[[1]]), dimnames(fitted(m)[[1]]))
	expect_identical(names(f$a), names(m$APM))
	expect_identical(rownames(f$U), rownames(m$U))
	# actors that enter and leave
	d$Y[[2]] = d$Y[[2]][-3, -3]
	d$Xd[[2]] = d$Xd[[2]][-3, -3, , drop = FALSE]
	f = suppressWarnings(lame_als(d$Y, Xdyad = d$Xd, R = 1, verbose = FALSE))
	expect_identical(rownames(fitted(f)[[1]]), sort(d$rn))
})

test_that("a covariate named with the actors in another order is matched by name", {
	d = ao_panel("directed", seed = 4)
	f = ao_fit(d, method = "als")
	k = rev(seq_along(d$rn))
	d$Xr = lapply(d$Xr, function(x) x[k, , drop = FALSE])
	d$Xd = lapply(d$Xd, function(x) x[k, k, , drop = FALSE])
	expect_equal(coef(ao_fit(d, method = "als")), coef(f))
})

test_that("a cross-section keeps the order of the data, as ame() does", {
	d = ao_panel("directed", seed = 5)
	f = suppressWarnings(ame_als(d$Y[[1]], Xdyad = d$Xd[[1]], R = 1, verbose = FALSE))
	expect_identical(rownames(fitted(f)), d$rn)
})

test_that("a refit or a bootstrap of a fit made from unsorted data is that of the sorted panel", {
	d = ao_panel("directed", seed = 5, x_named = FALSE)
	h = ao_by_hand(d)
	for (dyn in c(FALSE, TRUE)) {
		lab = paste("refit", dyn)
		f = ao_fit(d, method = "als", dynamic_uv = dyn)
		rf = suppressWarnings(ame_als_refit(f))
		ao_same(rf, suppressWarnings(ame_als_refit(ao_fit(h, method = "als", dynamic_uv = dyn))), lab)
		# the data the fit holds lead back to the fit
		expect_equal(unlist(fitted(rf)), unlist(fitted(f)), tolerance = 1e-3, label = lab)
	}
	bf = suppressWarnings(ame_als_bootstrap(ao_fit(d, method = "als"), R = 10, seed = 1, verbose = FALSE))
	bh = suppressWarnings(ame_als_bootstrap(ao_fit(h, method = "als"), R = 10, seed = 1, verbose = FALSE))
	expect_identical(bf$n_valid, 10L)
	expect_identical(bf$coefs, bh$coefs)
	expect_identical(bf$a_coefs, bh$a_coefs)
})

test_that("a panel given as an array is sorted like one given as a list", {
	d = ao_panel("directed", seed = 6)
	f = suppressWarnings(lame(simplify2array(d$Y), R = 1, method = "als", verbose = FALSE))
	ao_same(f, suppressWarnings(lame(d$Y, R = 1, method = "als", verbose = FALSE)), "array")
	expect_identical(rownames(fitted(f)[[1]]), sort(d$rn))
})
