# two chains of `m` independent stationary ar(1) series as the gate sees them (unit stationary sd),
# with chain 2 shifted by `shift` sds on the first series
ar1_chains = function(n, m, phi, shift = 0) {
	one = function() {
		x = matrix(0, n, m)
		x[1, ] = rnorm(m)
		for (i in 2:n) x[i, ] = phi * x[i - 1, ] + sqrt(1 - phi^2) * rnorm(m)
		colnames(x) = paste0("s", seq_len(m))
		x
	}
	a = one()
	b = one()
	b[, 1] = b[, 1] + shift
	lapply(list(a, b), function(s) list(
		BETA = cbind(intercept = s[, 1]),
		VC = cbind(va = 1 + 0.1 * s[, 2], rho = 0.1 * s[, 3]),
		chain_summary = s[, -(1:3), drop = FALSE],
		family = "normal",
		mode = "unipartite"
	))
}

# verdict of the gate on draws-by-chains matrices, one per series
verdict_of = function(...) {
	xs = list(...)
	tab = lame:::.lame_gate_verdict(vapply(xs, lame:::.lame_rhat_parts, numeric(5)))
	tab$parameter = names(xs)
	tab
}

test_that("the power of the threshold is the noncentral F tail and needs e effective draws", {
	a = 1 - 0.99^(1 / 60)
	thr = max(1.05, lame:::.lame_rhat_bound(a, 200, 4, 40))
	f = 39 * (thr^2 * 200 / 199 - 1)
	expect_equal(lame:::.lame_rhat_power(thr, 200, 4, 40), pf(f, 3, 156, ncp = 40, lower.tail = FALSE))
	# at ess_needed a one-sd separation is flagged with probability 0.95
	e = lame:::.lame_rhat_ess_needed(a, 200, 4)
	p = lame:::.lame_rhat_power(max(1.05, lame:::.lame_rhat_bound(a, 200, 4, e)), 200, 4, e)
	expect_equal(p, 0.95, tolerance = 1e-3)
	# more series to test or a stricter level need more effective draws
	expect_gt(lame:::.lame_rhat_ess_needed(1 - 0.99^(1 / 200), 200, 4), e)
	expect_lt(lame:::.lame_rhat_ess_needed(1 - 0.99^(1 / 2), 200, 4), e)
	# 30 series: about 38 effective draws per split chain, whatever the run length
	expect_equal(e, 37.6, tolerance = 0.01)
	expect_equal(lame:::.lame_rhat_ess_needed(a, 50, 4), e, tolerance = 0.01)
	# the threshold at that point is about 1.09, the detection level a one-sd separation clears
	expect_gt(max(1.05, lame:::.lame_rhat_bound(a, 200, 4, e)), 1.08)
	expect_lt(max(1.05, lame:::.lame_rhat_bound(a, 200, 4, e)), 1.10)
})

test_that("the threshold falls and the power rises with the effective draws and the run length", {
	skip_on_cran()
	for (m in c(1, 30, 1000)) {
		a = 1 - 0.99^(1 / (2 * m))
		for (n in c(30, 50, 200, 1000)) {
			e = seq(2, n, length.out = 500)
			thr = pmax(1.05, lame:::.lame_rhat_bound(a, n, 4, e))
			expect_true(all(diff(thr) <= 1e-12))
			expect_true(all(diff(lame:::.lame_rhat_power(thr, n, 4, e)) >= -1e-12))
		}
		# at a fixed autocorrelation time, longer runs
		for (tau in c(1, 4, 16)) {
			n = seq(30, 2000, by = 10)
			n = n[n / tau >= 2]
			thr = pmax(1.05, lame:::.lame_rhat_bound(a, n, 4, n / tau))
			expect_true(all(diff(thr) <= 1e-12))
			expect_true(all(diff(lame:::.lame_rhat_power(thr, n, 4, n / tau)) >= -1e-12))
		}
	}
	# a series is compared from ess_needed effective draws on (two series at the level of 30)
	a = 1 - 0.99^(1 / 60)
	e = lame:::.lame_rhat_ess_needed(a, 200, 4)
	parts = rbind(rhat = 1, ess = c(e - 1, e + 1), ess_chain = c(e - 1, e + 1), n = 200, chains = 4)
	expect_identical(lame:::.lame_gate_verdict(parts, alpha = 1 - 0.99^(4 / 60))$informative, c(FALSE, TRUE))
})

test_that("chains constant at one value pass and chains constant at different values fail", {
	set.seed(31)
	v = verdict_of(
		same = matrix(1, 100, 2),
		apart = cbind(rep(1, 100), rep(2, 100)),
		stuck = cbind(rep(5, 100), rnorm(100)),
		iid = matrix(rnorm(200), 100, 2)
	)
	expect_identical(lame:::.lame_rhat_chains(matrix(1, 100, 2)), NA_real_)
	expect_identical(lame:::.lame_rhat_chains(cbind(rep(1, 100), rep(2, 100))), Inf)
	bad = lame:::.lame_gate_flagged(v)
	expect_setequal(bad$parameter, c("apart", "stuck"))
	txt = lame:::.lame_gate_text(bad)
	expect_match(txt[bad$parameter == "apart"], "constant at different values", fixed = TRUE)
	expect_match(txt[bad$parameter == "stuck"], "R-hat [0-9.]+ above 1.0[0-9]+")
	# a non-finite series is not compared
	x = matrix(rnorm(200), 100, 2)
	x[3, 1] = NaN
	expect_identical(lame:::.lame_rhat_chains(x), NA_real_)
})

test_that("too few effective draws are flagged as such, with the nscan that would hold them", {
	set.seed(32)
	ar = function(n, phi) as.numeric(stats::filter(rnorm(n), phi, method = "recursive"))
	v = verdict_of(slow = cbind(ar(100, 0.95), ar(100, 0.95)), fast = matrix(rnorm(200), 100, 2))
	expect_identical(v$informative, c(FALSE, TRUE))
	bad = lame:::.lame_gate_flagged(v)
	expect_identical(bad$parameter, "slow")
	expect_match(lame:::.lame_gate_text(bad), "slow \\([0-9]+ effective draws per split chain, [0-9]+ needed to compare the chains\\)")
	# a quarter more than the effective draws needed, at the same autocorrelation
	nscan = lame:::.lame_gate_nscan(bad, 2000)
	expect_gte(nscan, 2000 * 1.25 * bad$ess_needed / bad$ess)
	expect_lt(nscan, 1.1 * 2000 * 1.25 * bad$ess_needed / bad$ess)
	expect_identical(lame:::.lame_gate_nscan(bad[0, ], 2000), 2000)
	# the block length is a multiple of odens, as lame() and ame() require
	for (od in c(7L, 30L, 25L)) {
		ns = lame:::.lame_gate_nscan(bad, 2100L, od)
		expect_equal(ns %% od, 0)
		expect_gte(ns, lame:::.lame_gate_nscan(bad, 2100L))
	}
	# a table without thresholds is read at 1.05
	plain = data.frame(parameter = c("a", "b"), rhat = c(1.04, 1.06))
	expect_identical(lame:::.lame_gate_flagged(plain)$parameter, "b")
	expect_identical(lame:::.lame_gate_text(lame:::.lame_gate_flagged(plain)), "b (R-hat 1.060 above 1.050)")
})

test_that("the drift threshold grows with the number of series tested and starts at 3", {
	expect_equal(lame:::.lame_chain_z_bound(1), 3)
	expect_gt(lame:::.lame_chain_z_bound(40), lame:::.lame_chain_z_bound(10))
	a = lame:::.lame_chain_alpha
	expect_equal(lame:::.lame_chain_z_bound(400), qnorm(1 - (1 - (1 - a)^(1 / 400)) / 2))
})

test_that("the R-hat bound is the null quantile of split R-hat for independent draws", {
	skip_on_cran()
	set.seed(21)
	# four split chains of 50 independent draws; R-hat is the larger of two tests (bulk, tail)
	r = replicate(4000, lame:::.lame_rhat(matrix(rnorm(200), 100, 2)))
	b = lame:::.lame_rhat_bound(1 - 0.95^(1 / 2), 50, 4, 50)
	expect_lt(abs(mean(r > b) - 0.05), 0.015)
	expect_gt(lame:::.lame_rhat_bound(0.05, 50, 4, 10), b)
	expect_lt(lame:::.lame_rhat_bound(0.05, 200, 4, 200), b)
})

test_that("settled chains are compared and pass at 400 draws while they mix, and are too slow to compare after", {
	skip_on_cran()
	set.seed(22)
	# 30 monitored series at one autocorrelation; per fit: any series flagged, any R-hat above threshold
	rate = function(n, phi, reps) {
		r = replicate(reps, {
			tab = lame:::.lame_gate_table(ar1_chains(n, 30, phi))
			thr = pmax(1.05, tab$threshold)
			c(any = nrow(lame:::.lame_gate_flagged(tab)) > 0L, over = any(tab$rhat > thr), few = mean(!tab$informative))
		})
		rowMeans(r)
	}
	r3 = rate(400, 0.3, 60)
	expect_lte(r3[["any"]], 0.05)
	# at phi 0.5 (67 effective draws per split chain) a series is rarely estimated below the 38 needed
	r5 = rate(400, 0.5, 60)
	expect_lte(r5[["few"]], 0.01)
	expect_lte(r5[["over"]], 0.05)
	r9 = rate(400, 0.9, 20)
	# a series at phi 0.9 holds about 11 effective draws per split chain of 38 needed
	expect_gt(r9[["few"]], 0.9)
	expect_lte(r9[["over"]], 0.1)
	# chance alone never flags R-hat above its threshold more than the family-wise level allows
	r7 = rate(100, 0.7, 30)
	expect_lte(r7[["over"]], 0.1)
	expect_gt(r7[["few"]], 0.99)
})

test_that("a one-sd separation on one series is flagged whenever the chains can be compared", {
	skip_on_cran()
	set.seed(23)
	for (phi in c(0.5, 0.7, 0.9)) {
		hit = replicate(30, {
			tab = lame:::.lame_gate_table(ar1_chains(400, 30, phi, shift = 1))
			row = tab[tab$parameter == "intercept", ]
			c(informative = row$informative, over = row$rhat > row$threshold, flagged = "intercept" %in% lame:::.lame_gate_flagged(tab)$parameter)
		})
		expect_true(all(hit["flagged", ] == 1))
		if (any(hit["informative", ] == 1)) expect_gte(mean(hit["over", hit["informative", ] == 1]), 0.9)
	}
	# two modes of a sign, as for mirrored latent positions
	fits = ar1_chains(100, 40, 0.3)
	fits[[1]]$chain_summary[, 1] = -0.88 + 0.02 * fits[[1]]$chain_summary[, 1]
	fits[[2]]$chain_summary[, 1] = 0.89 + 0.02 * fits[[2]]$chain_summary[, 1]
	tab = lame:::.lame_gate_table(fits)
	expect_true("s4" %in% lame:::.lame_gate_flagged(tab)$parameter)
	expect_true(all(tab$threshold >= 1.05, na.rm = TRUE))
	# one sticky chain beside a mixing one: the pooled effective draws are too few to compare
	set.seed(24)
	sticky = replicate(20, {
		fits = ar1_chains(400, 30, 0.3)
		x = numeric(400)
		for (i in 2:400) x[i] = 0.99 * x[i - 1] + sqrt(1 - 0.99^2) * rnorm(1)
		fits[[1]]$chain_summary[, 1] = x
		"s4" %in% lame:::.lame_gate_flagged(lame:::.lame_gate_table(fits))$parameter
	})
	expect_true(all(sticky))
})
