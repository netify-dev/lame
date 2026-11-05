# the translations of unipartite ame() fits (Poisson and ordinal), run from the hook in R/ame_unipartite.R
skip_on_cran()

ltmu_panel = function(family, sym, n = 8L) {
	set.seed(5)
	rn = sprintf('r%02d', 1:n)
	y = matrix(if (family == 'poisson') rpois(n * n, .4) else sample(1:3, n * n, TRUE), n, n, dimnames = list(rn, rn))
	if (family == 'poisson') { y[1:2, ] = 0; y[, 1:2] = 0 } else { y[1:2, ] = 1; y[, 1:2] = 1 }
	if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
	diag(y) = NA
	x = matrix(rnorm(n * n), n, n, dimnames = list(rn, rn))
	if (sym) x[lower.tri(x)] = t(x)[lower.tri(x)]
	list(Y = y, X = x, n = n)
}
ltmu_run = function(p, family, sym, prior = list(), R = 0L, burn = 10) suppressMessages(suppressWarnings(ame(
	p$Y,
	Xdyad = p$X,
	family = family,
	symmetric = sym,
	R = R,
	burn = burn,
	nscan = 40,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = modifyList(list(actor_move_every = 2), prior)
)))

test_that('unipartite ame() records the translations of Poisson and ordinal fits and honours the switches', {
	for (family in c('poisson', 'ordinal')) for (sym in c(FALSE, TRUE)) {
		p = ltmu_panel(family, sym)
		on = ltmu_run(p, family, sym); off = ltmu_run(p, family, sym, list(actor_move = FALSE)); never = ltmu_run(
			p,
			family,
			sym,
			list(actor_move_every = 1000)
		)
		am = on$actor_move
		expect_equal(unname(am['visits']), 25)
		expect_equal(unname(am['stuck']), 0)
		expect_equal(unname(am['level_calls'] + am['skipped']), 25 * (if (sym) p$n else 2 * p$n))
		if (family == 'poisson') {
			expect_equal(unname(am['coef_calls']), 25 * ncol(on$BETA))
			expect_gt(unname(am['scale_attempts']), 0)
			expect_equal(unname(am['rho_attempts'] > 0), !sym)
		} else expect_equal(unname(am['location_calls']), 25)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_null(off$actor_move); expect_null(never$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	}
})

test_that('unipartite ame() moves the static positions of a Poisson fit with the latent array', {
	for (sym in c(FALSE, TRUE)) {
		p = ltmu_panel('poisson', sym)
		pos = ltmu_run(p, 'poisson', sym, R = 1L, burn = 20)
		expect_gt(unname(pos$actor_move['uv_calls']), 0); expect_equal(unname(pos$actor_move['uv_stuck']), 0)
		expect_true(all(is.finite(pos$BETA)) && all(is.finite(pos$VC)))
	}
})
