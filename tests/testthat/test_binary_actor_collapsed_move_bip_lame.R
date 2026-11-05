skip_on_cran()

bip_panel = function(seed, nA = 8L, nB = 10L, nt = 3L, dens = .12) {
	set.seed(seed)
	rn = sprintf('r%02d', 1:nA); cn = sprintf('c%02d', 1:nB)
	Y = lapply(
		1:nt,
		function(t) { y = matrix(
			rbinom(nA * nB, 1, dens),
			nA,
			nB,
			dimnames = list(rn, cn)
		); y[1:2, ] = 0; y[, 1:3] = 0; y }
	)
	X = lapply(1:nt, function(t) array(rnorm(nA * nB), c(nA, nB, 1), dimnames = list(rn, cn, 'x')))
	list(Y = Y, X = X, nA = nA, nB = nB)
}
bip_run = function(p, prior, dyn, R, ...) suppressMessages(lame(
	p$Y,
	Xdyad = p$X,
	mode = 'bipartite',
	family = 'binary',
	R_row = R,
	R_col = R,
	dynamic_ab = dyn,
	burn = 20,
	nscan = 60,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = prior,
	...
))
# the update records its work, a visit interval beyond the run never visits (the draws are those of the switch off),
# and a run with visits differs from the switch off
bip_check = function(p, dyn, R) {
	expect_no_warning({ on = bip_run(p, list(actor_move_every = 2), dyn, R) })
	expect_equal(unname(on$actor_move['visits']), 40)
	expect_equal(unname(on$actor_move['level_calls']), 40 * (p$nA + p$nB))
	expect_equal(unname(on$actor_move['uv_calls']), if (R > 0) 40 * (p$nA + p$nB) else 0)
	expect_equal(unname(on$actor_move['stuck']), 0)
	# under dynamic effects the deviation paths of random actor pairs are drawn too, nA %/% 2 + nB %/% 2 pairs per visit
	if (dyn) expect_equal(
		unname(on$actor_move['dev_calls']),
		40 * (p$nA %/% 2 + p$nB %/% 2)
	) else expect_false('dev_calls' %in% names(on$actor_move))
	expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
	off = bip_run(p, list(actor_move = FALSE), dyn, R)
	expect_null(off$actor_move)
	never = bip_run(p, list(actor_move_every = 1000), dyn, R)
	expect_null(never$actor_move)
	expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
	expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
}

test_that('lame() runs the bipartite actor update as a prior switch and records its work', {
	p = bip_panel(7306L)
	bip_check(p, FALSE, 0L)
	bip_check(p, TRUE, 1L)
	expect_error(bip_run(p, list(actor_move = NA), FALSE, 0L), 'actor_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(
		bip_run(p, list(actor_move_every = 0), FALSE, 0L),
		'actor_move_every must be a positive whole number',
		fixed = TRUE
	)
})

test_that('the bipartite actor update covers static and dynamic effects, ranks 0 and 1, and moving positions', {
	skip_on_cran()
	p = bip_panel(7308L)
	bip_check(p, FALSE, 1L)
	bip_check(p, TRUE, 0L)
	# positions that move over time keep their own updates: only the levels are visited
	dy = bip_run(p, list(actor_move_every = 2), TRUE, 1L, dynamic_uv = TRUE)
	expect_equal(
		unname(dy$actor_move['level_calls']),
		40 * (p$nA + p$nB)
	); expect_equal(unname(dy$actor_move['uv_calls']), 0)
})

test_that('lame() runs the bipartite actor update on a panel with actors absent from a period', {
	skip_on_cran()
	p = bip_panel(7307L, dens = .15)
	p$Y[[2]] = p$Y[[2]][-3, ]; p$Y[[3]] = p$Y[[3]][, -4]
	p$X[[2]] = p$X[[2]][-3, , , drop = FALSE]; p$X[[3]] = p$X[[3]][, -4, , drop = FALSE]
	for (dyn in c(FALSE, TRUE)) {
		expect_no_warning({ on = bip_run(p, list(actor_move_every = 2), dyn, 1L) })
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['stuck']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		off = bip_run(p, list(actor_move = FALSE), dyn, 1L)
		expect_null(off$actor_move)
	}
})

test_that('lame() runs the bipartite actor update with dynamic coefficients and records its work', {
	skip_on_cran()
	p = bip_panel(7309L, nt = 5L)
	for (kind in c('rw1', 'ar1')) for (blocks in list(
		'intercept',
		c('intercept', 'x_dyad'),
		'x_dyad'
	)) for (dyn in c(FALSE, TRUE)) {
		on = bip_run(p, list(actor_move_every = 2), dyn, 1L, dynamic_beta = blocks, dynamic_beta_kind = kind)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['level_calls']), 40 * (p$nA + p$nB))
		expect_gt(unname(on$actor_move['uv_calls']), 0)
		expect_equal(unname(on$actor_move['stuck']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		off = bip_run(p, list(actor_move = FALSE), dyn, 1L, dynamic_beta = blocks, dynamic_beta_kind = kind)
		expect_null(off$actor_move)
		# a schedule that never reaches a visit consumes no random numbers: the draws are those of the switch-off
		never = bip_run(p, list(actor_move_every = 10^6), dyn, 1L, dynamic_beta = blocks, dynamic_beta_kind = kind)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	}
	# kinds and grids the update does not cover keep the sampler's own updates
	expect_null(
		bip_run(
			p,
			list(actor_move_every = 2),
			FALSE,
			0L,
			dynamic_beta = 'intercept',
			dynamic_beta_kind = 'rw2'
		)$actor_move
	)
	expect_null(
		bip_run(
			p,
			list(actor_move_every = 2),
			FALSE,
			0L,
			dynamic_beta = 'intercept',
			time_index = c(1, 2, 4, 5, 6)
		)$actor_move
	)
})
