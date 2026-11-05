# receiver-level move of the cbin and frn samplers in bipartite ame() fits (R/nomination_actor_move.R)
# the sweep itself is tested in test_nomination_receiver_move.R

test_that('bipartite ame() runs the receiver move for cbin and frn and records it', {
	skip_on_cran()
	set.seed(8801)
	nA = 14; nB = 9; cap = 3
	a = rnorm(nA, 0, .6); b = rnorm(nB, 0, .6)
	for (family in c('cbin', 'frn')) {
		z = -0.7 + outer(a, b, '+') + matrix(rnorm(nA * nB), nA, nB); z[, nB] = -Inf
		y = matrix(0, nA, nB)
		for (i in seq_len(nA)) {
			pos = which(z[i, ] > 0); k = min(cap, length(pos))
			if (k > 0) { top = pos[order(z[i, pos], decreasing = TRUE)[seq_len(k)]]; y[i, top] = if (family == 'cbin') 1 else rank(z[i, top]) }
		}
		y[sample(length(y), 8)] = NA
		if (family == 'frn') for (i in seq_len(nA)) { r = which(!is.na(y[i, ]) & y[i, ] > 0); if (length(r)) y[i, r] = rank(y[i, r]) }
		dimnames(y) = list(sprintf('a%02d', seq_len(nA)), sprintf('b%02d', seq_len(nB)))
		xd = array(rnorm(nA * nB), c(nA, nB, 1), dimnames = list(rownames(y), colnames(y), 'x1'))
		run = function(prior, ...) suppressMessages(suppressWarnings(ame(
			y,
			Xdyad = xd,
			family = family,
			odmax = cap,
			mode = 'bipartite',
			R = 0,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		)))
		on = run(list(actor_move_every = 2))
		expect_identical(
			names(on$actor_move),
			c('visits', 'col_calls', 'col_skips', 'level_calls', 'level_rows', 'ridge_calls')
		)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['col_calls'] + on$actor_move['col_skips']), 40 * nB)
		expect_gt(unname(on$actor_move['col_calls']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_null(run(list(actor_move = FALSE))$actor_move)
		# the default visits every sweep (20 burn-in and 60 production sweeps)
		expect_equal(unname(run(list())$actor_move['visits']), 80)
		# a schedule that never reaches a visit uses no random numbers, so its draws match the switched-off run
		never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_identical(never$a_samples, off$a_samples); expect_identical(never$b_samples, off$b_samples)
	}
})
