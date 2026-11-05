test_that('actor-specific predictors retain alignment when composition changes', {
	skip_on_cran()
	set.seed(927531)
	for (mode in c('unipartite', 'bipartite')) for (kind in c('row', 'col')) {
		rows = list(c('a', 'b', 'c'), c('b', 'c', 'd', 'e'), letters[1:5])
		cols = if (mode == 'bipartite') list(c('u', 'v'), c('v', 'w', 'x'), c('u', 'v', 'w', 'x')) else rows
		Y = Map(
			function(r, c) {x = matrix(
				rnorm(length(r) * length(c)),
				length(r),
				dimnames = list(r, c)
			);if (mode == 'unipartite') diag(x) = NA;x},
			rows,
			cols
		)
		X = Map(
			function(r, c) array(rnorm(length(r) * length(c)), c(length(r), length(c), 1), dimnames = list(r, c, 'x')),
			rows,
			cols
		)
		f = lame(
			Y,
			Xdyad = X,
			R = 0,
			family = 'normal',
			mode = mode,
			dynamic_beta_per_actor = kind,
			burn = 20,
			nscan = 80,
			odens = 4,
			seed = 8317,
			verbose = FALSE,
			gof = FALSE
		)
		expected_names = sort(unique(unlist(if (kind == 'row') rows else cols)))
		expect_identical(rownames(f$theta_actor_mean), expected_names)
		expect_identical(dimnames(f$THETA_ACTOR)[[2]], expected_names)
		predicted = predict(f, newdata = X, type = 'link')
		for (t in seq_along(Y)) {
			r = rownames(predicted[[t]]);c = colnames(predicted[[t]])
			x = X[[t]][r, c, 1, drop = FALSE];x = matrix(x, length(r), length(c), dimnames = list(r, c))
			b = colMeans(f$BETA);th = f$theta_actor_mean[if (kind == 'row') r else c, t]
			expected = b[1] + b[2] * x + outer(f$APM[r], f$BPM[c], '+') + sweep(x, if (kind == 'row') 1 else 2, th, '*')
			expect_equal(predicted[[t]], expected, tolerance = 1e-10)
			stored = f$EZ[[t]][r, c, drop = FALSE];keep = is.finite(stored)
			expect_equal(stored[keep], expected[keep], tolerance = 1e-10)
		}
		expect_true(all(unlist(f$tryErrorChecks) == 0))
	}
})
