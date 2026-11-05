.actor_prediction_fixture = function(kind = 'row') {
	rn = c('a', 'b', 'c');cn = c('d', 'e', 'f', 'g');nt = 2L
	X = lapply(seq_len(nt), function(t) array(seq_len(12) * t / 10, c(3, 4, 1), dimnames = list(rn, cn, 'x')))
	design = lapply(
		X,
		function(x) {a = array(1, c(3, 4, 2));a[, , 2] = x[, , 1];dimnames(a) = list(
			rn,
			cn,
			c('intercept', 'x_dyad')
		);a}
	)
	theta = if (kind == 'row') rbind(c(-1, .5), c(0, -1), c(1, .5)) else rbind(
		c(-1, 1),
		c(0, 0),
		c(.5, -.5),
		c(.5, -.5)
	)
	f = list(
		family = 'normal',
		mode = 'bipartite',
		dynamic_beta_per_actor = kind,
		per_actor_model_version = 2L,
		per_actor_design_idx = 2L,
		theta_actor_mean = theta,
		Xlist = design,
		BETA = matrix(rep(c(.2, .4), each = 5), 5, 2, dimnames = list(NULL, c('intercept', 'x_dyad'))),
		APM = setNames(c(.1, -.2, .1), rn),
		BPM = setNames(c(.2, -.1, .1, -.2), cn),
		EZ = lapply(X, function(x) matrix(0, 3, 4, dimnames = list(rn, cn))),
		YPM = lapply(X, function(x) matrix(0, 3, 4, dimnames = list(rn, cn)))
	)
	class(f) = c('lame', 'ame');list(fit = f, X = X, theta = theta)
}

test_that('counterfactuals include actor slopes with aligned rows and columns', {
	for (kind in c('row', 'col')) {
		z = .actor_prediction_fixture(kind);f = z$fit;X = z$X
		expected = lapply(
			seq_along(X),
			function(t) .2 + .4 * X[[t]][, , 1] + outer(f$APM, f$BPM, '+') +
				sweep(X[[t]][, , 1], if (kind == 'row') 1 else 2, z$theta[, t], '*')
		)
		expect_equal(predict(f, newdata = X, type = 'link'), expected)
		changed = lapply(X, function(x) x[c(3, 1, 2), c(4, 2, 1, 3), , drop = FALSE])
		expect_equal(predict(f, newdata = changed, type = 'link'), expected)
		subset = lapply(X, function(x) x[c(1, 3), c(2, 4), , drop = FALSE])
		expect_equal(
			predict(f, newdata = subset, type = 'link'),
			lapply(expected, function(x) x[c(1, 3), c(2, 4), drop = FALSE])
		)
		old = f;old$per_actor_model_version = NULL
		expect_error(predict(old, newdata = X), 'corrected actor-specific')
		none = f;none$theta_actor_mean = NULL
		expect_error(predict(none, newdata = X), 'corrected actor-specific')
	}
})

test_that('combined streaming actor-slope summaries match pooled draws', {
	set.seed(963164);d1 = array(rnorm(5 * 3 * 2), c(5, 3, 2));d2 = array(rnorm(9 * 3 * 2) + 2, c(9, 3, 2))
	dimnames(d1) = list(NULL, c('a', 'b', 'c'), c('1', '2'));dimnames(d2) = dimnames(d1)
	make = function(d) {f = .actor_prediction_fixture()$fit;f$BETA = f$BETA[rep(1, dim(d)[1]), , drop = FALSE];f$THETA_ACTOR = d;f$theta_actor_mean = apply(d, c(2, 3), mean);f$theta_actor_sd = apply(d, c(2, 3), sd);f}
	a = make(d1);b = make(d2);full = combine_ame_chains(list(a, b), diagnostics = FALSE)
	a$THETA_ACTOR = b$THETA_ACTOR = NULL;stream = combine_ame_chains(list(a, b), diagnostics = FALSE)
	pooled = abind::abind(d1, d2, along = 1)
	expect_equal(full$theta_actor_mean, apply(pooled, c(2, 3), mean))
	expect_equal(full$theta_actor_sd, apply(pooled, c(2, 3), sd))
	expect_equal(stream$theta_actor_mean, full$theta_actor_mean)
	expect_equal(stream$theta_actor_sd, full$theta_actor_sd)
})
