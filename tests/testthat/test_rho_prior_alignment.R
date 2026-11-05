# ame() and lame() share the dyadic correlation prior 2 Beta(2, 2) - 1, density proportional to 1 - rho^2
# (R/rho_prior.R): the shared definition, exactness of each ame() update of rho, prior-only runs in every family,
# and agreement of ame() with a single-period lame()

# paired differences of functions of rho before and after one update, as z values
rp_fun = function(rho) c(rho, rho^2, abs(rho), rho^3, rho^4)
rp_z = function(d) colMeans(d) / (apply(d, 2, sd) / sqrt(nrow(d)))

# n x n latent pairs: dyads (e_ij, e_ji) with variance s2 and correlation rho
rp_pairs = function(n, rho, s2 = 1) {
	up = which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	e1 = rnorm(nrow(up), 0, sqrt(s2))
	e2 = rho * e1 + sqrt(1 - rho^2) * rnorm(nrow(up), 0, sqrt(s2))
	z = matrix(0, n, n)
	z[up] = e1
	z[up[, 2:1, drop = FALSE]] = e2
	z
}

test_that('ame() and lame() read one prior of the dyadic correlation', {
	expect_identical(lame:::.lame_rho_prior_power(), 1)
	# the draw is 2 Beta(2, 2) - 1 and the log density is 1 - rho^2 up to a constant
	set.seed(5101)
	draws = replicate(40000, lame:::.lame_rho_prior_draw())
	expect_lt(abs(mean(draws)), 0.01)
	expect_lt(abs(var(draws) - 0.2), 0.005)
	expect_lt(abs(mean(draws^4) - 3 / 35), 0.003)
	expect_lt(mean(abs(draws) > 0.9), 0.02)
	r = c(-0.7, 0, 0.3, 0.95)
	expect_equal(diff(lame:::.lame_rho_prior_log(r)), diff(log(1 - r^2)))
	expect_equal(diff(lame:::.lame_rho_prior_log(r, power = -0.5)), diff(-0.5 * log(1 - r^2)))
	# the arc-sine and uniform priors of rrho_mh() are the exponents -1/2 and 0
	set.seed(5102)
	expect_equal(mean(abs(replicate(20000, lame:::.lame_rho_prior_draw(-0.5))) > 0.9), 0.295, tolerance = 0.05)
	# every update reads the same exponent: the nomination update by default, the compiled kernel through its argument
	expect_identical(eval(formals(lame:::.lame_nom_rho_collapsed)$pw), lame:::.lame_rho_prior_power())
	expect_null(formals(rrho_mh)$asp)
	# the collapsed binary updates of ame() and lame() propose from the same prior
	cx = lame:::.lame_binary_dynamic_collapsed_rho_context(array(NA_real_, c(3, 3, 1)))
	for (seed in 1:5) {
		set.seed(5103 + seed)
		a = lame:::.ame_binary_collapsed_rho_mh(0.3, array(0, c(3, 3, 1)), cx, tries = 3L)
		set.seed(5103 + seed)
		b = lame:::.lame_binary_dynamic_collapsed_rho_mh(0.3, array(0, c(3, 3, 1)), cx, tries = 3L)
		expect_identical(a, b)
	}
})

test_that('the Metropolis update of ame() leaves the joint law of rho and the latent pairs unchanged', {
	skip_on_cran()
	one_update = function(n, reps, seed, asp = NULL, prior_draw = lame:::.lame_rho_prior_draw) {
		set.seed(seed)
		d = NULL
		for (k in seq_len(reps)) {
			rho = prior_draw()
			s2 = 1.7
			z = rp_pairs(n, rho, s2)
			d = rbind(d, rp_fun(rrho_mh(z, rho, s2, asp = asp)) - rp_fun(rho))
		}
		rp_z(d)
	}
	for (n in 3:4) expect_lt(max(abs(one_update(n, 15000, 5200 + n))), 4.6)
	# the explicit arc-sine prior is exact under itself
	arcsine = function() sin(runif(1, -pi / 2, pi / 2))
	expect_lt(max(abs(one_update(4, 15000, 5210, asp = TRUE, prior_draw = arcsine))), 4.6)
	# power: the arc-sine kernel against the package prior, and the reverse, fail
	expect_gt(max(abs(one_update(4, 15000, 5211, asp = TRUE))), 8)
	expect_gt(max(abs(one_update(4, 15000, 5212, prior_draw = arcsine))), 8)
})

test_that('the compiled kernel and the binary update of lame() leave the joint law unchanged under the same prior', {
	skip_on_cran()
	set.seed(5251)
	d = NULL
	for (k in 1:12000) {
		rho = lame:::.lame_rho_prior_draw()
		s2 = 1.3
		z = array(c(rp_pairs(4L, rho, s2), rp_pairs(4L, rho, s2)), c(4, 4, 2))
		d = rbind(d, rp_fun(lame:::rrho_mh_rep_cpp(z, rho, s2, 0.1, lame:::.lame_rho_prior_power())) - rp_fun(rho))
	}
	expect_lt(max(abs(rp_z(d))), 4.6)
	# mutant: the arc-sine exponent against the package prior
	d = NULL
	for (k in 1:12000) {
		rho = lame:::.lame_rho_prior_draw()
		z = array(c(rp_pairs(4L, rho), rp_pairs(4L, rho)), c(4, 4, 2))
		d = rbind(d, rp_fun(lame:::rrho_mh_rep_cpp(z, rho, 1, 0.1, -0.5)) - rp_fun(rho))
	}
	expect_gt(max(abs(rp_z(d))), 8)
	set.seed(5252)
	d = NULL
	for (k in 1:10000) {
		rho = lame:::.lame_rho_prior_draw()
		eta = array(rnorm(16, 0, 0.8), c(4, 4, 1))
		Y = array(1 * (eta + array(rp_pairs(4L, rho), c(4, 4, 1)) > 0), c(4, 4, 1))
		Y[cbind(1:4, 1:4, 1)] = NA
		cx = lame:::.lame_binary_dynamic_collapsed_rho_context(Y)
		d = rbind(d, rp_fun(lame:::.lame_binary_dynamic_collapsed_rho_mh(rho, eta, cx)$rho) - rp_fun(rho))
	}
	expect_lt(max(abs(rp_z(d))), 4.6)
})

test_that('the binary update of ame() with the latent responses integrated out leaves the joint law unchanged', {
	skip_on_cran()
	one_update = function(n, reps, seed, tries = 4L, kernel = lame:::.ame_binary_collapsed_rho_mh) {
		set.seed(seed)
		d = NULL
		for (k in seq_len(reps)) {
			rho = lame:::.lame_rho_prior_draw()
			eta = array(rnorm(n * n, 0, 0.8), c(n, n, 1))
			Z = eta + array(rp_pairs(n, rho), c(n, n, 1))
			Y = array(1 * (Z > 0), c(n, n, 1))
			Y[cbind(1:n, 1:n, 1)] = NA
			cx = lame:::.lame_binary_dynamic_collapsed_rho_context(Y)
			d = rbind(d, rp_fun(kernel(rho, eta, cx, tries = tries)$rho) - rp_fun(rho))
		}
		rp_z(d)
	}
	expect_lt(max(abs(one_update(4, 12000, 5301))), 4.6)
	expect_lt(max(abs(one_update(3, 12000, 5302, tries = 1L))), 4.6)
	# power: proposals from the arc-sine prior fail
	arcsine_kernel = function(rho, eta, cx, tries = 4L) {
		old = lame:::.lame_binary_dynamic_collapsed_rho_loglik(rho, eta, cx)
		for (k in seq_len(tries)) {
			proposal = sin(runif(1L, -pi / 2, pi / 2))
			new = lame:::.lame_binary_dynamic_collapsed_rho_loglik(proposal, eta, cx)
			if (log(runif(1L)) < new - old) { rho = proposal; old = new }
		}
		list(rho = rho)
	}
	expect_gt(max(abs(one_update(4, 12000, 5303, kernel = arcsine_kernel))), 8)
})

test_that('ame() with next to no data returns the prior of the dyadic correlation, in every family', {
	skip_on_cran()
	n = 6L
	nm = sprintf('a%02d', 1:n)
	for (family in c('normal', 'binary', 'poisson', 'ordinal', 'cbin', 'frn')) {
		# two observed cells in different dyads carry no information on rho
		Y = matrix(NA_real_, n, n, dimnames = list(nm, nm))
		Y[1, 2] = switch(family, normal = 0.3, binary = 1, poisson = 2, ordinal = 2, 1)
		Y[3, 4] = switch(family, normal = -0.2, 0)
		fit = suppressMessages(suppressWarnings(ame(
			Y,
			family = family,
			R = 0,
			burn = 200,
			nscan = 3000,
			odens = 1,
			seed = 1,
			verbose = FALSE,
			gof = FALSE,
			odmax = if (family %in% c('cbin', 'frn')) 2 else NULL
		)))
		rho = fit$VC[, 'rho']
		# 2 Beta(2, 2) - 1 has variance 0.2 and little mass beyond 0.9; the arc-sine prior has variance 0.5
		expect_lt(abs(var(rho) - 0.2), 0.07)
		expect_lt(mean(abs(rho) > 0.9), 0.05)
		expect_lt(abs(mean(rho)), 0.2)
	}
})

test_that('the nomination update of ame() leaves the joint law unchanged under the package prior', {
	skip_on_cran()
	# with the exponent omitted the update reads the package prior
	n = 6L; nt = 2L
	set.seed(5401)
	draw = function(family) {
		rho = lame:::.lame_rho_prior_draw()
		a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
		Z = EZ = Y = array(NA_real_, c(n, n, nt))
		for (t in seq_len(nt)) {
			eta = -1.1 + outer(a, b, '+')
			E = rp_pairs(n, rho)
			diag(E) = rnorm(n, 0, sqrt(1 + rho))
			z = eta + E
			y = matrix(0, n, n)
			for (i in seq_len(n)) {
				zi = z[i, ]; zi[i] = -Inf
				pos = which(zi > 0); k = min(2, length(pos))
				if (k > 0) { top = pos[order(zi[pos], decreasing = TRUE)[seq_len(k)]]; y[i, top] = if (family == 'cbin') 1 else rank(zi[top]) }
			}
			diag(y) = NA
			Z[, , t] = z; EZ[, , t] = eta; Y[, , t] = y
		}
		YL = if (family == 'frn') lapply(seq_len(nt), function(t) {
			mx = max(c(0, Y[, , t]), na.rm = TRUE); m = matrix(NA_integer_, n, max(mx, 1))
			for (i in seq_len(n)) for (r in seq_len(mx)) { j = which(Y[i, , t] == r); if (length(j)) m[i, r] = j }
			m
		}) else NULL
		list(rho = rho, Z = Z, EZ = EZ, Y = Y, YL = YL)
	}
	for (family in c('cbin', 'frn')) {
		d = NULL
		for (r in 1:2000) {
			s = draw(family)
			o = lame:::.lame_nom_rho_collapsed(s$Z, s$EZ, s$rho, s$Y, family, s$YL, rep(2, n))
			d = rbind(d, rp_fun(o$rho) - rp_fun(s$rho))
		}
		expect_lt(max(abs(rp_z(d))), 4.6)
	}
})

test_that('ame() and a single-period lame() return the same posterior of rho', {
	skip_on_cran()
	n = 10L
	set.seed(5501)
	nm = sprintf('a%02d', 1:n)
	x = array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, 'x1'))
	a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
	E = rp_pairs(n, 0.5)
	diag(E) = rnorm(n)
	Z = -0.4 + 0.6 * x[, , 1] + outer(a, b, '+') + E
	for (family in c('binary', 'normal')) {
		Y = if (family == 'binary') 1 * (Z > 0) else Z
		diag(Y) = NA
		dimnames(Y) = list(nm, nm)
		means = sapply(c('ame', 'lame'), function(kind) sapply(1:4, function(s) {
			fit = suppressMessages(suppressWarnings(
				if (kind == 'ame') ame(
					Y,
					Xdyad = x,
					family = family,
					R = 0,
					burn = 500,
					nscan = 2000,
					odens = 5,
					seed = 100 + s,
					verbose = FALSE,
					gof = FALSE
				)
				else lame(
					list(Y),
					Xdyad = list(x),
					family = family,
					R = 0,
					burn = 500,
					nscan = 2000,
					odens = 5,
					seed = 100 + s,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE
				)
			))
			mean(fit$VC[, 'rho'])
		}))
		z = (mean(means[, 'ame']) - mean(means[, 'lame'])) / sqrt(var(means[, 'ame']) / 4 + var(means[, 'lame']) / 4)
		expect_lt(abs(z), 3.5)
	}
})
