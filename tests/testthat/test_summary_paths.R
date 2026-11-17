# summary, print, coef and confint of coefficient paths: per period, never a pooled z

skip_on_cran()

smp_panel = function(n = 10L, N = 4L, seed = 9) {
	set.seed(seed)
	ids = sprintf('a%02d', seq_len(n))
	Y = X = vector('list', N)
	for (t in seq_len(N)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(ids, ids, 'x1'))
		y = 1 * (-0.6 + 0.4 * t + 0.4 * x[, , 1] + matrix(rnorm(n * n), n, n) > 0)
		diag(y) = NA
		dimnames(y) = list(ids, ids)
		Y[[t]] = y
		X[[t]] = x
	}
	names(Y) = names(X) = paste0('t', seq_len(N))
	list(Y = Y, X = X)
}
smp_fit = function(d, ...) suppressWarnings(suppressMessages(lame(
	d$Y,
	Xdyad = d$X,
	R = 0,
	family = 'binary',
	burn = 20,
	nscan = 100,
	odens = 2,
	chains = 1,
	verbose = FALSE,
	gof = FALSE,
	plot = FALSE,
	prior = list(chain_extensions = 0),
	...
)))

test_that('a path has no pooled z or p-value in summary, and its periods are reported', {
	d = smp_panel()
	cases = list(
		list(dynamic_ab = TRUE),
		list(dynamic_beta = 'intercept', dynamic_beta_kind = 'ar1'),
		list(dynamic_beta = c('intercept', 'dyad'), dynamic_beta_kind = 'rw1'),
		list(dynamic_beta = 'dyad', dynamic_beta_kind = 'rw2'),
		list(dynamic_beta = 'intercept', dynamic_beta_kind = 'matern32')
	)
	for (cs in cases) {
		fit = do.call(smp_fit, c(list(d), cs))
		s = summary(fit)
		path = s$beta_path
		expect_identical(names(path), colnames(fit$BETA))
		expect_true(any(path))
		expect_true(all(is.na(s$beta[path, c('z_value', 'p_value')])))
		# the row of a path is its mean level over the periods
		expect_equal(unname(s$beta[, 'Estimate']), unname(apply(fit$BETA, 2, mean)), tolerance = 1e-12)
		# a static coefficient keeps its z and p-value
		if (any(!path)) expect_true(all(is.finite(s$beta[!path, 'z_value'])))
		# per period: the means of coef() and the intervals of confint()
		cf = coef(fit)
		expect_equal(s$beta_dynamic_per_t, cf, tolerance = 1e-12)
		for (k in names(path)[path]) {
			ci = confint(fit, parm = k)
			expect_identical(rownames(ci), paste0(k, '[', colnames(cf), ']'))
			expect_equal(unname(ci[, 1]), unname(s$beta_dynamic_per_t_lower[k, ]), tolerance = 1e-12)
			expect_equal(unname(ci[, 2]), unname(s$beta_dynamic_per_t_upper[k, ]), tolerance = 1e-12)
		}
		out = paste(capture.output(print(s)), collapse = '\n')
		for (k in names(path)[path]) {
			expect_match(out, paste0(k, ' by period'), fixed = TRUE)
		}
		expect_match(out, 'coefficient path', fixed = TRUE)
		# no significance star on a path row
		tab_line = grep(paste0('^', names(path)[path][1], ' '), capture.output(print(s)), value = TRUE)[1]
		expect_false(grepl('\\*', tab_line))
	}
	# one period of a path is still selected by its own name
	fit = smp_fit(d, dynamic_ab = TRUE)
	expect_identical(rownames(confint(fit, parm = 'intercept[t2]')), 'intercept[t2]')
	expect_error(confint(fit, parm = 'nope'), 'Unknown parameter')
})

test_that('a static fit keeps its z and p-values', {
	fit = smp_fit(smp_panel())
	s = summary(fit)
	expect_null(s$beta_path)
	expect_true(all(is.finite(s$beta[, 'z_value'])))
})
