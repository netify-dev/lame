skip_on_cran()

.gaussian_score_fixture = function(symmetric, rank, missing, nr = 12L, nc = nr) {
	nt = 3L;U = matrix(rnorm(nr * rank, sd = .4), nr, rank)
	V = if (symmetric) U else matrix(rnorm(nc * rank, sd = .4), nc, rank)
	a = rnorm(nr, sd = .3);b = if (symmetric) a else rnorm(nc, sd = .3)
	signs = rep(c(1, -1), length.out = rank)
	O = if (symmetric) tcrossprod(sweep(U, 2, signs, '*'), U) else tcrossprod(U, V)
	Y = array(0, c(nr, nc, nt))
	for (tt in seq_len(nt)) {
		noise = matrix(rnorm(nr * nc, sd = .4), nr, nc)
		if (symmetric) noise = (noise + t(noise)) / 2
		z = .2 + outer(a, b, '+') + O + noise
		if (nr == nc) diag(z) = NA_real_
		if (missing) {
			if (symmetric) {ij = which(upper.tri(z), arr.ind = TRUE);ij = ij[sample(
				nrow(ij),
				floor(nrow(ij) * .3)
			), , drop = FALSE];z[ij] = z[ij[, 2:1]] = NA_real_}
			else z[sample(which(is.finite(z)), floor(sum(is.finite(z)) * .3))] = NA_real_
		}
		Y[, , tt] = z
	}
	list(Z = Y, X = array(numeric(), c(nr, nc, 0L, nt)))
}

test_that('Gaussian convergence requires a score certificate on the fitted objective', {
	skip_on_cran()
	set.seed(911900L);refinements = 0L
	for (symmetric in c(TRUE, FALSE)) for (rank in 0:3) for (missing in c(FALSE, TRUE)) {
		dat = .gaussian_score_fixture(symmetric, rank, missing)
		old = .ae_bcd_fit(dat$Z, dat$X, rank, symmetric, 200L, 1e-6)
		fit = .ae_bcd_fit(dat$Z, dat$X, rank, symmetric, 200L, 1e-6, certify = TRUE)
		expect_true(fit$converged)
		expect_true(fit$gaussian_diagnostics$pass)
		expect_lte(fit$gaussian_diagnostics$max_score, 1e-3)
		expect_lte(fit$sse, old$sse + 1e-9)
		eta = array(rep(fit$mu + outer(fit$a, fit$b, '+') + fit$Omat, 3L), dim(dat$Z))
		expect_equal(sum((dat$Z - eta)^2, na.rm = TRUE), fit$sse, tolerance = 1e-9)
		expect_equal(tail(fit$dev_history, 1), fit$sse)
		if (!length(fit$gaussian_diagnostics$refinement)) {
			copy = fit;copy$gaussian_diagnostics = NULL;expect_identical(copy, old)
		} else refinements = refinements + 1L
	}
	expect_gt(refinements, 0L)
})

test_that('centered Gaussian gradients match independent finite differences', {
	set.seed(912001L)
	for (symmetric in c(FALSE, TRUE)) {
		dat = .gaussian_score_fixture(symmetric, 2L, TRUE)
		dat$X = array(rnorm(12 * 12 * 2 * 3), c(12, 12, 2, 3))
		if (symmetric) for (k in 1:2) for (tt in 1:3) dat$X[, , k, tt] = (dat$X[, , k, tt] + t(dat$X[, , k, tt])) / 2
		f = .ae_bcd_fit(dat$Z, dat$X, 2L, symmetric, 200L, 1e-6)
		cx = .ae_gaussian_context(dat$Z, dat$X, 2L, symmetric)
		if (symmetric) cx$sign = ifelse(f$L < 0, -1, 1)
		theta = .ae_gaussian_encode(f, cx)
		expect_equal(.ae_gaussian_eval(theta, cx)$value, f$sse, tolerance = 1e-9)
		theta = theta + rnorm(length(theta), sd = .04);ev = .ae_gaussian_eval(theta, cx)
		ng = vapply(
			seq_along(theta),
			function(k) {a = b = theta;a[k] = a[k] + 1e-6;b[k] = b[k] - 1e-6;(.ae_gaussian_eval(a, cx)$value - .ae_gaussian_eval(b, cx)$value) / 2e-6},
			numeric(1)
		)
		expect_lt(max(abs(ng - ev$gradient) / pmax(1, abs(ng))), 2e-7)
	}
})

test_that('complete bipartite Gaussian fits match an independent SVD solution', {
	set.seed(912002L)
	for (rank in 0:2) {
		dat = .gaussian_score_fixture(FALSE, rank, FALSE, 12L, 10L)
		f = .ae_bcd_fit(dat$Z, dat$X, rank, FALSE, 200L, 1e-6, certify = TRUE)
		M = apply(dat$Z, c(1, 2), mean);S = sweep(sweep(M, 1, rowMeans(M)), 2, colMeans(M)) + mean(M)
		sv = svd(S);O = if (rank) sv$u[, seq_len(rank), drop = FALSE] %*% (sv$d[seq_len(rank)] * t(
			sv$v[, seq_len(rank), drop = FALSE]
		)) else matrix(0, 12, 10)
		eta = outer(rowMeans(M), colMeans(M), '+') - mean(M) + O
		reference = sum((dat$Z - array(rep(eta, 3L), dim(dat$Z)))^2)
		expect_true(f$converged);expect_equal(f$sse, reference, tolerance = 1e-6)
	}
})

test_that('public Gaussian fits and warm refits retain their convergence evidence', {
	set.seed(912003L);d = .gaussian_score_fixture(FALSE, 2L, TRUE)
	f = lame_als(d$Z, R = 2L, family = 'normal', bias_correct = 'none', multistart = 'cheap', seed = 3, verbose = FALSE)
	expect_true(f$converged);expect_true(f$gaussian_diagnostics$pass)
	expect_length(f$multistart_gaussian_diagnostics, 4L)
	expect_true(all(vapply(f$multistart_gaussian_diagnostics, `[[`, logical(1), 'pass')))
	g = ame_als_refit(f, verbose = FALSE)
	expect_true(g$converged);expect_true(g$gaussian_diagnostics$pass)
	one = ame_als_refit(f, max_iter = 1L, verbose = FALSE)
	expect_false(one$converged)
	expect_error(.ae_bcd_fit(d$Z, d$X, 2L, FALSE, 2L, 1e-6, weights = 1 + 0 * d$Z, certify = TRUE), 'unweighted')
})


test_that('rank-zero public refinement preserves absent factors and the least-squares solution', {
	Y = array(NA_real_, c(3L, 3L, 3L))
	index = cbind(c(1L, 2L, 1L, 3L, 2L, 3L), c(1L, 1L, 2L, 2L, 3L, 3L), rep(1:3, each = 2))
	Y[index] = c(1, 0, 2, 1, 0, 3)
	X = array(0, c(3L, 3L, 1L, 3L))
	X[cbind(index[, 1:2], 1L, index[, 3])] = seq(.1, .6, .1)
	refinements = 0L
	for (with_x in c(FALSE, TRUE)) {
		args = list(
			Y = Y,
			R = 0L,
			mode = 'bipartite',
			family = 'normal',
			bias_correct = 'none',
			max_iter = 20L,
			verbose = FALSE
		)
		if (with_x) args$Xdyad = X
		fit = do.call(lame_als, args)
		expect_true(fit$converged)
		expect_true(fit$gaussian_diagnostics$pass)
		refinements = refinements + length(fit$gaussian_diagnostics$refinement)
		expect_null(fit$U)
		expect_null(fit[['V']])
		expect_null(fit$L)
		row = factor(index[, 1]); col = factor(index[, 2]); x = seq(.1, .6, .1)
		design = if (with_x) model.matrix(~ row + col + x) else model.matrix(~ row + col)
		reference = lm.fit(design, Y[index])
		expect_equal(fit$gaussian_diagnostics$final_sse, sum(reference$residuals^2), tolerance = 1e-6)
	}
	expect_gt(refinements, 0L)
})
