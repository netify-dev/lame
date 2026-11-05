skip_on_cran()

library(lame)

# the penalised dynamic route (bias_correct = "none") converges on ordinary panels, or is refused with the state:
# a fit that stops short of its criterion is an error of class lame_als_convergence whose field `fit` is the last iterate

dc_names = function(n) paste0("a", sprintf("%02d", seq_len(n)))

# a normal panel in which one actor is present but has no observed cell in any period, and another none in one period
dc_unobserved = function(n = 10L, Tt = 4L, seed = 7) {
	set.seed(seed)
	rn = dc_names(n)
	a = rnorm(n, 0, 0.6)
	b = rnorm(n, 0, 0.6)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x1"))
		y = 0.3 + 0.5 * x[, , 1] + outer(a, b, "+") * (1 + 0.2 * t) + matrix(rnorm(n * n), n)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
		X[[t]] = x
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	for (t in seq_len(Tt)) {
		Y[[t]][3, ] = NA
		Y[[t]][, 3] = NA
	}
	Y[[2]][5, ] = NA
	Y[[2]][, 5] = NA
	list(Y = Y, X = X)
}

test_that("an actor that is present but has no observed cell does not stall the penalised fit: it is not part of the gauge", {
	d = dc_unobserved()
	drop3 = function(L) lapply(L, function(m) if (length(dim(m)) == 3L) m[-3, -3, , drop = FALSE] else m[-3, -3])
	# (the observed-dyad graph of such a panel has an isolated actor: the fit says so, and it is not the subject here)
	quiet = function(expr) suppressWarnings(expr)
	fit = quiet(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "normal",
		dynamic_ab = TRUE,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_true(fit$converged)
	expect_true(fit$diagnostics$stationarity$pass)
	expect_lt(fit$iterations, 20L)
	# its path is held at zero, which is the minimiser of its penalty
	expect_true(all(fit$a_dynamic[3, ] == 0))
	expect_true(all(fit$b_dynamic[3, ] == 0))
	# and the others are what they are without it
	ref = quiet(lame:::lame_dynamic_als(
		drop3(d$Y),
		drop3(d$X),
		family = "normal",
		dynamic_ab = TRUE,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_true(ref$converged)
	expect_equal(fit$coef_path, ref$coef_path, tolerance = 1e-3)
	keep = match(rownames(ref$a_dynamic), rownames(fit$a_dynamic))
	expect_equal(fit$a_dynamic[keep, ], ref$a_dynamic, tolerance = 1e-3, ignore_attr = TRUE)
	expect_equal(fit$b_dynamic[keep, ], ref$b_dynamic, tolerance = 1e-3, ignore_attr = TRUE)
	# symmetric networks hold the same convention
	Ys = lapply(d$Y, function(y) { y[lower.tri(y)] = t(y)[lower.tri(y)]; y })
	# (a symmetric fit needs a symmetric covariate)
	Xs = lapply(d$X, function(x) {
		if (length(dim(x)) == 3L) { for (k in seq_len(dim(x)[3L])) x[, , k] = (x[, , k] + t(x[, , k])) / 2 } else x = (x + t(x)) / 2
		x })
	fs = quiet(lame:::lame_dynamic_als(
		Ys,
		Xs,
		family = "normal",
		dynamic_ab = TRUE,
		symmetric = TRUE,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_true(fs$converged)
	expect_true(all(fs$a_dynamic[3, ] == 0))
})

# a binary panel with strong rank-2 structure and a sparse symmetric count panel, both fit at rank 2
dc_binary = function(n = 14L, Tt = 3L, seed = 31) {
	set.seed(seed)
	rn = dc_names(n)
	U = matrix(rnorm(n * 2), n, 2)
	V = matrix(rnorm(n * 2), n, 2)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		y = 1 * (-0.2 + 2.5 * tcrossprod(U, V) + rnorm(n * n, 0, 0.3) > 0)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
	}
	names(Y) = paste0("t", seq_len(Tt))
	Y
}

dc_poisson_sym = function(n = 12L, Tt = 3L, seed = 41) {
	set.seed(seed)
	rn = dc_names(n)
	U = matrix(rnorm(n * 2), n, 2)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		y = matrix(rpois(n * n, exp(-1.2 + 0.9 * tcrossprod(U, U))), n)
		y[lower.tri(y)] = t(y)[lower.tri(y)]
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
	}
	names(Y) = paste0("t", seq_len(Tt))
	Y
}

test_that("the static factors of the binary and count families are penalised, so that their fit has a minimiser and converges", {
	skip_on_cran()
	fb = lame:::lame_dynamic_als(
		dc_binary(),
		family = "binary",
		R = 2,
		dynamic_ab = TRUE,
		bias_correct = "none",
		verbose = FALSE
	)
	expect_true(fb$converged)
	expect_true(fb$diagnostics$stationarity$pass)
	expect_false(fb$diagnostics$descent_guard$active)
	expect_lt(fb$iterations, 400L)
	# the penalised product stays bounded
	expect_lt(sqrt(sum(fb$UVPM^2, na.rm = TRUE)), 100)
	expect_true(all(is.finite(fb$EZ[[1]][!is.na(fb$Y[, , 1])])))
	fp = lame:::lame_dynamic_als(
		dc_poisson_sym(),
		family = "poisson",
		R = 2,
		symmetric = TRUE,
		dynamic_ab = TRUE,
		bias_correct = "none",
		verbose = FALSE
	)
	expect_true(fp$converged)
	expect_true(fp$diagnostics$stationarity$pass)
	expect_lt(sqrt(sum(fp$UVPM^2, na.rm = TRUE)), 100)
	# the ridge is the level prior of the dynamic and bipartite factors: the penalised objective carries it
	obj = fp$diagnostics$objective
	expect_true(
		all(
			diff(obj) <= 1e-8 * abs(utils::head(obj, -1L)) + 1e-8 | seq_along(diff(obj)) %in% fp$diagnostics$reweight_steps
		)
	)
})

test_that("the ridge solvers minimise their penalised objective and are descent steps", {
	set.seed(3)
	n = 9
	O0 = tcrossprod(matrix(rnorm(n * 2), n, 2), matrix(rnorm(n * 2), n, 2))
	Rbar = O0 + matrix(rnorm(n * n, 0, 0.3), n)
	C = matrix(runif(n * n, 0.2, 1), n)
	C[sample(n * n, 12)] = 0
	obs = C > 0
	for (ridge in c(0.1, 1, 4)) {
		for (method in c("mm", "als", "hybrid")) {
			Om = matrix(0, n, n)
			prev = Inf
			for (it in 1:8) {
				up = lame:::.lda_lowrank_update_ridge(Om, Rbar, C, obs, 2L, FALSE, method, ridge)
				val = lame:::.ae_lowrank_obj(up$Omat, Rbar, C, obs) + 2 * ridge * sum(svd(up$Omat, 0, 0)$d)
				expect_lte(val, prev + 1e-9, label = paste(method, ridge, it))
				prev = val
				Om = up$Omat
			}
		}
	}
	# the eigenmodel: the same, symmetric, with the signed eigenvalues soft-thresholded
	Rs = (Rbar + t(Rbar)) / 2
	Cs = (C + t(C)) / 2
	obs_s = Cs > 0
	Om = matrix(0, n, n)
	prev = Inf
	for (it in 1:8) {
		up = lame:::.lda_lowrank_update_ridge(Om, Rs, Cs, obs_s, 2L, TRUE, "mm", 1)
		val = lame:::.ae_lowrank_obj(up$Omat, Rs, Cs, obs_s) + 2 * 1 * sum(abs(up$L))
		expect_lte(val, prev + 1e-9)
		prev = val
		Om = up$Omat
	}
})

test_that("a fit that stops on its iteration limit is an error that carries the last iterate; a normal fit's limit is raised once", {
	d = dc_unobserved()
	quiet = function(expr) suppressWarnings(expr)
	err = expect_error(
		quiet(lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "normal",
			dynamic_ab = TRUE,
			R = 1,
			max_iter = 1,
			bias_correct = "none",
			verbose = FALSE
		)),
		class = "lame_als_convergence"
	)
	expect_false(inherits(err, "warning"))
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(msg, "did not reach its convergence criterion within `max_iter` = 1 (raised once to 4)", fixed = TRUE)
	expect_match(msg, "largest scaled first-order score")
	fit = err$fit
	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$converged)
	expect_identical(err$kind, "bcd")
	expect_identical(fit$iterations, 4L)
	# the limit of a normal fit is raised once, by the static estimators' factor: a fit that needs k sweeps is refused
	# at a max_iter whose fourfold falls short of k and fitted at one whose fourfold reaches it (k is read off a fit with room)
	fit_at = function(mi, tol) quiet(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "normal",
		dynamic_ab = TRUE,
		max_iter = mi,
		tol = tol,
		bias_correct = "none",
		verbose = FALSE
	))
	k = fit_at(200, 1e-12)$iterations
	expect_gt(k, 4L)
	m_short = (k - 1L) %/% 4L
	m_enough = (k + 3L) %/% 4L
	expect_lt(4L * m_short, k)
	expect_gte(4L * m_enough, k)
	expect_error(fit_at(m_short, 1e-12), class = "lame_als_convergence")
	ok = fit_at(m_enough, 1e-12)
	expect_true(ok$converged)
	expect_identical(ok$iterations, k)
	expect_gt(ok$iterations, m_enough)
	# the same panel at the default limit is a fit, with no condition of any kind
	seen = character()
	ok = withCallingHandlers(
		lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, bias_correct = "none", verbose = FALSE),
		warning = function(w) {
			seen <<- c(seen, class(w))
			invokeRestart("muffleWarning")
		}
	)
	expect_false("lame_als_convergence" %in% seen)
	expect_true(ok$converged)
})

test_that("a reweighted step that raises the penalised deviance is damped, not the end of the fit", {
	set.seed(12)
	n = 12L
	Tt = 3L
	rn = dc_names(n)
	a = rnorm(n, 0, 0.7)
	Y = lapply(seq_len(Tt), function(t) {
		y = 1 * (-0.3 + 0.8 * outer(a, a, "+") + matrix(rnorm(n * n), n) > 0)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		y
	})
	names(Y) = paste0("t", seq_len(Tt))
	fit_it = function() lame:::lame_dynamic_als(
		Y,
		family = "binary",
		dynamic_ab = TRUE,
		bias_correct = "none",
		verbose = FALSE
	)
	ref = fit_it()
	expect_true(ref$converged)
	orig = get(".lda_penalized_deviance", envir = asNamespace("lame"))
	# the deviance at the second and the third gate is raised by 1,000: the step to that iterate overshoots
	for (hit in 2:3) {
		calls = 0L
		fit = testthat::with_mocked_bindings(
			fit_it(),
			.lda_penalized_deviance = function(...) {
				calls <<- calls + 1L
				v = orig(...)
				if (calls == hit) v + 1e3 else v
			},
			.package = "lame"
		)
		expect_true(fit$converged, info = hit)
		expect_gt(fit$iterations, ref$iterations, label = paste("iterations after an overshoot at gate", hit))
		expect_equal(fit$coefficients, ref$coefficients, tolerance = 1e-6, info = hit)
	}
})

test_that("an MCMC chain starts from the last iterate of a penalised fit that was refused", {
	set.seed(3)
	n = 10L
	Tt = 3L
	rn = dc_names(n)
	u = rnorm(n)
	Y = lapply(seq_len(Tt), function(t) {
		y = matrix(rnorm(n * n), n) + 0.8 * tcrossprod(u)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		y
	})
	names(Y) = paste0("t", seq_len(Tt))
	# the premise: the penalised fit that starts a directed dynamic_uv chain is refused at this limit
	err = expect_error(
		suppressWarnings(lame:::lame_dynamic_als(
			Y,
			family = "normal",
			R = 1,
			dynamic_uv = TRUE,
			max_iter = 1,
			bias_correct = "none",
			verbose = FALSE
		)),
		class = "lame_als_convergence"
	)
	expect_s3_class(err$fit, "lame_dynamic_als")
	# and the chain still starts from it: the refusal carries the last iterate, which is a start, not a certified fit
	fit = suppressWarnings(suppressMessages(
		lame(
			Y,
			family = "normal",
			R = 1,
			dynamic_uv = TRUE,
			als_max_iter = 1,
			nscan = 60,
			burn = 20,
			odens = 4,
			seed = 11,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	))
	expect_s3_class(fit, "lame")
	expect_true(isTRUE(fit$als_start))
	expect_true(all(is.finite(fit$BETA)))
})

test_that("a jackknife half that is refused is retried from its last iterate", {
	skip_on_cran()
	set.seed(3112001)
	n = 14L
	Tt = 4L
	al = rnorm(n, 0, .3)
	be = rnorm(n, 0, .3)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n)
		y = matrix(rpois(n * n, exp(.5 + .5 * x + outer(al, be, "+") + rnorm(n * n, 0, sqrt(.4)))), n)
		diag(y) = NA_real_
		dimnames(y) = list(dc_names(n), dc_names(n))
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1L), dimnames = c(dimnames(y), list("x")))
	}
	# the halves get two iterations, which refuses them (the full fit has its own, ample limit); the retry at five times
	# that limit starts from the iterate the refusal carries
	orig = get(".lda_count_jackknife", envir = asNamespace("lame"))
	fit = testthat::with_mocked_bindings(
		suppressWarnings(lame(
			Y,
			Xdyad = X,
			R = 0,
			family = "poisson",
			method = "als",
			dynamic_ab = TRUE,
			als_max_iter = 200,
			als_tol = 1e-6,
			bias_correct = "jackknife",
			verbose = FALSE
		)),
		.lda_count_jackknife = function(object, call_args, max_iter, tol, ...) orig(
			object,
			call_args,
			max_iter = 2L,
			tol,
			...
		),
		.package = "lame"
	)
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_true(all(fit$bias_correct$converged))
	for (h in 1:2) {
		att = fit$bias_correct$attempts[[h]]
		expect_length(att, 2L)
		expect_false(is.null(att[[1L]]$error))
		expect_true(att[[2L]]$converged)
	}
})
