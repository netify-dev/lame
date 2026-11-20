# the posterior mode of a two-mode fit holds the factor of its smaller side at unit scale (the receivers' unit scale is
# the prior of its variational stage),
# the variational budget, and the lower-bound rule of column retirement (R/ame_als_lowrank.R, R/ame_als_retire.R)

kn = function(f) get(f, envir = asNamespace("lame"))

# the reference sweeps: sender scale from the incoming factors, then ridged row solves at 1 / su and column
# solves at unit precision until the sender factors stop moving
ref_fixed = function(z, w, u, v, su, suv0, kappa0, maxit = 15L, tol = 1e-6) {
	R = ncol(u)
	tr = 0
	for (i in seq_len(nrow(z))) tr = tr + diag(solve(crossprod(v, w[i, ] * v) + diag(1 / su, R)))
	su = (suv0 + colSums(u^2) + tr) / (nrow(z) + kappa0)
	for (it in seq_len(maxit)) {
		u0 = u
		for (i in seq_len(nrow(z))) u[i, ] = solve(crossprod(v, w[i, ] * v) + diag(1 / su, R), crossprod(v, w[i, ] * z[i, ]))
		for (j in seq_len(ncol(z))) v[j, ] = solve(crossprod(u, w[, j] * u) + diag(R), crossprod(u, w[, j] * z[, j]))
		if (max(abs(u - u0)) < tol * max(1, max(abs(u)))) break
	}
	list(U = u, V = v, sig2_u = su)
}

test_that("the factor update of a two-mode fit holds the receiver scale at 1 and updates the sender scale", {
	set.seed(2201)
	for (rank in c(1L, 2L, 7L)) {
		nr = 30L; nc = 12L
		u = matrix(rnorm(nr * rank), nr, rank); v = matrix(rnorm(nc * rank), nc, rank)
		z = matrix(rnorm(nr * nc), nr, nc); w = matrix(runif(nr * nc), nr, nc); w[w < 0.2] = 0
		su = runif(rank, 0.2, 2); sv = runif(rank, 0.2, 2)
		got = lame:::als_factor_balanced_cpp(z, w, u, v, su, sv, 1, 2 * rank + 2, 15L, 1e-6, TRUE)
		free = lame:::als_factor_update_cpp(z, w, u, v, su, sv, 1, 2 * rank + 2, 1L, 1e-6)
		expect_false(got$balanced)
		expect_identical(got$sig2_v, rep(1, rank))
		# the sender scale takes the same EM image as with a free receiver scale
		expect_equal(got$sig2_u, free$sig2_u, tolerance = 1e-12)
		ref = ref_fixed(z, w, u, v, su, 1, 2 * rank + 2)
		expect_equal(got$sig2_u, ref$sig2_u, tolerance = 1e-12)
		expect_equal(got$U, ref$U, tolerance = 1e-8)
		expect_equal(got$V, ref$V, tolerance = 1e-8)
		expect_equal(got$Omat, tcrossprod(got$U, got$V), tolerance = 1e-12)
	}
})

test_that("the penalized update holds the receiver scale only when asked, and a one-mode update is unchanged", {
	set.seed(2202)
	nr = 25L; nc = 10L; R = 2L
	z = matrix(rnorm(nr * nc), nr, nc); w = matrix(runif(nr * nc), nr, nc)
	u = matrix(rnorm(nr * R), nr, R); v = matrix(rnorm(nc * R), nc, R)
	obs = w > 0
	two = kn(".ae_lowrank_update_pen")(tcrossprod(u, v), z, w, obs, R, FALSE, u, v, c(0.7, 0.4), c(3, 2), 1, 6, fixed_receiver = TRUE)
	expect_identical(two$sig2_v, rep(1, R))
	ref = ref_fixed(z, w, u, v, c(0.7, 0.4), 1, 6)
	expect_equal(two$U, ref$U, tolerance = 1e-8)
	one = kn(".ae_lowrank_update_pen")(tcrossprod(u, v), z, w, obs, R, FALSE, u, v, c(0.7, 0.4), c(3, 2), 1, 6)
	old = kn(".ae_lowrank_em_als")(z, w, obs, R, u, v, c(0.7, 0.4), c(3, 2), 1, 6)
	expect_identical(one, old)
	expect_false(isTRUE(all(one$sig2_v == 1)))
})

test_that("the posterior mode of a two-mode fit holds the scale of its smaller side at 1 and certifies the other only", {
	skip_on_cran()
	set.seed(2203)
	nr = 60L; nc = 15L
	U = matrix(rnorm(nr * 2), nr); V = matrix(rnorm(nc * 2), nc)
	eta = -0.3 + rnorm(nr, 0, 0.5) + rep(rnorm(nc, 0, 0.5), each = nr) + 0.6 * tcrossprod(U, V)
	Y = matrix(rbinom(nr * nc, 1, pnorm(eta)), nr, dimnames = list(sprintf("r%02d", 1:nr), sprintf("c%02d", 1:nc)))
	modes = list()
	stage = kn(".ae_mf_stage")
	local_mocked_bindings(.ae_mf_stage = function(prep, ir, ...) {
		modes[[length(modes) + 1L]] <<- list(sv = ir$fit$sig2_v, su = ir$fit$sig2_u, pass = isTRUE(ir$stationarity$pass))
		stage(prep, ir, ...)
	}, .package = "lame")
	fit = suppressMessages(ame_als(Y, R = 2, family = "binary", mode = "bipartite"))
	expect_identical(fit$R, 2L)
	expect_null(fit$eb$retirement)
	expect_true(length(modes) >= 1L)
	for (m in modes) {
		expect_true(m$pass)
		expect_identical(m$sv, rep(1, 2))
		expect_true(all(is.finite(m$su) & m$su > 0))
	}
	# a table with more columns than rows holds the sender scale instead, the scale of the side with fewer actors
	modes = list()
	fit_t = suppressMessages(ame_als(t(Y), R = 2, family = "binary", mode = "bipartite"))
	expect_identical(fit_t$R, 2L)
	expect_true(length(modes) >= 1L)
	for (m in modes) {
		expect_true(m$pass)
		expect_identical(m$su, rep(1, 2))
		expect_true(all(is.finite(m$sv) & m$sv > 0))
	}
	# a directed one-mode fit keeps a receiver scale of its own
	n = 40L
	ed = -0.3 + outer(rnorm(n, 0, 0.5), rnorm(n, 0, 0.5), "+") + tcrossprod(rnorm(n), rnorm(n))
	Yd = matrix(rbinom(n * n, 1, pnorm(ed)), n, dimnames = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n))); diag(Yd) = NA
	modes = list()
	suppressMessages(ame_als(Yd, R = 1, family = "binary"))
	expect_true(length(modes) >= 1L)
	expect_false(identical(modes[[1L]]$sv, 1))
})

test_that("the variational stage runs at most 20 times the caller's budget under the cavity rule", {
	expect_identical(kn(".ae_mf_budget")(200L, "cavity"), 4000L)
	expect_identical(kn(".ae_mf_budget")(30L, "cavity"), 600L)
	expect_identical(kn(".ae_mf_budget")(200L, "plugin"), 400L)
})

test_that("an unsettled rank keeps its columns when its lower bound is above that of the next rank down", {
	mk = function(Rr, settled, elbo) {
		U = matrix(1, 6, Rr) %*% diag(seq(Rr, 1), Rr)
		list(
			fit = list(
				U = U,
				V = U,
				eb_mf = list(record = list(converged = settled, elbo = elbo,
					termination = if (settled) "converged" else "iteration limit"))
			),
			stationarity = list(pass = TRUE, termination = "converged"),
			retained = Rr
		)
	}
	init = list(U = matrix(seq_len(18), 6, 3), V = matrix(seq_len(18), 6, 3) + 100)
	go = function(elbos, settled_at = 1L, formed = identity) {
		seen = integer(0)
		run = function(init, Rr) { seen <<- c(seen, Rr); mk(Rr, Rr <= settled_at, elbos[Rr + 1L]) }
		ir = kn(".ae_retire_columns")(run, init, 3L, FALSE, TRUE, formed = formed)
		list(ir = ir, seen = seen)
	}
	# rank 3 has the higher bound: kept, unsettled, nothing retired (the caller returns it with a warning)
	k = go(c(-120, -110, -105, -100))
	expect_identical(k$seen, c(3L, 2L))
	expect_identical(k$ir$retained, 3L)
	expect_null(k$ir$retirement)
	expect_identical(k$ir$unsettled$rank, 3L)
	expect_identical(k$ir$unsettled$stage, "variational refinement")
	expect_identical(k$ir$unsettled$termination, "iteration limit")
	expect_identical(c(k$ir$unsettled$elbo, k$ir$unsettled$elbo_lower), c(-100, -105))
	# when its estimates cannot be formed it is retired as a rank that did not converge, and the rank below,
	# already fitted, goes on: rank 2 is above rank 1 but cannot be formed either, and rank 1 settles
	f = go(c(-120, -110, -105, -100), formed = function(ir) NULL)
	expect_identical(f$seen, c(3L, 2L, 1L))
	expect_identical(f$ir$retained, 1L)
	expect_null(f$ir$unsettled)
	expect_true(f$ir$retirement$settled)
	expect_identical(vapply(f$ir$retirement$steps, function(z) z$reason, ""), c("not formed", "not formed"))
	# rank 2 at least as high as rank 3 (margin 0): retired by the bound; rank 2 still unsettled and above rank 1: kept at 2
	e = go(c(-120, -110, -100, -100))
	expect_identical(e$seen, c(3L, 2L, 1L))
	expect_identical(e$ir$retained, 2L)
	expect_identical(e$ir$retirement$effective_rank, 2L)
	expect_length(e$ir$retirement$steps, 1L)
	expect_false(e$ir$retirement$settled)
	expect_identical(e$ir$retirement$steps[[1L]]$reason, "lower bound")
	expect_identical(e$ir$unsettled$rank, 2L)
	# every lower rank at least as high: retired down to the rank that settles
	d = go(c(-120, -90, -95, -100))
	expect_identical(d$seen, c(3L, 2L, 1L))
	expect_identical(d$ir$retained, 1L)
	expect_true(d$ir$retirement$settled)
	# a lower rank without a bound: rank 3 goes to rank 2 on the bound; a rank 1 that settled is preferred
	# to the unsettled rank 2
	n = go(c(-120, NA, -95, -100))
	expect_identical(n$seen, c(3L, 2L, 1L))
	expect_identical(n$ir$retained, 1L)
	expect_null(n$ir$unsettled)
	expect_true(n$ir$retirement$settled)
	# a rank 1 that did not settle and has no bound: rank 2 is never retired against it
	u = go(c(-120, NA, -95, -100), settled_at = 0L)
	expect_identical(u$seen, c(3L, 2L, 1L))
	expect_identical(u$ir$retained, 2L)
	expect_identical(u$ir$unsettled$rank, 2L)
	expect_true(is.na(u$ir$unsettled$elbo_lower))
	expect_identical(kn(".ae_ir_elbo")(list(fit = list())), NA_real_)
})
