# default prior of the innovation variance of dynamic latent positions: an inverse gamma on the stationary
# variance of a coordinate path (?lame, prior$uv_var_shape and prior$uv_var_scale; R/uv_var_prior.R)

uvp_default = lame:::.lame_uv_var_prior_resolve(list(), 1)
uvp_grid = seq(-0.9895, 0.9895, by = 0.001)

# the posterior of |rho| given known paths x (m x TT) with the variance integrated out in closed form,
# under rho ~ N(0.7, 0.3) and an inverse gamma (a, b) on the stationary variance or on sigma^2
uvp_posterior = function(x, a, b, stationary = TRUE, grid = uvp_grid) {
	TT = ncol(x)
	sxy = sum(x[, -1] * x[, -TT])
	syy = sum(x[, -TT]^2)
	stt = sum(x[, -1]^2)
	s1 = sum(x[, 1]^2)
	n1 = nrow(x)
	ntr = nrow(x) * (TT - 1)
	om = 1 - grid^2
	ss = stt - 2 * grid * sxy + grid^2 * syy
	lp = if (stationary) {
		-(ntr / 2) * log(om) - ((n1 + ntr) / 2 + a) * log((s1 + ss / om + 2 * b) / 2)
	} else {
		-(ntr / 2 + a) * log(om) - ((n1 + ntr) / 2 + a) * log((s1 + (ss + 2 * b) / om) / 2)
	}
	lp = lp + dnorm(grid, 0.7, 0.3, log = TRUE)
	w = exp(lp - max(lp))
	w_fold = w[grid > 0] + rev(w[grid < 0])
	w_fold / sum(w_fold)
}

uvp_calibration = function(rho, w0, m, TT, reps, ...) {
	a = uvp_grid[uvp_grid > 0]
	prior_args = list(...)
	replicate(reps, {
		x = matrix(0, m, TT)
		x[, 1] = rnorm(m, 0, sqrt(w0))
		for (t in 2:TT) x[, t] = rho * x[, t - 1] + sqrt(w0 * (1 - rho^2)) * rnorm(m)
		w = do.call(uvp_posterior, c(list(x), prior_args))
		cw = cumsum(w)
		mu = sum(w * a)
		s = sqrt(sum(w * (a - mu)^2))
		c(cover = a[which(cw >= 0.05)[1]] <= rho && rho <= a[which(cw >= 0.95)[1]], z = (mu - rho) / s)
	})
}

test_that("the default variance prior keeps the persistence calibrated on the latent scale", {
	skip_on_cran()
	# 40 coordinate paths over four periods with the stationary variance of the calibration study's
	# binary and poisson panels; the former inverse gamma on sigma^2 pulls a true 0.9 down
	set.seed(2027)
	for (rho in c(0.3, 0.6, 0.9)) {
		res = uvp_calibration(
			rho,
			w0 = 0.36,
			m = 40,
			TT = 4,
			reps = 300,
			a = uvp_default$uv_var_shape,
			b = uvp_default$uv_var_scale
		)
		# Monte Carlo sd of the coverage is 0.017 at 300 panels
		expect_gt(mean(res["cover", ]), 0.85)
		expect_lt(mean(res["cover", ]), 0.96)
		expect_lt(abs(mean(res["z", ])), 0.5)
	}
	former = uvp_calibration(0.9, w0 = 0.36, m = 40, TT = 4, reps = 300, a = 2, b = 1, stationary = FALSE)
	expect_lt(mean(former["z", ]), -0.6)
})

test_that("the stationary prior is the innovation prior with scale b (1 - rho^2) at each rho", {
	pr = lame:::.lame_uv_var_prior_resolve(list(), 0.5)
	expect_identical(pr$uv_var_scale, 0.25)
	expect_true(lame:::.lame_uv_var_stationary(lame:::.lame_uv_var_prior_resolve(list(), 1, "t")))
	expect_false(lame:::.lame_uv_var_stationary(
		lame:::.lame_uv_var_prior_resolve(list(sigma_uv_shape = 2), 1, "t")
	))
	expect_true(lame:::.lame_uv_var_stationary(pr))
	expect_null(pr$sigma_uv_shape)
	expect_null(pr$sigma_uv_scale)
	expect_equal(lame:::.lame_uv_ig_scale(pr, 0.6), pr$uv_var_scale * 0.64)
	expect_identical(lame:::.lame_uv_ig_shape(pr), pr$uv_var_shape)
	# the density of w at s2 / (1 - rho^2) on (rho, s2) is the inverse gamma of s2 with that scale
	ldig = function(x, a, b) a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x
	for (r in c(0.2, 0.9)) {
		s2 = c(0.05, 0.3, 1.2)
		expect_equal(
			ldig(s2 / (1 - r^2), pr$uv_var_shape, pr$uv_var_scale) - log(1 - r^2),
			ldig(s2, pr$uv_var_shape, lame:::.lame_uv_ig_scale(pr, r))
		)
	}
	innovation = lame:::.lame_uv_var_prior_resolve(list(sigma_uv_scale = 3), 0.5)
	expect_false(lame:::.lame_uv_var_stationary(innovation))
	expect_identical(innovation$sigma_uv_shape, 2)
	expect_identical(lame:::.lame_uv_ig_scale(innovation, 0.6), 3)
	expect_null(innovation$uv_var_shape)
	expect_error(
		lame:::.lame_uv_var_prior_resolve(list(sigma_uv_shape = 2, uv_var_scale = 1), 1),
		class = "lame_uv_var_prior"
	)
	expect_error(lame:::.lame_uv_var_prior_resolve(list(uv_var_shape = -1), 1), class = "lame_uv_var_prior")
	expect_error(lame:::.lame_uv_var_prior_resolve(list(uv_var_scale = NA_real_), 1), class = "lame_uv_var_prior")
})

skip_on_cran()

# known AR(1) coordinate paths on both sides of a bipartite panel
uvp_paths = function(seed, nA = 15L, nB = 12L, N = 5L, rho = 0.9, w0 = 0.4) {
	set.seed(seed)
	mk = function(n) {
		a = array(0, c(n, 1L, N))
		a[, , 1] = rnorm(n, 0, sqrt(w0))
		for (t in 2:N) a[, , t] = rho * a[, , t - 1] + rnorm(n, 0, sqrt(w0 * (1 - rho^2)))
		a
	}
	list(U = mk(nA), V = mk(nB))
}

# mean and sd of |rho| from a chain of draws against the exact posterior
uvp_exact = function(p, a, b) {
	x = rbind(matrix(p$U, ncol = dim(p$U)[3]), matrix(p$V, ncol = dim(p$V)[3]))
	w = uvp_posterior(x, a, b)
	g = uvp_grid[uvp_grid > 0]
	m = sum(w * g)
	c(mean = m, sd = sqrt(sum(w * (g - m)^2)))
}

test_that("the block draw of persistence and stationary variance targets the stationary prior", {
	p = uvp_paths(11)
	ex = uvp_exact(p, 1.5, 0.5)
	set.seed(3)
	s = c(0.5, 0.3)
	draws = numeric(6000)
	for (i in seq_along(draws)) {
		o = lame:::.lame_uv_rho_w_step(
			p$U,
			p$V,
			s[1],
			s[2],
			0.7,
			0.3,
			1.5,
			0.5,
			stationary = TRUE
		)
		s = c(o$rho, o$sigma)
		draws[i] = abs(s[1])
	}
	draws = draws[-(1:500)]
	expect_lt(abs(mean(draws) - ex[["mean"]]), 0.15 * ex[["sd"]])
	expect_lt(abs(sd(draws) / ex[[2]] - 1), 0.15)
})

test_that("the snap persistence and variance draws target the stationary prior without resets", {
	p = uvp_paths(12)
	ex = uvp_exact(p, 1.5, 0.5)
	pr = list(uv_var_shape = 1.5, uv_var_scale = 0.5)
	set.seed(4)
	rho = 0.5
	sigma = 0.3
	du = matrix(0, dim(p$U)[1], dim(p$U)[3])
	dv = matrix(0, dim(p$V)[1], dim(p$V)[3])
	draws = numeric(6000)
	for (i in seq_along(draws)) {
		rho = lame:::.lame_sample_rho_uv_drift_only(
			p$U,
			p$V,
			sigma,
			rho,
			FALSE,
			du,
			dv,
			prior_mean = 0.7,
			prior_sd = 0.3,
			stationary_prior = c(pr$uv_var_shape, pr$uv_var_scale)
		)
		sigma = lame:::.lame_sample_sigma_uv_drift_only(
			p$U,
			p$V,
			rho,
			sigma,
			FALSE,
			du,
			dv,
			lame:::.lame_uv_ig_shape(pr),
			lame:::.lame_uv_ig_scale(pr, rho)
		)
		draws[i] = abs(rho)
	}
	draws = draws[-(1:500)]
	expect_lt(abs(mean(draws) - ex[["mean"]]), 0.15 * ex[["sd"]])
	expect_lt(abs(sd(draws) / ex[[2]] - 1), 0.15)
})

test_that("the non-centred persistence move keeps the persistence prior at a fixed stationary variance", {
	# with a flat likelihood the move leaves rho at its prior given w, whatever the variance prior
	set.seed(5)
	paths = matrix(rnorm(40 * 4), 40, 4)
	rho = 0.6
	sigma = 0.5
	draws = numeric(4000)
	for (i in seq_along(draws)) {
		o = lame:::.lame_binary_dynamic_nc_ar_step(
			paths,
			rho,
			sigma,
			function(x) 0,
			0.7,
			0.3,
			2,
			1,
			stationary = TRUE
		)
		paths = o$paths
		rho = o$rho
		sigma = o$sigma
		draws[i] = rho
	}
	expect_equal(sigma^2 / (1 - rho^2), 0.25 / 0.64, tolerance = 1e-10)
	ref = rnorm(2e5, 0.7, 0.3)
	ref = ref[abs(ref) < 0.99]
	expect_lt(abs(mean(draws) - mean(ref)), 0.04)
	expect_lt(abs(sd(draws) - sd(ref)), 0.04)
})

test_that("the non-centred persistence and variance step leaves the joint prior under a flat likelihood", {
	set.seed(6)
	n = 6
	N = 4
	U = array(rnorm(n * 2 * N), c(n, 2, N))
	V = array(rnorm(n * 2 * N), c(n, 2, N))
	E = array(0, c(n, n, N))
	rho = 0.5
	sigma = 0.6
	for (stationary in c(TRUE, FALSE)) {
		rho = 0.5
		sigma = 0.6
		draws = matrix(NA_real_, 6000, 2)
		for (i in seq_len(nrow(draws))) {
			o = lame:::.lame_uv_nc_hyper_step(
				U,
				V,
				E,
				1e12,
				rho,
				sigma,
				0.7,
				0.3,
				2,
				1,
				stationary = stationary,
				step_rho = 0.6,
				step_w = 0.8
			)
			U = o$U
			V = o$V
			rho = o$rho
			sigma = o$sigma
			draws[i, ] = c(rho, if (stationary) sigma^2 / (1 - rho^2) else sigma^2)
		}
		draws = draws[-(1:500), ]
		ref = rnorm(2e5, 0.7, 0.3)
		ref = ref[abs(ref) < 0.99]
		expect_lt(abs(mean(draws[, 1]) - mean(ref)), 0.05)
		# the variance the prior is on has E[log x] = log(scale) - digamma(shape)
		expect_lt(abs(mean(log(draws[, 2])) + digamma(2)), 0.1)
	}
	# the paths follow the draws with their whitened innovations unchanged
	whiten = function(x, r, s) {
		e = x
		e[, , 1] = x[, , 1] * sqrt(1 - r^2) / s
		for (t in 2:N) e[, , t] = (x[, , t] - r * x[, , t - 1]) / s
		e
	}
	e0 = whiten(U, rho, sigma)
	o = lame:::.lame_uv_nc_hyper_step(U, V, E, 1e12, rho, sigma, 0.7, 0.3, 2, 1, step_rho = 0.6, step_w = 0.8)
	expect_gt(o$accepts, 0)
	expect_equal(whiten(o$U, o$rho, o$sigma), e0)
	# the carried latent array moves with the product and keeps every observed binary cell on its side;
	# a zero residual with a variance near zero leaves only the carried proposals to be accepted
	prod_of = function(U, V) {
		P = array(0, c(n, n, N))
		for (t in 1:N) P[, , t] = U[, , t] %*% t(V[, , t])
		P
	}
	Y = array(rbinom(n * n * N, 1, 0.5), c(n, n, N))
	Z = ifelse(Y > 0, 1, -1) * (abs(rnorm(n * n * N)) + 2)
	moved = 0
	for (i in 1:50) {
		o = lame:::.lame_uv_nc_hyper_step(
			U,
			V,
			prod_of(U, V),
			1e-12,
			rho,
			sigma,
			0.7,
			0.3,
			2,
			1,
			step_rho = 0.6,
			step_w = 0.8,
			carry = list(family = "binary", Z = Z, Y = Y)
		)
		expect_equal(o$Z - prod_of(o$U, o$V), Z - prod_of(U, V))
		expect_true(all((o$Z > 0) == (Y > 0)))
		moved = moved + o$accepts
		U = o$U
		V = o$V
		Z = o$Z
		rho = o$rho
		sigma = o$sigma
	}
	expect_gt(moved, 0)
	# with no observed poisson cell the carried proposals leave the joint prior too
	Yna = array(NA_real_, c(n, n, N))
	draws = matrix(NA_real_, 6000, 2)
	for (i in seq_len(nrow(draws))) {
		o = lame:::.lame_uv_nc_hyper_step(
			U,
			V,
			E,
			1e12,
			rho,
			sigma,
			0.7,
			0.3,
			2,
			1,
			step_rho = 0.6,
			step_w = 0.8,
			carry = list(family = "poisson", Z = Z, Y = Yna, log_exposure = 0)
		)
		U = o$U
		V = o$V
		Z = o$Z
		rho = o$rho
		sigma = o$sigma
		draws[i, ] = c(rho, sigma^2 / (1 - rho^2))
	}
	draws = draws[-(1:500), ]
	expect_lt(abs(mean(draws[, 1]) - mean(ref)), 0.05)
	expect_lt(abs(mean(log(draws[, 2])) + digamma(2)), 0.1)
})

test_that("the joint persistence HMC takes the stationary prior in its density and gradient", {
	set.seed(921865)
	n = 3L
	nt = 2L
	Y = array(NA_real_, c(n, n, nt))
	Y[cbind(c(1, 2, 2, 3, 1, 3), c(2, 1, 3, 2, 3, 1), c(1, 1, 1, 1, 2, 2))] = c(1, 0, 0, 0, 0, 1)
	X = lapply(seq_len(nt), function(t) {
		x = array(rnorm(n * n * 2), c(n, n, 2))
		x[, , 1] = 1
		x
	})
	cx = lame:::.lame_binary_dynamic_pair_hmc_context(Y, X, 2L)
	P = matrix(c(3, 0.2, 0.2, 2), 2)
	Sab = matrix(c(0.2, 0.03, 0.03, 0.15), 2)
	h_inn = lame:::.lame_binary_dynamic_ar_hmc_hyper(cx, P, Sab, 0.6, 0.25, -0.4, 1.5, 0.4, 0.7, 0.3)
	h_sta = lame:::.lame_binary_dynamic_ar_hmc_hyper(cx, P, Sab, 0.6, 0.25, -0.4, 1.5, 0.4, 0.7, 0.3, stationary = TRUE)
	for (rr in c(-0.85, 0, 0.7, 0.989)) {
		q = c(rnorm(cx$size, sd = 0.3), log(0.4), atanh(rr / 0.99))
		l = q[cx$size + 1]
		r = 0.99 * tanh(q[cx$size + 2])
		w = exp(2 * l)
		v = w * (1 - r^2)
		# the two densities differ by the variance priors in (l, alpha)
		expect_equal(
			lame:::.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h_sta)$logp -
				lame:::.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h_inn)$logp,
			(-2 * 1.5 * l - 0.4 / w) - (-1.5 * log(v) - 0.4 / v),
			tolerance = 1e-10
		)
		ev = lame:::.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h_sta)
		eps = 1e-6
		grad = vapply(
			seq_along(q),
			function(j) {
				up = q
				dn = q
				up[j] = up[j] + eps
				dn[j] = dn[j] - eps
				(lame:::.lame_binary_dynamic_ar_hmc_evaluate(up, cx, h_sta)$logp -
					lame:::.lame_binary_dynamic_ar_hmc_evaluate(dn, cx, h_sta)$logp) / (2 * eps)
			},
			numeric(1)
		)
		expect_lt(max(abs(grad - ev$gradient) / (1 + abs(grad))), 1e-6)
	}
})

# a directed panel with persistent positions
uvp_panel = function(family = "normal", n = 8, TT = 4, rho = 0.8, seed = 41) {
	set.seed(seed)
	U = V = array(0, c(n, 2, TT))
	U[, , 1] = rnorm(2 * n)
	V[, , 1] = rnorm(2 * n)
	for (t in 2:TT) {
		U[, , t] = rho * U[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n)
		V[, , t] = rho * V[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n)
	}
	nm = sprintf("a%02d", 1:n)
	Y = array(0, c(n, n, TT), dimnames = list(nm, nm, paste0("t", 1:TT)))
	for (t in 1:TT) {
		eta = U[, , t] %*% t(V[, , t]) + matrix(rnorm(n * n, 0, 0.5), n, n)
		Y[, , t] = switch(family, normal = eta, binary = 1 * (eta > 0), poisson = rpois(n * n, exp(0.3 * eta)))
		diag(Y[, , t]) = NA
	}
	Y
}

uvp_fit = function(args) suppressWarnings(suppressMessages(do.call(lame, args)))

uvp_args = function(Y, family = "normal", ...) {
	modifyList(
		list(
			Y = Y,
			R = 1,
			family = family,
			dynamic_uv = TRUE,
			burn = 0,
			nscan = 8,
			odens = 1,
			seed = 116,
			gof = FALSE,
			plot = FALSE,
			verbose = FALSE
		),
		list(...)
	)
}

test_that("every dynamic_uv route resolves the stationary-variance prior and keeps explicit settings", {
	Y = uvp_panel()
	Yb = Y[, 1:6, ]
	dimnames(Yb)[[2]] = paste0("c", 1:6)
	fits = list(
		normal = uvp_fit(uvp_args(Y)),
		binary = uvp_fit(uvp_args(uvp_panel("binary"), "binary")),
		poisson = uvp_fit(uvp_args(uvp_panel("poisson"), "poisson")),
		snap = uvp_fit(uvp_args(Y, dynamic_uv_kind = "snap")),
		symmetric = uvp_fit(uvp_args(
			array(apply(Y, 3, function(y) (y + t(y)) / 2), dim(Y), dimnames(Y)),
			symmetric = TRUE
		)),
		bipartite = uvp_fit(uvp_args(Yb, mode = "bipartite")),
		bipartite_binary = uvp_fit(uvp_args(1 * (Yb > 0), "binary", mode = "bipartite"))
	)
	for (nm in names(fits)) {
		pr = fits[[nm]]$prior
		expect_identical(pr$uv_var_shape, uvp_default$uv_var_shape, info = nm)
		expect_true(is.numeric(pr$uv_var_scale) && pr$uv_var_scale > 0, info = nm)
		expect_null(pr$sigma_uv_shape, info = nm)
		expect_null(pr$sigma_uv_scale, info = nm)
		expect_identical(fits[[nm]]$refit_options$prior, list(), info = nm)
	}
	# the scale is half the unit of the positions: the factor scale of a directed normal panel, one otherwise
	expect_equal(fits$normal$prior$uv_var_scale, 0.5 * lame:::.lame_sigma_uv_scale_default(Y, "normal"))
	expect_identical(fits$binary$prior$uv_var_scale, 0.5)
	expect_identical(fits$bipartite$prior$uv_var_scale, 0.5)
	# the "t" kind takes the stationary prior too, and keeps the innovation prior when it is set
	tk = uvp_fit(uvp_args(Y, dynamic_uv_kind = "t"))
	expect_identical(tk$prior$uv_var_shape, uvp_default$uv_var_shape)
	expect_equal(tk$prior$uv_var_scale, 0.5 * lame:::.lame_sigma_uv_scale_default(Y, "normal"))
	expect_null(tk$prior$sigma_uv_shape)
	expect_null(tk$prior$sigma_uv_scale)
	tk_innovation = uvp_fit(uvp_args(Y, dynamic_uv_kind = "t", prior = list(sigma_uv_shape = 2)))
	expect_identical(tk_innovation$prior$sigma_uv_shape, 2)
	expect_equal(tk_innovation$prior$sigma_uv_scale, 2 * tk$prior$uv_var_scale)
	expect_null(tk_innovation$prior$uv_var_shape)
	# an explicit innovation prior is used as given, with the other entry at its former default
	inn = uvp_fit(uvp_args(Y, prior = list(sigma_uv_shape = 3)))
	expect_identical(inn$prior$sigma_uv_shape, 3)
	expect_equal(inn$prior$sigma_uv_scale, lame:::.lame_sigma_uv_scale_default(Y, "normal"))
	expect_null(inn$prior$uv_var_shape)
	sta = uvp_fit(uvp_args(Y, prior = list(uv_var_scale = 2)))
	expect_identical(sta$prior$uv_var_scale, 2)
	expect_identical(sta$prior$uv_var_shape, uvp_default$uv_var_shape)
	expect_error(
		uvp_fit(uvp_args(Y, prior = list(sigma_uv_scale = 1, uv_var_scale = 1))),
		class = "lame_uv_var_prior"
	)
})

test_that("the default fit is the fit that states the default variance prior", {
	for (family in c("normal", "binary")) {
		Y = uvp_panel(family, n = 10, TT = 4)
		implicit = uvp_fit(uvp_args(Y, family, burn = 50, nscan = 100, R = 2))
		stated = uvp_fit(uvp_args(
			Y,
			family,
			burn = 50,
			nscan = 100,
			R = 2,
			prior = implicit$prior[c("uv_var_shape", "uv_var_scale")]
		))
		for (nm in c("BETA", "VC", "rho_uv", "sigma_uv", "UVPM", "U_draws", "EZ")) {
			expect_identical(implicit[[nm]], stated[[nm]], info = paste(family, nm))
		}
	}
})

test_that("checkpoints keep the variance prior they were written with", {
	Y = uvp_panel(n = 8, TT = 4)
	for (pr in list(list(), list(sigma_uv_shape = 2, sigma_uv_scale = 1))) {
		ck = tempfile(fileext = ".rds")
		withr::defer(unlink(ck))
		original = uvp_fit(uvp_args(
			Y,
			prior = pr,
			checkpoint_path = ck,
			max_seconds = 1e-12
		))
		resumed = suppressMessages(lame_resume(ck, nscan_more = 8))
		expect_identical(resumed$prior, original$prior)
		expect_identical(
			lame:::.lame_uv_var_stationary(resumed$prior),
			!length(pr)
		)
	}
})
