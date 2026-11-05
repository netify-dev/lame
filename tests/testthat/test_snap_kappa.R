# the snap kind's reset scale: by default a reset is a fresh draw from the stationary law of the position paths
# (kappa^2 = sigma_uv^2 / (1 - rho_uv^2)) and enters the rho_uv and sigma_uv conditionals; a user-set
# prior$snap_kappa is a fixed scale that keeps the reset positions out of them
skip_on_cran()

snap_state = function(seed, n = 40, R = 2, TT = 5, rho = 0.7, sigma = 0.4, pi = 0.2) {
	set.seed(seed)
	s_stat = sigma / sqrt(1 - rho^2)
	draw = function() {
		X = array(0, c(n, R, TT)); d = matrix(0, n, TT)
		X[, , 1] = rnorm(n * R, 0, s_stat)
		for (t in 2:TT) {
			d[, t] = runif(n) < pi
			X[, , t] = rho * X[, , t - 1] + rnorm(n * R, 0, sigma)
			X[d[, t] == 1, , t] = rnorm(sum(d[, t]) * R, 0, s_stat)
		}
		list(X = X, d = d)
	}
	u = draw(); v = draw()
	list(U = u$X, V = v$X, du = u$d, dv = v$d, rho = rho, sigma = sigma)
}

test_that("the default reset scale is the stationary sd of the paths", {
	expect_equal(lame:::.lame_snap_kappa(NULL, 0.8, 0.36), 0.36 / sqrt(1 - 0.64))
	expect_equal(lame:::.lame_snap_kappa(2, 0.8, 0.36), 2)
	# the same bound on rho_uv as the samplers
	expect_equal(lame:::.lame_snap_kappa(NULL, 0.999, 0.1), 0.1 / sqrt(1 - 0.9801))
})

test_that("sigma_uv is drawn from its exact conditional with stationary resets", {
	skip_on_cran()
	st = snap_state(1)
	one_m = 1 - st$rho^2
	ss_drift = 0; n_drift = 0; ss_stat = 0; n_stat = 0
	for (X in list(list(st$U, st$du), list(st$V, st$dv))) {
		M = X[[1]]; d = X[[2]]
		ss_stat = ss_stat + sum(M[, , 1]^2); n_stat = n_stat + length(M[, , 1])
		for (t in 2:dim(M)[3]) {
			k = d[, t] == 0
			ss_drift = ss_drift + sum((M[k, , t] - st$rho * M[k, , t - 1])^2); n_drift = n_drift + sum(k) * dim(M)[2]
			ss_stat = ss_stat + sum(M[!k, , t]^2); n_stat = n_stat + sum(!k) * dim(M)[2]
		}
	}
	shape = (n_drift + n_stat) / 2 + 2
	rate = (ss_drift + one_m * ss_stat) / 2 + 1
	set.seed(2)
	prec = replicate(
		4000,
		1 / lame:::.lame_sample_sigma_uv_drift_only(
			st$U,
			st$V,
			st$rho,
			0.4,
			FALSE,
			st$du,
			st$dv,
			prior_shape = 2,
			prior_scale = 1,
			reset_stationary = TRUE
		)^2
	)
	# gamma(shape, rate) for the precision
	expect_lt(abs(mean(prec) - shape / rate) / (sqrt(shape) / rate / sqrt(4000)), 4.5)
	expect_lt(abs(var(prec) / (shape / rate^2) - 1), 0.1)
	# with a fixed reset scale the reset positions stay out of the conditional
	set.seed(2)
	prec0 = replicate(
		4000,
		1 / lame:::.lame_sample_sigma_uv_drift_only(
			st$U,
			st$V,
			st$rho,
			0.4,
			FALSE,
			st$du,
			st$dv,
			prior_shape = 2,
			prior_scale = 1
		)^2
	)
	n_reset = sum(st$du[, -1]) * 2 + sum(st$dv[, -1]) * 2
	ss_reset = 0
	for (t in 2:5) ss_reset = ss_reset + sum(st$U[st$du[, t] == 1, , t]^2) + sum(st$V[st$dv[, t] == 1, , t]^2)
	shape0 = shape - n_reset / 2
	rate0 = rate - one_m * ss_reset / 2
	expect_lt(abs(mean(prec0) - shape0 / rate0) / (sqrt(shape0) / rate0 / sqrt(4000)), 4.5)
})

test_that("rho_uv is drawn from its exact conditional with stationary resets", {
	skip_on_cran()
	st = snap_state(3, n = 30, TT = 4)
	sig = st$sigma
	grid = seq(-0.985, 0.985, by = 0.0005)
	lp = dnorm(grid, 0.9, 0.1, log = TRUE)
	for (X in list(list(st$U, st$du), list(st$V, st$dv))) {
		M = X[[1]]; d = X[[2]]
		S = sum(M[, , 1]^2); m = length(M[, , 1])
		for (t in 2:dim(M)[3]) {
			k = d[, t] == 0
			a = sum(M[k, , t] * M[k, , t - 1]); b = sum(M[k, , t - 1]^2)
			lp = lp + (grid * a - 0.5 * grid^2 * b) / sig^2
			S = S + sum(M[!k, , t]^2); m = m + sum(!k) * dim(M)[2]
		}
		lp = lp + 0.5 * m * log(1 - grid^2) - 0.5 * (1 - grid^2) * S / sig^2
	}
	w = exp(lp - max(lp)); w = w / sum(w)
	m_ex = sum(w * grid); sd_ex = sqrt(sum(w * (grid - m_ex)^2))
	set.seed(4)
	N = 6000; r = numeric(N); cur = 0.8
	for (i in seq_len(N)) {
		cur = lame:::.lame_sample_rho_uv_drift_only(
			st$U,
			st$V,
			sig,
			cur,
			FALSE,
			st$du,
			st$dv,
			prior_mean = 0.9,
			prior_sd = 0.1,
			reset_stationary = TRUE
		)
		r[i] = cur
	}
	# an independence sampler with a close proposal: half the draws as the
	# effective size is conservative
	expect_lt(abs(mean(r) - m_ex) / (sd_ex / sqrt(N / 2)), 4.5)
	expect_lt(abs(sd(r) / sd_ex - 1), 0.1)
})

test_that("a snap fit takes the linked scale by default and refuses a bad fixed one", {
	set.seed(6)
	n = 12; TT = 4; nm = sprintf("a%02d", 1:n)
	Y = lapply(1:TT, function(t) { y = matrix(rnorm(n * n), n, n); diag(y) = NA; dimnames(y) = list(nm, nm); y })
	ctrl = list(nscan = 60, burn = 20, odens = 4, seed = 3, gof = FALSE, verbose = FALSE, plot = FALSE)
	f = suppressWarnings(suppressMessages(do.call(
		lame,
		c(
			list(
				Y,
				family = "normal",
				R = 1,
				dynamic_uv = TRUE,
				dynamic_uv_kind = "snap"
			),
			ctrl
		)
	)))
	expect_true(all(is.finite(f$rho_uv)) && all(is.finite(f$sigma_uv)))
	# at rank one the reset rate is reported (see test_dynamic_hyper_scale.R
	# for the rank-two rule)
	sm = summary(f)
	expect_false(isTRUE(sm$pi_snap_unidentified))
	expect_true(all(is.finite(f$pi_snap)))
	expect_equal(sm$pi_snap, mean(f$pi_snap))
	expect_match(paste(capture.output(print(sm)), collapse = " "), "pi_snap =")
	# an ar1 fit reports no such state
	f_ar1 = suppressWarnings(suppressMessages(do.call(
		lame,
		c(
			list(
				Y,
				family = "normal",
				R = 1,
				dynamic_uv = TRUE,
				prior = list(als_start = FALSE)
			),
			ctrl
		)
	)))
	expect_false(isTRUE(summary(f_ar1)$pi_snap_unidentified))
	expect_error(
		suppressMessages(do.call(
			lame,
			c(
				list(
					Y,
					family = "normal",
					R = 1,
					dynamic_uv = TRUE,
					dynamic_uv_kind = "snap",
					prior = list(snap_kappa = -1)
				),
				ctrl
			)
		)),
		class = "lame_snap_kappa"
	)
})
