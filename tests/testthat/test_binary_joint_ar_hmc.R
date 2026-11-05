skip_on_cran()

.binary_joint_ar_fixture = function() {
	set.seed(921865);n = 3L;nt = 2L;rank = 2L
	Y = array(NA_real_, c(n, n, nt));Y[cbind(
		c(1, 2, 2, 3, 1, 3),
		c(2, 1, 3, 2, 3, 1),
		c(1, 1, 1, 1, 2, 2)
	)] = c(1, 0, 0, 0, 0, 1);Y[2, 1, 2] = 0
	X = lapply(seq_len(nt), function(t) {x = array(rnorm(n * n * 2), c(n, n, 2));x[, , 1] = 1;x})
	cx = .lame_binary_dynamic_pair_hmc_context(Y, X, rank);P = matrix(c(3, .2, .2, 2), 2);Sab = matrix(c(.2, .03, .03, .15), 2)
	h = .lame_binary_dynamic_ar_hmc_hyper(cx, P, Sab, .6, .25, -.4, 5, .4, .7, .1)
	K = function(r, s) s * s / (1 - r * r) * r^abs(outer(seq_len(nt), seq_len(nt), '-'))
	ia = solve(K(.6, .25));il = solve(Sab)
	# independent density in original coordinates with every log determinant, prior and Jacobian term;
	# pair probabilities by conditional quadrature
	density = function(q, likelihood = TRUE) {z = .lame_binary_dynamic_ar_hmc_decode(q, cx, h);v = z$sigma^2;lp = 0
		if (likelihood) for (t in seq_len(nt)) {
			eta = matrix(matrix(X[[t]], n * n, 2) %*% z$beta, n, n) + outer(
				z$level[, 1] + z$adev[, t],
				z$level[, 2] + z$bdev[, t],
				'+'
			) + z$U[, , t] %*% t(z$V[, , t])
			for (i in 1:(n - 1)) for (j in (i + 1):n) {yy = c(Y[i, j, t], Y[j, i, t]);mu = c(eta[i, j], eta[j, i]);ok = is.finite(yy)
				if (all(ok)) {sg = 2 * yy - 1;rr = h$rho * prod(sg);aa = sg * mu;pr = integrate(
					function(x) dnorm(x) * pnorm((aa[2] - rr * x) / sqrt(1 - rr^2)),
					-Inf,
					aa[1],
					rel.tol = 1e-11
				)$value;lp = lp + log(pr)} else if (any(ok)) lp = lp + pnorm((2 * yy[ok] - 1) * mu[ok], log.p = TRUE)
			}
		}
		ku = K(z$rho_uv, z$sigma);iu = solve(ku);u = matrix(z$U, ncol = nt);vv = matrix(z$V, ncol = nt)
		gaussian = -.5 * (drop(z$beta %*% P %*% z$beta) + sum(z$level * (z$level %*% il)) + sum(z$adev * (z$adev %*% ia)) + sum(z$bdev * (z$bdev %*% ia)) + sum(u * (u %*% iu)) + sum(vv * (vv %*% iu))) - n * rank * as.numeric(
			determinant(ku, logarithm = TRUE)$modulus
		)
		jacobian = n * rank * as.numeric(determinant(ku, logarithm = TRUE)$modulus) + log(v)
		lp + gaussian - (h$shape + 1) * log(v) - h$scale / v + jacobian + dnorm(
			z$rho_uv,
			h$prior_mean,
			h$prior_sd,
			log = TRUE
		) + log(h$limit * (1 - (z$rho_uv / h$limit)^2))
	}
	list(cx = cx, h = h, density = density)
}

test_that('joint stationary-scale and persistence HMC matches original-coordinate density', {
	f = .binary_joint_ar_fixture();cx = f$cx;h = f$h;density = f$density
	set.seed(921865);q = c(
		rnorm(cx$size, sd = .3),
		log(.3 / sqrt(1 - .7^2)),
		atanh(.7 / .99)
	);q2 = q + rnorm(length(q), sd = .05)
	ev = .lame_binary_dynamic_ar_hmc_evaluate(q, cx, h)
	expect_equal(
		ev$logp - .lame_binary_dynamic_ar_hmc_evaluate(q2, cx, h)$logp,
		density(q) - density(q2),
		tolerance = 1e-10
	)
	for (rr in c(-.85, 0, .7, .989)) {
		q = c(
			rnorm(cx$size, sd = .3),
			log(.4),
			atanh(rr / .99)
		);ev = .lame_binary_dynamic_ar_hmc_evaluate(q, cx, h);eps = 1e-6
		grad = vapply(
			seq_along(q),
			function(j) {p = z = q;p[j] = p[j] + eps;z[j] = z[j] - eps;(density(p) - density(z)) / (2 * eps)},
			numeric(1)
		)
		expect_lt(max(abs(grad - ev$gradient) / (1 + abs(grad))), 1e-6)
		z = .lame_binary_dynamic_ar_hmc_decode(q, cx, h)
		expect_equal(
			.lame_binary_dynamic_ar_hmc_encode(z$beta, z$level, z$adev, z$bdev, z$U, z$V, z$sigma, z$rho_uv, cx, h),
			q,
			tolerance = 1e-10
		)
	}
	q = c(
		rnorm(cx$size, sd = .3),
		log(.4),
		atanh(.7 / .99)
	);mass = .lame_binary_dynamic_ar_hmc_mass(cx, h);p = rnorm(length(q)) * sqrt(mass)
	a = .lame_binary_dynamic_ar_hmc_integrate(q, p, cx, h, mass, .12, 9)
	b = .lame_binary_dynamic_ar_hmc_integrate(a$q, -a$p, cx, h, mass, .12, 9)
	expect_equal(b$q, q, tolerance = 1e-10);expect_equal(b$p, -p, tolerance = 1e-10)
})

test_that('joint AR trajectories preserve bounded and unbounded posterior starts', {
	skip_on_cran()
	f = .binary_joint_ar_fixture();cx = f$cx;h = f$h;density = f$density
	stat = function(q) {z = .lame_binary_dynamic_ar_hmc_decode(q, cx, h);c(
		z$beta,
		z$level[1, 1],
		z$adev[1, cx$nt],
		z$U[1, 1, cx$nt],
		atan(sum(z$U[1, , 2] * z$V[2, , 2])),
		z$sigma,
		z$rho_uv,
		q[cx$size + 1],
		sum(q[seq_len(cx$size)]^2) / cx$size
	)}
	set.seed(283716);N = 5000L
	for (bounded in c(FALSE, TRUE)) {
		h$variance_max = if (bounded) .2 else Inf;h$stationary_max = if (bounded) .4 else Inf;h$max_abs = if (bounded) .9 else Inf
		weighted = replicate(N, {
			repeat {rr = qnorm(
				runif(1, pnorm(-h$limit, h$prior_mean, h$prior_sd), pnorm(h$limit, h$prior_mean, h$prior_sd)),
				h$prior_mean,
				h$prior_sd
			);vv = 1 / rgamma(1, h$shape, rate = h$scale);q = c(
				rnorm(cx$size),
				.5 * log(vv / (1 - rr^2)),
				atanh(rr / h$limit)
			);if (is.finite(.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h)$logp)) break}
			lp = density(q) - density(q, FALSE);z = .lame_binary_dynamic_ar_hmc_step(q, cx, h, .18, L = 7)
			c(lp, stat(z$q) - stat(q), numerical = z$nonfinite)
		});wt = exp(
			weighted[1, ] - max(weighted[1, ])
		);d = weighted[2:11, , drop = FALSE];est = drop(d %*% wt) / sum(wt);se = sqrt(
			rowSums(sweep(d - est, 2, wt, '*')^2)
		) / sum(wt)
		expect_lt(max(abs(est / se)), 6);expect_gt(sum(wt)^2 / sum(wt^2), 100);expect_equal(sum(weighted[12, ]), 0)
	}
})

test_that('exact AR endpoints retain other conditional updates without a Fisher transform', {
	# the persistence goes through the block draw of persistence and stationary variance
	set.seed(316825);n = 6L;Y = lapply(
		1:2,
		function(t) {y = matrix(
			rbinom(n * n, 1, .3),
			n,
			n,
			dimnames = list(as.character(1:n), as.character(1:n))
		);diag(y) = NA;y[1, 2] = 1;y[2, 1] = 0;y}
	)
	for (endpoint in c(-.99, .99)) {
		local_mocked_bindings(
			.lame_uv_rho_w_step = function(U_cube, V_cube, rho, sigma, ...) list(rho = endpoint, sigma = sigma),
			.package = 'lame'
		)
		f = suppressMessages(lame(
			Y,
			family = 'binary',
			R = 2L,
			dynamic_uv = TRUE,
			dynamic_ab = TRUE,
			prior = list(uv_var_cap = Inf, uv_max_abs = 50),
			burn = 0L,
			nscan = 20L,
			odens = 5L,
			seed = 28361L,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		))
		expect_true(all(f$rho_uv == endpoint));expect_null(f$dynamic_marginal_move)
		expect_true(all(f$sampler_health$numerical_failures == 0));expect_gt(stats::sd(f$BETA[, 1]), 0)
	}
})
