test_that('tempered gradients and acceptance ratios match independent densities', {
	withr::local_seed(73136)
	for (rank in 1:2) for (kind in c('unit', 'hierarchical')) {
		n = 5L; nt = 3L
		Y = lapply(1:nt, function(t) {
			y = outer(1:n, 1:n, function(i, j) as.numeric((i + j + t) %% 3 == 0))
			diag(y) = NA_real_; y
		})
		input = lame:::.lame_tempered_input(Y, NULL, rank)
		prior = lame:::.lame_tempered_prior(
			list(
				uv_prior = kind,
				Sab0 = diag(.4, 2),
				etaab = 6,
				rho_ab_mean = .5,
				rho_ab_sd = .2,
				sigma_ab_shape = 3,
				sigma_ab_scale = .3
			),
			input$Y
		)
		design = lame:::.lame_tempered_design(input, NA)
		cx = lame:::.lame_tempered_context(design$Y, design$Xlist, rank, design$Q, prior)
		d = n * rank + rank + cx$p + (n - 1L) * (nt + 1L)
		q = sin(seq_len(d + 3L)) / 5
		q[d + 1:3] = c(-1, -1, atanh(.3 / .99))
		theta_index = n * rank + seq_len(rank)
		density = function(v, temperature) {
			z = lame:::.lame_tempered_decode(v, cx, design$Q)
			normal_index = if (kind == 'unit') seq_len(d) else setdiff(seq_len(d), theta_index)
			lp = -.5 * sum(v[normal_index]^2) - cx$av * log(z$va) - cx$bv / z$va -
				2 * cx$as * log(z$sigma) - cx$bs / z$sigma^2 -
				.5 * ((z$rho - cx$rm) / cx$rs)^2 + log(.99 * (1 - (z$rho / .99)^2))
			if (kind == 'hierarchical') for (theta in v[theta_index]) {
				val = integrate(
					function(x) x * exp(-.5 * theta^2 * x^2 - .5 * x),
					0,
					Inf,
					rel.tol = 1e-11
				)$value
				lp = lp + log(.5 * val / sqrt(2 * pi))
			}
			likelihood = 0
			for (t in seq_len(nt)) {
				eta = z$beta[1L] + outer(z$path[, t], z$path[, t], '+') + tcrossprod(z$U, z$V)
				mask = upper.tri(Y[[t]]) & is.finite(Y[[t]])
				likelihood = likelihood + sum(pnorm((2 * Y[[t]][mask] - 1) * eta[mask], log.p = TRUE))
			}
			lp + temperature * likelihood
		}
		temperature = .37
		ev = lame:::tempered_binary_evaluate(q, cx, temperature)
		expect_true(ev$valid)
		expect_equal(ev$logdensity, density(q, temperature), tolerance = 1e-9)
		gradient = vapply(seq_along(q), function(k) {
			delta = 1e-5 * (1 + abs(q[k])); plus = minus = q
			plus[k] = plus[k] + delta; minus[k] = minus[k] - delta
			(density(plus, temperature) - density(minus, temperature)) / (2 * delta)
		}, 0)
		expect_lt(max(abs(drop(ev$gradient) - gradient) / (1 + abs(gradient))), 1e-6)
		proposal = lame:::tempered_binary_transition(q, cx, temperature, .3)
		expect_true(proposal$valid)
		candidate = drop(proposal$candidate)
		reverse = sum(dnorm(q, candidate + drop(proposal$new_drift), drop(proposal$new_sd), log = TRUE))
		forward = sum(dnorm(candidate, q + drop(proposal$drift), drop(proposal$sd), log = TRUE))
		ratio = density(candidate, temperature) - density(q, temperature) + reverse - forward
		expect_equal(proposal$log_ratio, ratio, tolerance = 1e-9)
		reflection = lame:::tempered_binary_reflection(q, cx, temperature, 0L)
		reflected = drop(reflection$proposal)
		expect_equal(reflection$log_ratio, density(reflected, temperature) - density(q, temperature), tolerance = 1e-9)
		inverse = lame:::tempered_binary_reflection(reflected, cx, temperature, 0L)
		expect_identical(drop(inverse$proposal), q)
	}
})
