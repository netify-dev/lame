skip_on_cran()

test_that("circular conditional draws match independent Bessel moments", {
	skip_on_cran()
	set.seed(641177)
	for (k in c(0, 0.5, 3, 50)) {
		x = replicate(6000, lame:::.lame_uv_von_mises(k))
		truth = if (k == 0) 0 else besselI(k, 1, TRUE) / besselI(k, 0, TRUE)
		expect_lt(abs(mean(cos(x)) - truth), 6 * sd(cos(x)) / sqrt(length(x)))
		expect_lt(abs(mean(sin(x))), 6 * sd(sin(x)) / sqrt(length(x)))
	}
})

test_that("orthogonal conditional agrees with numerical Haar integration", {
	M = matrix(c(0.7, -0.4, 1.2, -0.3), 2)
	angles = seq(-pi, pi, length.out = 8193)[-8193]
	grid = lapply(c(1, -1), function(d) {
		q = vapply(angles, function(a) c(cos(a), sin(a), -d * sin(a), d * cos(a)), numeric(4))
		list(q = q, weight = exp(colSums(q * as.numeric(M))))
	})
	den = sum(vapply(grid, function(x) sum(x$weight), numeric(1)))
	truth = Reduce(`+`, lapply(grid, function(x) rowSums(sweep(x$q, 2, x$weight, "*")))) / den
	p_positive = sum(grid[[1]]$weight) / den
	set.seed(921567)
	draws = replicate(7000, as.numeric(lame:::.lame_uv_o2(M)))
	se = apply(draws, 1, sd) / sqrt(ncol(draws))
	expect_true(all(abs(rowMeans(draws) - truth) < 6 * se))
	determinants = draws[1, ] * draws[4, ] - draws[2, ] * draws[3, ]
	expect_equal(abs(determinants), rep(1, ncol(draws)), tolerance = 1e-12)
	expect_lt(
		abs(mean(determinants > 0) - p_positive),
		6 * sqrt(p_positive * (1 - p_positive) / ncol(draws))
	)
})

test_that("orientation fields reproduce the full AR prior density", {
	set.seed(8224)
	U = array(rnorm(40), c(5, 2, 4))
	V = array(rnorm(40), c(5, 2, 4))
	rho = -0.6
	sigma = 0.7
	prior = function(U, V) {
		energy = function(A) (1 - rho^2) * sum(A[, , 1]^2) +
			sum((A[, , -1] - rho * A[, , -dim(A)[3]])^2)
		-(energy(U) + energy(V)) / (2 * sigma^2)
	}
	for (t in 1:4) {
		a = b = matrix(0, 5, 2)
		if (t > 1) { a = a + U[, , t - 1]; b = b + V[, , t - 1] }
		if (t < 4) { a = a + U[, , t + 1]; b = b + V[, , t + 1] }
		M = rho / sigma^2 * (crossprod(U[, , t], a) + crossprod(V[, , t], b))
		Q = lame:::.lame_uv_o2(matrix(rnorm(4), 2))
		A = U; B = V
		A[, , t] = U[, , t] %*% Q; B[, , t] = V[, , t] %*% Q
		expect_equal(prior(A, B) - prior(U, V), sum((Q - diag(2)) * M), tolerance = 1e-11)
	}
	aligned = lame:::.lame_uv_common_orientation(list(U = U, V = V))
	expect_equal(prior(aligned$U, aligned$V), prior(U, V), tolerance = 1e-11)
})

test_that("orientation updates preserve every product and singleton dimensions", {
	set.seed(66127)
	for (n in c(1L, 6L)) for (nt in c(1L, 4L)) {
		U = array(rnorm(n * 2 * nt), c(n, 2, nt))
		V = array(rnorm(n * 2 * nt), c(n, 2, nt))
		z = lame:::.lame_uv_orientation_ar1(list(U = U, V = V), 0.6, 0.7)
		expect_identical(dim(z$U), dim(U)); expect_identical(dim(z$V), dim(V))
		for (t in seq_len(nt)) {
			expect_equal(
				tcrossprod(matrix(z$U[, , t], n, 2), matrix(z$V[, , t], n, 2)),
				tcrossprod(matrix(U[, , t], n, 2), matrix(V[, , t], n, 2)),
				tolerance = 1e-11
			)
		}
	}
})

test_that("orientation updates retain prior cross-period moments", {
	skip_on_cran()
	set.seed(647191)
	for (rho in c(-0.7, 0, 0.7)) {
		differences = replicate(1000, {
			U = V = array(0, c(2, 2, 3))
			U[, , 1] = rnorm(4, 0, 0.6 / sqrt(1 - rho^2))
			V[, , 1] = rnorm(4, 0, 0.6 / sqrt(1 - rho^2))
			for (t in 2:3) {
				U[, , t] = rho * U[, , t - 1] + rnorm(4, 0, 0.6)
				V[, , t] = rho * V[, , t - 1] + rnorm(4, 0, 0.6)
			}
			z = lame:::.lame_uv_orientation_ar1(list(U = U, V = V), rho, 0.6)
			lag = function(U, V) mean(c(U[, , -1] * U[, , -3], V[, , -1] * V[, , -3]))
			lag(z$U, z$V) - lag(U, V)
		})
		expect_lt(abs(mean(differences)), 6 * sd(differences) / sqrt(length(differences)))
	}
})

test_that("bypassed and extreme-concentration moves do not consume randomness", {
	UV = list(U = array(1, c(2, 2, 2)), V = array(2, c(2, 2, 2)))
	set.seed(97516)
	before = .Random.seed
	expect_identical(lame:::.lame_uv_orientation_ar1(UV, 0.5, 1, enabled = FALSE), UV)
	expect_identical(.Random.seed, before)
	M = matrix(c(1e9, 2e9, -3e9, 4e9), 2)
	expect_identical(lame:::.lame_uv_o2(M), diag(2))
	expect_identical(.Random.seed, before)
})


test_that("the extra update respects the existing coordinate alignment option", {
	set.seed(561144)
	UV = list(U = array(rnorm(48), c(6, 2, 4)), V = array(rnorm(48), c(6, 2, 4)))
	set.seed(921861)
	plain = lame:::.lame_uv_orientation_ar1(UV, 0.6, 0.7, align = FALSE)
	rng = .Random.seed
	set.seed(921861)
	aligned = lame:::.lame_uv_orientation_ar1(UV, 0.6, 0.7, align = TRUE)
	expect_identical(.Random.seed, rng)
	expect_equal(aligned, lame:::.lame_uv_common_orientation(plain), tolerance = 1e-12)
	expect_false(isTRUE(all.equal(plain$U, aligned$U)))
})


test_that("global shear draws match an independent conditional density", {
	set.seed(419812)
	U = array(rnorm(24), c(4, 2, 3))
	V = array(rnorm(24), c(4, 2, 3))
	rho = -0.4
	sigma = 0.6
	log_density = function(q) {
		u = U; v = V
		u[, 2, ] = u[, 2, ] + q * u[, 1, ]
		v[, 1, ] = v[, 1, ] - q * v[, 2, ]
		ans = 0
		for (x in list(u, v)) {
			ans = ans + sum(dnorm(x[, , 1], sd = sigma / sqrt(1 - rho^2), log = TRUE))
			for (t in 2:3) ans = ans + sum(dnorm(x[, , t], rho * x[, , t - 1], sigma, log = TRUE))
		}
		ans
	}
	precision = 2 * log_density(0) - log_density(1) - log_density(-1)
	truth_mean = (log_density(1) - log_density(-1)) / (2 * precision)
	truth_var = 1 / precision
	# the second shear leaves column two of U unchanged, so its first displacement is measured directly
	draws = replicate(4000, {
		z = lame:::.lame_uv_shear_ar1(list(U = U, V = V), rho, sigma)
		sum((z$U[, 2, ] - U[, 2, ]) * U[, 1, ]) / sum(U[, 1, ]^2)
	})
	expect_lt(abs(mean(draws) - truth_mean), 6 * sqrt(truth_var / length(draws)))
	expect_lt(abs(var(draws) - truth_var), 6 * truth_var * sqrt(2 / (length(draws) - 1)))
})

test_that("global shears preserve centered factor constraints", {
	set.seed(324183)
	U = array(rnorm(48), c(6, 2, 4))
	V = array(rnorm(48), c(6, 2, 4))
	for (t in 1:4) {
		U[, , t] = scale(U[, , t], center = TRUE, scale = FALSE)
		V[, , t] = scale(V[, , t], center = TRUE, scale = FALSE)
	}
	z = lame:::.lame_uv_shear_ar1(list(U = U, V = V), 0.7, 0.5)
	expect_lt(max(abs(apply(z$U, c(2, 3), sum))), 1e-12)
	expect_lt(max(abs(apply(z$V, c(2, 3), sum))), 1e-12)
	for (t in 1:4) expect_equal(
		tcrossprod(z$U[, , t], z$V[, , t]),
		tcrossprod(U[, , t], V[, , t]),
		tolerance = 1e-12
	)
})
