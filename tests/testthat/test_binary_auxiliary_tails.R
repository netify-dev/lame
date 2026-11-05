test_that("shared binary factors retain the exact conditional near both correlation boundaries", {
	set.seed(462193)
	Z = array(0, c(2, 2, 1)); Z[1, 2, 1] = 1; Z[2, 1, 1] = -.8
	EZ = array(0, dim(Z)); upper = upper.tri(matrix(0, 2, 2)); n = 6000
	for (rho in c(-1 + 1e-8, 1 - 1e-8)) {
		x = replicate(n, lame:::.lame_dyad_factor_draw(Z, EZ, rho, upper)[1, 2, 1])
		a = abs(rho); target_mean = a * (1 + sign(rho) * (-.8)) / (1 + a)
		target_variance = a * (1 - a) / (1 + a)
		expect_lt(abs(mean(x) - target_mean) / sqrt(target_variance / n), 6)
		expect_lt(abs(var(x) / target_variance - 1) / sqrt(2 / (n - 1)), 6)
	}
})

test_that("binary latent refresh uses the sampled conditional variance without clipping", {
	set.seed(648237)
	Z = EZ = off = array(0, c(2, 2, 1));n = 6000
	for (rho in c(-1 + 1e-8, 1 - 1e-8)) {
		x = replicate(
			n,
			lame:::.lame_binary_z_given_factor(Z, EZ, off, rho, integer(), integer(), c(2L, 3L), integer())[2]
		)
		variance = 1 - abs(rho)
		expect_lt(abs(mean(x)) / sqrt(variance / n), 6)
		expect_lt(abs(var(x) / variance - 1) / sqrt(2 / (n - 1)), 6)
	}
	EZ[2] = Inf
	expect_error(
		lame:::.lame_binary_z_given_factor(Z, EZ, off, .5, integer(), integer(), 2L, integer()),
		'Non-finite binary latent-response refresh',
		fixed = TRUE
	)
})
