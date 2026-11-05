# tempered transition of an eigen-direction for a symmetric count fit with dynamic positions (R/uv_tempered_move.R)

# stationarity of the joint law of (weights, paths, latent array, counts) on 6 actors, four periods, a missing pair
# and an absent actor; returns the largest z score of paired differences and the acceptance rate
utm_exact = function(mutant = 'none', ndraw = 4000L, seed = 61L) {
	set.seed(seed)
	n = 6L; R = 2L; N = 4L; rho = .85; s2 = .4; lvar = 1; off = c(0, .3, -.2, .1); r = 1L
	bs = c(1, .75, .5, .25)
	draw = function() {
		L = rnorm(R, 0, sqrt(lvar)); U = lame:::.lame_uv_unit_paths(array(rnorm(n * R * N), c(n, R, N)), rho)
		Z = Y = array(0, c(n, n, N))
		for (t in 1:N) {
			Ut = matrix(U[, , t], n, R); z = Ut %*% diag(L) %*% t(Ut) + matrix(rnorm(n * n, 0, sqrt(s2)), n); z[lower.tri(z)] = t(z)[lower.tri(z)]; Z[, , t] = z
			y = matrix(
				rpois(n * n, exp(z + off[t])),
				n
			); y[lower.tri(y)] = t(y)[lower.tri(y)]; diag(y) = NA; Y[, , t] = y
		}
		Y[1, 2, 1] = Y[2, 1, 1] = NA; Y[3, , 2] = Y[, 3, 2] = NA
		list(L = L, U = U, Z = Z, Y = Y)
	}
	prod_all = function(U, L) { P = array(0, c(n, n, N)); for (t in 1:N) { Ut = matrix(U[, , t], n, R); P[, , t] = Ut %*% diag(L) %*% t(Ut) }; P }
	stat = function(s) { e = s$Z - prod_all(s$U, s$L); c(
		s$L,
		s$L[1]^2,
		s$L[1] * s$L[2],
		sum(s$U[, 1, ]^2),
		sum(s$U[, 1, 1] * s$U[, 1, 2]),
		sum(s$U[, 2, 1] * s$U[, 2, 3]),
		e[1, 2, 1],
		e[3, 4, 3]^2,
		s$Z[2, 3, 2],
		sum(s$Z[, , 4]),
		sum(s$Y[, , 1] * s$Z[, , 1], na.rm = TRUE),
		sign(s$L[1]),
		s$Z[5, 5, 2]
	) }
	S0 = S1 = NULL; acc = 0
	for (k in seq_len(ndraw)) {
		s = draw(); yy = s$Y; for (t in 1:N) yy[, , t][lower.tri(yy[, , t], TRUE)] = NA
		dt = lame:::.lame_temper_data(lame:::.lame_ltm_context(yy, 'poisson', 'sym', offset = off), n, N)
		if (mutant == 'wrongoff') dt$off = dt$off + .5
		o = lame:::.lame_uv_sym_temper(s$U, s$L, s$Z, dt, r, rho, lvar, bs, cycles = 1L)
		s1 = s
		if (o$accepted) { acc = acc + 1; s1$U = o$U; s1$L = o$L; s1$Z = o$Z }
		S0 = rbind(S0, stat(s)); S1 = rbind(S1, stat(s1))
	}
	d = S1 - S0
	sdv = apply(d, 2L, sd)
	list(maxz = max(abs(colMeans(d)[sdv > 0] / (sdv[sdv > 0] / sqrt(ndraw)))), accept = acc / ndraw)
}

test_that('a tempered transition of one direction leaves the joint law of the weights, paths, latent array and counts unchanged', {
	skip_on_cran()
	o = utm_exact()
	expect_lt(o$maxz, 4.2)
	expect_gt(o$accept, .05)
})

test_that('a wrong exposure in the factor is caught', {
	skip_on_cran()
	expect_gt(utm_exact('wrongoff', 3000L)$maxz, 5)
})

test_that('lame() makes the tempered transitions of a symmetric Poisson fit with dynamic positions, and the switch omits them', {
	skip_on_cran()
	set.seed(3)
	n = 10L; N = 4L; nm = sprintf('a%02d', 1:n)
	Y = setNames(
		lapply(
			1:N,
			function(t) { y = matrix(rpois(n * n, .3), n); y[lower.tri(y)] = t(y)[lower.tri(y)]; diag(y) = NA; dimnames(y) = list(nm, nm); y }
		),
		paste0('p', 1:N)
	)
	f = suppressWarnings(suppressMessages(lame(
		Y,
		family = 'poisson',
		symmetric = TRUE,
		R = 2L,
		dynamic_uv = TRUE,
		burn = 10,
		nscan = 40,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 1,
		prior = list(uv_temper_every = 10)
	)))
	expect_equal(unname(f$uv_temper_move['calls']), 5)
	f0 = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = 'poisson',
				symmetric = TRUE,
				R = 2L,
				dynamic_uv = TRUE,
				burn = 10,
				nscan = 40,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(uv_temper_move = FALSE)
			)
		)
	)
	expect_null(f0$uv_temper_move)
	expect_error(
		lame(
			Y,
			family = 'poisson',
			symmetric = TRUE,
			R = 2L,
			dynamic_uv = TRUE,
			burn = 5,
			nscan = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(uv_temper_every = 0)
		),
		'positive integer'
	)
})
