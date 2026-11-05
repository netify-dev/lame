# the coefficient/additive-effect block and the factor columns against independent dense conditionals built cell by cell
skip_on_cran()

# dense conditional of (beta, a, b): each unordered pair is bivariate with covariance s2 [1 rho; rho 1], each diagonal cell
# has variance s2 (1 + rho), (a_i, b_i) ~ N(0, Sab), beta ~ N(m0, P0^-1); centering makes both actor margins sum to zero
.dense_block = function(Z, mX, Sab, rho, s2, P0, m0, center) {
	n = dim(Z)[1]; nt = dim(Z)[3]; p = dim(mX)[2]; rho = rep_len(rho, nt)
	q = p + 2 * n
	L = matrix(0, q, q); h = numeric(q)
	for (t in seq_len(nt)) {
		iS = solve(s2 * matrix(c(1, rho[t], rho[t], 1), 2))
		X = matrix(mX[, , t], n * n, p)
		row = function(i, j) {
			r = numeric(q); if (p) r[seq_len(p)] = X[i + (j - 1) * n, ]
			r[p + i] = 1; r[p + n + j] = 1; r
		}
		for (j in seq_len(n)) for (i in seq_len(j)) {
			if (i < j) {
				D = rbind(row(i, j), row(j, i)); y = c(Z[i, j, t], Z[j, i, t])
				L = L + t(D) %*% iS %*% D; h = h + drop(t(D) %*% iS %*% y)
			} else {
				d = row(i, i); v = s2 * (1 + rho[t])
				L = L + tcrossprod(d) / v; h = h + d * Z[i, i, t] / v
			}
		}
	}
	if (p) {
		L[seq_len(p), seq_len(p)] = L[seq_len(p), seq_len(p)] + P0
		h[seq_len(p)] = h[seq_len(p)] + drop(P0 %*% m0)
	}
	iSab = solve(Sab)
	for (i in seq_len(n)) { k = p + c(i, n + i); L[k, k] = L[k, k] + iSab }
	S = solve(L); mu = drop(S %*% h)
	if (center) {
		C = rbind(c(rep(0, p), rep(1, n), rep(0, n)), c(rep(0, p), rep(0, n), rep(1, n)))
		K = S %*% t(C) %*% solve(C %*% S %*% t(C))
		mu = mu - drop(K %*% C %*% mu); S = S - K %*% C %*% S
	}
	list(mean = mu, covariance = (S + t(S)) / 2)
}

.block_fixture = function(n = 5L, nt = 2L, p = 2L, seed = 1L) {
	set.seed(seed)
	Z = array(rnorm(n * n * nt), c(n, n, nt))
	mX = array(rnorm(n * n * p * nt), c(n * n, p, nt))
	if (p) mX[, 1L, ] = 1
	list(Z = Z, mX = mX, n = n, nt = nt, p = p)
}

.block_moments = function(f, Sab, rho, s2, g, center) {
	lame:::.lame_replicated_gaussian_moments(f$Z, f$mX, Sab, rho, s2, g, center)
}

.block_prior = function(f, g) {
	if (!f$p) return(matrix(0, 0L, 0L))
	lame:::.lame_replicated_prior(f$mX, NULL, g, f$p, f$nt)
}

test_that("block moments match the dense conditional with and without centering", {
	Sab = matrix(c(1, .3, .3, .7), 2)
	for (cfg in list(
		list(p = 2L, nt = 3L, center = FALSE, rho = .4),
		list(p = 2L, nt = 3L, center = TRUE, rho = c(-.3, .5, .1)),
		list(p = 0L, nt = 2L, center = FALSE, rho = .6),
		list(p = 0L, nt = 2L, center = TRUE, rho = -.5),
		list(p = 3L, nt = 1L, center = FALSE, rho = 0)
	)) {
		f = .block_fixture(nt = cfg$nt, p = cfg$p)
		got = .block_moments(f, Sab, cfg$rho, 1.3, 40, cfg$center)
		ref = .dense_block(
			f$Z,
			f$mX,
			Sab,
			cfg$rho,
			1.3,
			.block_prior(f, 40),
			numeric(f$p),
			cfg$center
		)
		scale = max(abs(ref$covariance))
		expect_lt(max(abs(drop(got$mean) - ref$mean)), 1e-10 * max(1, abs(ref$mean)))
		expect_lt(max(abs(got$covariance - ref$covariance)), 1e-10 * scale)
	}
})

test_that("block moments stay exact close to both reciprocity boundaries", {
	f = .block_fixture(n = 4L, nt = 1L, p = 2L, seed = 2L)
	Sab = matrix(c(.8, .25, .25, 1.1), 2)
	for (rho in c(-1 + 1e-6, 1 - 1e-6)) {
		got = .block_moments(f, Sab, rho, 1.7, 16, FALSE)
		ref = .dense_block(f$Z, f$mX, Sab, rho, 1.7, .block_prior(f, 16), numeric(2), FALSE)
		sd = sqrt(diag(ref$covariance))
		expect_lt(max(abs(drop(got$mean) - ref$mean) / sd), 1e-6)
		expect_lt(max(abs(got$covariance - ref$covariance) / outer(sd, sd)), 1e-6)
	}
})

test_that("singular actor covariances keep their supported Gaussian law", {
	f = .block_fixture(n = 5L, nt = 2L, p = 1L, seed = 3L)
	prior = .block_prior(f, 25)
	for (S in list(matrix(1, 2, 2), diag(c(1, 0)), matrix(0, 2, 2))) {
		got = .block_moments(f, S, .3, 1, 25, FALSE)
		expect_equal(got$rank, as.integer(sum(eigen(S)$values > 1e-12)))
		# the singular law is the limit of the regular one
		ref = .dense_block(f$Z, f$mX, S + diag(1e-10, 2), .3, 1, prior, 0, FALSE)
		expect_lt(max(abs(drop(got$mean) - ref$mean)), 1e-6)
		expect_lt(max(abs(got$covariance - ref$covariance)), 1e-6)
	}
	# a disabled sender margin is a structural zero
	set.seed(4)
	x = lame:::.lame_static_actor_draw(f$Z, diag(c(.7, 1.2)), .2, 1, rvar = FALSE)
	expect_identical(x$a, numeric(5))
	expect_true(all(x$b != 0))
})

test_that("coefficient draws are invariant to covariate units", {
	f = .block_fixture(n = 6L, nt = 2L, p = 2L, seed = 5L)
	Sab = matrix(c(.5, .1, .1, .4), 2)
	base = .block_moments(f, Sab, .2, .9, 36, TRUE)
	f$mX[, 2L, ] = 1000 * f$mX[, 2L, ]
	scaled = .block_moments(f, Sab, .2, .9, 36, TRUE)
	u = c(1, 1000, rep(1, 12))
	expect_equal(drop(scaled$mean) * u, drop(base$mean), tolerance = 1e-9)
	expect_equal(scaled$covariance * outer(u, u), base$covariance, tolerance = 1e-9)
})

test_that("block draws reproduce the dense conditional moments", {
	f = .block_fixture(n = 4L, nt = 2L, p = 2L, seed = 6L)
	Sab = matrix(c(1, -.2, -.2, .6), 2)
	ref = .dense_block(
		f$Z,
		f$mX,
		Sab,
		c(.3, -.4),
		1.1,
		.block_prior(f, 16),
		numeric(2),
		TRUE
	)
	set.seed(7)
	N = 4000L
	x = replicate(
		N,
		unlist(lame:::.lame_replicated_gaussian_draw(
			f$Z,
			f$mX,
			Sab,
			c(.3, -.4),
			1.1,
			16,
			TRUE
		))
	)
	sd = sqrt(diag(ref$covariance))
	expect_lt(max(abs(rowMeans(x) - ref$mean) / (sd / sqrt(N))), 4.5)
	expect_lt(max(abs(cov(t(x)) - ref$covariance) / outer(sd, sd)), 0.1)
	# centered margins sum to zero in every draw
	expect_lt(max(abs(colSums(x[2 + 1:4, ])), abs(colSums(x[2 + 5:8, ]))), 1e-12)
})

test_that("rbeta_ab_fc uses the supplied coefficient prior and the dense law", {
	f = .block_fixture(n = 5L, nt = 1L, p = 2L, seed = 8L)
	X = precomputeX(array(f$mX[, , 1L], c(5L, 5L, 2L)))
	Z = f$Z[, , 1L]
	Sab = matrix(c(.6, .2, .2, .9), 2)
	P0 = matrix(c(2, .3, .3, 1), 2); m0 = c(.5, -1)
	ref = .dense_block(f$Z, f$mX, Sab, .35, 1.4, P0, m0, FALSE)
	set.seed(9)
	N = 4000L
	x = replicate(N, unlist(rbeta_ab_fc(Z, Sab, .35, X, 1.4, iV0 = P0, m0 = m0)))
	sd = sqrt(diag(ref$covariance))
	expect_lt(max(abs(rowMeans(x) - ref$mean) / (sd / sqrt(N))), 4.5)
	expect_lt(max(abs(cov(t(x)) - ref$covariance) / outer(sd, sd)), 0.1)
	# the public replicated wrapper returns plain vectors
	z = rbeta_ab_rep_fc(f$Z, Sab, .35, array(f$mX, c(5L, 5L, 2L, 1L)), center = TRUE)
	expect_length(z$beta, 2L); expect_length(z$a, 5L); expect_length(z$b, 5L)
	expect_lt(abs(sum(z$a)) + abs(sum(z$b)), 1e-12)
})

test_that("invalid block inputs fail explicitly", {
	f = .block_fixture(n = 4L, nt = 1L, p = 1L)
	Z = f$Z; Z[1, 2, 1] = NaN
	expect_error(lame:::.lame_replicated_gaussian_draw(Z, f$mX, diag(2), 0, 1, 16), "finite")
	expect_error(lame:::.lame_replicated_gaussian_draw(f$Z, f$mX, diag(2), 1, 1, 16), "between")
	expect_error(lame:::.lame_replicated_gaussian_draw(f$Z, f$mX, diag(c(1, NA)), 0, 1, 16), "finite")
})

# dense conditional of one factor column x given its partner w and the
# partial residual P: transpose = FALSE has mean x_i w_j, TRUE has w_i x_j
.dense_factor = function(P, w, rho, s2, mu0, tau, transpose) {
	n = length(w); nt = dim(P)[3]; L = diag(tau, n); h = tau * mu0
	iS = solve(s2 * matrix(c(1, rho, rho, 1), 2))
	row = function(i, j) {
		r = numeric(n); if (!transpose) r[i] = r[i] + w[j] else r[j] = r[j] + w[i]; r
	}
	for (t in seq_len(nt)) for (j in seq_len(n)) for (i in seq_len(j)) {
		if (i < j) {
			D = rbind(row(i, j), row(j, i)); y = c(P[i, j, t], P[j, i, t])
			L = L + t(D) %*% iS %*% D; h = h + drop(t(D) %*% iS %*% y)
		} else {
			d = row(i, i); v = s2 * (1 + rho)
			L = L + tcrossprod(d) / v; h = h + d * P[i, i, t] / v
		}
	}
	S = solve(L); list(mean = drop(S %*% h), root = t(chol(S)))
}

.standardized = function(x, m) drop(forwardsolve(m$root, x - m$mean))

.expect_standard_normal = function(z) {
	z = do.call(rbind, z)
	expect_lt(max(abs(colMeans(z))), 4.5 / sqrt(nrow(z)))
	expect_lt(max(abs(cov(z) - diag(ncol(z)))), 0.12)
}

test_that("compiled factor columns follow their exact Gaussian conditionals", {
	skip_on_cran()
	set.seed(11)
	n = 5L; nt = 2L; R = 2L; s2 = .8; rho = .45; N = 2500L
	E = array(rnorm(n * n * nt, sd = 2), c(n, n, nt))
	U = matrix(rnorm(n * R), n); V = matrix(rnorm(n * R), n)
	for (r in seq_len(R)) {
		P = E
		for (t in seq_len(nt)) P[, , t] = E[, , t] - U[, -r, drop = FALSE] %*% t(V[, -r, drop = FALSE])
		# without shrinkage the factor prior is N(0, n I)
		mu = .dense_factor(P, V[, r], rho, s2, numeric(n), 1 / n, FALSE)
		zu = zv = vector("list", N)
		for (b in seq_len(N)) {
			x = lame:::rUV_rep_fc_cpp(
				E,
				U,
				V,
				rho,
				s2,
				diag(2),
				1e-6,
				FALSE,
				r - 1L,
				lame:::.lame_uv_prior_df(R),
				lame:::.lame_uv_prior_df(R)
			)
			zu[[b]] = .standardized(x$U[, r], mu)
			mv = .dense_factor(P, x$U[, r], rho, s2, numeric(n), 1 / n, TRUE)
			zv[[b]] = .standardized(x$V[, r], mv)
		}
		.expect_standard_normal(zu)
		.expect_standard_normal(zv)
	}
})

test_that("R factor updates follow their exact conditionals under a correlated prior", {
	skip_on_cran()
	set.seed(12)
	n = 5L; s2 = .8; rho = -.35; N = 2500L
	Suv = matrix(c(1.5, .6, .6, .9), 2); Puv = solve(Suv)
	E = matrix(rnorm(n * n, sd = 2), n); U = matrix(rnorm(n), n); V = matrix(rnorm(n), n)
	mu = .dense_factor(
		array(E, c(n, n, 1L)),
		V[, 1],
		rho,
		s2,
		-Puv[1, 2] / Puv[1, 1] * V[, 1],
		Puv[1, 1],
		FALSE
	)
	zu = zv = vector("list", N)
	for (b in seq_len(N)) {
		x = rUV_fc(E, U, V, Suv, rho, s2)
		zu[[b]] = .standardized(x$U[, 1], mu)
		mv = .dense_factor(
			array(E, c(n, n, 1L)),
			x$U[, 1],
			rho,
			s2,
			-Puv[2, 1] / Puv[2, 2] * x$U[, 1],
			Puv[2, 2],
			TRUE
		)
		zv[[b]] = .standardized(x$V[, 1], mv)
	}
	.expect_standard_normal(zu)
	.expect_standard_normal(zv)
})

test_that("factor prior precision tolerates rounding-level asymmetry", {
	set.seed(13)
	W = matrix(rnorm(40), 10)
	S = solve(rwish(solve(diag(4) + crossprod(W)), 14))
	S[1, 2] = S[1, 2] * (1 + 4 * .Machine$double.eps)
	P = lame:::.lame_factor_prior_precision(S)
	expect_true(isSymmetric(P))
	expect_equal(P %*% ((S + t(S)) / 2), diag(4), tolerance = 1e-10)
	expect_true(isSymmetric(rSuv_fc(W[, 1:2], W[, 3:4])))
})

test_that("directed normal ame fits with two latent dimensions run", {
	set.seed(2)
	n = 12; U = matrix(rnorm(n * 2), n); V = U + matrix(rnorm(n * 2, 0, .3), n)
	Y = U %*% t(V) + matrix(rnorm(n * n), n); diag(Y) = NA
	fit = ame(
		Y * 0.01,
		R = 2,
		nscan = 300,
		burn = 50,
		odens = 10,
		seed = 2,
		verbose = FALSE,
		gof = FALSE
	)
	expect_true(all(is.finite(fit$BETA)))
	expect_true(all(is.finite(fit$U)) && all(is.finite(fit$V)))
})
