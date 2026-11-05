# independent stationary AR and dyad-pair Gaussian design, including actor constraints
skip_on_cran()

.actor_reciprocity_private291 = function(name) get(name, environment(lame::lame), inherits = TRUE)

.actor_reciprocity_reference291 = function(n, nt, rho, s2, phi, sigma, center, rvar, cvar, partial = FALSE) {
	H = if (center) qr.Q(qr(contr.helmert(n))) else diag(n);k = ncol(H)
	active = c(if (rvar) 'a', if (cvar) 'b');dimension = k * nt * length(active)
	E = diag(nt);if (nt > 1) E[cbind(2:nt, 1:(nt - 1))] = -phi
	weight = rep(1 / sigma^2, nt);weight[1] = weight[1] * (1 - phi^2)
	Q = kronecker(diag(length(active)), kronecker(crossprod(E, weight * E), diag(k)))
	rhs = numeric(dimension);Z = array(rnorm(n * n * nt), c(n, n, nt))
	for (t in 1:nt) {
		diag(Z[, , t]) = NA_real_
		if (partial) {Z[1, 2, t] = NA_real_;if (t == 2L) Z[2, 3, t] = Z[3, 2, t] = NA_real_}
		ij = which(is.finite(Z[, , t]), arr.ind = TRUE);m = nrow(ij)
		covariance = diag(s2, m);D = matrix(0, m, dimension)
		for (j in seq_len(m)) {
			i = ij[j, 1];b = ij[j, 2];reverse = which(ij[, 1] == b & ij[, 2] == i)
			if (length(reverse)) covariance[j, reverse] = s2 * rho[t]
			for (a in seq_along(active)) D[j, (a - 1L) * k * nt + (t - 1L) * k + seq_len(k)] = H[if (active[a] == 'a') i else b, ]
		}
		W = solve(covariance);Q = Q + crossprod(D, W %*% D);rhs = rhs + drop(crossprod(D, W %*% Z[, , t][ij]))
	}
	V = solve(Q);mu = drop(V %*% rhs)
	unpack = function(theta) {
		out = list(a = matrix(0, n, nt), b = matrix(0, n, nt))
		for (j in seq_along(active)) out[[active[j]]] = H %*% matrix(theta[(j - 1L) * k * nt + seq_len(k * nt)], k, nt)
		out
	}
	pack = function(out) unlist(lapply(active, function(a) as.numeric(crossprod(H, out[[a]]))))
	list(Z = Z, EZ = array(0, dim(Z)), Q = Q, V = V, mu = mu, H = H, unpack = unpack, pack = pack)
}

test_that("period reciprocity actor sweep preserves its full stationary Gaussian law", {
	set.seed(47329129)
	cx = .actor_reciprocity_reference291(3L, 3L, c(-.8, .1, .7), .9, .6, .7, TRUE, TRUE, TRUE)
	M = 3000L;p = length(cx$mu);draws = matrix(NA_real_, M, p);rootV = chol(cx$V);center_error = 0
	sampler = .actor_reciprocity_private291(".lame_sample_dynamic_ab_period")
	for (j in seq_len(M)) {
		state = cx$unpack(cx$mu + drop(t(rootV) %*% rnorm(p)))
		out = sampler(state$a, state$b, cx$Z, cx$EZ, .6, .7, .9, FALSE, c(-.8, .1, .7), TRUE)
		draws[j, ] = cx$pack(out)
		center_error = max(center_error, abs(colSums(out$a)), abs(colSums(out$b)))
	}
	expect_lt(center_error, 1e-10)
	expect_lt(max(abs((colMeans(draws) - cx$mu) / sqrt(diag(cx$V) / M))), 5.5)
	expect_lt(max(abs((cov(draws) - cx$V) / sqrt((cx$V^2 + outer(diag(cx$V), diag(cx$V))) / (M - 1L)))), 5.5)
})

test_that("constant reciprocity vectors preserve scalar actor objects and RNG", {
	sampler = .actor_reciprocity_private291(".lame_sample_dynamic_ab_period")
	old = .actor_reciprocity_private291("sample_dynamic_ab_cpp")
	set.seed(47329141);a = matrix(rnorm(9), 3);b = matrix(rnorm(9), 3);Z = array(rnorm(27), c(3, 3, 3));EZ = Z * 0
	for (center in c(FALSE, TRUE)) for (rvar in c(FALSE, TRUE)) for (cvar in c(FALSE, TRUE)) for (rho in c(0, .4, -.6)) {
		set.seed(47329142);x = old(a, b, Z, EZ, .6, .7, .9, FALSE, rho, center, rvar, cvar);rng = .Random.seed
		set.seed(47329142);y = sampler(a, b, Z, EZ, .6, .7, .9, FALSE, rep(rho, 3), center, rvar, cvar)
		expect_identical(x, y);expect_identical(rng, .Random.seed)
	}
})

test_that("dynamic actor period correlations fail explicitly before consuming RNG", {
	sampler = .actor_reciprocity_private291(".lame_sample_dynamic_ab_period")
	set.seed(47329151);a = b = matrix(0, 3, 3);Z = EZ = array(0, c(3, 3, 3))
	for (rho in list(c(0, 1, .2), c(0, NA, .2), c(0, .2))) {
		before = .Random.seed
		expect_error(sampler(a, b, Z, EZ, .6, .7, .9, FALSE, rho, TRUE), "Residual correlations")
		expect_identical(before, .Random.seed)
	}
})


test_that("finite residual subtraction overflow refuses before every period draw", {
	a = b = matrix(0, 2, 2);Z = EZ = array(0, c(2, 2, 2))
	for (name in c("sample_dynamic_ab_cpp", "sample_dynamic_ab_rho_cpp")) {
		native = .actor_reciprocity_private291(name)
		rho = if (name == "sample_dynamic_ab_cpp") -.2 else c(-.2, .3)
		for (period in c(1L, 2L)) {
			z = Z;ez = EZ;z[1, 2, period] = 1e308;ez[1, 2, period] = -1e308
			set.seed(47329171);before = .Random.seed
			expect_error(
				native(a, b, z, ez, 0, .2, 1e308, FALSE, rho, FALSE, TRUE, FALSE),
				"Finite dynamic actor residual"
			)
			expect_identical(before, .Random.seed)
		}
		z = Z;z[1, 2, 2] = NA_real_
		set.seed(47329172);out = native(a, b, z, EZ, 0, .2, 1, FALSE, rho, FALSE, TRUE, FALSE)
		expect_true(all(is.finite(c(out$a, out$b))))
	}
})
