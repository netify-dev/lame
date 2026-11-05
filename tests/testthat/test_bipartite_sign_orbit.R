test_that("bipartite sign heat bath matches the independently enumerated AR prior", {
	skip_on_cran()
	set.seed(627149)
	nt = 4L; rho = .43; sigma = .8
	U = array(rnorm(3 * 2 * nt), c(3, 2, nt))
	V = array(rnorm(4 * 3 * nt), c(4, 3, nt))
	G = matrix(rnorm(6), 2, 3)
	energy = function(x) {
		x = matrix(x, ncol = nt)
		(1 - rho^2) * sum(x[, 1]^2) + sum((x[, -1, drop = FALSE] - rho * x[, -nt, drop = FALSE])^2)
	}
	paths = cbind(1, as.matrix(expand.grid(rep(list(c(-1, 1)), nt - 1L))))
	lp = apply(
		paths,
		1,
		function(s) -(energy(sweep(U, 3, s, "*")) +
				energy(sweep(V, 3, s, "*"))) / (2 * sigma^2)
	)
	prob = exp(lp - max(lp)); prob = prob / sum(prob)
	N = 6000L; counts = integer(nrow(paths)); product_error = 0
	for (j in seq_len(N)) {
		z = lame:::.lame_uv_sign_ar1(list(U = U, V = V, bound_rejections = 2L), rho, sigma)
		signs = sign(z$U[1, 1, ] / U[1, 1, ])
		k = which(rowSums(abs(sweep(paths, 2, signs, "-"))) == 0)
		counts[k] = counts[k] + 1L
		if (j <= 5L) for (tt in seq_len(nt)) product_error = max(
			product_error,
			max(abs(U[, , tt] %*% G %*% t(V[, , tt]) - z$U[, , tt] %*% G %*% t(z$V[, , tt])))
		)
	}
	expect_lt(max(abs(counts / N - prob) / sqrt(pmax(prob * (1 - prob) / N, 1 / N^2))), 6)
	expect_equal(product_error, 0)
	expect_identical(z$bound_rejections, 2L)
	expect_equal(abs(z$U), abs(U))
	expect_equal(abs(z$V), abs(V))
})

test_that("the bipartite dynamic sampler invokes the sign conditional", {
	set.seed(82819)
	Y = lapply(seq_len(3), function(tt) matrix(rnorm(30), 5, 6))
	original = lame:::.lame_uv_sign_ar1
	calls = 0L
	wrapped = function(...) { calls <<- calls + 1L; original(...) }
	testthat::with_mocked_bindings({
		fit = lame(
			Y,
			mode = "bipartite",
			family = "normal",
			R = 1,
			dynamic_uv = TRUE,
			burn = 2,
			nscan = 4,
			odens = 1,
			seed = 994,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
		expect_equal(calls, 6L)
		expect_true(all(is.finite(fit$rho_uv)))
	}, .lame_uv_sign_ar1 = wrapped, .package = "lame")
})

test_that("weighted sign draws match reset and Student orbit probabilities", {
	skip_on_cran()
	set.seed(70231)
	nt = 4L; rho = .48; sigma = .8
	U = array(rnorm(2 * 2 * nt, sd = .5), c(2, 2, nt))
	V = array(rnorm(3 * 2 * nt, sd = .5), c(3, 2, nt))
	paths = cbind(1, as.matrix(expand.grid(rep(list(c(-1, 1)), nt - 1L))))
	for (kind in c("snap", "t")) {
		wu = matrix(if (kind == "snap") rbinom(2 * nt, 1, .7) else rgamma(2 * nt, 2, 2), 2, nt)
		wv = matrix(if (kind == "snap") rbinom(3 * nt, 1, .7) else rgamma(3 * nt, 2, 2), 3, nt)
		energy = function(X, weights) {
			total = 0
			for (tt in 2:nt) for (actor in seq_len(dim(X)[1])) {
				total = total + weights[actor, tt] * sum((X[actor, , tt] - rho * X[actor, , tt - 1L])^2)
			}
			total
		}
		lp = apply(
			paths,
			1,
			function(s) -(energy(sweep(U, 3, s, "*"), wu) +
					energy(sweep(V, 3, s, "*"), wv)) / (2 * sigma^2)
		)
		prob = exp(lp - max(lp)); prob = prob / sum(prob)
		N = 5000L; counts = integer(nrow(paths))
		for (j in seq_len(N)) {
			z = lame:::.lame_uv_sign_ar1(list(U = U, V = V), rho, sigma, wu, wv)
			s = sign(z$U[1, 1, ] / U[1, 1, ])
			k = which(rowSums(abs(sweep(paths, 2, s, "-"))) == 0)
			counts[k] = counts[k] + 1L
		}
		expect_lt(max(abs(counts / N - prob) / sqrt(pmax(prob * (1 - prob) / N, 1 / N^2))), 6)
		expect_equal(abs(z$U), abs(U))
		expect_equal(abs(z$V), abs(V))
		for (tt in seq_len(nt)) expect_equal(tcrossprod(z$U[, , tt], z$V[, , tt]), tcrossprod(U[, , tt], V[, , tt]))
	}
})

test_that("singular-direction signs preserve the exact bounded and unbounded target", {
	skip_on_cran()
	set.seed(638917)
	nt = 3L; rho = .4; sigma = .8
	U = array(rnorm(2 * 2 * nt, sd = .5), c(2, 2, nt))
	V = array(rnorm(3 * 2 * nt, sd = .5), c(3, 2, nt))
	U[, , 1] = U[, , 1] * .1; V[, , 1] = V[, , 1] * .1
	G = matrix(c(1, .7, -.3, .5), 2)
	svd_g = svd(G)
	paths = as.matrix(expand.grid(rep(list(c(-1, 1)), 2 * (nt - 1))))
	states = lapply(seq_len(nrow(paths)), function(i) {
		u = U; v = V; signs = cbind(1, matrix(paths[i, ], 2, nt - 1))
		for (tt in seq_len(nt)) {
			u[, , tt] = U[, , tt] %*% svd_g$u %*% diag(signs[, tt]) %*% t(svd_g$u)
			v[, , tt] = V[, , tt] %*% svd_g$v %*% diag(signs[, tt]) %*% t(svd_g$v)
		}
		list(U = u, V = v)
	})
	energy = function(x) {
		x = matrix(x, ncol = nt)
		(1 - rho^2) * sum(x[, 1]^2) + sum((x[, -1] - rho * x[, -nt])^2)
	}
	lp = vapply(states, function(z) -(energy(z$U) + energy(z$V)) / (2 * sigma^2), numeric(1))
	limits = vapply(states, function(z) max(abs(z$U), abs(z$V)), numeric(1))
	bound = mean(range(limits))
	expect_true(any(limits > bound) && any(limits < bound))
	base = lapply(1:nt, function(tt) U[, , tt] %*% G %*% t(V[, , tt]))
	for (cap in c(Inf, bound)) {
		p = exp(lp - max(lp)); p[limits > cap] = 0; p = p / sum(p)
		N = 6000L; counts = integer(length(states))
		for (i in seq_len(N)) {
			z = lame:::.lame_uv_svd_sign_ar1(states[[sample.int(length(states), 1, prob = p)]], G, rho, sigma, cap)
			ss = sapply(2:nt, function(tt) sign((z$U[, , tt] %*% svd_g$u)[1, ] / (U[, , tt] %*% svd_g$u)[1, ]))
			index = which(rowSums(abs(sweep(paths, 2, c(ss), "-"))) == 0)
			counts[index] = counts[index] + 1L
			if (i <= 3L) {
				for (tt in seq_len(nt)) expect_equal(
					z$U[, , tt] %*% G %*% t(z$V[, , tt]),
					base[[tt]],
					tolerance = 1e-12
				)
				expect_lte(max(abs(z$U), abs(z$V)), cap)
			}
		}
		expect_lt(max(abs(counts / N - p) / sqrt(pmax(p * (1 - p) / N, 1 / N^2))), 6)
	}
})
