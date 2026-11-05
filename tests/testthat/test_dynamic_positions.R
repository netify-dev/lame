# reported positions of dynamic_uv fits factor each period's posterior mean product
skip_on_cran()

rank_part = function(P, R) {
	s = svd(P)
	s$u[, seq_len(R)] %*% diag(s$d[seq_len(R)], R) %*% t(s$v[, seq_len(R)])
}
rel_diff = function(A, B) sqrt(sum((A - B)^2)) / sqrt(sum(B^2))

dyn_panel = function(n, Tn, signal, seed, symmetric = FALSE) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	U0 = matrix(rnorm(n * 2), n)
	truth = if (symmetric) {
		U0 %*% diag(c(1.2, -0.8)) %*% t(U0)
	} else {
		signal * U0 %*% t(matrix(rnorm(n * 2), n))
	}
	dimnames(truth) = list(nm, nm)
	Y = lapply(seq_len(Tn), function(t) {
		e = matrix(rnorm(n * n), n)
		if (symmetric) e = (e + t(e)) / sqrt(2)
		m = truth + e
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	names(Y) = paste0("t", seq_len(Tn))
	list(Y = Y, truth = truth)
}

test_that("directed dynamic positions reproduce the product and its scale", {
	d = dyn_panel(15, 3, 1.5, 11)
	fit = lame(
		d$Y,
		R = 2,
		family = "normal",
		dynamic_uv = TRUE,
		burn = 150,
		nscan = 500,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		plot = FALSE
	)
	off = row(d$truth) != col(d$truth)
	for (t in 1:3) {
		P = fit$U[, , t] %*% t(fit[["V"]][, , t])
		expect_lt(rel_diff(P, rank_part(fit$UVPM[, , t], 2)), 1e-8)
		slope = unname(coef(lm(P[off] ~ d$truth[off]))[2])
		expect_gt(slope, 0.8)
	}
	expect_equal(dim(fit$U_sd), dim(fit$U))
	expect_true(all(fit$U_sd > 0))
	expect_lt(max(fit$U_sd), sd(fit$U))
	lp = latent_positions(fit)
	expect_false(anyNA(lp$posterior_sd))
})

test_that("symmetric dynamic positions keep U L U' and its signs", {
	d = dyn_panel(15, 3, 1, 12, symmetric = TRUE)
	fit = lame(
		d$Y,
		R = 2,
		family = "normal",
		symmetric = TRUE,
		dynamic_uv = TRUE,
		burn = 150,
		nscan = 500,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		plot = FALSE
	)
	expect_equal(sort(sign(fit$L)), c(-1, 1))
	off = row(d$truth) != col(d$truth)
	for (t in 1:3) {
		P = fit$U[, , t] %*% diag(fit$L) %*% t(fit$U[, , t])
		expect_lt(rel_diff(P, fit$ULUPM[, , t]), 0.05)
		expect_gt(cor(P[off], d$truth[off]), 0.9)
	}
	aligned = procrustes_align(fit)
	expect_equal(
		aligned$U[, , 3] %*% diag(fit$L) %*% t(aligned$U[, , 3]),
		fit$U[, , 3] %*% diag(fit$L) %*% t(fit$U[, , 3])
	)
	expect_equal(dim(fit$U_sd), dim(fit$U))
})

test_that("bipartite dynamic positions reproduce U G V' with a fixed diagonal G", {
	set.seed(13)
	ra = sprintf("r%02d", 1:12)
	cb = sprintf("c%02d", 1:8)
	truth = 1.5 * matrix(rnorm(24), 12) %*% t(matrix(rnorm(16), 8))
	Y = lapply(1:3, function(t) {
		m = truth + matrix(rnorm(96), 12)
		dimnames(m) = list(ra, cb)
		m
	})
	names(Y) = paste0("t", 1:3)
	fit = lame(
		Y,
		R = 2,
		family = "normal",
		mode = "bipartite",
		dynamic_uv = TRUE,
		burn = 150,
		nscan = 500,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		plot = FALSE
	)
	expect_equal(fit$G[1, 2], 0)
	expect_equal(fit$G[2, 1], 0)
	for (t in 1:3) {
		P = fit$U[, , t] %*% fit$G %*% t(fit[["V"]][, , t])
		expect_lt(rel_diff(P, rank_part(fit$UVPM[, , t], 2)), 1e-8)
		expect_gt(cor(c(P), c(truth)), cor(c(fit$UVPM[, , t]), c(truth)) - 0.02)
	}
	expect_equal(dim(fit$V_sd), dim(fit[["V"]]))
})

test_that("actors absent from a period have no position there", {
	d = dyn_panel(14, 3, 1.5, 14)
	nm = rownames(d$truth)
	d$Y[[2]] = d$Y[[2]][nm[-(1:2)], nm[-(1:2)]]
	fit = lame(
		d$Y,
		R = 2,
		family = "normal",
		dynamic_uv = TRUE,
		burn = 100,
		nscan = 300,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		plot = FALSE
	)
	absent = c("a01", "a02")
	expect_true(all(is.na(fit$U[absent, , 2])))
	expect_true(all(is.na(fit$U_sd[absent, , 2])))
	expect_false(anyNA(fit$U[, , 1]))
	sims = simulate(fit, nsim = 1, type = "approximate")
	expect_length(sims$Y, 1)
})

test_that("pooled chains refactor the pooled product", {
	d = dyn_panel(12, 3, 1.5, 15)
	fits = lapply(c(1, 2), function(k) {
		lame(
			d$Y,
			R = 2,
			family = "normal",
			dynamic_uv = TRUE,
			burn = 100,
			nscan = 300,
			odens = 5,
			chains = 1,
			seed = 100 * k,
			verbose = FALSE,
			plot = FALSE
		)
	})
	pooled = suppressWarnings(suppressMessages(combine_ame_chains(fits, diagnostics = FALSE)))
	for (t in 1:3) {
		P = pooled$U[, , t] %*% t(pooled[["V"]][, , t])
		expect_lt(rel_diff(P, rank_part(pooled$UVPM[, , t], 2)), 1e-8)
	}
	expect_equal(dim(pooled$U_sd), dim(pooled$U))
})

test_that("an early stop trims the bipartite G and static position draws", {
	set.seed(16)
	Y = lapply(1:2, function(t) {
		m = matrix(rnorm(80), 10)
		dimnames(m) = list(sprintf("r%02d", 1:10), sprintf("c%02d", 1:8))
		m
	})
	names(Y) = c("t1", "t2")
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		R = 1,
		family = "normal",
		mode = "bipartite",
		burn = 10,
		nscan = 200000,
		odens = 1,
		chains = 1,
		max_seconds = 1,
		posterior_opts = list(save_UV = TRUE),
		verbose = FALSE,
		plot = FALSE
	)))
	expect_true(isTRUE(fit$terminated_early))
	S = nrow(fit$BETA)
	expect_lt(S, 200000)
	expect_equal(dim(fit$G_samples)[3], S)
	expect_equal(dim(fit$U_samples)[3], S)
	expect_false(anyNA(fit$G_samples))
})
