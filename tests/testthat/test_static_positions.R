# reported positions of static fits: the factors of the posterior mean product, with
# posterior sds from draws matched onto them (R/static_positions.R)

bip_data = function(nA = 60, nB = 10, R = 2, family = "binary", seed = 1) {
	set.seed(seed)
	U = matrix(rnorm(nA * R), nA, R)
	V = sweep(matrix(rnorm(nB * R), nB, R), 2, seq(1.4, 0.8, length = R), "*")
	eta = -0.3 + U %*% t(V)
	Y = if (family == "binary") (eta + rnorm(nA * nB) > 0) * 1 else eta + rnorm(nA * nB)
	dimnames(Y) = list(paste0("p", seq_len(nA)), paste0("s", seq_len(nB)))
	Y
}

# the best rank-k part of a product
trunc_k = function(P, k) {
	s = svd(P, nu = k, nv = k)
	s$u %*% diag(s$d[seq_len(k)], k) %*% t(s$v)
}

fit_bip = function(Y, family = "binary", seed = 1, ...) {
	ame(
		Y,
		mode = "bipartite",
		family = family,
		burn = 50,
		nscan = 300,
		odens = 5,
		verbose = FALSE,
		chains = 1,
		seed = seed,
		posterior_opts = posterior_options(save_UV = TRUE, thin_UV = 1),
		...
	)
}

test_that("static bipartite positions factor UVPM at the fitted rank", {
	skip_on_cran()
	for (rk in list(c(2, 2), c(3, 2))) {
		fit = suppressWarnings(fit_bip(bip_data(R = 3), R_row = rk[1], R_col = rk[2]))
		k = min(rk)
		expect_equal(dim(fit$U), c(60L, rk[1]))
		expect_equal(dim(fit$V), c(10L, rk[2]))
		expect_equal(unname(fit$G), diag(1, rk[1], rk[2]))
		P = fit$U %*% fit$G %*% t(fit$V)
		expect_lt(max(abs(P - trunc_k(fit$UVPM, k))) / max(abs(fit$UVPM)), 1e-8)
		expect_equal(rownames(fit$U), rownames(fit$UVPM))
		# columns in decreasing order of the singular values, the scale split evenly
		d = colSums(fit$U[, seq_len(k), drop = FALSE]^2)
		expect_equal(d, colSums(fit$V[, seq_len(k), drop = FALSE]^2))
		expect_true(all(diff(d) <= 0))
		if (rk[1] > k) expect_true(all(fit$U[, (k + 1):rk[1]] == 0))
		# sign rule: the positive entries of each column carry more of its squared length
		x = sweep(rbind(fit$U[, seq_len(k), drop = FALSE], fit$V[, seq_len(k), drop = FALSE]), 2, sqrt(d), "/")
		expect_true(all(colSums(x * abs(x)) >= 0))
		# the sampler's own state stays the start of a continuation
		expect_equal(dim(fit$start_vals$G), as.integer(rk))
		expect_false(isTRUE(all.equal(fit$start_vals$U, fit$U)))
	}
})

test_that("two seeds of a static bipartite fit report the same positions", {
	skip_on_cran()
	Y = bip_data(nA = 80, family = "normal")
	f1 = fit_bip(Y, family = "normal", seed = 1, R = 2)
	f2 = fit_bip(Y, family = "normal", seed = 2, R = 2)
	expect_gt(cor(c(f1$U), c(f2$U)), 0.98)
	expect_gt(cor(c(f1$V), c(f2$V)), 0.98)
})

test_that("posterior sds of static positions come from draws matched onto them", {
	skip_on_cran()
	fit = fit_bip(bip_data(nA = 80, family = "normal"), family = "normal", R = 2)
	lp = latent_positions(fit)
	expect_true(all(is.finite(lp$posterior_sd)))
	expect_equal(dim(fit$U_sd), dim(fit$U))
	expect_equal(dim(fit$V_sd), dim(fit$V))
	# the raw draws carry an arbitrary rotation each; their spread is far larger
	raw_sd = apply(fit$U_samples, c(1, 2), sd)
	expect_lt(median(fit$U_sd), median(raw_sd) / 2)
	# the matching depends on each draw's product only: re-gauged draws U A, A^-1 G B^-T, V B
	# give the same sds
	set.seed(9)
	gauged = fit
	for (s in seq_len(dim(fit$U_samples)[3])) {
		A = matrix(rnorm(4), 2) + diag(2)
		B = matrix(rnorm(4), 2) + diag(2)
		gauged$U_samples[, , s] = fit$U_samples[, , s] %*% A
		gauged$V_samples[, , s] = fit$V_samples[, , s] %*% B
		gauged$G_samples[, , s] = solve(A) %*% fit$G_samples[, , s] %*% t(solve(B))
	}
	sds = lame:::.uv_static_sd(gauged, fit$U, fit$V, NULL, "bipartite")
	expect_lt(max(abs(sds$U_sd - fit$U_sd)), 1e-8)
	expect_lt(max(abs(sds$V_sd - fit$V_sd)), 1e-8)
})

test_that("unipartite and symmetric static positions follow the same rule", {
	skip_on_cran()
	set.seed(3)
	n = 25
	U = matrix(rnorm(n * 2), n, 2)
	E = matrix(rnorm(n * n), n, n)
	Y = U %*% t(U) + (E + t(E)) / 2
	diag(Y) = NA
	Yd = U %*% t(U) + E
	diag(Yd) = NA
	po = posterior_options(save_UV = TRUE, thin_UV = 1)
	fd = ame(Yd, R = 2, burn = 50, nscan = 300, odens = 5, verbose = FALSE, chains = 1, seed = 1, posterior_opts = po)
	fs = ame(Y, R = 2, symmetric = TRUE, burn = 50, nscan = 300, odens = 5, verbose = FALSE, chains = 1, seed = 1, posterior_opts = po)
	expect_lt(max(abs(fd$U %*% t(fd$V) - trunc_k(fd$UVPM, 2))) / max(abs(fd$UVPM)), 1e-8)
	expect_true(all(diff(abs(diag(fs$L))) <= 0))
	expect_lt(max(abs(fs$U %*% fs$L %*% t(fs$U) - trunc_k(fs$ULUPM, 2))) / max(abs(fs$ULUPM)), 1e-8)
	expect_true(all(is.finite(fd$U_sd)) && all(is.finite(fs$U_sd)))
	expect_lt(median(fd$U_sd), median(apply(fd$U_samples, c(1, 2), sd)))
	# the symmetric chain state keeps V = U L, so a continuation keeps its scale
	expect_false(isTRUE(all.equal(fs$start_vals$V, fs$start_vals$U)))
})

test_that("a fit saved with last-draw positions reads like a new one", {
	skip_on_cran()
	fit = fit_bip(bip_data(), R = 2)
	old = fit
	n_s = dim(fit$U_samples)[3]
	old$U = fit$U_samples[, , n_s]
	old$V = fit$V_samples[, , n_s]
	old$G = fit$G_samples[, , n_s]
	old$U_sd = old$V_sd = old$uv_positions = NULL
	expect_equal(latent_positions(old), latent_positions(fit))
	expect_equal(suppressMessages(procrustes_align(old))$U, fit$U)
})

test_that("static bipartite lame() reports one matrix of positions", {
	skip_on_cran()
	Yl = lapply(1:3, function(t) bip_data(nA = 30, nB = 8, seed = t))
	names(Yl) = paste0("t", 1:3)
	fit = lame(
		Yl,
		mode = "bipartite",
		R = 2,
		family = "binary",
		burn = 50,
		nscan = 200,
		odens = 5,
		verbose = FALSE,
		chains = 1,
		seed = 1,
		posterior_opts = posterior_options(save_UV = TRUE, thin_UV = 1)
	)
	expect_true(is.matrix(fit$U))
	expect_equal(rownames(fit$U), rownames(fit$UVPM))
	expect_lt(max(abs(fit$U %*% fit$G %*% t(fit$V) - trunc_k(fit$UVPM, 2))) / max(abs(fit$UVPM)), 1e-8)
	lp = latent_positions(fit)
	expect_equal(unique(lp$time), "1")
	expect_true(all(is.finite(lp$posterior_sd)))
})

# small fits of every static type, for the format and strength checks
static_fits = function() {
	set.seed(11)
	n = 16
	U = matrix(rnorm(n * 2), n, 2)
	E = matrix(rnorm(n * n), n, n)
	Yd = 2 * U %*% t(matrix(rnorm(n * 2), n, 2)) + E
	Ys = 2 * U %*% diag(c(1, -1)) %*% t(U) + (E + t(E)) / 2
	diag(Yd) = diag(Ys) = NA
	dimnames(Yd) = dimnames(Ys) = list(paste0("a", 1:n), paste0("a", 1:n))
	Yb = bip_data(nA = 20, nB = 8, family = "normal")
	panel = function(Y) {
		out = lapply(1:2, function(t) {
			E = matrix(rnorm(length(Y), 0, 0.3), nrow(Y))
			if (isSymmetric(unname(Y))) E = (E + t(E)) / 2
			Y + E
		})
		names(out) = c("t1", "t2")
		out
	}
	args = list(burn = 20, nscan = 100, odens = 5, verbose = FALSE, gof = FALSE, chains = 1, seed = 1,
		posterior_opts = posterior_options(save_UV = TRUE, thin_UV = 1))
	run = function(f, ...) suppressWarnings(suppressMessages(do.call(f, c(list(...), args), quote = TRUE)))
	list(
		ame_directed = run("ame", Yd, R = 2),
		ame_symmetric = run("ame", Ys, R = 2, symmetric = TRUE),
		ame_bipartite = run("ame", Yb, R = 2, mode = "bipartite"),
		lame_directed = run("lame", panel(Yd), R = 2),
		lame_symmetric = run("lame", panel(Ys), R = 2, symmetric = TRUE),
		lame_bipartite = run("lame", panel(Yb), R = 2, mode = "bipartite"),
		ame_bipartite_32 = run("ame", Yb, R_row = 3, R_col = 2, mode = "bipartite")
	)
}

shape = function(x) if (is.null(x)) "NULL" else paste0(class(x)[1], ":", paste(if (is.null(dim(x))) length(x) else dim(x), collapse = "x"))

test_that("static fits keep the U, V, G and L formats of lame 1.4.2", {
	skip_on_cran()
	fits = static_fits()
	got = vapply(fits[1:6], function(f) paste(shape(f[["U"]]), shape(f[["V"]]), shape(f[["G"]]), shape(f[["L"]])), character(1))
	# as 1.4.2 returned them, except static bipartite lame U and V (then [n, R, T] cubes)
	expect_equal(unname(got), c(
		"matrix:16x2 matrix:16x2 NULL NULL",
		"matrix:16x2 NULL NULL matrix:2x2",
		"matrix:20x2 matrix:8x2 matrix:2x2 NULL",
		"matrix:16x2 matrix:16x2 NULL NULL",
		"matrix:16x2 NULL NULL matrix:2x2",
		"matrix:20x2 matrix:8x2 matrix:2x2 NULL"
	))
	# a fit saved with L as a vector reads back with the diagonal matrix
	old = fits$lame_symmetric
	old$L = diag(old$L)
	old$uv_positions = NULL
	expect_equal(lame:::.lame_positions(old)$L, fits$lame_symmetric$L)
	expect_true(all(diff(abs(diag(fits$lame_symmetric$L))) <= 0))
})

test_that("unequal bipartite ranks list and print the rank of the product", {
	skip_on_cran()
	f = static_fits()$ame_bipartite_32
	expect_equal(dim(f$U), c(20L, 3L))
	expect_equal(dim(f$V), c(8L, 2L))
	lp = latent_positions(f)
	expect_equal(sort(unique(lp$dimension)), 1:2)
	expect_equal(nrow(lp), (20 + 8) * 2)
	expect_true(any(grepl("rank 2 (R_row = 3, R_col = 2)", capture.output(print(f)), fixed = TRUE)))
	expect_true(any(grepl("rank 2 (R_row = 3, R_col = 2)", capture.output(print(summary(f))), fixed = TRUE)))
	expect_equal(nrow(lame:::.lame_dim_strength(f)), 2L)
})

test_that("dimension strengths are the singular values of the reported product", {
	skip_on_cran()
	fits = static_fits()
	for (nm in names(fits)) {
		f = fits[[nm]]
		st = lame:::.lame_dim_strength(f)
		P = if (isTRUE(f$symmetric)) f$U %*% f$L %*% t(f$U) else if (!is.null(f[["G"]])) f$U %*% f$G %*% t(f$V) else f$U %*% t(f$V)
		d = if (isTRUE(f$symmetric)) sort(abs(eigen(P, symmetric = TRUE)$values), decreasing = TRUE)[1:2] else svd(P)$d[1:2]
		expect_equal(st$strength, d, tolerance = 1e-10, info = nm)
		expect_equal(st$share, d^2 / sum(d^2), tolerance = 1e-10, info = nm)
		expect_named(st, c("dimension", "strength", "share", "tied_with_next"))
		expect_true(is.logical(st$tied_with_next) && !anyNA(st$tied_with_next), info = nm)
		expect_identical(attr(suppressMessages(latent_positions(f)), "strength"), st, info = nm)
	}
	out = capture.output(print(summary(fits$ame_directed)))
	pct = paste0(formatC(100 * lame:::.lame_dim_strength(fits$ame_directed)$share, format = "f", digits = 0), "%", collapse = ", ")
	expect_true(any(grepl(paste0("Latent dimensions: ", pct, " of the latent product"), out, fixed = TRUE)))
	expect_true(any(grepl("Latent dimensions:", capture.output(print(summary(fits$lame_bipartite))), fixed = TRUE)))
})

test_that("static ALS fits report strengths of their factors without a tie flag", {
	skip_on_cran()
	set.seed(2)
	n = 20
	Y = 2 * matrix(rnorm(n * 2), n) %*% t(matrix(rnorm(n * 2), n)) + matrix(rnorm(n * n), n)
	diag(Y) = NA
	f = ame_als(Y, R = 2, verbose = FALSE, seed = 1)
	st = lame:::.lame_dim_strength(f)
	d = svd(f$U %*% t(f$V))$d[1:2]
	expect_equal(sort(st$strength, decreasing = TRUE), d, tolerance = 1e-10)
	expect_true(all(is.na(st$tied_with_next)))
	expect_identical(attr(latent_positions(f), "strength"), st)
	expect_true(any(grepl("Latent dimensions:", capture.output(print(summary(f))), fixed = TRUE)))
})

test_that("near-equal strengths are flagged as tied", {
	# the gap against the draw-to-draw spread of the gap
	set.seed(4)
	D = cbind(rnorm(200, 50, 2), rnorm(200, 48, 2), rnorm(200, 10, 1))
	st = lame:::.uv_strength_table(c(50, 49.5, 10), D)
	expect_equal(st$tied_with_next, c(TRUE, FALSE, FALSE))
	expect_equal(sum(st$share), 1)
	expect_false(lame:::.uv_strength_table(c(60, 40, 10), D)$tied_with_next[1])
	# symmetric dimensions of opposite sign keep identified eigenvectors
	expect_false(lame:::.uv_strength_table(c(50, 49.5), D[, 1:2], L = c(50, -49.5))$tied_with_next[1])
	# no draws, no flag
	expect_true(all(is.na(lame:::.uv_strength_table(c(50, 49.5))$tied_with_next)))
	expect_match(lame:::.lame_tie_sentence(st), "Dimensions 1 and 2 have nearly equal strength, so their two axes are not well separated; read them together as a plane.", fixed = TRUE)
	expect_null(lame:::.lame_tie_sentence(lame:::.uv_strength_table(c(60, 40, 10), D)))
})

test_that("a tied fit says so in summary() and latent_positions()", {
	skip_on_cran()
	f = static_fits()$ame_bipartite
	f$latent_strength$tied_with_next[1] = TRUE
	expect_message(lp <- latent_positions(f), "read them together as a plane")
	expect_true(attr(lp, "strength")$tied_with_next[1])
	expect_true(any(grepl("nearly equal strength", capture.output(print(summary(f))), fixed = TRUE)))
})
