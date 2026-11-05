library(lame)

# unit equivariance of the dynamic estimators: fitting c * y gives c times the fit of y, and a covariate x -> k x
# rescales its coefficient by 1 / k; the response unit enters priors and starting values only through R/normal_scale.R

ueq_net = function(mode = "dir", n = 20, nB = 15, Tt = 3, R = 0, seed = 21) {
	set.seed(seed)
	bip = mode == "bip"; sym = mode == "sym"
	m = if (bip) nB else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("c%02d", seq_len(m)) else rn
	a = rnorm(n, 0, 0.5); b = if (sym) a else rnorm(m, 0, 0.5)
	U = matrix(rnorm(n * max(R, 1)), n) * 0.7; V = if (sym) U else matrix(rnorm(m * max(R, 1)), m) * 0.7
	Y = X = Xr = Xc = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * m), n, m)
		if (sym) x = (x + t(x)) / 2
		mu = 0.5 + 0.8 * x + outer(a, b, "+") + if (R > 0) tcrossprod(U[, seq_len(R)], V[, seq_len(R)]) else 0
		e = matrix(rnorm(n * m), n, m)
		if (sym) e = (e + t(e)) / sqrt(2)
		y = mu + e
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, m, 1), dimnames = list(rn, cn, "xd"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, bip = bip, sym = sym, R = R)
}
ueq_fit = function(d, sc = 1, kx = 1, ...) {
	args = list(
		lapply(d$Y, function(y) y * sc),
		Xdyad = lapply(d$X, function(x) x * kx),
		family = "normal",
		method = "als",
		R = d$R,
		verbose = FALSE,
		...
	)
	if (d$bip) args$mode = "bipartite"
	if (d$sym) args$symmetric = TRUE
	suppressWarnings(suppressMessages(do.call(lame, args)))
}
ueq_rel = function(got, ref) max(abs(got - ref)) / max(abs(ref))

test_that("the empirical-Bayes dynamic fit of c * y is c times the fit of y, in every kind of network", {
	skip_on_cran()
	for (mode in c("dir", "bip", "sym")) for (spec in list(
		list(dynamic_ab = TRUE),
		list(dynamic_beta = TRUE),
		list(dynamic_uv = TRUE)
	)) {
		R = if (!is.null(spec$dynamic_uv)) 1L else 0L
		d = ueq_net(mode, R = R)
		f1 = do.call(ueq_fit, c(list(d), spec))
		for (cc in c(1e-3, 1e3)) {
			fc = do.call(ueq_fit, c(list(d, cc), spec))
			info = paste(mode, names(spec), cc)
			expect_lt(ueq_rel(unname(fc$coefficients), cc * unname(f1$coefficients)), 1e-5, label = info)
			expect_lt(ueq_rel(unname(fc$VC[["ve"]]), cc^2 * unname(f1$VC[["ve"]])), 1e-5, label = info)
			if (R > 0L) expect_lt(ueq_rel(unlist(fc$UVPM), cc * unlist(f1$UVPM)), 1e-4, label = info)
			if (R == 0L) expect_lt(ueq_rel(unlist(fc$a), cc * unlist(f1$a)), 1e-4, label = info)
		}
	}
})

test_that("the penalised dynamic fit of c * y is c times the fit of y, and the sampler started from it scales the same way", {
	skip_on_cran()
	d = ueq_net("dir", R = 2)
	f1 = ueq_fit(d, bias_correct = "none", dynamic_uv = TRUE, dynamic_ab = TRUE)
	for (cc in c(1e-6, 1e3)) {
		fc = ueq_fit(d, cc, bias_correct = "none", dynamic_uv = TRUE, dynamic_ab = TRUE)
		expect_identical(fc$iterations, f1$iterations)
		expect_lt(ueq_rel(unname(fc$coefficients), cc * unname(f1$coefficients)), 1e-8)
		expect_lt(ueq_rel(unlist(fc$UVPM), cc * unlist(f1$UVPM)), 1e-8)
	}
	mc = function(sc) suppressWarnings(suppressMessages(lame(
		lapply(d$Y, function(y) y * sc),
		Xdyad = d$X,
		family = "normal",
		R = 2,
		dynamic_uv = TRUE,
		nscan = 40,
		burn = 10,
		odens = 2,
		seed = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	m1 = mc(1)
	expect_true(isTRUE(m1$als_start))
	for (cc in c(1e-3, 1e3)) {
		mcc = mc(cc)
		expect_true(isTRUE(mcc$als_start))
		expect_lt(ueq_rel(mcc$BETA, cc * m1$BETA), 1e-6)
		expect_lt(ueq_rel(mcc$VC[, "ve"], cc^2 * m1$VC[, "ve"]), 1e-6)
		expect_lt(ueq_rel(unlist(mcc$UVPM), cc * unlist(m1$UVPM)), 1e-6)
	}
})

test_that("the coefficient priors follow the unit of the covariate: the g-prior of the static coefficients and the initial-state prior of the moving ones", {
	skip_on_cran()
	d = ueq_net("dir")
	for (spec in list(
		list(dynamic_ab = TRUE),
		list(dynamic_beta = TRUE),
		list(dynamic_beta = TRUE, bias_correct = "none")
	)) {
		f1 = do.call(ueq_fit, c(list(d), spec))
		fk = do.call(ueq_fit, c(list(d, 1, 1e3), spec))
		b1 = f1$coefficients
		bk = fk$coefficients
		sc = ifelse(grepl("xd", names(b1)), 1e-3, 1)
		expect_lt(max(abs(bk - b1 * sc) / abs(b1 * sc)), 1e-4, label = paste(names(spec), collapse = "+"))
	}
	# the static estimator and the sampler: the g-prior is invariant
	s1 = suppressWarnings(ame_als(d$Y[[1]], Xdyad = d$X[[1]], family = "normal", R = 0, verbose = FALSE))
	sk = suppressWarnings(ame_als(d$Y[[1]], Xdyad = d$X[[1]] * 1e3, family = "normal", R = 0, verbose = FALSE))
	expect_equal(unname(sk$coefficients[["xd_dyad"]]) * 1e3, unname(s1$coefficients[["xd_dyad"]]), tolerance = 1e-8)
})

test_that("Suv0 and kappa0 are refused for a bipartite network, which has no covariance of the rows of its factors", {
	d = ueq_net("bip", R = 1)
	for (arg in list(list(Suv0 = diag(2)), list(kappa0 = 6))) {
		err = expect_error(
			suppressWarnings(suppressMessages(lame(
				d$Y,
				Xdyad = d$X,
				family = "normal",
				mode = "bipartite",
				R_row = 1,
				R_col = 1,
				prior = arg,
				nscan = 20,
				burn = 10,
				odens = 2,
				verbose = FALSE,
				plot = FALSE,
				gof = FALSE
			))),
			class = "lame_prior_unused"
		)
		expect_identical(err$given, names(arg))
	}
	err = expect_error(
		suppressWarnings(suppressMessages(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = "normal",
			mode = "bipartite",
			R_row = 1,
			R_col = 1,
			prior = list(Suv0 = diag(2)),
			nscan = 20,
			burn = 10,
			odens = 2,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))),
		class = "lame_prior_unused"
	)
	expect_match(conditionMessage(err), "unit-free")
	# a unipartite network keeps its prior
	du = ueq_net("dir", R = 1)
	expect_s3_class(
		suppressWarnings(suppressMessages(ame(
			du$Y[[1]],
			Xdyad = du$X[[1]],
			family = "normal",
			R = 1,
			prior = list(Suv0 = diag(2)),
			nscan = 20,
			burn = 10,
			odens = 2,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))),
		"ame"
	)
})

test_that("the penalised fit of bipartite and undirected positions and static factors is unit-equivariant", {
	skip_on_cran()
	for (spec in list(
		list(mode = "bip", R = 2L, dynamic_uv = TRUE),
		list(mode = "bip", R = 1L, dynamic_ab = TRUE),
		list(mode = "sym", R = 2L, dynamic_uv = TRUE),
		list(mode = "sym", R = 1L, dynamic_ab = TRUE)
	)) {
		d = ueq_net(spec$mode, R = spec$R)
		extra = spec[setdiff(names(spec), c("mode", "R"))]
		f1 = do.call(ueq_fit, c(list(d, bias_correct = "none"), extra))
		for (cc in c(1e-3, 1e3)) {
			fc = do.call(ueq_fit, c(list(d, cc, bias_correct = "none"), extra))
			info = paste(spec$mode, spec$R, names(extra), cc)
			expect_lt(ueq_rel(unname(fc$coefficients), cc * unname(f1$coefficients)), 1e-6, label = info)
			expect_lt(ueq_rel(unlist(fc$UVPM), cc * unlist(f1$UVPM)), 1e-6, label = info)
			expect_lt(ueq_rel(unname(fc$VC[["ve"]]), cc^2 * unname(f1$VC[["ve"]])), 1e-6, label = info)
		}
	}
})

test_that("the sampler with dynamic additive effects scales exactly in the unit of a covariate and of the response", {
	skip_on_cran()
	d = ueq_net("dir", n = 12, nB = 12)
	mc = function(sc = 1, kx = 1) suppressWarnings(suppressMessages(lame(
		lapply(d$Y, function(y) y * sc),
		Xdyad = lapply(d$X, function(x) x * kx),
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 20,
		burn = 0,
		odens = 1,
		seed = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	m1 = mc()
	for (kx in c(1e-3, 1e3)) {
		mk = mc(1, kx)
		slope = grepl("dyad", colnames(m1$BETA))
		expect_lt(ueq_rel(mk$BETA[, slope], m1$BETA[, slope] / kx), 1e-8, label = kx)
		expect_lt(ueq_rel(mk$BETA[, !slope], m1$BETA[, !slope]), 1e-8, label = kx)
		expect_lt(ueq_rel(mk$VC[, "ve"], m1$VC[, "ve"]), 1e-8, label = kx)
	}
	mc3 = mc(1e3)
	expect_lt(ueq_rel(mc3$BETA, 1e3 * m1$BETA), 1e-8)
})
