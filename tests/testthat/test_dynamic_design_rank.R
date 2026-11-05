skip_on_cran()

library(lame)

# the dynamic ALS estimators (empirical-Bayes for the three families, and penalised with bias_correct = "none") refuse an
# exactly rank-deficient design over the observed cells with the classed error of ame() and lame() (R/design_rank.R,
# tested in test_design_rank.R); near-collinearity is fitted

dr_panel = function(family = "normal", n = 12, Tt = 3, seed = 4) {
	set.seed(seed)
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	a = rnorm(n, 0, 0.5)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = array(rnorm(n * n * 2), c(n, n, 2), dimnames = list(rn, rn, c("x1", "x2")))
		eta = -0.2 + outer(a, a, "+") + 0.4 * x[, , 1]
		y = switch(
			family,
			normal = eta + rnorm(n * n),
			poisson = matrix(rpois(n * n, exp(eta)), n),
			binary = 1 * (eta + rnorm(n * n) > 0)
		)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
		X[[t]] = x
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}
squash = function(e) gsub("\\s+", " ", conditionMessage(e))

test_that("the dynamic estimators refuse an aliased covariate, empirical-Bayes and penalised, for each family", {
	for (fam in c("normal", "binary", "poisson")) {
		d = dr_panel(fam)
		# the second covariate is a multiple of the first
		X = lapply(d$X, function(x) { x[, , 2] = -3 * x[, , 1]; x })
		for (bc in list(NULL, "none")) {
			args = list(d$Y, X, family = fam, dynamic_ab = TRUE, verbose = FALSE)
			if (!is.null(bc)) args$bias_correct = bc
			info = paste(fam, if (is.null(bc)) "eb" else bc)
			err = expect_error(do.call(lame:::lame_dynamic_als, args), class = "lame_design_rank", info = info)
			expect_identical(err$aliased, "x2_dyad", info = info)
			expect_identical(err$rank, 2L, info = info)
			expect_identical(err$n_col, 3L, info = info)
			expect_match(squash(err), "`x2_dyad` is aliased", fixed = TRUE, info = info)
			expect_match(squash(err), "x2_dyad = -3 x1_dyad", fixed = TRUE, info = info)
			# nothing precedes the refusal: no warning, no fit
			expect_no_warning(try(do.call(lame:::lame_dynamic_als, args), silent = TRUE))
		}
	}
	# the front end refuses it as well
	d = dr_panel("normal")
	X = lapply(d$X, function(x) { x[, , 2] = x[, , 1] + 1; x })
	err = expect_error(
		lame(d$Y, Xdyad = X, family = "normal", method = "als", dynamic_ab = TRUE, verbose = FALSE),
		class = "lame_design_rank"
	)
	expect_identical(err$aliased, "x2_dyad")
	expect_match(squash(err), "x2_dyad = 1 intercept + 1 x1_dyad", fixed = TRUE)
	# a design of full rank, with a nearly collinear covariate, is fitted
	X = lapply(d$X, function(x) { x[, , 2] = x[, , 1] + 1e-3 * rnorm(length(x[, , 1])); x })
	fit = lame:::lame_dynamic_als(d$Y, X, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
	expect_true(fit$converged)
})

test_that("a constant covariate beside the intercept is aliased with it", {
	d = dr_panel("poisson")
	X = lapply(d$X, function(x) { x[, , 2] = 5; x })
	for (bc in list(NULL, "none")) {
		args = list(d$Y, X, family = "poisson", dynamic_ab = TRUE, verbose = FALSE)
		if (!is.null(bc)) args$bias_correct = bc
		err = expect_error(do.call(lame:::lame_dynamic_als, args), class = "lame_design_rank")
		expect_identical(err$aliased, "x2_dyad")
		expect_match(squash(err), "x2_dyad = 5 intercept", fixed = TRUE)
	}
})

test_that("aliased node covariates are refused as well, whether or not their coefficients move", {
	d = dr_panel("normal", n = 10, Tt = 3, seed = 6)
	rn = rownames(d$Y[[1]])
	set.seed(8)
	w = rnorm(10)
	Xr = lapply(seq_along(d$Y), function(t) { m = cbind(w1 = w, w2 = 2 * w); rownames(m) = rn; m })
	names(Xr) = names(d$Y)
	for (bc in list(NULL, "none")) {
		args = list(d$Y, d$X[seq_along(d$Y)], Xrow = Xr, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
		if (!is.null(bc)) args$bias_correct = bc
		err = expect_error(do.call(lame:::lame_dynamic_als, args), class = "lame_design_rank")
		expect_identical(err$aliased, "w2_row")
		expect_match(squash(err), "w2_row = 2 w1_row", fixed = TRUE)
	}
	# independent node covariates are fitted
	Xr2 = lapply(Xr, function(m) { m[, "w2"] = rnorm(10); m })
	fit = lame:::lame_dynamic_als(
		d$Y,
		d$X[seq_along(d$Y)],
		Xrow = Xr2,
		family = "normal",
		dynamic_ab = TRUE,
		verbose = FALSE
	)
	expect_true(fit$converged)
})

test_that("the rank is that of the observed cells: columns that differ only on missing cells are aliased", {
	d = dr_panel("normal", n = 12, Tt = 3, seed = 9)
	hide = lapply(d$Y, function(y) { h = is.na(y); h[2, 3] = h[5, 1] = TRUE; h })
	Y = Map(function(y, h) { y[h] = NA; y }, d$Y, hide)
	# the second covariate repeats the first wherever the response is observed, and differs only on the hidden cells
	X = Map(function(x, h) { x[, , 2] = x[, , 1]; x[, , 2][h] = x[, , 2][h] + 7; x }, d$X, hide)
	for (bc in list(NULL, "none")) {
		args = list(Y, X, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
		if (!is.null(bc)) args$bias_correct = bc
		err = expect_error(do.call(lame:::lame_dynamic_als, args), class = "lame_design_rank")
		expect_identical(err$aliased, "x2_dyad")
	}
	# one observed cell where they differ separates them
	X[[2]][4, 6, 2] = X[[2]][4, 6, 1] + 7
	fit = lame:::lame_dynamic_als(Y, X, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
	expect_true(fit$converged)
})

test_that("bipartite and symmetric designs are checked on what the fit uses", {
	set.seed(11)
	nA = 9; nB = 6; Tt = 3
	ra = paste0("r", 1:nA); cb = paste0("c", 1:nB)
	Y = lapply(1:Tt, function(t) matrix(rnorm(nA * nB), nA, nB, dimnames = list(ra, cb)))
	names(Y) = paste0("t", 1:Tt)
	w = rnorm(nA); v = rnorm(nB)
	Xr = lapply(1:Tt, function(t) { m = cbind(w1 = w, w2 = 3 * w); rownames(m) = ra; m }); names(Xr) = names(Y)
	Xc = lapply(1:Tt, function(t) { m = cbind(v1 = v); rownames(m) = cb; m }); names(Xc) = names(Y)
	for (bc in list(NULL, "none")) {
		args = list(
			Y,
			Xrow = Xr,
			Xcol = Xc,
			family = "normal",
			mode = "bipartite",
			R_row = 0,
			R_col = 0,
			dynamic_ab = TRUE,
			verbose = FALSE
		)
		if (!is.null(bc)) args$bias_correct = bc
		err = expect_error(do.call(lame:::lame_dynamic_als, args), class = "lame_design_rank")
		expect_identical(err$aliased, "w2_row")
		expect_match(squash(err), "On the 162 observed cells of 3 periods", fixed = TRUE)
	}
	Xr[] = lapply(Xr, function(m) { m[, "w2"] = rnorm(nA); m })
	fit = lame:::lame_dynamic_als(
		Y,
		Xrow = Xr,
		Xcol = Xc,
		family = "normal",
		mode = "bipartite",
		R_row = 0,
		R_col = 0,
		dynamic_ab = TRUE,
		verbose = FALSE
	)
	expect_true(fit$converged)
	# symmetric: the node covariates enter as the sum over the two actors of a pair
	n = 10; rn = paste0("a", 1:n)
	Ys = lapply(1:Tt, function(t) {
		m = matrix(rnorm(n * n), n, n, dimnames = list(rn, rn)); m[lower.tri(m)] = t(m)[lower.tri(m)]; diag(m) = NA; m
	})
	names(Ys) = paste0("t", 1:Tt)
	ws = rnorm(n)
	Xs = lapply(1:Tt, function(t) { m = cbind(s1 = ws, s2 = -2 * ws); rownames(m) = rn; m }); names(Xs) = names(Ys)
	err = expect_error(
		lame:::lame_dynamic_als(
			Ys,
			Xrow = Xs,
			family = "normal",
			symmetric = TRUE,
			dynamic_ab = TRUE,
			verbose = FALSE
		),
		class = "lame_design_rank"
	)
	expect_identical(err$aliased, "s2_row")
	expect_match(squash(err), "s2_row = -2 s1_row", fixed = TRUE)
})
