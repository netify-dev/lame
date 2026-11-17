library(lame)

# every estimate of a dynamic empirical-Bayes ALS fit comes with a standard error: coefficients and their
# paths, the additive effects, their paths and levels, the variance components, the innovation scales of
# moving coefficients and the cells of the latent product; the AR(1) smoothing settings (rho_ab, sigma_ab,
# rho_beta, rho_uv, sigma_uv) stay on the fit and are not printed; a withheld static binary intercept is
# fitted as a path instead

dse_panel = function(family = "normal", n = 14, Tt = 3, strength = 0.6, seed = 7, mode = "unipartite") {
	set.seed(seed)
	nc = if (mode == "bipartite") 11 else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (mode == "bipartite") sprintf("c%02d", seq_len(nc)) else rn
	U = matrix(rnorm(n, 0, strength), n, 1); V = matrix(rnorm(nc, 0, strength), nc, 1)
	al = rnorm(n, 0, 0.5); be = rnorm(nc, 0, 0.5)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		al = al + rnorm(n, 0, 0.1); be = be + rnorm(nc, 0, 0.1)
		x = matrix(rnorm(n * nc), n, nc)
		eta = 0.3 + outer(al, be, "+") + 0.6 * x + U %*% t(V)
		y = switch(family,
			normal = eta + matrix(rnorm(n * nc), n, nc),
			binary = 1 * (eta + matrix(rnorm(n * nc), n, nc) > 0),
			poisson = matrix(rpois(n * nc, exp(0.5 * eta - 0.5)), n, nc))
		if (mode == "unipartite") diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, nc, 1), dimnames = list(rn, cn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

dse_fit = function(d, family = "normal", ...) {
	suppressWarnings(suppressMessages(lame(d$Y, Xdyad = d$X, family = family, method = "als", verbose = FALSE, ...)))
}

fin_pos = function(x) length(x) > 0 && all(is.finite(x)) && all(x > 0)

test_that("a dynamic_ab fit reports errors for its effects, paths and variance components", {
	skip_on_cran()
	d = dse_panel("normal")
	fit = dse_fit(d, R = 0, dynamic_ab = TRUE)
	expect_true(fin_pos(fit$eb$se))
	expect_identical(dim(fit$a_dynamic_sd), dim(fit$a_dynamic))
	expect_identical(names(fit$a_sd), names(fit$a))
	for (nm in c("a_sd", "b_sd", "a_dynamic_sd", "b_dynamic_sd", "a_level_sd", "b_level_sd")) expect_true(fin_pos(fit[[nm]]), info = nm)
	# the AR(1) law of the additive paths is a smoothing setting: no error, no table row
	expect_false(any(c("rho_ab", "sigma_ab") %in% names(fit$eb$hyper_se)))
	s = summary(fit)
	expect_true(all(c("Std.Error", "CI 2.5%", "CI 97.5%") %in% colnames(s$coefficients)))
	expect_true(all(c("va", "cab", "vb", "ve") %in% rownames(s$vc_table)))
	expect_true(fin_pos(s$vc_table[, "Std.Error"]))
	# the intervals of the variances stay positive
	expect_true(all(s$vc_table[c("va", "vb", "ve"), "CI 2.5%"] > 0))
	txt = c(capture.output(print(fit)), capture.output(print(s)))
	expect_false(any(grepl("Point estimate only", txt)))
})

test_that("a dynamic_beta fit reports its coefficient paths with errors, and the innovation scale of each block", {
	skip_on_cran()
	d = dse_panel("normal", seed = 9)
	fit = dse_fit(d, R = 0, dynamic_beta = "dyad")
	s = summary(fit)
	expect_true(fin_pos(fit$eb$se_path["x_dyad", ]))
	expect_true(all(paste0("x_dyad[", fit$time_names, "]") %in% rownames(s$coef_path)))
	expect_true(fin_pos(s$coef_path[, "Std.Error"]))
	expect_true("sigma_beta.dyad" %in% names(fit$eb$hyper_se))
	expect_true(fin_pos(fit$eb$hyper_se))
})

test_that("moving positions come with cell errors of the product; their AR(1) law has no error", {
	skip_on_cran()
	d = dse_panel("normal", strength = 1.2, seed = 11)
	fit = dse_fit(d, R = 1, dynamic_uv = TRUE)
	skip_if(!is.null(fit$eb$product_withheld))
	expect_identical(dim(fit$UVPM_sd), dim(fit$UVPM))
	expect_true(fin_pos(fit$UVPM_sd))
	expect_false(any(c("rho_uv", "sigma_uv") %in% names(fit$eb$hyper_se)))
	expect_true(is.finite(fit$rho_uv) && is.finite(fit$sigma_uv))
	expect_true(fin_pos(fit$U_sd[, 1, ]))
})

test_that("a bipartite count panel reports errors for its level variances and noise variance", {
	skip_on_cran()
	d = dse_panel("poisson", mode = "bipartite", seed = 13)
	fit = dse_fit(d, mode = "bipartite", R = 0, dynamic_ab = TRUE)
	vt = summary(fit)$vc_table
	expect_true(fin_pos(vt[c("va", "vb", "ve"), "Std.Error"]))
	# no dyad pairs in a bipartite panel: rho is not a row
	expect_false("rho" %in% rownames(vt))
	expect_true(fin_pos(fit$a_dynamic_sd) && fin_pos(fit$b_dynamic_sd))
})

test_that("a binary dynamic_beta fit with strong static factors fits its intercept as a path", {
	skip_on_cran()
	d = dse_panel("binary", n = 28, Tt = 4, strength = 1.1, seed = 31)
	fit = dse_fit(d, family = "binary", R = 1, dynamic_beta = "dyad")
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.list(fit$eb$intercept_path))
	expect_true(isTRUE(fit$beta_dynamic_mask[["intercept"]]))
	expect_true(all(is.finite(fit$coef_path["intercept", ])))
	expect_true(fin_pos(fit$eb$se_path["intercept", ]))
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
	expect_false(any(grepl("withheld", capture.output(print(fit)))))
	# the static factors come with cell errors of the product
	expect_identical(dim(fit$UVPM_sd), dim(fit$UVPM))
	expect_true(fin_pos(fit$UVPM_sd))
})

test_that("static factors of a dynamic fit report position errors in the reported gauge", {
	skip_on_cran()
	d = dse_panel("normal", n = 20, Tt = 3, strength = 1, seed = 17)
	fit = dse_fit(d, R = 2, dynamic_ab = TRUE)
	expect_identical(dim(fit$U_sd), dim(fit$U))
	expect_identical(dim(fit$V_sd), dim(fit$V))
	expect_true(fin_pos(fit$U_sd) && fin_pos(fit$V_sd))
	lp = latent_positions(fit)
	expect_true(fin_pos(lp$posterior_sd))
	# undirected: errors of the eigenvectors and the eigenvalues
	d$Y = lapply(d$Y, function(y) { y[lower.tri(y)] = t(y)[lower.tri(y)]; y })
	d$X = lapply(d$X, function(x) { x[, , 1] = (x[, , 1] + t(x[, , 1])) / 2; x })
	fs = dse_fit(d, R = 2, dynamic_ab = TRUE, symmetric = TRUE)
	expect_true(fin_pos(fs$U_sd) && fin_pos(fs$L_sd))
	# the eigenvalue intervals are centred on the posterior mean of each eigenvalue
	expect_identical(dim(fs$L_ci), c(2L, 2L))
	expect_true(all(fs$L_ci[, 1] < fs$L_ci[, 2]))
	expect_identical(fs$V_sd, fs$U_sd)
	expect_true(fin_pos(latent_positions(fs)$posterior_sd))
})

test_that("the eigenvector and eigenvalue errors of a symmetric product match a numerical derivative", {
	set.seed(3)
	n = 8; R = 2; sg = c(1, -1)
	fm = matrix(rnorm(n * R), n, R)
	es = lame:::.ae_eb_eigen_sd(fm, sg, diag(n * R), diag(n))
	at = function(x) {
		fx = matrix(x, n, R)
		eg = eigen(fx %*% (sg * t(fx)), symmetric = TRUE)
		idx = order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
		v = eg$vectors[, idx, drop = FALSE]
		for (r in seq_len(R)) if (sum(v[, r] * es$vectors[, r]) < 0) v[, r] = -v[, r]
		c(v, eg$values[idx])
	}
	h = 1e-6; x0 = as.vector(fm)
	jac = sapply(seq_along(x0), function(k) { e = numeric(length(x0)); e[k] = h; (at(x0 + e) - at(x0 - e)) / (2 * h) })
	sd_num = sqrt(rowSums(jac^2))
	expect_equal(as.vector(es$U), sd_num[seq_len(n * R)], tolerance = 1e-5)
	expect_equal(es$L, sd_num[n * R + seq_len(R)], tolerance = 1e-5)
})

test_that("Louis' identity widens the complete-data covariance by the variance of the statistics", {
	# an AR(1) law with statistics S and no uncertainty in them has the complete-data covariance
	lp = function(th, S) {
		r = th[1]; ls = th[2]
		if (abs(r) >= 1) return(NA_real_)
		a = (1 - r^2) * S[[1]] + S[[2]] - 2 * r * S[[3]] + r^2 * S[[4]]
		-(40 / 2 + 2) * ls - (a / 2 + 0.1) / exp(ls) + 10 / 2 * log1p(-r^2) - (r - 0.8)^2 / (2 * 0.15^2)
	}
	S = c(1.2, 4.1, 3.0, 4.0)
	th = c(0.7, log(0.12))
	V0 = lame:::.lde_louis_cov(lp, th, S, matrix(0, 4, 4))
	Ic = -lame:::.lde_num_hess(function(x) lp(x, S), th)
	expect_equal(V0, solve(Ic), tolerance = 1e-6)
	V1 = lame:::.lde_louis_cov(lp, th, S, diag(0.05, 4))
	expect_true(all(diag(V1) > diag(V0)))
	# the covariance of Gaussian quadratic forms matches its closed form for one coordinate
	m = 0.4; v = 0.3
	expect_equal(lame:::.lde_qf_cov(m, matrix(v), list(matrix(1)))[1, 1], 2 * v^2 + 4 * m^2 * v)
})

test_that("a dynamic fast-route summary prints no smoothing setting and reports an error for every row", {
	skip_on_cran()
	d = dse_panel("normal", strength = 1.2, seed = 11)
	fit = dse_fit(d, R = 1, dynamic_uv = TRUE, dynamic_ab = TRUE, dynamic_beta = "dyad")
	s = summary(fit)
	smooth = c("rho_ab", "sigma_ab", "rho_beta", "rho_uv", "sigma_uv")
	tabs = list(s$coefficients, s$coef_path, s$vc_table, s$hyper_table)
	tabs = tabs[!vapply(tabs, is.null, NA)]
	rows = unlist(lapply(tabs, rownames))
	expect_false(any(grepl(paste0("^(", paste(smooth, collapse = "|"), ")"), rows)))
	expect_false(any(grepl(paste0("^(", paste(smooth, collapse = "|"), ")"), names(fit$eb$hyper_se))))
	for (tb in tabs) expect_true(fin_pos(tb[, "Std.Error"]))
	expect_true(all(grepl("^sigma_beta\\.", rownames(s$hyper_table))) && "sigma_beta.dyad" %in% rownames(s$hyper_table))
	# the values stay on the fit, where prediction, simulation and warm starts read them
	for (k in c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv")) expect_true(is.numeric(fit[[k]]) && all(is.finite(fit[[k]])), info = k)
	expect_true(is.finite(fit$rho_beta[["dyad"]]))
	pr = unlist(predict(fit))
	expect_true(all(is.finite(pr[!is.na(pr)])) && sum(is.finite(pr)) > 0)
	sv = als_start_vals(fit)
	expect_equal(unname(c(sv$rho_uv, sv$sigma_uv, sv$rho_ab, sv$sigma_ab)), unname(c(fit$rho_uv, fit$sigma_uv, fit$rho_ab, fit$sigma_ab)))
	expect_s3_class(simulate(fit, nsim = 1, seed = 1), "ame_als.sim")
	txt = paste(capture.output(print(s)), collapse = " ")
	# neither print shows a smoothing setting, by name or by value
	expect_null(s$smoothing)
	vals = c(fit$rho_ab, fit$sigma_ab, fit$rho_beta[["dyad"]], fit$rho_uv, fit$sigma_uv)
	for (tx in list(txt, paste(capture.output(print(fit)), collapse = " "))) {
		expect_false(grepl("Smoothing chosen|rho_ab|sigma_ab|rho_beta|rho_uv|sigma_uv", tx))
		for (v in vals) expect_false(grepl(format(signif(v, 3)), tx, fixed = TRUE), info = format(v))
	}
	expect_false(grepl("point estimate only|no standard error", txt, ignore.case = TRUE))
})

test_that("the innovation scale of a moving coefficient has a standard error when its persistence is at the boundary", {
	# statistics of a short path whose fitted persistence sits just under one: the joint curvature in
	# (persistence, log scale) does not exist there, and the scale used to be reported without an error
	stats = list(cnt = 3, n_stat = 3, S = c(0.09, 0.085, 0.085, 0.03), VS = diag(1e-4, 4))
	local_mocked_bindings(.lde_beta_group_stats = function(...) stats, .package = "lame")
	ctx = list(
		beta = list(any = TRUE, ng = 1L, kind_g = "ar1", groups = "dyad"),
		betaprior = list(shape = 2, scale = 0.04, rho_mean = 0.8, rho_sd = 0.15),
		pool = list(sigma = FALSE)
	)
	se_at = function(rho) lame:::.lde_hyper_se(ctx, list(hyp = list(beta = list(s2 = 0.0225, rho = rho))), NULL, NULL)
	edge = se_at(1 - 1e-9)
	expect_named(edge, "sigma_beta.dyad")
	expect_true(is.finite(edge) && edge > 0)
	# away from the boundary the joint curvature is used, and the two agree to the order of the error
	inner = se_at(0.8)
	expect_true(is.finite(inner) && inner > 0)
	expect_lt(abs(log(edge / inner)), log(2))
})
