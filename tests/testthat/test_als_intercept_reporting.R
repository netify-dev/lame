library(lame)

# an intercept path, the rw1 path of dynamic_ab = TRUE or one from dynamic_beta, is reported by every ALS
# family in coef_path, coefficients, mu, coef(), summary() and predictions, binary fits with strong static factors too

air_panel = function(family, mode = "unipartite", n = 16, Tt = 4, seed = 31) {
	set.seed(seed)
	nr = n
	nc = if (mode == "bipartite") 12 else n
	rn = sprintf("r%02d", seq_len(nr))
	cn = if (mode == "bipartite") sprintf("c%02d", seq_len(nc)) else rn
	U = matrix(rnorm(nr, 0, 1.2), nr, 1); V = matrix(rnorm(nc, 0, 1.2), nc, 1)
	al = rnorm(nr, 0, 0.5); be = rnorm(nc, 0, 0.5)
	shift = seq(-0.5, 0.5, length.out = Tt)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		al = al + rnorm(nr, 0, 0.1); be = be + rnorm(nc, 0, 0.1)
		x = matrix(rnorm(nr * nc), nr, nc)
		eta = 0.3 + shift[t] + outer(al, be, "+") + 0.6 * x + U %*% t(V)
		y = switch(family,
			normal = eta + matrix(rnorm(nr * nc), nr, nc),
			binary = 1 * (eta + matrix(rnorm(nr * nc), nr, nc) > 0),
			poisson = matrix(rpois(nr * nc, exp(0.5 * eta - 0.5)), nr, nc))
		if (mode == "unipartite") diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(nr, nc, 1), dimnames = list(rn, cn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

air_fit = function(family, mode = "unipartite", ...) {
	d = air_panel(family, mode)
	args = list(d$Y, Xdyad = d$X, family = family, method = "als", verbose = FALSE, ...)
	if (mode == "bipartite") {
		args$mode = "bipartite"; args$R_row = 1; args$R_col = 1
	} else args$R = 1
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

air_specs = list(
	dynamic_ab = list(dynamic_ab = TRUE),
	dynamic_beta = list(dynamic_beta = "intercept"),
	dynamic_beta_rw1 = list(dynamic_beta = "intercept", dynamic_beta_kind = "rw1")
)

for (spec in names(air_specs)) for (fam in c("binary", "normal", "poisson")) for (mode in c("unipartite", "bipartite")) {
	if (spec != "dynamic_ab" && fam != "binary") next
	test_that(paste("a", mode, fam, "ALS fit with", spec, "reports its intercept path"), {
		skip_on_cran()
		fit = do.call(air_fit, c(list(fam, mode), air_specs[[spec]]))
		expect_s3_class(fit, "lame_dynamic_als")
		expect_true(isTRUE(fit$beta_dynamic_mask[["intercept"]]))
		# above the signal from which a static binary intercept is withheld
		if (fam == "binary") expect_gt(fit$eb$information[["signal_per_loading"]], 1.5)
		expect_null(fit$eb$intercept_withheld)
		path = fit$coef_path["intercept", ]
		expect_true(all(is.finite(path)))
		expect_identical(unname(path), unname(fit$BETA[1, "intercept", ]))
		expect_identical(fit$beta_dynamic_per_t, fit$coef_path)
		expect_true(is.finite(fit$coefficients[["intercept"]]))
		expect_identical(fit$mu[["intercept"]], fit$coefficients[["intercept"]])
		# beta holds the slopes, as for every family; the intercept is mu
		expect_identical(fit$beta, fit$coefficients[names(fit$coefficients) != "intercept"])
		expect_true(all(is.finite(fit$eb$se)))
		expect_true(all(is.finite(coef(fit)["intercept", ])))
		expect_true(all(is.finite(summary(fit)$coefficients["intercept", ])))
		expect_false(any(grepl("withheld", capture.output(print(fit)))))
		expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	})
}
