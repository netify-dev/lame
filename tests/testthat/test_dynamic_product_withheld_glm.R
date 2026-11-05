library(lame)

# a dynamic binary or count fit withholds the latent product the data do not determine: the statistic is
# the weighted squared double-centred product per loading (fit$eb$information), against .lde_glm_product_threshold

dpg_panel = function(family, suv, cfg = "dab", n = 30, Tt = 4L, seed = 71) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n))
	ar = function(k, v, rho) { x = matrix(0, k, Tt); x[, 1] = rnorm(k, 0, sqrt(v)); for (t in 2:Tt) x[, t] = rho * x[, t - 1] + rnorm(k, 0, sqrt(v * (1 - rho^2))); x }
	A = ar(n, 0.5, 0.8); B = ar(n, 0.5, 0.8)
	U = if (cfg == "duv") ar(n, suv, 0.9) else matrix(rnorm(n, 0, sqrt(suv)), n, Tt)
	V = if (cfg == "duv") ar(n, suv, 0.9) else matrix(rnorm(n, 0, sqrt(suv)), n, Tt)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n)
		eta = (if (family == "binary") 0.3 else -0.2) + outer(A[, t], B[, t], "+") + tcrossprod(U[, t], V[, t]) + 0.6 * x
		y = if (family == "binary") 1 * (eta + matrix(rnorm(n * n), n) > 0) else matrix(rpois(n * n, exp(eta)), n)
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(rn, rn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, cfg = cfg)
}
dpg_fit = function(d, family) {
	args = list(d$Y, Xdyad = d$X, family = family, method = "als", R = 1, dynamic_ab = TRUE, verbose = FALSE)
	if (d$cfg == "duv") args$dynamic_uv = TRUE
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

test_that("a binary or count fit with a weakly identified product withholds it, and one with a strong product keeps it", {
	skip_on_cran()
	thr = lame:::.lde_glm_product_threshold
	for (family in c("binary", "poisson")) {
		weak = dpg_fit(dpg_panel(family, 0.15, "duv"), family)
		strong = dpg_fit(dpg_panel(family, 0.9, "dab"), family)
		info = paste(family)
		expect_lt(weak$eb$information[["signal_per_loading"]], thr)
		expect_gte(strong$eb$information[["signal_per_loading"]], thr)
		expect_false(is.null(weak$eb$product_withheld), info = info)
		expect_null(strong$eb$product_withheld, info = info)
		expect_identical(weak$eb$product_withheld$family, family)
		expect_true(all(is.na(weak$UVPM)), info = info)
		expect_true(all(is.finite(unlist(strong$UVPM))), info = info)
		# withholding leaves the coefficients, fitted values and predictions alone
		expect_true(all(is.finite(weak$coefficients)))
		expect_true(all(is.finite(unlist(fitted(weak)))[!is.na(unlist(weak$EZ))]))
		expect_error(uv_plot(weak), class = "lame_product_withheld")
		# the option turns withholding off and reports the product
		old = options(lame.eb_withhold_product = FALSE)
		on.exit(options(old), add = TRUE)
		kept = dpg_fit(dpg_panel(family, 0.15, "duv"), family)
		options(old)
		expect_null(kept$eb$product_withheld)
		expect_equal(unname(kept$coefficients), unname(weak$coefficients))
		expect_equal(kept$eb$information[["signal_per_loading"]], weak$eb$information[["signal_per_loading"]])
	}
})
