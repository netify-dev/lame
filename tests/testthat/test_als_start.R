library(lame)
skip_on_cran()

set.seed(5)
n = 12; TT = 4; nm = sprintf("a%02d", 1:n)
mk = function(fam) lapply(1:TT, function(t) {
	z = matrix(rnorm(n * n), n, n)
	y = switch(
		fam,
		poisson = matrix(rpois(n * n, exp(0.3 + 0.5 * z)), n, n),
		normal = z,
		ordinal = matrix(cut(z, c(-Inf, -0.5, 0.5, Inf), labels = FALSE), n, n)
	)
	diag(y) = NA; dimnames(y) = list(nm, nm); y
})
quiet = function(expr) suppressWarnings(suppressMessages(expr))
ctrl = list(nscan = 120, burn = 60, odens = 4, seed = 3, gof = FALSE, verbose = FALSE, plot = FALSE)

test_that("ar1 dynamic-position fits start from a dynamic ALS fit by default", {
	f = quiet(do.call(lame, c(list(mk("poisson"), family = "poisson", R = 2, dynamic_uv = TRUE), ctrl)))
	expect_true(isTRUE(f$als_start))
	expect_true(all(is.finite(f$rho_uv)))
	fb = quiet(do.call(
		lame,
		c(
			list(
				lapply(
					1:TT,
					function(t) { m = matrix(rpois(10 * 7, 1.2), 10, 7); dimnames(m) = list(
						paste0("r", 1:10),
						paste0("c", 1:7)
					); m }
				),
				family = "poisson",
				mode = "bipartite",
				R_row = 1,
				R_col = 1,
				dynamic_uv = TRUE
			),
			ctrl
		)
	))
	# bipartite fits keep the cold start
	expect_false(isTRUE(fb$als_start))
})

test_that("the warm start is skipped where it does not apply", {
	f_off = quiet(
		do.call(
			lame,
			c(list(mk("normal"), family = "normal", R = 1, dynamic_uv = TRUE, prior = list(als_start = FALSE)), ctrl)
		)
	)
	expect_false(isTRUE(f_off$als_start))
	f_snap = quiet(
		do.call(
			lame,
			c(list(mk("normal"), family = "normal", R = 1, dynamic_uv = TRUE, dynamic_uv_kind = "snap"), ctrl)
		)
	)
	expect_false(isTRUE(f_snap$als_start))
	f_ord = quiet(do.call(lame, c(list(mk("ordinal"), family = "ordinal", R = 1, dynamic_uv = TRUE), ctrl)))
	expect_false(isTRUE(f_ord$als_start))
	# a binary fit starts from the dynamic ALS fit too, and the switch turns it off
	yb = lapply(mk("normal"), function(m) { y = 1 * (m > 0.3); diag(y) = NA; y })
	f_bin = quiet(do.call(lame, c(list(yb, family = "binary", R = 1, dynamic_uv = TRUE), ctrl)))
	expect_true(isTRUE(f_bin$als_start))
	expect_true(all(is.finite(f_bin$rho_uv)))
	f_bin_off = quiet(do.call(
		lame,
		c(
			list(
				yb,
				family = "binary",
				R = 1,
				dynamic_uv = TRUE,
				prior = list(als_start = FALSE)
			),
			ctrl
		)
	))
	expect_false(isTRUE(f_bin_off$als_start))
	Ys = lapply(mk("normal"), function(m) (m + t(m)) / 2)
	f_sym = quiet(do.call(lame, c(list(Ys, family = "normal", R = 1, symmetric = TRUE, dynamic_uv = TRUE), ctrl)))
	expect_false(isTRUE(f_sym$als_start))
	f0 = quiet(do.call(lame, c(list(mk("normal"), family = "normal", R = 1, dynamic_uv = TRUE), ctrl)))
	f_sv = quiet(
		do.call(
			lame,
			c(list(mk("normal"), family = "normal", R = 1, dynamic_uv = TRUE, start_vals = f0$start_vals), ctrl)
		)
	)
	expect_false(isTRUE(f_sv$als_start))
	f_static = quiet(do.call(lame, c(list(mk("normal"), family = "normal", R = 1), ctrl)))
	expect_null(f_static$als_start)
})

test_that("the warm start works when lame() is called from another function", {
	wrap = function(dat, ...) lame(dat, family = "poisson", R = 1, dynamic_uv = TRUE, ...)
	f = quiet(do.call(wrap, c(list(mk("poisson")), ctrl)))
	expect_true(isTRUE(f$als_start))
})
