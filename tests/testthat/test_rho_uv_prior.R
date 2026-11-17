# default prior of the persistence rho_uv of dynamic latent positions and the gauge of its draws
# (?lame, prior$rho_uv_mean and prior$rho_uv_sd; R/uv_rho_flip.R)

rup_default = c(mean = 0.7, sd = 0.3)

# the posterior of rho given known AR(1) paths with innovation variance s2, on a grid over the
# support, folded to |rho| as the stored draws are
rup_posterior = function(x, s2, prior_mean, prior_sd, grid) {
	TT = ncol(x)
	sxy = sum(x[, -1] * x[, -TT])
	syy = sum(x[, -TT]^2)
	stt = sum(x[, -1]^2)
	s1 = sum(x[, 1]^2)
	ll = 0.5 * nrow(x) * log(1 - grid^2) - (1 - grid^2) * s1 / (2 * s2) -
		(stt - 2 * grid * sxy + grid^2 * syy) / (2 * s2)
	lp = ll + dnorm(grid, prior_mean, prior_sd, log = TRUE)
	w = exp(lp - max(lp))
	w_fold = w[grid > 0] + rev(w[grid < 0])
	w_fold / sum(w_fold)
}

test_that("the default persistence prior is calibrated where the paths are short", {
	skip_on_cran()
	# 40 coordinate paths (20 actors at rank two) over four periods: the weakest information of the
	# calibration study's designs; 90 percent intervals of |rho| from the folded posterior
	grid = seq(-0.9895, 0.9895, by = 0.001)
	a = grid[grid > 0]
	set.seed(2026)
	for (rho in c(0.3, 0.6, 0.9)) {
		res = replicate(300, {
			x = matrix(0, 40, 4)
			x[, 1] = rnorm(40)
			for (t in 2:4) x[, t] = rho * x[, t - 1] + sqrt(1 - rho^2) * rnorm(40)
			w = rup_posterior(x, 1 - rho^2, rup_default[["mean"]], rup_default[["sd"]], grid)
			cw = cumsum(w)
			m = sum(w * a)
			s = sqrt(sum(w * (a - m)^2))
			c(cover = a[which(cw >= 0.05)[1]] <= rho && rho <= a[which(cw >= 0.95)[1]], z = (m - rho) / s)
		})
		# Monte Carlo sd of the coverage is 0.017 at 300 panels
		expect_gt(mean(res["cover", ]), 0.85)
		expect_lt(mean(res["cover", ]), 0.95)
		expect_lt(abs(mean(res["z", ])), 0.5)
	}
})

test_that("the period signs of symmetric paths keep each product and follow their conditional", {
	skip_on_cran()
	set.seed(5)
	n = 6
	N = 4
	rho = 0.6
	U = array(0, c(n, 2, N))
	U[, , 1] = rnorm(2 * n)
	for (t in 2:N) U[, , t] = rho * U[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n)
	L = diag(c(1.5, -0.8))
	draw = lame:::.lame_uv_sym_sign_draw(U, rho)
	for (t in 1:N) expect_equal(draw$U[, , t] %*% L %*% t(draw$U[, , t]), U[, , t] %*% L %*% t(U[, , t]))
	expect_equal(abs(draw$U), abs(U))
	# the share of draws keeping the relative sign of coordinate 1 between periods 1 and 2
	k = rho / (1 - rho^2) * sum(U[, 1, 1] * U[, 1, 2])
	keep = replicate(4000, {
		d = lame:::.lame_uv_sym_sign_draw(U, rho)$U
		sign(sum(d[, 1, 1] * d[, 1, 2])) == sign(sum(U[, 1, 1] * U[, 1, 2]))
	})
	p_keep = plogis(2 * k)
	expect_lt(abs(mean(keep) - p_keep), 4 * sqrt(p_keep * (1 - p_keep) / 4000))
	# without persistence or with a single period nothing is drawn from the paths
	expect_identical(lame:::.lame_uv_sym_sign_draw(U[, , 1, drop = FALSE], rho)$U, U[, , 1, drop = FALSE])
	expect_identical(lame:::.lame_uv_sym_sign_draw(U, NA_real_)$U, U)
})

skip_on_cran()

# a directed normal panel with persistent positions
rup_panel = function(n = 12, TT = 5, rho = 0.8, seed = 41) {
	set.seed(seed)
	U = V = array(0, c(n, 2, TT))
	U[, , 1] = rnorm(2 * n)
	V[, , 1] = rnorm(2 * n)
	for (t in 2:TT) {
		U[, , t] = rho * U[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n)
		V[, , t] = rho * V[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n)
	}
	nm = sprintf("a%02d", 1:n)
	Y = array(0, c(n, n, TT), dimnames = list(nm, nm, paste0("t", 1:TT)))
	for (t in 1:TT) {
		Y[, , t] = U[, , t] %*% t(V[, , t]) + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y[, , t]) = NA
	}
	Y
}

rup_fit = function(Y, prior = list(), nscan = 1000, ...) {
	suppressWarnings(suppressMessages(lame(
		Y,
		R = 2,
		family = "normal",
		dynamic_uv = TRUE,
		prior = prior,
		burn = 200,
		nscan = nscan,
		odens = 5,
		seed = 17,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		...
	)))
}

test_that("every dynamic_uv route takes the same default and explicit settings are kept", {
	Y = rup_panel(n = 8, TT = 4)
	args = list(
		Y = Y,
		R = 1,
		family = "normal",
		dynamic_uv = TRUE,
		burn = 0,
		nscan = 8,
		odens = 1,
		seed = 116,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
	fits = list(
		directed = do.call(lame, args),
		binary = do.call(lame, modifyList(args, list(Y = 1 * (Y > 0), family = "binary"))),
		symmetric = do.call(lame, modifyList(args, list(
			Y = array(apply(Y, 3, function(y) (y + t(y)) / 2), dim(Y), dimnames(Y)),
			symmetric = TRUE
		))),
		bipartite = do.call(lame, modifyList(args, list(Y = Y[, 1:6, ], mode = "bipartite")))
	)
	for (nm in names(fits)) {
		expect_equal(fits[[nm]]$prior$rho_uv_mean, rup_default[["mean"]], info = nm)
		expect_equal(fits[[nm]]$prior$rho_uv_sd, rup_default[["sd"]], info = nm)
		expect_identical(fits[[nm]]$refit_options$prior, list(), info = nm)
	}
	expect_equal(fits$symmetric$uv_sym_sign_move[["visits"]], 8)
	expect_null(fits$directed$uv_sym_sign_move)
	explicit = do.call(lame, modifyList(args, list(prior = list(rho_uv_mean = 0.9, rho_uv_sd = 0.1))))
	expect_equal(explicit$prior$rho_uv_mean, 0.9)
	expect_equal(explicit$prior$rho_uv_sd, 0.1)
	# a mean set alone keeps the default sd
	partial = do.call(lame, modifyList(args, list(prior = list(rho_uv_mean = 0.95))))
	expect_equal(partial$prior$rho_uv_sd, rup_default[["sd"]])
})

test_that("the default fit is the fit that states the default prior", {
	Y = rup_panel()
	implicit = rup_fit(Y, nscan = 200)
	stated = rup_fit(Y, prior = as.list(rup_default), nscan = 200)
	for (nm in c("BETA", "VC", "rho_uv", "sigma_uv", "UVPM", "U", "V")) {
		expect_identical(implicit[[nm]], stated[[nm]], info = nm)
	}
})

test_that("stored draws carry rho_uv >= 0 with the positions of the same mirror", {
	Y = rup_panel()
	ref = rup_fit(Y)
	# a prior symmetric about zero takes the mirror move at every visit, so the sampler state sits at
	# negative rho_uv half the time; the stored draws are each draw's mirror with rho_uv >= 0
	sym = rup_fit(Y, prior = list(rho_uv_mean = 0, rho_uv_sd = 0.5))
	expect_gt(sym$uv_rho_flip_move[["moved"]], 0.4 * sym$uv_rho_flip_move[["visits"]])
	expect_true(all(sym$rho_uv >= 0))
	expect_true(all(ref$rho_uv >= 0))
	# each stored path cube turns with its draw of rho_uv: the lag-one products are positive
	lag_one = apply(sym$U_draws, 4, function(u) sum(u[, , -1] * u[, , -dim(u)[3]]))
	expect_true(all(lag_one > 0))
	# the data see |rho_uv| only, so the two priors report the same persistence here
	se = sqrt(var(sym$rho_uv) / 40 + var(ref$rho_uv) / 40)
	expect_lt(abs(mean(sym$rho_uv) - mean(ref$rho_uv)), 4 * se)
	# the mirror leaves every product, so the posterior mean products agree as closely
	expect_lt(max(abs(sym$UVPM - ref$UVPM)), 0.25 * max(abs(ref$UVPM)))
})
